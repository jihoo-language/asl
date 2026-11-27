use std::rc::Rc;

use crate::error::{Error, Result};
use crate::integer_smith::Track;
use crate::model::{Model, Value};
use crate::rv_lens::gcd;
use crate::rv_qp_quad::Wide;
use crate::syntax::{Label, ParsedExpr};

use super::{matrix, Basis};
use linear::{self as lin, add, mul, Matrix};

mod congruence;
mod linear;
mod polynomial;

pub(crate) enum Command {
    Define(Label, Label, i128, Vec<Vec<i128>>),
    Describe(Label),
    Congruence(Label, usize, usize, usize),
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<Option<String>> {
        match self {
            Self::Define(name, base, coefficient, factors) => {
                let basis = model.fricke_ss(&base)?;
                let order = Order::build(model, base, &basis, coefficient, factors)
                    .map_err(|error| error.at(name.span))?;
                model.define_fricke_order(name, order)?;
                Ok(None)
            }
            Self::Describe(name) => Ok(Some(model.fricke_order(&name)?.describe())),
            Self::Congruence(name, scalar, block, ell) => {
                let order = model.fricke_order(&name)?;
                order.congruence(model, scalar, block, ell).map(Some)
            }
        }
    }
}

pub(crate) enum Expression {
    Basis(Label, usize, usize),
    Project(Label, usize, Box<ParsedExpr>),
    Extract(bool, Box<ParsedExpr>),
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Fraction {
    pub(crate) numerator: Vec<i128>,
    pub(crate) denominator: i128,
}

impl Fraction {
    pub(super) fn new(mut numerator: Vec<i128>, mut denominator: i128) -> Self {
        let divisor = numerator
            .iter()
            .fold(denominator as u128, |d, n| gcd(d, n.unsigned_abs()))
            as i128;
        for coefficient in &mut numerator {
            *coefficient /= divisor;
        }
        denominator /= divisor;
        Self {
            numerator,
            denominator,
        }
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(5 + 2 * self.numerator.len())?;
        Ok(Value::FrickeSSQDiv(Rc::new(self)))
    }
}

pub(crate) struct Order {
    pub(crate) base: Label,
    coefficient: i128,
    blocks: Vec<Block>,
    smith: Vec<i128>,
    index: i128,
    operator_names: Vec<&'static str>,
}

#[derive(Clone, Copy, Eq, PartialEq)]
enum FactorKind {
    Irreducible,
    Reducible,
    Uncertified,
}

struct Block {
    factor: Vec<i128>,
    kind: FactorKind,
    projector: Matrix,
    denominator: i128,
    kernel: Vec<Vec<i128>>,
    restrictions: Vec<Matrix>,
    images: Option<Vec<Fraction>>,
}

impl Order {
    fn build(
        model: &mut Model,
        base: Label,
        basis: &Basis,
        coefficient: i128,
        factors: Vec<Vec<i128>>,
    ) -> Result<Self> {
        let n = basis.rank();
        if n > 16 {
            return Err(Error::plain(
                "Brandt Hecke orders support divisor rank at most 16",
            ));
        }
        if coefficient.unsigned_abs() > 2 || (basis.p == 5 && coefficient != 0) {
            return Err(Error::plain(
                "Brandt generator must be B3 + c B5 with |c| ≤ 2 and c=0 in characteristic five",
            ));
        }
        if factors.is_empty()
            || factors.len() > n
            || factors.iter().any(|f| {
                f.len() < 2
                    || f.len() > n + 1
                    || f.last() != Some(&1)
                    || f.iter().any(|v| v.unsigned_abs() > 1_000_000_000_000)
            })
        {
            return Err(Error::plain(
                "Brandt factors must be positive-degree monic integral polynomials with coefficients of absolute value at most 10^12",
            ));
        }
        let degree: usize = factors.iter().map(|f| f.len() - 1).sum();
        if degree > n {
            return Err(Error::plain(
                "Brandt factor degrees exceed the divisor rank",
            ));
        }
        let three = basis.operator(model, 3)?;
        let five = if basis.p == 5 {
            None
        } else {
            Some(basis.operator(model, 5)?)
        };
        let mut operator_refs = vec![three];
        operator_refs.extend(five);
        let (weights, _) = matrix::balance(model, &operator_refs, &basis.frobenius)?;
        let mut generator = lin::transpose(&three.0);
        if coefficient != 0 {
            for (a, &b) in generator
                .iter_mut()
                .flatten()
                .zip(lin::transpose(&five.unwrap().0).iter().flatten())
            {
                *a = add(*a, mul(coefficient, b)?)?;
            }
        }
        let mut frobenius = vec![vec![0; n]; n];
        for (j, &i) in basis.frobenius.iter().enumerate() {
            frobenius[i][j] = 1;
        }
        let mut operator_names = vec!["B3"];
        let mut operators = vec![lin::transpose(&three.0)];
        if let Some(five) = five {
            operator_names.push("B5");
            operators.push(lin::transpose(&five.0));
        }
        operator_names.push("Frob");
        operators.push(frobenius);
        certify_operators(model, &operators, &weights)?;
        let powers = lin::powers(model, &generator, degree)?;
        let total = factors
            .iter()
            .try_fold(vec![1], |f, g| polynomial::multiply(&f, g))?;
        if lin::combine(model, &total, &powers)?
            .iter()
            .flatten()
            .any(|&a| a != 0)
        {
            return Err(Error::plain(
                "supplied Brandt factors do not annihilate the chosen operator",
            ));
        }
        let lower = powers[..degree].iter().map(lin::flat).collect::<Vec<_>>();
        if lin::solve(model, &lower, &vec![0; n * n])?.is_none() {
            return Err(Error::plain(
                "the supplied product is not the minimal polynomial",
            ));
        }
        if !polynomial::squarefree(model, &total)? {
            return Err(Error::plain(
                "Brandt factors must be squarefree and pairwise coprime",
            ));
        }
        let mut blocks = Vec::with_capacity(factors.len());
        for factor in factors {
            let (projector, denominator) = polynomial::projector(model, &total, &factor, &powers)?;
            let annihilator = lin::combine(model, &factor, &powers)?;
            let frame = lin::kernel_frame(model, &annihilator)
                .map_err(|error| Error::plain(format!("Brandt block kernel: {}", error.message)))?;
            let kernel = lin::kernel(&frame, n);
            if kernel.is_empty() {
                return Err(certificate("nonempty factor kernel"));
            }
            for v in &kernel {
                if lin::apply(model, &annihilator, v)?.iter().any(|&c| c != 0)
                    || lin::apply(model, &projector, v)?
                        != v.iter()
                            .map(|&c| mul(c, denominator))
                            .collect::<Result<Vec<_>>>()?
                {
                    return Err(certificate(
                        "saturated eigenbasis and image of CRT projector",
                    ));
                }
            }
            let restricted_generator = lin::restrict(model, &generator, &frame, &kernel)?;
            let restrictions = operators
                .iter()
                .map(|op| lin::restrict(model, op, &frame, &kernel))
                .collect::<Result<Vec<_>>>()?;
            let kind = polynomial::factor_kind(model, &factor, 4 + 6 * coefficient.abs())?;
            let images = if kernel.len() + 1 == factor.len() {
                Some(simple_images(model, &restricted_generator, &restrictions)?)
            } else {
                None
            };
            blocks.push(Block {
                factor,
                kind,
                projector,
                denominator,
                kernel,
                restrictions,
                images,
            });
        }
        certify_projectors(model, &blocks, &weights, &operators)?;
        let joined: Matrix = blocks
            .iter()
            .flat_map(|block| block.kernel.iter().cloned())
            .collect();
        if joined.len() != n {
            return Err(certificate("full-rank direct sum of saturated blocks"));
        }
        let smith = lin::smith(model, &lin::transpose(&joined), Track::None)
            .map_err(|error| {
                Error::plain(format!("Brandt congruence cokernel: {}", error.message))
            })?
            .factors;
        if smith.len() != n {
            return Err(certificate("finite gluing cokernel"));
        }
        let index = smith.iter().try_fold(1, |index, &d| mul(index, d))?;
        let cells = blocks
            .iter()
            .map(|b| {
                let m = b.kernel.len();
                30 + 4 * b.factor.len()
                    + 2 * n * n
                    + 2 * n * m
                    + operators.len() * (3 * m + 2 * m * m)
            })
            .sum::<usize>();
        model.allocate(20 + cells)?;
        Ok(Self {
            base,
            coefficient,
            blocks,
            smith,
            index,
            operator_names,
        })
    }

    fn block(&self, i: usize) -> Result<&Block> {
        self.blocks
            .get(i)
            .ok_or_else(|| Error::plain("Brandt factor index is out of range"))
    }

    pub(crate) fn basis(&self, block: usize, column: usize) -> Result<Vec<i128>> {
        self.block(block)?
            .kernel
            .get(column)
            .cloned()
            .ok_or_else(|| Error::plain("Brandt block basis index is out of range"))
    }

    pub(crate) fn check_block(&self, block: usize) -> Result<()> {
        self.block(block).map(|_| ())
    }

    pub(crate) fn project(
        &self,
        model: &mut Model,
        block: usize,
        input: &[i128],
    ) -> Result<Fraction> {
        let block = self.block(block)?;
        model.spend(9 * input.len() * input.len())?;
        let divisor = block.denominator as u128;
        let mut result = Vec::with_capacity(input.len());
        // Reduce the common denominator before enforcing the output range.
        let mut common = block.denominator as u128;
        for row in &block.projector {
            let amount = row
                .iter()
                .zip(input)
                .try_fold(Wide::<9>::from_i128(0), |amount, (&a, &b)| {
                    amount.add(Wide::product(a, b)).map_err(|_| arithmetic())
                })?;
            let (_, remainder) = amount.magnitude.divide(model, divisor)?;
            common = gcd(common, remainder);
            result.push(amount);
        }
        let common = i128::try_from(common).map_err(|_| arithmetic())?;
        let mut numerator = Vec::with_capacity(result.len());
        for amount in result {
            numerator.push(
                amount
                    .divide_exact(model, Wide::from_i128(common))?
                    .and_then(Wide::signed)
                    .ok_or_else(arithmetic)?,
            );
        }
        Ok(Fraction {
            numerator,
            denominator: block.denominator / common,
        })
    }

    fn describe(&self) -> String {
        let generator = if self.coefficient == 0 {
            "B3".to_string()
        } else {
            format!("B3 + ({}) B5", self.coefficient)
        };
        let mut lines = vec![format!(
            "T = {generator}; saturated-block cokernel Smith = {:?}; index = {}",
            self.smith, self.index
        )];
        for (i, block) in self.blocks.iter().enumerate() {
            let kind = match block.kind {
                FactorKind::Irreducible => "certified irreducible",
                FactorKind::Reducible => "certified reducible",
                FactorKind::Uncertified => "irreducibility unproved",
            };
            lines.push(format!(
                "block {i}: factor (ascending) = {:?}; {kind}; rank = {}; projector denominator = {}; basis columns = {:?}",
                block.factor, block.kernel.len(), block.denominator, block.kernel,
            ));
            if let Some(images) = &block.images {
                for (name, fraction) in self.operator_names.iter().zip(images) {
                    lines.push(format!(
                        "{name} on block {i} = {:?} / {}",
                        fraction.numerator, fraction.denominator
                    ));
                }
            } else {
                lines.push(format!(
                    "block {i} has multiplicity; a single polynomial image is uncertified"
                ));
            }
            for (name, matrix) in self.operator_names.iter().zip(&block.restrictions) {
                lines.push(format!("{name} integral matrix on block {i} = {matrix:?}"));
            }
        }
        lines.join("\n")
    }
}

fn certify_operators(model: &mut Model, operators: &[Matrix], weights: &[i128]) -> Result<()> {
    for (index, operator) in operators.iter().enumerate() {
        certify_adjoint(operator, weights)?;
        for other in &operators[..index] {
            if lin::multiply(model, operator, other)? != lin::multiply(model, other, operator)? {
                return Err(certificate("the specified operators commute"));
            }
        }
    }
    Ok(())
}

fn certify_adjoint(matrix: &Matrix, weights: &[i128]) -> Result<()> {
    for (i, row) in matrix.iter().enumerate() {
        for (j, &value) in row.iter().enumerate() {
            if mul(weights[j], value)? != mul(weights[i], matrix[j][i])? {
                return Err(certificate("weighted self-adjointness on divisors"));
            }
        }
    }
    Ok(())
}

fn certify_projectors(
    model: &mut Model,
    blocks: &[Block],
    weights: &[i128],
    operators: &[Matrix],
) -> Result<()> {
    let common = blocks
        .iter()
        .try_fold(1, |d, b| lin::lcm(d, b.denominator))?;
    let n = weights.len();
    let mut sum = vec![0i128; n * n];
    for (i, block) in blocks.iter().enumerate() {
        let numerator = &block.projector;
        certify_adjoint(numerator, weights)?;
        let square = lin::flat(&lin::multiply(model, numerator, numerator)?);
        for ((sum, square), entry) in sum.iter_mut().zip(square).zip(lin::flat(numerator)) {
            if square != mul(block.denominator, entry)? {
                return Err(certificate("CRT idempotents"));
            }
            *sum = add(*sum, mul(common / block.denominator, entry)?)?;
        }
        for other in &blocks[..i] {
            if lin::multiply(model, numerator, &other.projector)?
                .iter()
                .flatten()
                .any(|&a| a != 0)
            {
                return Err(certificate("orthogonal CRT idempotents"));
            }
        }
        for operator in operators {
            if lin::multiply(model, numerator, operator)?
                != lin::multiply(model, operator, numerator)?
            {
                return Err(certificate(
                    "Hecke and Frobenius commute with each projector",
                ));
            }
        }
    }
    if sum
        != lin::flat(&lin::identity(n))
            .into_iter()
            .map(|e| e * common)
            .collect::<Vec<_>>()
    {
        return Err(certificate("CRT idempotents sum to identity"));
    }
    Ok(())
}

fn simple_images(
    model: &mut Model,
    generator: &Matrix,
    operators: &[Matrix],
) -> Result<Vec<Fraction>> {
    let powers = lin::powers(model, generator, generator.len() - 1)?;
    let columns: Vec<_> = powers.iter().map(lin::flat).collect();
    operators
        .iter()
        .map(|operator| {
            let coefficients = lin::solve(model, &columns, &lin::flat(operator))?
                .ok_or_else(|| certificate("simple block is cyclic, with commutant Q[T]/f"))?;
            let (numerator, denominator) = lin::clear(&coefficients)?;
            if lin::flat(&lin::combine(model, &numerator, &powers)?)
                != lin::flat(operator)
                    .into_iter()
                    .map(|c| mul(c, denominator))
                    .collect::<Result<Vec<_>>>()?
            {
                return Err(certificate(
                    "exact matrix evaluation of the polynomial image",
                ));
            }
            Ok(Fraction::new(numerator, denominator))
        })
        .collect()
}

fn arithmetic() -> Error {
    Error::plain("exact Brandt Hecke order arithmetic exceeds the signed 128-bit range")
}

fn certificate(identity: &str) -> Error {
    Error::plain(format!("Brandt Hecke order certificate failed: {identity}"))
}
