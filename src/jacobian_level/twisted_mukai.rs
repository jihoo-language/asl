use std::rc::Rc;

use super::brauer::{identity, Exterior};
use super::cm_cycles::{Cycle, Cycles};
use super::linking::{add_vectors, image, intersection, product, scaled, Fixed};
use super::transfers::{action, extension_ratio, Direction, Quotient, Transfer};
use crate::bp_smith;
use crate::error::{Error, Result};
use crate::hw_algebra::{self, Matrix};
use crate::hyperelliptic::Residues;
use crate::model::Model;
use crate::prime_linear;
use crate::rv_unram::{add_mod, mul_mod, sub_mod};
use crate::syntax::{Label, ParsedExpr};

pub(crate) type Vector = [u64; 8];

pub(crate) enum Statement {
    Define(Label, Label, usize, Box<ParsedExpr>),
    Shift(Label, Label, [i128; 6]),
    Algebraic(Label, Label, Box<ParsedExpr>),
    Extend(Label, Label, Label),
    Marking(Label, Label, Label),
    Data(Label, usize),
    Pair(Box<[ParsedExpr; 2]>),
}

pub(crate) enum Expression {
    Vector(Label, [i128; 8]),
    Witness(Label),
    Frobenius(usize, Box<ParsedExpr>),
    Transport(Label, Box<ParsedExpr>),
}

pub(crate) struct Chart {
    quotient: Rc<Quotient>,
    precision: usize,
    modulus: u64,
    denominator: u64,
    numerator: Fixed,
    frobenius: Operator,
    period: u64,
    relation: Option<Box<Relation>>,
}

struct Relation {
    parent: Rc<Chart>,
    forward: Operator,
    reverse: Option<Operator>,
}

#[derive(Clone, Copy, Eq, PartialEq)]
pub(crate) struct Operator {
    linear: Exterior,
    shift: Fixed,
}

impl Chart {
    pub(crate) fn new(
        model: &mut Model,
        quotient: Rc<Quotient>,
        precision: usize,
        numerator: Fixed,
    ) -> Result<Self> {
        let pairing = &quotient.pairing;
        let precision_needed = quotient.level.checked_add(precision)
            .filter(|&level| precision > 0 && level <= pairing.ring.precision)
            .ok_or_else(|| Error::plain("Brauer Mukai chart needs positive precision m with s+m at most the marking height"))?;
        let modulus = pairing.modulus(precision)?;
        let denominator = pairing.modulus(quotient.level)?;
        let working_modulus = pairing.modulus(precision_needed)?;
        let period = quotient.order(quotient.project(model, numerator.map(|x| x % denominator))?);
        model.spend(72)?;
        let difference = image(&pairing.difference, numerator, working_modulus);
        let mut shift = [0; 6];
        for (value, &entry) in shift.iter_mut().zip(&difference) {
            if !entry.is_multiple_of(denominator) {
                return Err(inconsistent());
            }
            *value = sub_mod(0, entry / denominator, modulus);
        }
        let frobenius = Operator {
            linear: action(pairing, modulus),
            shift,
        };
        frobenius.verify_isometry(model, modulus)?;
        model.allocate(70)?;
        Ok(Self {
            quotient,
            precision,
            modulus,
            denominator,
            numerator,
            frobenius,
            period,
            relation: None,
        })
    }

    pub(crate) fn modulus(&self) -> u64 {
        self.modulus
    }

    pub(crate) fn shift(model: &mut Model, parent: Rc<Self>, change: [i128; 6]) -> Result<Self> {
        let h = parent.quotient.pairing.ring.modulus;
        let numerator = std::array::from_fn(|i| {
            let coordinate = change[i].rem_euclid(i128::from(h)) as u64;
            add_mod(
                parent.numerator[i],
                mul_mod(parent.denominator, coordinate, h),
                h,
            )
        });
        let change = change.map(|entry| entry.rem_euclid(i128::from(parent.modulus)) as u64);
        let forward = Operator::translation(scaled(change, parent.modulus - 1, parent.modulus));
        let reverse = Operator::translation(change);
        let next = Self::new(model, parent.quotient.clone(), parent.precision, numerator)?;
        next.equivalent(model, parent, forward, reverse)
    }

    pub(crate) fn algebraic_shift(
        model: &mut Model,
        parent: Rc<Self>,
        cycles: &Cycles,
        cycle: Cycle,
    ) -> Result<Self> {
        let pairing = &parent.quotient.pairing;
        if !Rc::ptr_eq(&cycles.pairing, pairing) {
            return Err(Error::plain("Brauer Mukai algebraic shift requires a geometric cycle in the chart's named pairing"));
        }
        let addition = cycles.lift(model, pairing.ring.precision, cycle)?;
        let numerator = add_vectors(parent.numerator, addition, pairing.ring.modulus);
        let next = Self::new(model, parent.quotient.clone(), parent.precision, numerator)?;
        // An exact rational Frobenius-fixed shift does not change the chart operator.
        next.equivalent(model, parent, Operator::identity(), Operator::identity())
    }

    pub(crate) fn extend(
        model: &mut Model,
        parent: Rc<Self>,
        quotient: Rc<Quotient>,
    ) -> Result<Self> {
        compatible_levels(&parent, &quotient)?;
        let base = &parent.quotient.pairing;
        let extension = &quotient.pairing;
        let ratio = extension_ratio(base, extension)?;
        let transfer = Transfer::new(model, base, extension, Direction::Restrict)?;
        transfer.apply(
            model,
            extension,
            quotient.level,
            parent.numerator.map(|x| x % parent.denominator),
        )?;
        let mut next = Self::new(model, quotient, parent.precision, parent.numerator)?;
        if parent.frobenius.power(model, ratio, parent.modulus)? != next.frobenius
            || !parent.period.is_multiple_of(next.period)
        {
            return Err(inconsistent());
        }
        model.allocate(50)?;
        next.relation = Some(Box::new(Relation {
            parent,
            forward: Operator::identity(),
            reverse: None,
        }));
        Ok(next)
    }

    pub(crate) fn marking(
        model: &mut Model,
        parent: Rc<Self>,
        quotient: Rc<Quotient>,
    ) -> Result<Self> {
        compatible_levels(&parent, &quotient)?;
        let change = quotient.pairing.change(model, &parent.quotient.pairing)?;
        let height_modulus = quotient.pairing.ring.modulus;
        model.spend(72)?;
        let numerator = image(&change.matrix, parent.numerator, height_modulus);
        let next = Self::new(model, quotient, parent.precision, numerator)?;
        let matrix = change
            .matrix
            .map(|row| row.map(|entry| entry % next.modulus));
        let matrix_ring = Residues::new(next.quotient.pairing.ring.prime, next.precision)?;
        let square = Matrix::from_entries(6, matrix.into_iter().flatten().collect())?;
        let inverse =
            hw_algebra::inverse(model, &matrix_ring, &square)?.ok_or_else(inconsistent)?;
        let inverse = std::array::from_fn(|i| std::array::from_fn(|j| *inverse.entry(i, j)));
        let forward = Operator {
            linear: matrix,
            shift: [0; 6],
        };
        let reverse = Operator {
            linear: inverse,
            shift: [0; 6],
        };
        next.equivalent(model, parent, forward, reverse)
    }

    fn equivalent(
        mut self,
        model: &mut Model,
        parent: Rc<Self>,
        forward: Operator,
        reverse: Operator,
    ) -> Result<Self> {
        let n = parent.modulus;
        if parent.period != self.period
            || forward.compose(model, &parent.frobenius, n)?
                != self.frobenius.compose(model, &forward, n)?
            || forward.compose(model, &reverse, n)? != Operator::identity()
        {
            return Err(inconsistent());
        }
        forward.verify_isometry(model, n)?;
        model.allocate(94)?;
        self.relation = Some(Box::new(Relation {
            parent,
            forward,
            reverse: Some(reverse),
        }));
        Ok(self)
    }

    pub(crate) fn transport(from: &Rc<Self>, to: &Rc<Self>) -> Result<Operator> {
        if Rc::ptr_eq(from, to) {
            return Ok(Operator::identity());
        }
        if let Some(relation) = to
            .relation
            .as_ref()
            .filter(|relation| Rc::ptr_eq(from, &relation.parent))
        {
            return Ok(relation.forward);
        }
        if let Some(relation) = from
            .relation
            .as_ref()
            .filter(|relation| Rc::ptr_eq(to, &relation.parent))
        {
            return relation.reverse.ok_or_else(|| {
                Error::plain("Brauer Mukai base-extension restriction has no inverse transport")
            });
        }
        Err(Error::plain(
            "Brauer Mukai transport requires a direct registered chart relation",
        ))
    }

    pub(crate) fn power(&self, model: &mut Model, exponent: usize) -> Result<Operator> {
        self.frobenius.power(model, exponent, self.modulus)
    }

    pub(crate) fn vector(&self, entries: [i128; 8]) -> Vector {
        entries.map(|x| x.rem_euclid(i128::from(self.modulus)) as u64)
    }

    pub(crate) fn witness(&self, model: &mut Model) -> Result<Vector> {
        // Solve F_B(r,c,0) = (r,c,0) with the independently certified rank r.
        // The last row is necessary over a truncated ring when r is a nonunit.
        let n = self.modulus;
        if self.period == self.denominator {
            // Here c = -rB is an exact integral choice, even if r reduces to zero.
            let mut vector = [0; 8];
            vector[0] = self.period % n;
            for (entry, numerator) in vector[1..7].iter_mut().zip(self.numerator) {
                *entry = sub_mod(0, numerator % n, n);
            }
            return self.verify_witness(model, vector);
        }
        let prime = self.quotient.pairing.ring.prime.get();
        let columns = self.frobenius.columns(model, n)?;
        let mut equations = bp_smith::Matrix::new(7, 6)?;
        model.spend(64)?;
        for (j, column) in columns.iter().skip(1).take(6).enumerate() {
            for (i, &entry) in column.iter().skip(1).enumerate() {
                equations.set(i, j, i128::from(sub_mod(entry, u64::from(i == j), n)));
            }
        }
        let target: Vec<_> = columns[0]
            .iter()
            .skip(1)
            .map(|&x| i128::from(sub_mod(0, x, n)))
            .collect();
        let mut rank = 6;
        while bp_smith::factors_at_precision(model, &equations, prime, n, rank)?.is_none() {
            rank -= 1;
        }
        let solution = bp_smith::modular_preimage(
            model,
            &equations,
            prime,
            n,
            &target,
            self.period % n,
            rank,
        )?;
        let mut vector = [0; 8];
        vector[0] = self.period % n;
        for (entry, coeff) in vector[1..7].iter_mut().zip(solution) {
            *entry = coeff.rem_euclid(i128::from(n)) as u64;
        }
        self.verify_witness(model, vector)
    }

    fn verify_witness(&self, model: &mut Model, vector: Vector) -> Result<Vector> {
        if self.frobenius.apply(model, vector, self.modulus)? != vector {
            return Err(inconsistent());
        }
        Ok(vector)
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str, power: usize) -> Result<String> {
        let operator = self.power(model, power)?;
        let columns = operator.columns(model, self.modulus)?;
        let ranks = defect_ranks(model, &columns, self.quotient.pairing.ring.prime)?;
        let denominator = self.denominator;
        let numerator = self.numerator;
        let period = self.period;
        let delta = operator.shift;
        let degree = self.quotient.pairing.degree;
        Ok(format!("rvBrauerMukaiData[{name}, {power}] = base_degree={degree}, precision={}, modulus={}, denominator={denominator}, numerator={numerator:?}, cohomological_index={period}, delta={delta:?}, frobenius_columns={columns:?}, defect_power_ranks_mod_prime={ranks:?}", self.precision, self.modulus))
    }
}

impl Operator {
    fn identity() -> Self {
        Self::translation([0; 6])
    }

    fn translation(shift: Fixed) -> Self {
        Self {
            linear: identity(),
            shift,
        }
    }

    fn compose(&self, model: &mut Model, rhs: &Self, n: u64) -> Result<Self> {
        model.spend(258)?;
        Ok(Self {
            linear: product(&self.linear, &rhs.linear, n),
            shift: add_vectors(self.shift, image(&self.linear, rhs.shift, n), n),
        })
    }

    fn power(mut self, model: &mut Model, mut exponent: usize, n: u64) -> Result<Self> {
        let mut result = Self::identity();
        while exponent != 0 {
            if exponent & 1 != 0 {
                result = result.compose(model, &self, n)?;
            }
            exponent >>= 1;
            if exponent != 0 {
                self = self.compose(model, &self, n)?;
            }
        }
        Ok(result)
    }

    pub(crate) fn apply(&self, model: &mut Model, value: Vector, n: u64) -> Result<Vector> {
        model.spend(65)?;
        let (rank, center, last) = components(value);
        let rotated = image(&self.linear, center, n);
        let corrected = add_vectors(rotated, scaled(self.shift, rank, n), n);
        let half_square = mul_mod(intersection(self.shift, self.shift, n), n / 2 + 1, n);
        let last = add_mod(
            last,
            add_mod(
                intersection(rotated, self.shift, n),
                mul_mod(rank, half_square, n),
                n,
            ),
            n,
        );
        let mut result = [rank; 8];
        result[1..7].copy_from_slice(&corrected);
        result[7] = last;
        Ok(result)
    }

    fn columns(&self, model: &mut Model, n: u64) -> Result<[Vector; 8]> {
        let mut columns = [[0; 8]; 8];
        for (i, column) in columns.iter_mut().enumerate() {
            *column = self.apply(model, unit(i), n)?;
        }
        Ok(columns)
    }

    fn verify_isometry(&self, model: &mut Model, n: u64) -> Result<()> {
        let columns = self.columns(model, n)?;
        model.spend(720)?;
        for i in 0..8 {
            for j in i..8 {
                if pair(columns[i], columns[j], n) != pair(unit(i), unit(j), n) {
                    return Err(inconsistent());
                }
            }
        }
        Ok(())
    }
}

pub(crate) fn pair(v: Vector, w: Vector, n: u64) -> u64 {
    let (a, x, b) = components(v);
    let (r, y, s) = components(w);
    sub_mod(
        intersection(x, y, n),
        add_mod(mul_mod(a, s, n), mul_mod(b, r, n), n),
        n,
    )
}

fn components(v: Vector) -> (u64, Fixed, u64) {
    let [rank, a, b, c, d, e, f, last] = v;
    (rank, [a, b, c, d, e, f], last)
}

fn unit(i: usize) -> Vector {
    std::array::from_fn(|j| u64::from(i == j))
}

fn defect_ranks(
    model: &mut Model,
    columns: &[Vector; 8],
    prime: crate::odd_cobar::Prime,
) -> Result<Vec<usize>> {
    let p = u64::from(prime.get());
    let defect = columns.map(|column| column.map(|x| x % p));
    let defect: [Vector; 8] = std::array::from_fn(|i| {
        std::array::from_fn(|j| sub_mod(defect[i][j], u64::from(i == j), p))
    });
    let mut powers = defect;
    let mut ranks = Vec::new();
    for _ in 0..8 {
        let mut basis = prime_linear::Basis::default();
        for column in &powers {
            let entries = column
                .iter()
                .enumerate()
                .filter(|(_, x)| **x != 0)
                .map(|(i, &x)| (i, x as u32))
                .collect();
            basis.insert(model, prime, entries)?;
        }
        let rank = basis.rank();
        let previous = ranks.last().copied();
        ranks.push(rank);
        if rank == 0 || previous == Some(rank) {
            break;
        }
        model.spend(512)?;
        powers = powers.map(|column| {
            std::array::from_fn(|i| {
                column
                    .iter()
                    .zip(&defect)
                    .fold(0, |acc, (&x, v)| add_mod(acc, mul_mod(x, v[i], p), p))
            })
        });
    }
    Ok(ranks)
}

fn compatible_levels(chart: &Chart, target: &Quotient) -> Result<()> {
    if chart.quotient.level != target.level || chart.quotient.pairing.ring != target.pairing.ring {
        return Err(Error::plain(
            "Brauer Mukai chart change requires the same torsion exponent and marking level",
        ));
    }
    Ok(())
}

fn inconsistent() -> Error {
    Error::plain(
        "Brauer Mukai computation disagrees with its Frobenius or intersection certificate",
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn higher_precision_uses_both_halves_of_the_mukai_square_formula() {
        let n = 9;
        let mut model = Model::new();
        let mut diagonal = identity();
        diagonal[0][0] = 8;
        diagonal[5][5] = 8;
        let operator = Operator {
            linear: diagonal,
            shift: [1, 0, 0, 0, 0, 1],
        };
        assert!(operator.verify_isometry(&mut model, n).is_ok());
        let vector = unit(0);
        let image = operator
            .apply(&mut model, vector, n)
            .ok()
            .expect("modulus nine");
        assert_eq!(image, [1, 1, 0, 0, 0, 0, 1, 4]);
        assert_eq!(operator.apply(&mut model, image, n).ok(), Some(vector));
        let even_power = operator.power(&mut model, 8, n).ok().expect("even power");
        assert!(even_power == Operator::identity());
    }
}
