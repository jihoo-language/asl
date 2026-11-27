use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use super::{Kernel, Point};
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_cyclotomic;
use crate::rv_lens::{self, Rational};
use crate::rv_unram::{add_mod, mul_mod, sub_mod};
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Statement {
    Verify(Label),
}

pub(crate) enum Expression {
    Named(Named, Label),
    Basis(Label, usize),
    Lift(u64, Box<ParsedExpr>),
    Automorphism(Label, Box<[ParsedExpr; 4]>),
    Projector(Label, Vec<ParsedExpr>),
    Vector(Label, Vec<ParsedExpr>),
    Unary(Unary, Box<ParsedExpr>),
    Binary(Binary, Box<[ParsedExpr; 2]>),
    Entry([usize; 2], Box<ParsedExpr>),
    Coefficient(usize, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Named {
    Fourier(bool),
    Chirp(bool),
    Identity,
    Parity,
    Gauss(bool),
    FourierAuto,
    ChirpAuto,
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Heisenberg,
    Intertwine,
    Trace,
}

#[derive(Clone, Copy)]
pub(crate) enum Binary {
    Compose,
    Apply,
    Scale,
    AutoCompose,
    Cocycle,
    LiftCompose,
    Image,
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Lift {
    pub(crate) point: Point,
    pub(crate) phase: u64,
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Auto(pub(crate) [Lift; 4]);

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Matrix(pub(crate) Rc<[Rational]>);

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Vector(pub(crate) Rc<[Rational]>);

#[derive(Clone, Copy)]
pub(crate) struct Space {
    divisors: [u64; 2],
    pub(crate) conductor: usize,
    size: usize,
}

type Position = [u64; 2];

const ZERO: Lift = Lift {
    point: [0; 4],
    phase: 0,
};
const MAX_DIMENSION: usize = 32;

impl Space {
    pub(crate) fn new(kernel: &Kernel) -> Result<Self> {
        let divisors = kernel.divisors;
        let size = usize::try_from(kernel.dimension)
            .ok()
            .filter(|&n| n <= MAX_DIMENSION)
            .ok_or_else(|| {
                Error::plain("CM metaplectic dimension exceeds the 32-dimensional matrix limit")
            })?;
        let conductor = 2 * divisors[1] as usize;
        Ok(Self {
            divisors,
            conductor,
            size,
        })
    }

    fn modulus(self) -> u64 {
        self.conductor as u64
    }

    fn moduli(self) -> Point {
        let [a, b] = self.divisors;
        [a, a, b, b]
    }

    fn position(self, i: usize) -> Position {
        [i as u64 / self.divisors[1], i as u64 % self.divisors[1]]
    }

    fn index(self, [a, b]: Position) -> usize {
        (a * self.divisors[1] + b) as usize
    }

    fn checked_index(self, index: usize) -> Result<usize> {
        if index < self.size {
            Ok(index)
        } else {
            Err(Error::plain("CM metaplectic basis index is out of range"))
        }
    }

    fn root(self, model: &mut Model, exponent: u64) -> Result<Rational> {
        let mut coefficients = vec![0; exponent as usize + 1];
        coefficients[exponent as usize] = 1;
        let numerator = rv_cyclotomic::reduced(model, self.conductor, coefficients)?;
        Ok(Rational {
            numerator: numerator.into(),
            denominator: 1,
        })
    }

    fn b(self, left: Position, right: Position) -> u64 {
        let n = self.modulus();
        let terms = self.divisors.map(|d| n / d);
        left.into_iter()
            .zip(right)
            .zip(terms)
            .fold(0, |sum, ((x, y), weight)| {
                add_mod(sum, mul_mod(mul_mod(x, y, n), weight, n), n)
            })
    }

    fn q(self, position: Position) -> u64 {
        let n = self.modulus();
        position
            .into_iter()
            .zip(self.divisors)
            .fold(0, |sum, (x, d)| {
                let numerator = if d.is_multiple_of(2) { 1 } else { d + 1 };
                let weight = n / (2 * d) * numerator;
                add_mod(sum, mul_mod(mul_mod(x, x, n), weight, n), n)
            })
    }

    fn position_add(self, a: Position, b: Position) -> Position {
        std::array::from_fn(|i| add_mod(a[i], b[i], self.divisors[i]))
    }

    fn positions(point: Point) -> (Position, Position) {
        ([point[0], point[2]], [point[1], point[3]])
    }

    pub(crate) fn lift(self, point: Point, phase: u64) -> Lift {
        Lift {
            point,
            phase: phase % self.modulus(),
        }
    }

    pub(crate) fn product(self, left: Lift, right: Lift) -> Lift {
        let n = self.modulus();
        let (_, b) = Self::positions(left.point);
        let (a, _) = Self::positions(right.point);
        Lift {
            point: std::array::from_fn(|i| {
                add_mod(left.point[i], right.point[i], self.moduli()[i])
            }),
            phase: sub_mod(add_mod(left.phase, right.phase, n), self.b(b, a), n),
        }
    }

    fn power(self, mut x: Lift, mut k: u64) -> Lift {
        let mut answer = ZERO;
        while k > 0 {
            if k % 2 == 1 {
                answer = self.product(answer, x);
            }
            k /= 2;
            x = self.product(x, x);
        }
        answer
    }

    fn commutator(self, a: Lift, b: Lift) -> u64 {
        sub_mod(
            self.product(a, b).phase,
            self.product(b, a).phase,
            self.modulus(),
        )
    }

    fn action(self, lift: Lift, source: usize) -> (usize, u64) {
        let (a, b) = Self::positions(lift.point);
        let position = self.position(source);
        let target = self.index(self.position_add(a, position));
        let phase = sub_mod(lift.phase, self.b(b, position), self.modulus());
        (target, phase)
    }

    fn basis_lifts(self) -> [Lift; 4] {
        let moduli = self.moduli();
        std::array::from_fn(|i| Lift {
            point: std::array::from_fn(|j| u64::from(i == j) % moduli[j]),
            phase: 0,
        })
    }

    pub(crate) fn fourier_auto(self) -> Auto {
        let mut target = self.basis_lifts();
        for (i, d) in self.divisors.into_iter().enumerate() {
            target[2 * i].point.swap(2 * i, 2 * i + 1);
            target[2 * i + 1].point[2 * i + 1] = 0;
            target[2 * i + 1].point[2 * i] = (d - 1) % d;
        }
        Auto(target)
    }

    pub(crate) fn chirp_auto(self) -> Auto {
        let mut target = self.basis_lifts();
        for (i, d) in self.divisors.into_iter().enumerate() {
            if d > 1 {
                target[2 * i].point[2 * i + 1] = d - 1;
                target[2 * i].phase = self.q(Self::positions(target[2 * i].point).0);
            }
        }
        Auto(target)
    }

    pub(crate) fn automorphism(self, model: &mut Model, target: [Lift; 4]) -> Result<Auto> {
        let source = self.basis_lifts();
        for (i, image) in target.into_iter().enumerate() {
            if self.power(image, self.moduli()[i]) != ZERO {
                return Err(Error::plain(
                    "CM metaplectic generator image does not satisfy its order relation",
                ));
            }
            for j in 0..i {
                if self.commutator(image, target[j]) != self.commutator(source[i], source[j]) {
                    return Err(Error::plain(
                        "CM metaplectic generator images do not preserve the commutator",
                    ));
                }
            }
        }
        let auto = Auto(target);
        let mut seen = BTreeSet::new();
        let moduli = self.moduli();
        for i in 0..self.size.pow(2) {
            model.spend(16)?;
            let mut x = i as u64;
            let point = moduli.map(|d| {
                let out = x % d;
                x /= d;
                out
            });
            if !seen.insert(self.image(auto, self.lift(point, 0)).point) {
                return Err(Error::plain(
                    "CM metaplectic generator images do not induce a bijection on the kernel",
                ));
            }
        }
        Ok(auto)
    }

    pub(crate) fn compose_auto(self, model: &mut Model, a: Auto, b: Auto) -> Result<Auto> {
        let images = b.0.map(|x| self.image(a, x));
        self.automorphism(model, images)
    }

    pub(crate) fn image(self, auto: Auto, lift: Lift) -> Lift {
        // The normal order X1^a1 Z1^b1 X2^a2 Z2^b2 has phase zero.
        auto.0.into_iter().zip(lift.point).fold(
            self.lift([0; 4], lift.phase),
            |result, (image, exponent)| self.product(result, self.power(image, exponent)),
        )
    }

    pub(crate) fn identity(self) -> Matrix {
        self.diagonal(|_| Rational::integer(1))
    }

    fn diagonal(self, value: impl Fn(usize) -> Rational) -> Matrix {
        let mut matrix = vec![Rational::zero(); self.size.pow(2)];
        for i in 0..self.size {
            matrix[i * (self.size + 1)] = value(i);
        }
        Matrix(matrix.into())
    }

    pub(crate) fn parity(self) -> Matrix {
        let mut matrix = vec![Rational::zero(); self.size.pow(2)];
        for x in 0..self.size {
            let position = self.position(x);
            let y = self.index(std::array::from_fn(|i| {
                sub_mod(0, position[i], self.divisors[i])
            }));
            matrix[y * self.size + x] = Rational::integer(1);
        }
        Matrix(matrix.into())
    }

    pub(crate) fn fourier(self, model: &mut Model, inverse: bool) -> Result<Matrix> {
        let mut matrix = Vec::with_capacity(self.size.pow(2));
        let n = self.modulus();
        for i in 0..self.size {
            for j in 0..self.size {
                let phase = self.b(self.position(i), self.position(j));
                let mut value =
                    self.root(model, if inverse { phase } else { sub_mod(0, phase, n) })?;
                if inverse {
                    value.denominator = self.size as u128;
                }
                matrix.push(value);
            }
        }
        Ok(Matrix(matrix.into()))
    }

    pub(crate) fn chirp(self, model: &mut Model, inverse: bool) -> Result<Matrix> {
        let n = self.modulus();
        let mut coefficients = Vec::with_capacity(self.size);
        for i in 0..self.size {
            let phase = self.q(self.position(i));
            coefficients
                .push(self.root(model, if inverse { sub_mod(0, phase, n) } else { phase })?);
        }
        Ok(self.diagonal(|i| coefficients[i].clone()))
    }

    pub(crate) fn gauss(self, model: &mut Model, conjugate: bool) -> Result<Rational> {
        let n = self.modulus();
        let mut coefficients = vec![0; self.conductor];
        for i in 0..self.size {
            let phase = self.q(self.position(i));
            coefficients[if conjugate {
                sub_mod(0, phase, n)
            } else {
                phase
            } as usize] += 1;
        }
        let numerator = rv_cyclotomic::reduced(model, self.conductor, coefficients)?;
        Ok(Rational {
            numerator: numerator.into(),
            denominator: 1,
        })
    }

    pub(crate) fn heisenberg(self, model: &mut Model, lift: Lift) -> Result<Matrix> {
        let mut matrix = vec![Rational::zero(); self.size.pow(2)];
        for i in 0..self.size {
            let (j, phase) = self.action(lift, i);
            matrix[j * self.size + i] = self.root(model, phase)?;
        }
        Ok(Matrix(matrix.into()))
    }

    pub(crate) fn basis(self, index: usize) -> Result<Vector> {
        let i = self.checked_index(index)?;
        let mut entries = vec![Rational::zero(); self.size];
        entries[i] = Rational::integer(1);
        Ok(Vector(entries.into()))
    }

    pub(crate) fn vector(self, entries: Vec<Rational>) -> Result<Vector> {
        if entries.len() != self.size {
            return Err(Error::plain(
                "CM metaplectic vector needs exactly the dimension many coefficients",
            ));
        }
        Ok(Vector(entries.into()))
    }

    pub(crate) fn entry(self, m: &Matrix, [row, column]: [usize; 2]) -> Result<Rational> {
        let i = self.checked_index(row)?;
        let j = self.checked_index(column)?;
        Ok(m.0[i * self.size + j].clone())
    }

    pub(crate) fn coefficient(self, v: &Vector, i: usize) -> Result<Rational> {
        Ok(v.0[self.checked_index(i)?].clone())
    }

    fn dot<'a>(
        self,
        model: &mut Model,
        a: &[Rational],
        b: impl Iterator<Item = &'a Rational>,
    ) -> Result<Rational> {
        let mut value = Rational::zero();
        for (x, y) in a.iter().zip(b) {
            let product = rv_lens::multiply(model, self.conductor, x, y)?;
            value = rv_lens::sum(model, &value, &product)?;
        }
        Ok(value)
    }

    pub(crate) fn multiply(self, model: &mut Model, a: &Matrix, b: &Matrix) -> Result<Matrix> {
        model.spend(self.size.pow(3))?;
        let mut entries = Vec::with_capacity(self.size.pow(2));
        for row in a.0.chunks_exact(self.size) {
            for j in 0..self.size {
                let column = b.0[j..].iter().step_by(self.size);
                entries.push(self.dot(model, row, column)?);
            }
        }
        Ok(Matrix(entries.into()))
    }

    pub(crate) fn apply(self, model: &mut Model, a: &Matrix, b: &Vector) -> Result<Vector> {
        model.spend(self.size.pow(2))?;
        let mut entries = Vec::with_capacity(self.size);
        for row in a.0.chunks_exact(self.size) {
            entries.push(self.dot(model, row, b.0.iter())?);
        }
        Ok(Vector(entries.into()))
    }

    pub(crate) fn scale(self, model: &mut Model, a: &Matrix, scalar: &Rational) -> Result<Matrix> {
        let entries =
            a.0.iter()
                .map(|x| rv_lens::multiply(model, self.conductor, x, scalar))
                .collect::<Result<Vec<_>>>()?;
        Ok(Matrix(entries.into()))
    }

    pub(crate) fn trace(self, model: &mut Model, a: &Matrix) -> Result<Rational> {
        let mut value = Rational::zero();
        for x in a.0.iter().step_by(self.size + 1) {
            value = rv_lens::sum(model, &value, x)?;
        }
        Ok(value)
    }

    fn shift_vector(
        self,
        model: &mut Model,
        a: Lift,
        vector: &[Rational],
    ) -> Result<Vec<Rational>> {
        let mut output = vec![Rational::zero(); self.size];
        for (i, value) in vector.iter().enumerate() {
            let (j, phase) = self.action(a, i);
            output[j] = rv_lens::shifted(model, self.conductor, value, i128::from(phase))?;
        }
        Ok(output)
    }

    fn projector_column(
        self,
        model: &mut Model,
        subgroup: &[Lift],
        i: usize,
    ) -> Result<Vec<Rational>> {
        model.spend(subgroup.len() * self.conductor)?;
        let mut entries = vec![vec![0; self.conductor]; self.size];
        for &element in subgroup {
            let (j, phase) = self.action(element, i);
            entries[j][phase as usize] += 1;
        }
        entries
            .into_iter()
            .map(|entry| {
                let entry = rv_cyclotomic::reduced(model, self.conductor, entry)?;
                rv_lens::fraction(model, &entry, subgroup.len() as u128)
            })
            .collect()
    }

    pub(crate) fn projector(self, model: &mut Model, generators: &[Lift]) -> Result<Matrix> {
        let mut subgroup = BTreeMap::from([(ZERO.point, ZERO)]);
        let mut pending = vec![ZERO];
        while let Some(a) = pending.pop() {
            for &b in generators {
                model.spend(1)?;
                let next = self.product(a, b);
                match subgroup.get(&next.point) {
                    Some(&old) if old != next => {
                        return Err(Error::plain(
                            "CM metaplectic projector lifts generate a nontrivial central phase",
                        ))
                    }
                    Some(_) => {}
                    None => {
                        if subgroup.len() == self.size {
                            return Err(Error::plain(
                                "CM metaplectic projector requires a lifted maximal isotropic subgroup",
                            ));
                        }
                        subgroup.insert(next.point, next);
                        pending.push(next);
                    }
                }
            }
        }
        if subgroup.len() != self.size {
            return Err(Error::plain(
                "CM metaplectic projector requires a lifted maximal isotropic subgroup",
            ));
        }
        let subgroup = subgroup.into_values().collect::<Vec<_>>();
        let mut entries = vec![Rational::zero(); self.size.pow(2)];
        for i in 0..self.size {
            for (j, value) in self
                .projector_column(model, &subgroup, i)?
                .into_iter()
                .enumerate()
            {
                entries[j * self.size + i] = value;
            }
        }
        let result = Matrix(entries.into());
        if self.trace(model, &result)? != Rational::integer(1) {
            return Err(inconsistent());
        }
        Ok(result)
    }

    fn intertwiner(self, model: &mut Model, auto: Auto) -> Result<(Matrix, usize)> {
        let clocks = (0..self.size)
            .map(|i| {
                let [a, b] = self.position(i);
                self.product(self.power(auto.0[1], a), self.power(auto.0[3], b))
            })
            .collect::<Vec<_>>();
        let mut data = None;
        for i in 0..self.size {
            let v = self.projector_column(model, &clocks, i)?;
            if !v[i].is_zero() {
                data = Some((i, v));
                break;
            }
        }
        let (pivot, mut zero_column) = data.ok_or_else(inconsistent)?;
        // The diagonal is the proportion of clock lifts stabilizing the
        // chosen position with trivial total character. It is rational.
        let diagonal = &zero_column[pivot];
        let [value] = diagonal.numerator.as_ref() else {
            return Err(inconsistent());
        };
        if *value <= 0 {
            return Err(inconsistent());
        }
        let numerator = i128::try_from(diagonal.denominator).map_err(|_| inconsistent())?;
        let inverse = rv_lens::fraction(model, &[numerator], *value as u128)?;
        for coefficient in &mut zero_column {
            *coefficient = rv_lens::multiply(model, self.conductor, coefficient, &inverse)?;
        }
        let mut columns = Vec::with_capacity(self.size);
        for i in 0..self.size {
            let [a, b] = self.position(i);
            let target_shift = self.product(self.power(auto.0[0], a), self.power(auto.0[2], b));
            columns.push(self.shift_vector(model, target_shift, &zero_column)?);
        }
        if columns[0][pivot] != Rational::integer(1) {
            return Err(inconsistent());
        }
        for (source, target) in self.basis_lifts().into_iter().zip(auto.0) {
            for (i, column) in columns.iter().enumerate() {
                let (j, phase) = self.action(source, i);
                let lhs = self.shift_vector(model, self.lift([0; 4], phase), &columns[j])?;
                let rhs = self.shift_vector(model, target, column)?;
                if lhs != rhs {
                    return Err(inconsistent());
                }
            }
        }
        let mut entries = Vec::with_capacity(self.size.pow(2));
        for i in 0..self.size {
            entries.extend(columns.iter().map(|column| column[i].clone()));
        }
        Ok((Matrix(entries.into()), pivot))
    }

    pub(crate) fn normalized(self, model: &mut Model, auto: Auto) -> Result<Matrix> {
        Ok(self.intertwiner(model, auto)?.0)
    }

    pub(crate) fn cocycle(self, model: &mut Model, a: Auto, b: Auto) -> Result<Rational> {
        let product = self.compose_auto(model, a, b)?;
        let (normalized, pivot) = self.intertwiner(model, product)?;
        let a = self.normalized(model, a)?;
        let b = self.normalized(model, b)?;
        let composed = self.multiply(model, &a, &b)?;
        let scalar = self.entry(&composed, [pivot, 0])?;
        if scalar.is_zero() || self.scale(model, &normalized, &scalar)? != composed {
            return Err(inconsistent());
        }
        Ok(scalar)
    }

    fn verify_full_lifts(self, model: &mut Model, fourier: Auto, chirp: Auto) -> Result<()> {
        let n = self.modulus();
        for i in 0..self.size {
            let a = self.position(i);
            for j in 0..self.size {
                model.spend(32)?;
                let b = self.position(j);
                let input = self.lift([a[0], b[0], a[1], b[1]], 0);
                let rotation = self.lift(
                    [
                        sub_mod(0, b[0], self.divisors[0]),
                        a[0],
                        sub_mod(0, b[1], self.divisors[1]),
                        a[1],
                    ],
                    self.b(a, b),
                );
                let shear = self.lift(
                    [
                        a[0],
                        sub_mod(b[0], a[0], self.divisors[0]),
                        a[1],
                        sub_mod(b[1], a[1], self.divisors[1]),
                    ],
                    self.q(a),
                );
                if self.image(fourier, input) != rotation
                    || self.image(chirp, input) != shear
                    || self.q(self.position_add(a, b))
                        != add_mod(add_mod(self.q(a), self.q(b), n), self.b(a, b), n)
                {
                    return Err(inconsistent());
                }
            }
        }
        Ok(())
    }

    pub(crate) fn verify(self, model: &mut Model) -> Result<String> {
        let n = Rational::integer(self.size as i128);
        let f = self.fourier(model, false)?;
        let fi = self.fourier(model, true)?;
        let t = self.chirp(model, false)?;
        let ti = self.chirp(model, true)?;
        let identity = self.identity();
        if self.multiply(model, &f, &fi)? != identity
            || self.multiply(model, &fi, &f)? != identity
            || self.multiply(model, &t, &ti)? != identity
            || self.multiply(model, &ti, &t)? != identity
        {
            return Err(inconsistent());
        }
        let f2 = self.multiply(model, &f, &f)?;
        let ft = self.multiply(model, &f, &t)?;
        let ft2 = self.multiply(model, &ft, &ft)?;
        let ft3 = self.multiply(model, &ft2, &ft)?;
        let g = self.gauss(model, false)?;
        let conjugate = self.gauss(model, true)?;
        if f2 != self.scale(model, &self.parity(), &n)?
            || ft3 != self.scale(model, &f2, &g)?
            || rv_lens::multiply(model, self.conductor, &g, &conjugate)? != n
        {
            return Err(inconsistent());
        }
        let fa = self.automorphism(model, self.fourier_auto().0)?;
        let ta = self.automorphism(model, self.chirp_auto().0)?;
        self.verify_full_lifts(model, fa, ta)?;
        for (matrix, auto) in [(&f, fa), (&t, ta)] {
            for (from, to) in self.basis_lifts().into_iter().zip(auto.0) {
                let from = self.heisenberg(model, from)?;
                let to = self.heisenberg(model, to)?;
                if self.multiply(model, matrix, &from)? != self.multiply(model, &to, matrix)? {
                    return Err(inconsistent());
                }
            }
        }
        let first = self.basis_lifts();
        let clocks = self.projector(model, &[first[1], first[3]])?;
        let shifts = self.projector(model, &[first[0], first[2]])?;
        let chirped = self.projector(model, &[ta.0[0], ta.0[2]])?;
        let product = self.multiply(model, &clocks, &shifts)?;
        let product = self.multiply(model, &product, &chirped)?;
        let triple = self.trace(model, &product)?;
        let expected = rv_lens::fraction(model, &g.numerator, self.size.pow(2) as u128)?;
        if triple != expected {
            return Err(inconsistent());
        }
        let phases: Vec<_> = (0..self.size)
            .map(|i| mul_mod(self.divisors[1], self.q(self.position(i)), self.modulus()))
            .collect();
        let clock = self.divisors.map(|d| {
            if d.is_multiple_of(2) && !(self.divisors[1] / d).is_multiple_of(2) {
                d / 2
            } else {
                0
            }
        });
        let clock_lift = self.lift([0, clock[0], 0, clock[1]], 0);
        for (i, &phase) in phases.iter().enumerate() {
            if self.action(clock_lift, i) != (i, phase) {
                return Err(inconsistent());
            }
        }
        Ok(format!(
            "rvCMMetaVerify = dimension={}, conductor={}, gauss={:?}/{}, triple_trace={:?}/{}, top_chirp_power_phases={phases:?}, top_chirp_power_clock={clock:?}",
            self.size, self.conductor, g.numerator, g.denominator, triple.numerator, triple.denominator,
        ))
    }
}

pub(crate) fn persist(model: &mut Model, values: &[Rational]) -> Result<()> {
    model.allocate(values.len() + 2 * values.iter().map(|x| x.numerator.len()).sum::<usize>())
}

fn inconsistent() -> Error {
    Error::plain("CM metaplectic exact intertwining certificate is inconsistent")
}
