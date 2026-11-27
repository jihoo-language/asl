use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::hyperelliptic::{self, Residues};
use crate::model::{Model, TypeId, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::invert;
use crate::rv_unram::Ring;

use super::{connection, Operator};

mod fiber;

const MAX_PRIME: u32 = 257;
type Poly = Vec<u32>;
type PolyMatrix = Matrix<Poly>;

pub(crate) struct PrimePoint {
    pub(crate) prime: usize,
    pub(crate) rank: usize,
    pub(crate) point: i128,
    pub(crate) exponent: usize,
}

impl PrimePoint {
    pub(crate) fn value(self, model: &mut Model) -> Result<(TypeId, Value)> {
        let rank = Operator::with_rank(self.rank)?.rank();
        let p = checked_prime(model, self.prime)?;
        let ring = Residues::new(p, 1)?;
        let ty = model.hw_matrix_type(ring, rank)?;
        let matrix = Data::new(model, p)?.specialize(
            model,
            rank,
            &ring,
            &ring.coefficient(self.point),
            u64::from,
        )?;
        let powered = fiber::iterate(model, &ring, matrix, self.exponent)?;
        Ok((ty, hyperelliptic::persist_matrix(model, powered)?))
    }

    pub(crate) fn report(self, model: &mut Model) -> Result<String> {
        let rank = Operator::with_rank(self.rank)?.rank();
        let p = checked_prime(model, self.prime)?;
        let ring = Residues::new(p, 1)?;
        let z = ring.coefficient(self.point);
        let description = Data::new(model, p)?.fiber(model, rank, &ring, &z, u64::from)?;
        Ok(format!(
            "rvDworkPCurvatureFiber[{}, {rank}, {}] = {description}",
            p.get(),
            self.point
        ))
    }
}

pub(crate) fn require_closed(rank: usize, ring: &Ring) -> Result<()> {
    Operator::with_rank(rank)?;
    if ring.precision() != 1 {
        return Err(Error::plain(
            "Dwork p-curvature specialization requires a residue field of precision one",
        ));
    }
    require_prime(ring.prime())
}

pub(crate) fn closed_value(
    model: &mut Model,
    rank: usize,
    ring: &Ring,
    z: &[u64],
    exponent: usize,
) -> Result<Value> {
    require_closed(rank, ring)?;
    let p = Prime::validated(ring.prime() as u32);
    let data = Data::new(model, p)?;
    let matrix = data.specialize(model, rank, ring, &z.to_vec(), |c| {
        ring.scalar(u64::from(c))
    })?;
    let powered = fiber::iterate(model, ring, matrix, exponent)?;
    crate::hw_closed::persist(model, ring, powered)
}

pub(crate) fn closed_report(
    model: &mut Model,
    rank: usize,
    ring: &Ring,
    z: &[u64],
) -> Result<String> {
    require_closed(rank, ring)?;
    let p = Prime::validated(ring.prime() as u32);
    let data = Data::new(model, p)?;
    let description = data.fiber(model, rank, ring, &z.to_vec(), |c| {
        ring.scalar(u64::from(c))
    })?;
    Ok(format!(
        "rvDworkPCurvatureClosedFiber[{rank}] = p={}, z={z:?}, {description}",
        p.get()
    ))
}

pub(crate) fn report(model: &mut Model, prime: usize, rank: usize) -> Result<String> {
    let rank = Operator::with_rank(rank)?.rank();
    let p = checked_prime(model, prime)?;
    let data = Data::new(model, p)?;
    let entries = data
        .numerator(rank)
        .entries()
        .chunks(rank)
        .collect::<Vec<_>>();
    let ranks = (0..=rank).rev().collect::<Vec<_>>();
    Ok(format!(
        "rvDworkPCurvature[{prime}, {rank}] = denominator=(1-256z)^{prime}, numerator={entries:?}, generic power ranks={ranks:?}, nilpotency index={rank}, Hasse={:?}, certified nilpotence, horizontality, skew-adjointness, symmetric square and Hasse section",
        data.hasse[0],
    ))
}

fn require_prime(prime: u64) -> Result<()> {
    if prime < 3 || prime.is_multiple_of(2) {
        return Err(Error::plain("Dwork p-curvature requires an odd prime"));
    }
    if prime > u64::from(MAX_PRIME) {
        return Err(Error::plain(format!(
            "Dwork p-curvature requires p at most {MAX_PRIME}"
        )));
    }
    Ok(())
}

fn checked_prime(model: &mut Model, prime: usize) -> Result<Prime> {
    require_prime(prime as u64)?;
    Prime::new(model, prime)
}

struct Data {
    prime: Prime,
    rank_two: PolyMatrix,
    rank_three: PolyMatrix,
    hasse: [Poly; 3],
}

impl Data {
    fn new(model: &mut Model, prime: Prime) -> Result<Self> {
        let ring = Polynomials(prime);
        let d = ring.linear(connection::D);
        let two = ring.matrix(connection::companion_two())?;
        let three = ring.matrix(connection::companion_three())?;
        let rank_two = numerator(model, ring, &d, &two)?;
        let rank_three = numerator(model, ring, &d, &three)?;
        let second = hw_algebra::multiply(model, &ring, &rank_three, &rank_three)?;
        let third = hw_algebra::multiply(model, &ring, &second, &rank_three)?;
        let square = hw_algebra::multiply(model, &ring, &rank_two, &rank_two)?;
        // A(0) is the upper Jordan block, so P_p(0)=0 and N(0)=A(0).
        // The constant coefficients give uniform maximal-index witnesses.
        if !matrix_zero(&square)
            || !matrix_zero(&third)
            || rank_two.entry(0, 1).first() != Some(&1)
            || second.entry(0, 2).first() != Some(&1)
        {
            return Err(integrity());
        }
        check_horizontal(model, ring, &d, &two, &rank_two)?;
        check_horizontal(model, ring, &d, &three, &rank_three)?;
        check_skew(
            model,
            ring,
            &rank_two,
            &ring.matrix([[[0; 3], [1, 0, 0]], [[-1, 0, 0], [0; 3]]])?,
        )?;
        check_skew(
            model,
            ring,
            &rank_three,
            &ring.matrix(connection::polarization())?,
        )?;
        check_symmetric_square(model, ring, &d, &rank_two, &rank_three)?;
        let hasse = hasse(model, ring, &d, &three, &rank_three)?;
        Ok(Self {
            prime,
            rank_two,
            rank_three,
            hasse,
        })
    }

    fn numerator(&self, rank: usize) -> &PolyMatrix {
        if rank == 2 {
            &self.rank_two
        } else {
            &self.rank_three
        }
    }
}

// P_j represents (theta-A)^j on constant basis columns, with denominator D^j.
fn numerator(model: &mut Model, ring: Polynomials, d: &Poly, a: &PolyMatrix) -> Result<PolyMatrix> {
    let p = ring.0;
    let rank = a.size();
    let d_derivative = ring.theta(d);
    let mut current = hw_algebra::identity(model, &ring, rank)?;
    for j in 0..p.get() {
        let left = hw_algebra::multiply(model, &ring, a, &current)?;
        let j_derivative = ring.scale(&d_derivative, j);
        let mut next = Vec::with_capacity(rank * rank);
        for (entry, left) in current.entries().iter().zip(left.entries()) {
            let first = ring.multiply(model, d, &ring.theta(entry))?;
            let second = ring.multiply(model, &j_derivative, entry)?;
            let difference = ring.subtract(model, &first, &second)?;
            next.push(ring.subtract(model, &difference, left)?);
        }
        current = Matrix::from_entries(rank, next)?;
    }
    // D^(p-1) = (1-(256z)^p)/(1-256z) in F_p[z].
    let mut power = 1;
    let z_coefficient = p.negate(d[1]);
    let d_power = (0..p.get())
        .map(|_| {
            let coefficient = power;
            power = p.mul(power, z_coefficient);
            coefficient
        })
        .collect::<Vec<_>>();
    let entries = current
        .entries()
        .iter()
        .zip(a.entries())
        .map(|(p, a)| {
            let term = ring.multiply(model, a, &d_power)?;
            ring.add(model, p, &term)
        })
        .collect::<Result<_>>()?;
    Matrix::from_entries(rank, entries)
}

fn check_horizontal(
    model: &mut Model,
    ring: Polynomials,
    d: &Poly,
    a: &PolyMatrix,
    n: &PolyMatrix,
) -> Result<()> {
    let an = hw_algebra::multiply(model, &ring, a, n)?;
    let na = hw_algebra::multiply(model, &ring, n, a)?;
    // theta(D^p)=0, so both sides have denominator D^(p+1).
    for ((entry, an), na) in n.entries().iter().zip(an.entries()).zip(na.entries()) {
        if ring.multiply(model, d, &ring.theta(entry))? != ring.subtract(model, an, na)? {
            return Err(integrity());
        }
    }
    Ok(())
}

fn check_skew(model: &mut Model, ring: Polynomials, n: &PolyMatrix, b: &PolyMatrix) -> Result<()> {
    let product = hw_algebra::multiply(model, &ring, b, n)?;
    let transpose = transpose(n)?;
    let second = hw_algebra::multiply(model, &ring, &transpose, b)?;
    for (a, b) in product.entries().iter().zip(second.entries()) {
        if !ring.add(model, a, b)?.is_empty() {
            return Err(integrity());
        }
    }
    Ok(())
}

fn transpose(m: &PolyMatrix) -> Result<PolyMatrix> {
    let rank = m.size();
    Matrix::from_entries(
        rank,
        (0..rank)
            .flat_map(|i| (0..rank).map(move |j| m.entry(j, i).clone()))
            .collect(),
    )
}

fn check_symmetric_square(
    model: &mut Model,
    ring: Polynomials,
    d: &Poly,
    two: &PolyMatrix,
    three: &PolyMatrix,
) -> Result<()> {
    let [v, w, x, y] = two.entries() else {
        unreachable!()
    };
    let trace = ring.add(model, v, y)?;
    let action = Matrix::from_entries(
        3,
        vec![
            ring.scale(v, 2),
            ring.scale(w, 2),
            Vec::new(),
            x.clone(),
            trace,
            w.clone(),
            Vec::new(),
            ring.scale(x, 2),
            ring.scale(y, 2),
        ],
    )?;
    let zero = Vec::new();
    let companion = ring.matrix(connection::companion_two())?;
    // G maps (UV,(U'V+UV')/2,U'V') to the three-jet frame.
    let g = Matrix::from_entries(
        3,
        vec![
            d.clone(),
            zero.clone(),
            zero.clone(),
            zero,
            ring.scale(d, 2),
            Vec::new(),
            ring.scale(companion.entry(1, 0), 2),
            ring.scale(companion.entry(1, 1), 2),
            ring.scale(d, 2),
        ],
    )?;
    if hw_algebra::multiply(model, &ring, three, &g)?
        != hw_algebra::multiply(model, &ring, &g, &action)?
    {
        return Err(integrity());
    }
    Ok(())
}

fn hasse(
    model: &mut Model,
    ring: Polynomials,
    d: &Poly,
    a: &PolyMatrix,
    n: &PolyMatrix,
) -> Result<[Poly; 3]> {
    let p = ring.0;
    let mut h = vec![1];
    for i in 1..=(p.get() - 1) / 4 {
        let inverse = invert(model, p, i)?;
        let dot = [
            4,
            4 * i - 3,
            4 * i - 2,
            4 * i - 1,
            inverse,
            inverse,
            inverse,
        ];
        let value = dot
            .into_iter()
            .fold(*h.last().unwrap(), |v, x| p.mul(v, x % p.get()));
        h.push(value);
    }
    let first = ring.theta(&h);
    let second = ring.theta(&first);
    let section = [h, first, second];
    for i in 0..3 {
        let mut horizontal = Vec::new();
        let mut curvature = Vec::new();
        for (j, value) in section.iter().enumerate() {
            let term = ring.multiply(model, a.entry(i, j), value)?;
            horizontal = ring.add(model, &horizontal, &term)?;
            let term = ring.multiply(model, n.entry(i, j), value)?;
            curvature = ring.add(model, &curvature, &term)?;
        }
        if horizontal != ring.multiply(model, d, &ring.theta(&section[i]))? || !curvature.is_empty()
        {
            return Err(integrity());
        }
    }
    Ok(section)
}

fn matrix_zero(matrix: &PolyMatrix) -> bool {
    matrix.entries().iter().all(Vec::is_empty)
}

fn integrity() -> Error {
    Error::plain("Dwork p-curvature failed an exact differential or linear-algebra identity")
}

#[derive(Clone, Copy)]
struct Polynomials(Prime);

impl Polynomials {
    fn linear(self, coefficients: connection::Polynomial) -> Poly {
        trim(
            coefficients
                .map(|c| c.rem_euclid(i128::from(self.0.get())) as u32)
                .to_vec(),
        )
    }

    fn matrix<const N: usize>(
        self,
        entries: [[connection::Polynomial; N]; N],
    ) -> Result<PolyMatrix> {
        Matrix::from_entries(
            N,
            entries
                .into_iter()
                .flatten()
                .map(|p| self.linear(p))
                .collect(),
        )
    }

    fn theta(self, a: &[u32]) -> Poly {
        trim(
            a.iter()
                .enumerate()
                .map(|(i, &c)| self.0.mul(self.0.reduce(i), c))
                .collect(),
        )
    }

    fn scale(self, a: &[u32], b: u32) -> Poly {
        trim(a.iter().map(|&c| self.0.mul(c, b)).collect())
    }

    fn sum(self, model: &mut Model, a: &[u32], b: &[u32], subtract: bool) -> Result<Poly> {
        let width = a.len().max(b.len());
        model.spend(width)?;
        Ok(trim(
            (0..width)
                .map(|i| {
                    let x = a.get(i).copied().unwrap_or(0);
                    let mut y = b.get(i).copied().unwrap_or(0);
                    if subtract {
                        y = self.0.negate(y);
                    }
                    self.0.add(x, y)
                })
                .collect(),
        ))
    }
}

fn trim(mut poly: Poly) -> Poly {
    poly.truncate(poly.iter().rposition(|&c| c != 0).map_or(0, |i| i + 1));
    poly
}

impl CoefficientRing for Polynomials {
    type Element = Poly;

    fn cells(&self) -> usize {
        1
    }

    fn zero(&self) -> Poly {
        Vec::new()
    }

    fn one(&self) -> Poly {
        vec![1]
    }

    fn is_zero(&self, value: &Poly) -> bool {
        value.is_empty()
    }

    fn is_unit(&self, value: &Poly) -> bool {
        value.len() == 1
    }

    fn add(&self, model: &mut Model, a: &Poly, b: &Poly) -> Result<Poly> {
        self.sum(model, a, b, false)
    }

    fn subtract(&self, model: &mut Model, a: &Poly, b: &Poly) -> Result<Poly> {
        self.sum(model, a, b, true)
    }

    fn multiply(&self, model: &mut Model, a: &Poly, b: &Poly) -> Result<Poly> {
        if a.is_empty() || b.is_empty() {
            return Ok(self.zero());
        }
        model.spend(a.len() * b.len())?;
        let mut result = vec![0; a.len() + b.len() - 1];
        for (i, &x) in a.iter().enumerate().filter(|&(_, &x)| x != 0) {
            for (out, &y) in result[i..].iter_mut().zip(b) {
                *out = self.0.add(*out, self.0.mul(x, y));
            }
        }
        Ok(trim(result))
    }

    fn inverse_unit(&self, model: &mut Model, value: &Poly) -> Result<Poly> {
        match value.as_slice() {
            &[constant] => Ok(vec![invert(model, self.0, constant)?]),
            _ => Err(integrity()),
        }
    }
}
