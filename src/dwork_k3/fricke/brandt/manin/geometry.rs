use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::lee_linear::Rational as Q;
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear;
use crate::rv_lens::gcd;

use super::linear::{self as lin, Matrix, Vector, ONE, ZERO};
use super::{certificate, Operator};

const MAX_PRIME: usize = 43;
const MAX_HECKE: usize = 31;
pub(super) const MAX_CUSP_ENTRY: u128 = 1_000_000_000_000;

type Transform = [i128; 4];

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
struct Cusp(i128, i128);

impl Cusp {
    const INFINITY: Self = Self(1, 0);

    fn new(a: i128, b: i128) -> Result<Self> {
        let common = gcd(a.unsigned_abs(), b.unsigned_abs());
        if common == 0 {
            return Err(Error::plain("a Manin cusp cannot be [0,0]"));
        }
        if b == 0 {
            return Ok(Self::INFINITY);
        }
        let sign = if b < 0 { -1 } else { 1 };
        let common = sign * common as i128;
        Ok(Self(a / common, b / common))
    }

    fn transform(self, [a, b, c, d]: Transform) -> Result<Self> {
        Self::new(a * self.0 + b * self.1, c * self.0 + d * self.1)
    }
}

pub(super) struct Geometry {
    pub(super) prime: usize,
    pub(super) level: usize,
    pub(super) boundary: Matrix,
    free: Vec<usize>,
    projection: Vec<Vector>,
    lifts: Vec<[Cusp; 2]>,
    inverses: Vec<u32>,
}

impl Geometry {
    pub(super) fn new(model: &mut Model, p: usize) -> Result<Self> {
        if !(5..=MAX_PRIME).contains(&p) {
            return Err(Error::plain(format!(
                "Manin symbols require prime 5 ≤ p ≤ {MAX_PRIME}"
            )));
        }
        let prime = Prime::new(model, p)?;
        let mut inverses = vec![0];
        for c in 1..p {
            inverses.push(prime_linear::invert(model, prime, c as u32)?);
        }
        let level = 2 * p;
        let n = 3 * (p + 1);
        let mut geometry = Self {
            prime: p,
            level,
            boundary: Vec::new(),
            free: Vec::new(),
            projection: Vec::new(),
            lifts: Vec::new(),
            inverses,
        };
        let representatives: Vec<_> = [(1, 0), (0, 1), (1, 1)]
            .into_iter()
            .flat_map(|(c, d)| {
                (0..=p).map(move |s| {
                    let pair = if s == p { (1, 0) } else { (s, 1) };
                    (crt(c, pair.0, p) as i128, crt(d, pair.1, p) as i128)
                })
            })
            .collect();
        let mut relations = vec![vec![ZERO; n]; 2 * n];
        for (i, &(c, d)) in representatives.iter().enumerate() {
            if geometry.index(c, d)? != i {
                return Err(certificate("CRT enumeration of the Manin projective line"));
            }
            let s = geometry.index(d, -c)?;
            let r = geometry.index(d, -c - d)?;
            let r_squared = geometry.index(-c - d, c)?;
            for j in [i, s] {
                relations[i][j] = relations[i][j].add(ONE)?;
            }
            for j in [i, r, r_squared] {
                relations[n + i][j] = relations[n + i][j].add(ONE)?;
            }
        }
        let normal = lin::Frame::nullspace(model, n, relations)?;
        geometry.projection = lin::columns(&normal.vectors, n);
        geometry.free = normal.free;
        geometry.lifts = representatives
            .iter()
            .map(|&(c, d)| geometry.lift(c, d))
            .collect::<Result<_>>()?;
        geometry.boundary = lin::columns(
            &geometry
                .free
                .iter()
                .map(|&j| geometry.raw_boundary(geometry.lifts[j]))
                .collect::<Vec<_>>(),
            4,
        );
        for (i, lift) in geometry.lifts.iter().enumerate() {
            if lin::apply(model, &geometry.boundary, &geometry.projection[i])?
                != geometry.raw_boundary(*lift)
            {
                return Err(certificate("the cusp boundary respects Manin relations"));
            }
        }
        model.allocate(25 + 10 * n + 4 * n * geometry.dimension())?;
        Ok(geometry)
    }

    pub(super) fn dimension(&self) -> usize {
        self.free.len()
    }

    pub(super) fn count(&self) -> usize {
        self.projection.len()
    }

    pub(super) fn symbol(&self, c: i128, d: i128) -> Result<Vector> {
        Ok(self.projection[self.index(c, d)?].clone())
    }

    pub(super) fn path(
        &self,
        model: &mut Model,
        start: [i128; 2],
        end: [i128; 2],
    ) -> Result<Vector> {
        let start = Cusp::new(start[0], start[1])?;
        let end = Cusp::new(end[0], end[1])?;
        let mut path = self.partial(model, end)?;
        let previous = self.partial(model, start)?;
        lin::add_scaled(model, &mut path, &previous, Q::integer(-1))?;
        if lin::apply(model, &self.boundary, &path)? != self.raw_boundary([end, start]) {
            return Err(certificate(
                "the continued-fraction path has the prescribed boundary",
            ));
        }
        Ok(path)
    }

    fn index(&self, c: i128, d: i128) -> Result<usize> {
        let i = match (c.rem_euclid(2), d.rem_euclid(2)) {
            (1, 0) => 0,
            (0, 1) => 1,
            (1, 1) => 2,
            _ => return Err(primitive()),
        };
        let p = self.prime as i128;
        let c = c.rem_euclid(p) as usize;
        let d = d.rem_euclid(p) as usize;
        let j = if d == 0 {
            if c == 0 {
                return Err(primitive());
            }
            self.prime
        } else {
            c * self.inverses[d] as usize % self.prime
        };
        Ok(i * (self.prime + 1) + j)
    }

    fn lift(&self, mut c: i128, mut d: i128) -> Result<[Cusp; 2]> {
        let index = self.index(c, d)?;
        if d == 0 {
            c = 1;
        } else if c == 0 {
            d = 1;
        } else {
            for _ in 1..d {
                if gcd(c as u128, d as u128) == 1 {
                    break;
                }
                c += self.level as i128;
            }
        }
        let (a, b, divisor) = bezout(d, c);
        if divisor != 1 || a * d + b * c != 1 || self.index(c, d)? != index {
            return Err(certificate("SL2 lift of a primitive projective row"));
        }
        Ok([Cusp::new(a, c)?, Cusp::new(-b, d)?])
    }

    fn partial(&self, model: &mut Model, Cusp(mut a, mut b): Cusp) -> Result<Vector> {
        let mut vector = vec![ZERO; self.dimension()];
        let (mut u, mut v) = (Cusp::INFINITY, Cusp(0, 1));
        while b != 0 {
            let (quotient, remainder) = (a.div_euclid(b), a.rem_euclid(b));
            let t = Cusp(quotient * u.0 + v.0, quotient * u.1 + v.1);
            let det = t.0 * u.1 - u.0 * t.1;
            if det.unsigned_abs() != 1 {
                return Err(certificate("unimodular continued-fraction edge"));
            }
            let index = self.index(t.1, det * u.1)?;
            lin::add_scaled(model, &mut vector, &self.projection[index], ONE)?;
            (v, u, a, b) = (u, t, b, remainder);
        }
        Ok(vector)
    }

    fn cusp_index(&self, cusp: Cusp) -> usize {
        match gcd(cusp.1.unsigned_abs(), self.level as u128) as usize {
            1 => 0,
            2 => 1,
            p if p == self.prime => 2,
            _ => 3,
        }
    }

    fn raw_boundary(&self, [end, start]: [Cusp; 2]) -> Vector {
        let mut result = vec![ZERO; 4];
        let (i, j) = (self.cusp_index(end), self.cusp_index(start));
        if i != j {
            result[i] = ONE;
            result[j] = Q::integer(-1);
        }
        result
    }

    pub(super) fn operator(&self, model: &mut Model, operation: Operator) -> Result<Matrix> {
        let transforms = self.transforms(operation)?;
        let mut partials: BTreeMap<Cusp, Vector> = BTreeMap::new();
        let mut images = Vec::with_capacity(self.count());
        let n = self.level as i128;
        for &lift in &self.lifts {
            let image = self.correspondence(model, lift, &transforms, &mut partials)?;
            // Change both rows by a matrix in Gamma_0(N). The second lift
            // checks the actual cusp algorithm independently of the first.
            let alternate = [
                lift[0].transform([1, 1, n, n + 1])?,
                lift[1].transform([1, 1, n, n + 1])?,
            ];
            let witness = self.correspondence(model, alternate, &transforms, &mut partials)?;
            if image != witness {
                return Err(certificate("correspondence is independent of the SL2 lift"));
            }
            images.push(image);
        }
        let matrix = lin::columns(
            &self
                .free
                .iter()
                .map(|&j| images[j].clone())
                .collect::<Vec<_>>(),
            self.dimension(),
        );
        for (image, symbol) in images.iter().zip(&self.projection) {
            if *image != lin::apply(model, &matrix, symbol)? {
                return Err(certificate(
                    "the cusp correspondence respects Manin relations",
                ));
            }
        }
        let mut cusp_action = vec![vec![ZERO; 4]; 4];
        for (j, denominator) in [1, 2, self.prime, self.level].into_iter().enumerate() {
            let cusp = Cusp(1, denominator as i128);
            for &transform in &transforms {
                let i = self.cusp_index(cusp.transform(transform)?);
                cusp_action[i][j] = cusp_action[i][j].add(ONE)?;
            }
        }
        if lin::multiply(model, &self.boundary, &matrix)?
            != lin::multiply(model, &cusp_action, &self.boundary)?
        {
            return Err(certificate(
                "the correspondence commutes with the cusp boundary",
            ));
        }
        Ok(matrix)
    }

    fn correspondence(
        &self,
        model: &mut Model,
        [end, start]: [Cusp; 2],
        transforms: &[Transform],
        cache: &mut BTreeMap<Cusp, Vector>,
    ) -> Result<Vector> {
        let mut output = vec![ZERO; self.dimension()];
        for &transform in transforms {
            for (cusp, sign) in [(end, 1), (start, -1)] {
                let image = cusp.transform(transform)?;
                let path = match cache.entry(image) {
                    std::collections::btree_map::Entry::Occupied(path) => path.into_mut(),
                    std::collections::btree_map::Entry::Vacant(path) => {
                        path.insert(self.partial(model, image)?)
                    }
                };
                lin::add_scaled(model, &mut output, path, Q::integer(sign))?;
            }
        }
        Ok(output)
    }

    pub(super) fn transforms(&self, op: Operator) -> Result<Vec<Transform>> {
        match op {
            Operator::Scale(_) => Err(certificate("scalar action is handled directly")),
            Operator::J => Ok(vec![[-1, 0, 0, 1]]),
            Operator::Atkin(q) => {
                if q != 2 && q != self.prime {
                    return Err(Error::plain("Manin Atkin–Lehner index must be 2 or p"));
                }
                let q = q as i128;
                let n = self.level as i128;
                let r = n / q;
                let y = bezout(q, r).0.rem_euclid(r);
                let b = (q * y - 1) / r;
                Ok(vec![[q, b, n, q * y]])
            }
            Operator::Hecke(m) => {
                if !(1..=MAX_HECKE).contains(&m) || gcd(m as u128, self.level as u128) != 1 {
                    return Err(Error::plain(
                        "Manin Hecke index must be ≤ 31 and coprime to the level",
                    ));
                }
                let mut matrices = Vec::new();
                for d in (1..=m).filter(|d| m.is_multiple_of(*d)) {
                    for b in 0..d {
                        matrices.push([(m / d) as i128, b as i128, 0, d as i128]);
                    }
                }
                Ok(matrices)
            }
        }
    }
}

fn crt(parity: usize, residue: usize, p: usize) -> usize {
    if residue % 2 == parity {
        residue
    } else {
        residue + p
    }
}

fn bezout(mut a: i128, mut b: i128) -> (i128, i128, i128) {
    let (mut x, mut y, mut u, mut v) = (1, 0, 0, 1);
    while b != 0 {
        let quotient = a / b;
        (a, b, x, y, u, v) = (
            b,
            a - quotient * b,
            u,
            v,
            x - quotient * u,
            y - quotient * v,
        );
    }
    (x, y, a)
}

fn primitive() -> Error {
    Error::plain("a Manin bottom row must be primitive modulo 2p")
}
