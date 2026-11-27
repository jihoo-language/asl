use std::collections::BTreeMap;
use std::rc::Rc;

use crate::ec_poly::Scalar;
use crate::elliptic::Curve;
use crate::elliptic_local::{self, Test};

use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_unram::Ring;

use super::{certificate, require_good, Space};

const MAX_FIELD: u64 = 121;

#[derive(Clone, Copy)]
struct Lines {
    trace: i32,
    two: usize,
    level_prime: usize,
    scalar: bool,
}

struct Tables {
    ring_name: Rc<str>,
    ring: Rc<Ring>,
    elements: Vec<Scalar>,
    scales: Vec<(usize, usize)>,
    scalar_orbits: BTreeMap<(usize, usize), bool>,
    prime: usize,
    q: usize,
    add: Vec<Vec<usize>>,
    multiply: Vec<Vec<usize>>,
    cubes: Vec<usize>,
    character: Vec<i32>,
}

impl Tables {
    fn new(model: &mut Model, ring_name: Rc<str>, ring: Rc<Ring>) -> Result<Self> {
        if ring.precision() != 1 || ring.prime() <= 3 {
            return Err(Error::plain("geometric Manin counts require a residue field of characteristic greater than three at precision one"));
        }
        let q = ring.size_of_residue_field()?;
        if q > MAX_FIELD {
            return Err(Error::plain(format!(
                "geometric Manin enumeration requires q ≤ {MAX_FIELD}"
            )));
        }
        let q = q as usize;
        let prime = ring.prime() as usize;
        let elements: Vec<_> = (0..q)
            .map(|mut i| {
                (0..ring.degree())
                    .map(|_| {
                        let coefficient = (i % prime) as u64;
                        i /= prime;
                        coefficient
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        let mut add = Vec::with_capacity(q);
        let mut multiply = Vec::with_capacity(q);
        for a in &elements {
            let mut sum = Vec::with_capacity(q);
            let mut product = Vec::with_capacity(q);
            for b in &elements {
                sum.push(index(&ring.add(model, a, b)?, prime));
                product.push(index(&ring.multiply(model, a, b)?, prime));
            }
            add.push(sum);
            multiply.push(product);
        }
        let mut character = vec![-1; q];
        let mut cubes = vec![0; q];
        for x in 0..q {
            let square = multiply[x][x];
            cubes[x] = multiply[x][square];
            character[square] = 1;
        }
        character[0] = 0;
        let scales = (1..q)
            .map(|u| {
                let square = multiply[u][u];
                let fourth = multiply[square][square];
                (fourth, multiply[fourth][square])
            })
            .collect();
        Ok(Self {
            ring_name,
            ring,
            elements,
            scales,
            scalar_orbits: BTreeMap::new(),
            prime,
            q,
            add,
            multiply,
            cubes,
            character,
        })
    }

    fn nonsingular(&self, a: usize, b: usize) -> bool {
        let u = self.multiply[4 % self.prime][self.cubes[a]];
        let v = self.multiply[27 % self.prime][self.multiply[b][b]];
        self.add[u][v] != 0
    }

    fn offset(&self, a: usize) -> Vec<usize> {
        self.cubes
            .iter()
            .zip(&self.multiply[a])
            .map(|(&c, &ax)| self.add[c][ax])
            .collect()
    }

    fn lines(
        &mut self,
        model: &mut Model,
        a: usize,
        offset: &[usize],
        b: usize,
        p: usize,
    ) -> Result<Lines> {
        model.spend(self.q)?;
        let mut trace = 0;
        let mut two = 0;
        for &c in offset {
            let value = self.character[self.add[c][b]];
            trace -= value;
            two += usize::from(value == 0);
        }
        let delta = trace * trace - 4 * self.q as i32;
        if delta > 0 {
            return Err(certificate("the enumerated elliptic trace satisfies Hasse"));
        }
        // If delta=0, the difference of Frobenius from its repeated
        // eigenvalue has square zero in the division algebra End(E)⊗Q.
        // Otherwise a scalar integral Tate matrix mod p forces p²|delta.
        let scalar = delta == 0
            || (delta % (p * p) as i32 == 0 && self.scalar_frobenius(model, a, b, p, trace)?);
        let level_prime = if scalar {
            p + 1
        } else {
            let discriminant = delta.rem_euclid(p as i32) as usize;
            (1 + quadratic_character(discriminant, p)) as usize
        };
        Ok(Lines {
            trace,
            two,
            level_prime,
            scalar,
        })
    }

    fn scalar_frobenius(
        &mut self,
        model: &mut Model,
        a: usize,
        b: usize,
        p: usize,
        trace: i32,
    ) -> Result<bool> {
        // Short models over an odd residue field of characteristic >3
        // are isomorphic over F_q exactly by these simultaneous scalings.
        model.spend(self.scales.len())?;
        let (a, b) = self
            .scales
            .iter()
            .map(|&(u4, u6)| (self.multiply[u4][a], self.multiply[u6][b]))
            .min()
            .ok_or_else(|| certificate("the isomorphism orbit includes the identity scaling"))?;
        if let Some(&result) = self.scalar_orbits.get(&(a, b)) {
            return Ok(result);
        }
        let curve = Curve::new(
            model,
            self.ring_name.clone(),
            self.ring.clone(),
            self.elements[a].clone(),
            self.elements[b].clone(),
        )?;
        let candidate = (i128::from(trace) * p.div_ceil(2) as i128).rem_euclid(p as i128);
        let proof = elliptic_local::certificate(
            model,
            &curve,
            Test {
                prime: p,
                level: 1,
                candidate: Some(candidate),
            },
        )?;
        let result = proof.checks() == (true, true);
        self.scalar_orbits.insert((a, b), result);
        Ok(result)
    }
}

impl Space {
    pub(super) fn geometric(
        &self,
        model: &mut Model,
        ring_name: Rc<str>,
        ring: Rc<Ring>,
    ) -> Result<String> {
        let ell = ring.prime() as usize;
        require_good(model, self.geometry.prime, ell)?;
        let mut table = Tables::new(model, ring_name, ring.clone())?;
        let mut histogram = BTreeMap::new();
        let mut total = 0usize;
        let mut nonsingular = 0;
        let mut scalar = 0;
        for a in 0..table.q {
            let offset = table.offset(a);
            for b in 0..table.q {
                if !table.nonsingular(a, b) {
                    continue;
                }
                let lines = table.lines(model, a, &offset, b, self.geometry.prime)?;
                *histogram
                    .entry((lines.trace, lines.two, lines.level_prime))
                    .or_insert(0usize) += 1;
                total += lines.two * lines.level_prime;
                nonsingular += 1;
                scalar += usize::from(lines.scalar);
            }
        }
        if nonsingular != table.q * (table.q - 1) || !total.is_multiple_of(table.q - 1) {
            return Err(certificate(
                "the nonsingular model count and weighted quotient are exact",
            ));
        }
        let count = 4 + total / (table.q - 1);
        let prediction = self.frob_data(model, ell, ring.degree())?;
        if prediction.points.signed() != Some(count as i128) {
            return Err(certificate(
                "geometric enumeration equals the Eichler–Shimura Newton count",
            ));
        }
        let certificates = table.scalar_orbits.len();
        let histogram: Vec<_> = histogram
            .into_iter()
            .map(|((t, a, b), n)| (t, a, b, n))
            .collect();
        Ok(format!(
            "X0({}) over F{}: {} points (4 cusps, {} noncuspidal); Eichler–Shimura count agrees\n\
             nonsingular models = {nonsingular}; scalar on E[{}] models = {scalar}; division-scheme certificates = {certificates}; weighted line count = {total}\n\
             (trace, 2-lines, {}-lines, models) = {histogram:?}",
            self.geometry.level,table.q,count,count-4,self.geometry.prime,self.geometry.prime
        ))
    }

    pub(crate) fn geometric_lines(
        &self,
        model: &mut Model,
        ring_name: Rc<str>,
        ring: Rc<Ring>,
        a: &[u64],
        b: &[u64],
    ) -> Result<String> {
        require_good(model, self.geometry.prime, ring.prime() as usize)?;
        let mut table = Tables::new(model, ring_name, ring.clone())?;
        let (a, b) = (index(a, table.prime), index(b, table.prime));
        if !table.nonsingular(a, b) {
            return Err(Error::plain(
                "geometric Manin lines require a nonsingular short Weierstrass model",
            ));
        }
        let offset = table.offset(a);
        let lines = table.lines(model, a, &offset, b, self.geometry.prime)?;
        Ok(format!(
            "F{} elliptic model: trace = {}; #E = {}; rational 2-lines = {}; rational {}-lines = {}",
            table.q,lines.trace,table.q as i32+1-lines.trace,lines.two,self.geometry.prime,lines.level_prime
        ))
    }
}

fn index(coords: &[u64], p: usize) -> usize {
    coords.iter().rev().fold(0, |acc, &c| acc * p + c as usize)
}

fn quadratic_character(a: usize, p: usize) -> i32 {
    if a == 0 {
        0
    } else if (1..p).any(|x| x * x % p == a) {
        1
    } else {
        -1
    }
}
