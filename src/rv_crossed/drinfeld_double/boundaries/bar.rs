use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::qz::Scalar;
use crate::integer_smith::{self, Matrix, Track, Vector};
use crate::model::Model;

use super::super::Double;

pub(super) struct Subgroup {
    pub(super) members: Vec<usize>,
    indices: Vec<Option<usize>>,
    products: Vec<Vec<usize>>,
}

pub(super) enum Solution {
    Cochain(Vec<Scalar>),
    Obstruction(Witness),
}

pub(super) struct Witness {
    pub(super) terms: Vec<([usize; 3], i128)>,
    pub(super) value: Scalar,
}

impl Subgroup {
    pub(super) fn new(model: &mut Model, double: &Double, mut members: Vec<usize>) -> Result<Self> {
        if members.len() > double.group.order() {
            return Err(Error::plain(
                "double boundary subgroup has too many elements",
            ));
        }
        members.sort_unstable();
        if members.first() != Some(&0) || members.windows(2).any(|w| w[0] == w[1]) {
            return Err(Error::plain(
                "double boundary subgroup must contain identity and distinct elements",
            ));
        }
        let mut indices = vec![None; double.group.order()];
        for (i, &g) in members.iter().enumerate() {
            *indices.get_mut(g).ok_or_else(|| {
                Error::plain("double boundary subgroup element is out of range")
            })? = Some(i);
        }
        model.spend(members.len() * members.len())?;
        let products = members
            .iter()
            .map(|&g| {
                members
                    .iter()
                    .map(|&h| {
                        indices[double.group.mul(g, h)].ok_or_else(|| {
                            Error::plain(
                                "double boundary subgroup is not closed under multiplication",
                            )
                        })
                    })
                    .collect::<Result<_>>()
            })
            .collect::<Result<_>>()?;
        Ok(Self {
            members,
            indices,
            products,
        })
    }

    pub(super) fn cells(&self) -> usize {
        self.members.len() * (self.members.len() + 4) + 2 * self.indices.len() + 6
    }

    pub(super) fn order(&self) -> usize {
        self.members.len()
    }

    pub(super) fn index(&self, g: usize) -> Result<usize> {
        self.indices
            .get(g)
            .copied()
            .flatten()
            .ok_or_else(|| Error::plain("double boundary argument is outside the subgroup"))
    }

    pub(super) fn contains(&self, g: usize) -> bool {
        self.indices.get(g).is_some_and(Option::is_some)
    }

    fn pair(&self, a: usize, b: usize) -> Option<usize> {
        Some(a.checked_sub(1)? * (self.order() - 1) + b.checked_sub(1)?)
    }

    pub(super) fn mu(&self, cochain: &[Scalar], g: usize, h: usize) -> Result<Scalar> {
        let pair = self.pair(self.index(g)?, self.index(h)?);
        Ok(pair.map_or_else(Scalar::zero, |i| cochain[i]))
    }

    pub(super) fn cochain(
        &self,
        model: &mut Model,
        dense: Vec<(i128, i128)>,
    ) -> Result<Vec<Scalar>> {
        let count = (self.order() - 1).pow(2);
        if dense.len() != count {
            return Err(Error::plain(format!(
                "double boundary requires {count} normalized cochain coefficients"
            )));
        }
        dense
            .into_iter()
            .map(|(num, den)| {
                let den = u128::try_from(den).map_err(|_| {
                    Error::plain("double boundary cochain denominator must be positive")
                })?;
                Scalar::literal(model, num, den)
            })
            .collect()
    }

    fn triple(&self, mut index: usize) -> [usize; 3] {
        let base = self.order() - 1;
        let mut result = [0; 3];
        for i in result.iter_mut().rev() {
            *i = index % base + 1;
            index /= base;
        }
        result
    }

    fn triple_index(&self, global: [usize; 3]) -> Result<usize> {
        let mut out = 0;
        for x in global {
            let digit = self.index(x)?.checked_sub(1).ok_or_else(|| {
                Error::plain("normalized boundary bar cycles use nonidentity triples")
            })?;
            out = out * (self.order() - 1) + digit;
        }
        Ok(out)
    }

    fn row(&self, [a, b, c]: [usize; 3]) -> Vector {
        let ab = self.products[a][b];
        let bc = self.products[b][c];
        let mut out = Vector::new();
        for (key, sign) in [
            (self.pair(b, c), 1),
            (self.pair(ab, c), -1),
            (self.pair(a, bc), 1),
            (self.pair(a, b), -1),
        ] {
            if let Some(i) = key {
                *out.entry(i).or_default() += sign;
            }
        }
        out.retain(|_, value| *value != 0);
        out
    }

    fn rhs(&self, double: &Double) -> Vec<usize> {
        (0..(self.order() - 1).pow(3))
            .map(|i| double.cocycle_phase(self.triple(i).map(|j| self.members[j])))
            .collect()
    }

    pub(super) fn verify(
        &self,
        model: &mut Model,
        double: &Double,
        cochain: &[Scalar],
    ) -> Result<()> {
        let rhs = self.rhs(double);
        model.spend(rhs.len())?;
        for (i, &phase) in rhs.iter().enumerate() {
            let mut delta = Scalar::zero();
            let triple = self.triple(i);
            for (j, c) in self.row(triple) {
                let a = cochain[j];
                if a == Scalar::zero() {
                    continue;
                }
                let term = a.scale(model, c)?;
                delta = delta.add(model, term)?;
            }
            if delta.root_index(double.conductor) != Some(phase) {
                let [a, b, c] = triple.map(|j| self.members[j]);
                return Err(Error::plain(format!(
                    "double boundary cochain fails dµ = k on ({a},{b},{c})"
                )));
            }
        }
        Ok(())
    }

    pub(super) fn solve(&self, model: &mut Model, double: &Double) -> Result<Solution> {
        let rhs = self.rhs(double);
        let count = (self.order() - 1).pow(2);
        model.spend(rhs.len() + count)?;
        let mut cochain = vec![Scalar::zero(); count];
        if rhs.iter().all(|&x| x == 0) {
            return Ok(Solution::Cochain(cochain));
        }
        let mut matrix = Matrix::new(model, rhs.len(), count)?;
        for i in 0..rhs.len() {
            model.spend(1)?;
            for (j, c) in self.row(self.triple(i)) {
                matrix.write(model, i, j, c)?;
            }
        }
        let smith = integer_smith::reduce(model, matrix, Track::Coordinate, Track::Basis)?;
        let order = double.conductor;
        // U A V = D. Positive diagonal entries are surjective on Q/Z;
        // a nonzero obstruction must lie in the rows beyond Smith rank.
        for i in smith.factors.len()..rhs.len() {
            let row = smith.rows.coordinate(i);
            if dot(model, row, &rhs, order)? != 0 {
                let proof = self.witness(model, double, row)?;
                return Ok(Solution::Obstruction(proof));
            }
        }
        for (i, &d) in smith.factors.iter().enumerate() {
            let value = dot(model, smith.rows.coordinate(i), &rhs, order)?;
            if value == 0 {
                continue;
            }
            let denominator = (d as u128).checked_mul(order as u128).ok_or_else(|| {
                Error::plain("double boundary Smith lift denominator exceeds its range")
            })?;
            let lifted = Scalar::literal(model, value as i128, denominator)?;
            for (&j, &coefficient) in smith.columns.basis(i) {
                let term = lifted.scale(model, coefficient)?;
                cochain[j] = cochain[j].add(model, term)?;
            }
        }
        Ok(Solution::Cochain(cochain))
    }

    pub(super) fn explicit_cycle(
        &self,
        model: &mut Model,
        double: &Double,
        terms: Vec<([usize; 3], i128)>,
    ) -> Result<Witness> {
        let mut chain = Vector::new();
        for (triple, coefficient) in terms {
            let i = self.triple_index(triple)?;
            let current = chain.entry(i).or_default();
            *current = current
                .checked_add(coefficient)
                .ok_or_else(integer_smith::arithmetic)?;
        }
        chain.retain(|_, value| *value != 0);
        self.witness(model, double, &chain)
    }

    fn witness(&self, model: &mut Model, double: &Double, chain: &Vector) -> Result<Witness> {
        let mut boundary = Vector::new();
        let mut terms = Vec::with_capacity(chain.len());
        for (&i, &coefficient) in chain {
            let triple = self.triple(i);
            integer_smith::add_scaled(model, &mut boundary, &self.row(triple), coefficient)?;
            terms.push((triple.map(|j| self.members[j]), coefficient));
        }
        if !boundary.is_empty() {
            return Err(Error::plain(
                "supplied integral subgroup bar chain is not a cycle",
            ));
        }
        let rhs = self.rhs(double);
        model.spend(rhs.len())?;
        let value = dot(model, chain, &rhs, double.conductor)?;
        Ok(Witness {
            terms,
            value: Scalar::literal(model, value as i128, double.conductor as u128)?,
        })
    }
}

fn dot(model: &mut Model, row: &Vector, rhs: &[usize], order: usize) -> Result<usize> {
    model.spend(row.len())?;
    Ok(row.iter().fold(0, |sum, (&i, &c)| {
        (sum + c.rem_euclid(order as i128) as usize * rhs[i]) % order
    }))
}
