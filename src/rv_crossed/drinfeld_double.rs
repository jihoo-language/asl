use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_cyclotomic;
use crate::rv_lens::{self, Rational};
use crate::syntax::{Label, ParsedExpr};

use super::path_integral::Quotient;
use super::postnikov::{Cochain, Postnikov};

pub(crate) mod boundaries;
pub(crate) mod braid;
mod characters;
pub(crate) mod defects;
pub(crate) mod fusion_symbols;
mod matrix;
mod representations;

const MAX_SIMPLES: usize = 256;

pub(crate) enum Command {
    Declare(Label, Label, Box<ParsedExpr>),
    Data(Label),
    Sector(Label, usize),
    Phase(Label, [usize; 2]),
    Coefficient(usize, Box<ParsedExpr>),
}

pub(crate) enum Expression {
    Basis(Label, usize),
    Simple(Label, usize, Vec<usize>),
    Linear(Label, Vec<i128>),
    S(Label, [usize; 2]),
    T(Label, usize),
    Character(Label, [usize; 2]),
    Lens(Label, i128),
    Scale(i128, Box<ParsedExpr>),
    Dual(Box<ParsedExpr>),
    Fuse(Box<[ParsedExpr; 2]>),
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Sector(pub(crate) Vec<(usize, i128)>);

impl Sector {
    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(1 + 3 * self.0.len())?;
        Ok(Value::DoubleSector(Rc::new(self)))
    }

    pub(crate) fn coefficient(&self, i: usize) -> i128 {
        self.0
            .binary_search_by_key(&i, |&(key, _)| key)
            .map_or(0, |j| self.0[j].1)
    }

    pub(crate) fn scale(&self, model: &mut Model, factor: i128) -> Result<Self> {
        if factor == 0 {
            return Ok(Self(Vec::new()));
        }
        model.spend(self.0.len())?;
        let data = self
            .0
            .iter()
            .map(|&(i, x)| {
                let y = x.checked_mul(factor).ok_or_else(arithmetic)?;
                Ok((i, y))
            })
            .collect::<Result<_>>()?;
        Ok(Self(data))
    }

    pub(crate) fn add(&self, model: &mut Model, other: &Self) -> Result<Self> {
        model.spend(self.0.len() + other.0.len())?;
        let mut answer = BTreeMap::new();
        for &(i, coefficient) in self.0.iter().chain(&other.0) {
            accumulate(&mut answer, i, coefficient)?;
        }
        Ok(Self(answer.into_iter().filter(|&(_, x)| x != 0).collect()))
    }
}

fn accumulate(out: &mut BTreeMap<usize, i128>, i: usize, coefficient: i128) -> Result<()> {
    let entry = out.entry(i).or_default();
    *entry = entry.checked_add(coefficient).ok_or_else(arithmetic)?;
    Ok(())
}

struct Simple {
    flux: usize,
    dimension: usize,
    phases: Vec<Option<usize>>,
}

#[derive(Clone, Copy)]
struct Monomial {
    phase: usize,
    weight: i128,
}

pub(crate) struct Double {
    pub(crate) conductor: usize,
    group: Quotient,
    simples: Vec<Simple>,
    representations: Vec<representations::Irrep>,
    cocycle: Vec<usize>,
    s: Vec<Vec<Monomial>>,
    twists: Vec<usize>,
    duals: Vec<usize>,
    flux_offsets: Vec<usize>,
}

impl Double {
    pub(crate) fn new(model: &mut Model, post: &Postnikov, k: &Cochain) -> Result<Self> {
        if post.action.iter().any(|&x| x != 1) {
            return Err(Error::plain(
                "twisted double requires trivial coefficient action",
            ));
        }
        post.require_cocycle(model, k)?;
        let group = Quotient::new(model, &post.crossed)?;
        let n = group.order();
        model.spend(n * n)?;
        if (0..n).any(|x| (0..x).any(|y| group.mul(x, y) != group.mul(y, x))) {
            return Err(Error::plain(
                "twisted double currently requires an abelian group",
            ));
        }
        let exponent = characters::exponent(model, &group)?;
        let conductor = post.prime.get() as usize * exponent;
        rv_cyclotomic::check_order(conductor)?;
        let mut simples = Vec::new();
        let mut representations = Vec::new();
        let mut flux_offsets = Vec::with_capacity(n + 1);
        for g in 0..n {
            flux_offsets.push(simples.len());
            let alpha = characters::transgress(model, &group, post, k, g, exponent)?;
            let radical = characters::radical(model, &alpha)?;
            if !radical.contains(&g) {
                return Err(Error::plain("flux is missing from its commutator radical"));
            }
            if simples.len() + radical.len() > MAX_SIMPLES {
                return Err(Error::plain(
                    "twisted double exceeds its simple-sector limit of 256",
                ));
            }
            let sectors = characters::extend(model, &group, g, &alpha, &radical, conductor)?;
            representations.extend(representations::induce_all(
                model, &group, &alpha, &radical, &sectors, conductor,
            )?);
            simples.extend(sectors);
        }
        let mut cocycle = Vec::with_capacity(n * n * n);
        for g in 0..n {
            for h in 0..n {
                for j in 0..n {
                    cocycle.push(exponent * post.coefficient(model, k, &[g, h, j])? as usize);
                }
            }
        }
        flux_offsets.push(simples.len());
        let total: usize = simples.iter().map(|s| s.dimension * s.dimension).sum();
        if total != n * n || simples[0].phases.iter().any(|&x| x != Some(0)) {
            return Err(Error::plain(
                "twisted double fails the dimension or vacuum check",
            ));
        }
        let twists = simples
            .iter()
            .map(|s| s.phases[s.flux].expect("flux is radical"))
            .collect();
        let mut s = Vec::with_capacity(simples.len());
        model.spend(simples.len() * simples.len())?;
        for a in &simples {
            s.push(
                simples
                    .iter()
                    .map(|b| {
                        if let (Some(x), Some(y)) = (a.phases[b.flux], b.phases[a.flux]) {
                            Monomial {
                                phase: negate_phase(conductor, x + y),
                                weight: (a.dimension * b.dimension) as i128,
                            }
                        } else {
                            Monomial {
                                phase: 0,
                                weight: 0,
                            }
                        }
                    })
                    .collect(),
            );
        }
        let mut double = Self {
            conductor,
            group,
            simples,
            representations,
            cocycle,
            s,
            twists,
            duals: Vec::new(),
            flux_offsets,
        };
        double.duals = double.certify(model)?;
        let count = double.simples.len();
        model.allocate(n * n * n + 4 * n * n + count * (9 + 2 * n + 3 * count) + 12)?;
        Ok(double)
    }

    fn cocycle_phase(&self, [g, h, j]: [usize; 3]) -> usize {
        let n = self.group.order();
        self.cocycle[(g * n + h) * n + j]
    }

    fn tensor_phase(&self, x: usize, g: usize, h: usize) -> usize {
        (self.cocycle_phase([x, g, h])
            + self.cocycle_phase([g, h, x])
            + negate_phase(self.conductor, self.cocycle_phase([g, x, h])))
            % self.conductor
    }

    pub(crate) fn count(&self) -> usize {
        self.simples.len()
    }

    pub(crate) fn check_index(&self, index: usize) -> Result<()> {
        if index >= self.count() {
            return Err(Error::plain("twisted double simple index is out of range"));
        }
        Ok(())
    }

    fn simple(&self, index: usize) -> Result<&Simple> {
        self.check_index(index)?;
        Ok(&self.simples[index])
    }

    pub(crate) fn lookup(&self, model: &mut Model, flux: usize, phases: &[usize]) -> Result<usize> {
        if flux >= self.group.order() {
            return Err(Error::plain("twisted double flux index is out of range"));
        }
        let candidates = self.flux_offsets[flux]..self.flux_offsets[flux + 1];
        model.spend(candidates.len() * self.group.order())?;
        for i in candidates {
            if self.simples[i]
                .phases
                .iter()
                .flatten()
                .copied()
                .eq(phases.iter().copied())
            {
                return Ok(i);
            }
        }
        Err(Error::plain(
            "no simple at this flux has the supplied radical phases",
        ))
    }

    pub(crate) fn linear(&self, model: &mut Model, dense: Vec<i128>) -> Result<Sector> {
        if dense.len() != self.count() {
            return Err(Error::plain(format!(
                "twisted double linear value requires {} coefficients",
                self.count()
            )));
        }
        model.spend(dense.len())?;
        Ok(Sector(
            dense
                .into_iter()
                .enumerate()
                .filter(|&(_, x)| x != 0)
                .collect(),
        ))
    }

    pub(crate) fn data(&self, name: &str) -> String {
        let mut dimensions = BTreeMap::new();
        let mut twists = BTreeMap::new();
        for (s, &phase) in self.simples.iter().zip(&self.twists) {
            *dimensions.entry(s.dimension).or_insert(0) += 1;
            *twists.entry(phase).or_insert(0) += 1;
        }
        let order = self.group.order();
        format!(
            "rvDoubleData[{name}] = group_order={order}, conductor={}, simples={}, dimension_counts={dimensions:?}, twist_counts={twists:?}, total_squared_dimension={}, modular_relations=true",
            self.conductor, self.count(), order*order,
        )
    }

    pub(crate) fn sector_report(&self, name: &str, i: usize) -> Result<String> {
        let simple = self.simple(i)?;
        let phases: Vec<_> = simple
            .phases
            .iter()
            .enumerate()
            .filter_map(|(x, phase)| phase.map(|p| (x, p)))
            .collect();
        Ok(format!("rvDoubleSector[{name},{i}] = flux={}, dimension={}, radical_phases={phases:?}, dual={}, twist={}", simple.flux, simple.dimension, self.duals[i], self.twists[i]))
    }

    pub(crate) fn phase(&self, i: usize, x: usize) -> Result<Option<usize>> {
        self.simple(i)?
            .phases
            .get(x)
            .copied()
            .ok_or_else(|| Error::plain("twisted double evaluation element is out of range"))
    }

    pub(crate) fn character(&self, model: &mut Model, i: usize, x: usize) -> Result<Value> {
        let Some(phase) = self.phase(i, x)? else {
            return Ok(rv_lens::zero());
        };
        self.value(
            model,
            Monomial {
                phase,
                weight: self.simple(i)?.dimension as i128,
            },
            1,
        )
    }

    pub(crate) fn s_value(&self, model: &mut Model, i: usize, j: usize) -> Result<Value> {
        self.check_index(i)?;
        self.check_index(j)?;
        self.value(model, self.s[i][j], self.group.order() as u128)
    }

    pub(crate) fn t_value(&self, model: &mut Model, i: usize) -> Result<Value> {
        self.check_index(i)?;
        self.value(
            model,
            Monomial {
                phase: self.twists[i],
                weight: 1,
            },
            1,
        )
    }

    pub(crate) fn lens(&self, model: &mut Model, power: i128) -> Result<Value> {
        let mut numerator = vec![0; self.conductor];
        let power = power.rem_euclid(self.conductor as i128) as usize;
        model.spend(self.count())?;
        for (simple, &phase) in self.simples.iter().zip(&self.twists) {
            numerator[(phase * power) % self.conductor] +=
                (simple.dimension * simple.dimension) as i128;
        }
        self.rational(
            model,
            numerator,
            (self.group.order() * self.group.order()) as u128,
        )
    }

    fn value(&self, model: &mut Model, mono: Monomial, denominator: u128) -> Result<Value> {
        if mono.weight == 0 {
            return Ok(rv_lens::zero());
        }
        let mut numerator = vec![0; mono.phase + 1];
        numerator[mono.phase] = mono.weight;
        self.rational(model, numerator, denominator)
    }

    fn rational(
        &self,
        model: &mut Model,
        numerator: Vec<i128>,
        denominator: u128,
    ) -> Result<Value> {
        let numerator = rv_cyclotomic::reduced(model, self.conductor, numerator)?;
        let result: Rational = rv_lens::fraction(model, &numerator, denominator)?;
        rv_lens::rational_value(model, result)
    }

    pub(crate) fn dual(&self, model: &mut Model, value: &Sector) -> Result<Sector> {
        model.spend(value.0.len())?;
        let mut answer: Vec<_> = value.0.iter().map(|&(i, x)| (self.duals[i], x)).collect();
        answer.sort_unstable_by_key(|&(i, _)| i);
        Ok(Sector(answer))
    }

    pub(crate) fn fuse(&self, model: &mut Model, a: &Sector, b: &Sector) -> Result<Sector> {
        let mut answer = BTreeMap::new();
        for &(i, a) in &a.0 {
            for &(j, b) in &b.0 {
                model.spend(1)?;
                let factor = a.checked_mul(b).ok_or_else(arithmetic)?;
                let g = self.group.mul(self.simples[i].flux, self.simples[j].flux);
                let candidates = self.flux_offsets[g]..self.flux_offsets[g + 1];
                let mut weighted_dimension = 0i128;
                for k in candidates {
                    let multiplicity = self.verlinde(model, i, j, k)?;
                    if multiplicity != 0 {
                        let coeff = factor.checked_mul(multiplicity).ok_or_else(arithmetic)?;
                        accumulate(&mut answer, k, coeff)?;
                        weighted_dimension = multiplicity
                            .checked_mul(self.simples[k].dimension as i128)
                            .and_then(|term| weighted_dimension.checked_add(term))
                            .ok_or_else(arithmetic)?;
                    }
                }
                if weighted_dimension
                    != (self.simples[i].dimension * self.simples[j].dimension) as i128
                {
                    return Err(Error::plain(
                        "twisted double fusion fails quantum-dimension conservation",
                    ));
                }
            }
        }
        Ok(Sector(
            answer.into_iter().filter(|&(_, x)| x != 0).collect(),
        ))
    }

    fn verlinde(&self, model: &mut Model, i: usize, j: usize, k: usize) -> Result<i128> {
        let n2 = (self.group.order() * self.group.order()) as i128;
        let mut numerator = vec![0; self.conductor];
        model.spend(self.count() + self.conductor)?;
        for (l, simple) in self.simples.iter().enumerate() {
            let (a, b, c) = (self.s[i][l], self.s[j][l], self.s[k][l]);
            let phase = (a.phase + b.phase + self.conductor - c.phase) % self.conductor;
            numerator[phase] += a.weight * b.weight * c.weight / simple.dimension as i128;
        }
        let result = rv_cyclotomic::reduced(model, self.conductor, numerator)?;
        if result.is_empty() {
            return Ok(0);
        }
        if let &[value] = result.as_slice() {
            if value >= 0 && value % n2 == 0 {
                return Ok(value / n2);
            }
        }
        Err(Error::plain(
            "Verlinde coefficient is not a nonnegative integer",
        ))
    }

    fn certify(&self, model: &mut Model) -> Result<Vec<usize>> {
        let n = self.group.order() as i128;
        let mut duals = vec![None; self.count()];
        for (i, row) in self.s.iter().enumerate() {
            for (j, column) in self.s.iter().enumerate().skip(i) {
                model.spend(self.count() + 2 * self.conductor)?;
                let mut square = vec![0; self.conductor];
                let mut sts = vec![0; self.conductor];
                for (l, (a, b)) in row.iter().zip(column).enumerate() {
                    let phase = (a.phase + b.phase) % self.conductor;
                    let coeff = a.weight * b.weight;
                    square[phase] += coeff;
                    sts[(phase + self.twists[l]) % self.conductor] += coeff;
                }
                let reduced = rv_cyclotomic::reduced(model, self.conductor, square)?;
                if !reduced.is_empty() {
                    if reduced != [n * n]
                        || duals[i].is_some()
                        || duals[j].is_some()
                        || self.group.mul(self.simples[i].flux, self.simples[j].flux) != 0
                        || self.simples[i].dimension != self.simples[j].dimension
                        || self.twists[i] != self.twists[j]
                    {
                        return Err(Error::plain(
                            "twisted double fails the S² charge-conjugation check",
                        ));
                    }
                    duals[i] = Some(j);
                    duals[j] = Some(i);
                }
                // After S²=C, the modular relation is equivalent to
                // S T S = T^-1 S T^-1. Compare with denominator |G|².
                let s = row[j];
                let expected_phase = (s.phase
                    + negate_phase(self.conductor, self.twists[i] + self.twists[j]))
                    % self.conductor;
                sts[expected_phase] -= n * s.weight;
                if !rv_cyclotomic::reduced(model, self.conductor, sts)?.is_empty() {
                    return Err(Error::plain(
                        "twisted double fails the (ST)³ charge-conjugation check",
                    ));
                }
            }
        }
        if duals[0] != Some(0) {
            return Err(Error::plain("twisted double vacuum is not self-dual"));
        }
        duals
            .into_iter()
            .map(|dual| {
                dual.ok_or_else(|| Error::plain("twisted double has no dual for one simple"))
            })
            .collect()
    }
}

fn negate_phase(conductor: usize, phase: usize) -> usize {
    (conductor - phase % conductor) % conductor
}

fn arithmetic() -> Error {
    Error::plain("twisted double fusion arithmetic exceeds the signed 128-bit coefficient range")
}
