use std::ops::ControlFlow;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Poset, TypeId, Value};
use crate::set_family::{self, Diagram};

type Sieve = Vec<u64>;

fn contains(sieve: &[u64], point: usize) -> bool {
    sieve[point / 64] & (1 << (point % 64)) != 0
}

fn insert(sieve: &mut [u64], point: usize) {
    sieve[point / 64] |= 1 << (point % 64);
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Top,
    Bottom,
    And,
    Or,
    Implies,
    Not,
    Exists,
    Forall,
    Classify,
    Inclusion,
}

impl Operation {
    pub(crate) fn detect(name: &str) -> Option<Self> {
        Some(match name {
            "setTop" => Self::Top,
            "setBottom" => Self::Bottom,
            "setAnd" => Self::And,
            "setOr" => Self::Or,
            "setImplies" => Self::Implies,
            "setNot" => Self::Not,
            "setExists" => Self::Exists,
            "setForall" => Self::Forall,
            "setClassify" => Self::Classify,
            "setInclusion" => Self::Inclusion,
            _ => return None,
        })
    }

    pub(crate) fn label_count(self) -> usize {
        match self {
            Self::Top | Self::Bottom => 2,
            Self::Classify | Self::Inclusion => 1,
            _ => 0,
        }
    }

    pub(crate) fn arity(self) -> usize {
        match self {
            Self::Top | Self::Bottom | Self::Inclusion => 0,
            Self::Classify | Self::Not => 1,
            _ => 2,
        }
    }
}

#[derive(Clone)]
pub(crate) enum Origin {
    Omega(Rc<Omega>),
    Subobject(Rc<Subobject>),
}

pub(crate) struct Omega {
    pub(crate) base: TypeId,
    site: Rc<Poset>,
    futures: Vec<Vec<usize>>,
    masks: Vec<Sieve>,
    sieves: Vec<Vec<Sieve>>,
    words: usize,
}

impl Omega {
    pub(crate) fn build(model: &mut Model, base: TypeId) -> Result<(Diagram, Self)> {
        let site = model.poset(base).expect("named space");
        let n = site.names.len();
        let words = n.div_ceil(64);
        let mut futures = Vec::with_capacity(n);
        let mut masks = Vec::with_capacity(n);
        model.spend(n * n)?;
        for x in 0..n {
            let mut mask = vec![0; words];
            let points = (0..n)
                .filter(|&y| x == y || site.strict_less(x, y))
                .collect::<Vec<_>>();
            for &y in &points {
                insert(&mut mask, y);
            }
            futures.push(points);
            masks.push(mask);
        }

        // All strict successors occur later in this order.
        let mut order: Vec<_> = (0..n).collect();
        order.sort_unstable_by_key(|&x| (std::cmp::Reverse(futures[x].len()), x));
        let mut sieves = Vec::with_capacity(n);
        let mut count = 0;
        let mut cells = 3 * n + n * words + futures.iter().map(Vec::len).sum::<usize>();
        for points in &futures {
            let mut positions = vec![None; n];
            for (index, &point) in points.iter().enumerate() {
                positions[point] = Some(index);
            }
            let order: Vec<_> = order.iter().filter_map(|&x| positions[x]).collect();
            let domains = vec![2; points.len()];
            let mut rows = Vec::new();
            set_family::solve(
                model,
                &domains,
                &order,
                |model, x, value, pending| {
                    if value == 1 {
                        let upper = &futures[points[x]];
                        model.spend(upper.len())?;
                        pending.extend(
                            upper
                                .iter()
                                .map(|&p| (positions[p].expect("upward scope"), 1)),
                        );
                    }
                    Ok(())
                },
                |row| {
                    count += 1;
                    if count > set_family::MAX_OCCURRENCES {
                        return Err(Error::plain(
                            "truth classifier exceeds the limit of 4096 cosieves",
                        ));
                    }
                    cells += words + 1;
                    set_family::check_cells(cells)?;
                    let mut sieve = vec![0; words];
                    for (&point, &value) in points.iter().zip(row) {
                        if value == 1 {
                            insert(&mut sieve, point);
                        }
                    }
                    rows.push(sieve);
                    Ok(ControlFlow::Continue(()))
                },
            )?;
            let steps = rows.len().max(1).ilog2() as usize + 1;
            model.spend(rows.len() * (words + 1) * steps)?;
            rows.sort_unstable();
            sieves.push(rows);
        }

        let omega = Self {
            base,
            site,
            futures,
            masks,
            sieves,
            words,
        };
        let sizes = omega.sieves.iter().map(Vec::len).collect();
        let diagram = Diagram::generated(model, base, sizes, |model, x, y| {
            let mut restricted = vec![0; words];
            omega.sieves[x]
                .iter()
                .map(|sieve| {
                    model.spend(words)?;
                    for ((out, &u), &mask) in restricted.iter_mut().zip(sieve).zip(&omega.masks[y])
                    {
                        *out = u & mask;
                    }
                    omega.index(model, y, &restricted)
                })
                .collect()
        })?;
        model.allocate(cells)?;
        Ok((diagram, omega))
    }

    fn index(&self, model: &mut Model, x: usize, sieve: &[u64]) -> Result<usize> {
        let steps = self.sieves[x].len().ilog2() as usize + 1;
        model.spend((self.words + 1) * steps)?;
        self.sieves[x]
            .binary_search_by(|row| row.as_slice().cmp(sieve))
            .map_err(|_| Error::plain("computed truth value is not an upward subset"))
    }

    fn sieve(&self, point: usize, index: usize) -> Result<&[u64]> {
        self.sieves[point]
            .get(index)
            .map(Vec::as_slice)
            .ok_or_else(|| Error::plain("cosieve index is outside its fiber"))
    }

    pub(crate) fn members(
        &self,
        model: &mut Model,
        point: usize,
        index: usize,
    ) -> Result<Vec<&str>> {
        let sieve = self.sieve(point, index)?;
        model.spend(self.futures[point].len())?;
        Ok(self.futures[point]
            .iter()
            .filter(|&&y| contains(sieve, y))
            .map(|&y| self.site.names[y].as_str())
            .collect())
    }

    fn holds(&self, point: usize, ordinal: usize) -> bool {
        contains(&self.sieves[point][ordinal], point)
    }

    pub(crate) fn constant(&self, model: &mut Model, family: &Diagram, top: bool) -> Result<Value> {
        model.spend(family.count() + self.futures.len())?;
        model.allocate(family.count())?;
        // The empty and entire future are respectively the least and greatest masks.
        Ok(Value::SetNat(
            (0..self.futures.len())
                .flat_map(|point| {
                    let ordinal = if top { self.sieves[point].len() - 1 } else { 0 };
                    std::iter::repeat_n(ordinal, family.size(point))
                })
                .collect(),
        ))
    }

    pub(crate) fn connect(
        &self,
        model: &mut Model,
        family: &Diagram,
        operation: Operation,
        left: &[usize],
        right: Option<&[usize]>,
    ) -> Result<Value> {
        let mut sieve = vec![0; self.words];
        let mut rows = Vec::with_capacity(family.count());
        for (x, future) in self.futures.iter().enumerate() {
            let right = right.map(|right| family.components(x, right));
            for (i, &u) in family.components(x, left).iter().enumerate() {
                let u = &self.sieves[x][u];
                let v = &self.sieves[x][right.map_or(0, |r| r[i])];
                match operation {
                    Operation::And | Operation::Or => {
                        model.spend(self.words + 1)?;
                        for ((out, &u), &v) in sieve.iter_mut().zip(u).zip(v) {
                            *out = if matches!(operation, Operation::And) {
                                u & v
                            } else {
                                u | v
                            };
                        }
                    }
                    Operation::Implies | Operation::Not => {
                        model.spend((self.words + 1) * (future.len() + 1))?;
                        sieve.fill(0);
                        for &y in future {
                            if u.iter()
                                .zip(v)
                                .zip(&self.masks[y])
                                .all(|((&u, &v), &mask)| (u & !v & mask) == 0)
                            {
                                insert(&mut sieve, y);
                            }
                        }
                    }
                    _ => unreachable!("propositional connective"),
                }
                rows.push(self.index(model, x, &sieve)?);
            }
        }
        model.allocate(rows.len())?;
        Ok(Value::SetNat(rows.into()))
    }

    fn predicate_from_members(
        &self,
        model: &mut Model,
        family: &Diagram,
        members: &[bool],
    ) -> Result<Value> {
        let mut rows = Vec::with_capacity(family.count());
        let mut sieve = vec![0; self.words];
        for (x, future) in self.futures.iter().enumerate() {
            for i in 0..family.size(x) {
                model.spend(future.len() + self.words)?;
                sieve.fill(0);
                for &y in future {
                    let j = family.apply(x, y, i);
                    if members[family.occurrence(y, j)] {
                        insert(&mut sieve, y);
                    }
                }
                rows.push(self.index(model, x, &sieve)?);
            }
        }
        model.allocate(rows.len())?;
        Ok(Value::SetNat(rows.into()))
    }

    pub(crate) fn classify(
        &self,
        model: &mut Model,
        sub: &Diagram,
        family: &Diagram,
        mono: &[usize],
    ) -> Result<Value> {
        let mut members = vec![false; family.count()];
        model.spend(sub.count() + self.futures.len())?;
        for (x, name) in self.site.names.iter().enumerate() {
            for &i in sub.components(x, mono) {
                let position = family.occurrence(x, i);
                if std::mem::replace(&mut members[position], true) {
                    return Err(Error::plain(format!(
                        "setClassify requires an injective component at '{}'",
                        name
                    )));
                }
            }
        }
        self.predicate_from_members(model, family, &members)
    }

    pub(crate) fn quantify(
        &self,
        model: &mut Model,
        source: &Diagram,
        target: &Diagram,
        map: &[usize],
        predicate: &[usize],
        universal: bool,
    ) -> Result<Value> {
        let mut local = vec![universal; target.count()];
        model.spend(source.count() + self.futures.len())?;
        for x in 0..self.futures.len() {
            for (&b, &p) in source
                .components(x, map)
                .iter()
                .zip(source.components(x, predicate))
            {
                if self.holds(x, p) != universal {
                    local[target.occurrence(x, b)] = !universal;
                }
            }
        }
        if !universal {
            return self.predicate_from_members(model, target, &local);
        }

        // A true component must remain true after every future transport.
        let mut members = vec![false; target.count()];
        for (x, future) in self.futures.iter().enumerate() {
            for i in 0..target.size(x) {
                model.spend(future.len())?;
                members[target.occurrence(x, i)] = future.iter().all(|&y| {
                    let j = target.apply(x, y, i);
                    local[target.occurrence(y, j)]
                });
            }
        }
        self.predicate_from_members(model, target, &members)
    }
}

pub(crate) struct Subobject {
    pub(crate) family: Rc<str>,
    inclusion: Rc<[usize]>,
}

impl Subobject {
    pub(crate) fn build(
        model: &mut Model,
        family: Rc<str>,
        domain: &Diagram,
        omega: &Omega,
        predicate: &[usize],
    ) -> Result<(Diagram, Self)> {
        let mut fibers = Vec::new();
        let mut positions = vec![None; domain.count()];
        model.spend(domain.count() + domain.site.names.len())?;
        for x in 0..domain.site.names.len() {
            let mut members = Vec::new();
            for (i, &truth) in domain.components(x, predicate).iter().enumerate() {
                if omega.holds(x, truth) {
                    positions[domain.occurrence(x, i)] = Some(members.len());
                    members.push(i);
                }
            }
            fibers.push(members);
        }
        let sizes = fibers.iter().map(Vec::len).collect();
        let diagram = Diagram::generated(model, domain.base, sizes, |model, x, y| {
            model.spend(fibers[x].len())?;
            fibers[x]
                .iter()
                .map(|&i| {
                    let j = domain.apply(x, y, i);
                    positions[domain.occurrence(y, j)].ok_or_else(|| {
                        Error::plain("predicate comprehension is not stable under transport")
                    })
                })
                .collect()
        })?;
        model.allocate(diagram.count() + 1)?;
        let inclusion = fibers.into_iter().flatten().collect();
        Ok((diagram, Self { family, inclusion }))
    }

    pub(crate) fn inclusion(&self) -> Value {
        Value::SetNat(self.inclusion.clone())
    }
}
