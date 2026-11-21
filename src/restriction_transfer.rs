use std::collections::HashMap;
use std::hash::Hash;
use std::rc::Rc;

use crate::adams;
use crate::cobar::{self, Word};
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::linear;
use crate::model::{Model, Value};
use crate::profile::Profile;
use crate::transfer::{self, odd_stem, Bidegree, Input, Intervals};

pub(crate) fn bidegree(
    inputs: impl IntoIterator<Item = Bidegree>,
    arity: usize,
) -> Result<Bidegree> {
    let shift = arity
        .checked_sub(1)
        .ok_or_else(|| Error::plain("resHigher requires at least one Ext class"))?;
    transfer::sum_degrees(inputs, shift, "resHigher")
}

pub(crate) fn evaluate(
    model: &mut Model,
    source: Profile,
    target: Profile,
    inputs: &[Input],
) -> Result<Value> {
    let arity = inputs.len();
    let (s, t) = bidegree(inputs.iter().map(|&(degree, _)| degree), arity)?;
    model.spend(arity)?;
    if inputs.iter().any(|(_, coordinates)| coordinates.is_empty()) {
        return Ok(Value::Class([].into()));
    }
    if let [(_, class)] = inputs {
        if source == target {
            return Ok(Value::Class(class.clone()));
        }
        let from = model.profile_ext_type(source, s, t)?;
        let to = model.profile_ext_type(target, s, t)?;
        return adams::restrict(model, from, class, to);
    }
    if source == target || inputs.iter().any(|((s, _), _)| *s == 0) {
        return Ok(Value::Class([].into()));
    }

    check_partitions(model, arity)?;
    let lifts = Intervals::build(model, source, inputs, arity)?;
    let mut projection = TensorProjection::new(target);
    let coordinates = projection.project_blocks(model, arity, |model, start, end| {
        let words = target.restrict(model, lifts.lift(start, end))?;
        Ok((lifts.grade(start, end)?, words))
    })?;
    model.allocate(coordinates.len())?;
    Ok(Value::Class(coordinates.into()))
}

pub(crate) fn check_partitions(model: &mut Model, arity: usize) -> Result<()> {
    // Check the 2^(k-1) consecutive partitions before allocating the
    // polynomial source lift table. Zero branches can be skipped later.
    let partitions = arity
        .checked_sub(1)
        .and_then(|exponent| u32::try_from(exponent).ok())
        .and_then(|exponent| 1usize.checked_shl(exponent))
        .ok_or_else(partition_limit)?;
    model.spend(partitions).map_err(|_| partition_limit())
}

fn partition_limit() -> Error {
    Error::plain("resHigher exceeds the consecutive partition work limit")
}

type ChainId = usize;
type Coordinates<T> = Rc<[T]>;
type Splitting = (Option<ChainId>, Option<ChainId>);
type CochainProjection<T, C> = (Vec<T>, Vec<C>);
type CochainSplit<T> = (Vec<T>, Vec<T>);

pub(crate) trait TensorContraction {
    type Term: Eq + Hash;
    type Coordinate;
    const COORDINATE_CELLS: usize = 1;

    fn cells(term: &Self::Term) -> usize;
    fn project(
        &self,
        model: &mut Model,
        grade: Bidegree,
        cochain: &[Self::Term],
    ) -> Result<CochainProjection<Self::Term, Self::Coordinate>>;
    fn split_closed(
        &self,
        model: &mut Model,
        grade: Bidegree,
        closed: &[Self::Term],
        coordinates: &[Self::Coordinate],
    ) -> Result<CochainSplit<Self::Term>>;
    fn multiply(
        &self,
        model: &mut Model,
        left: &[Self::Term],
        right: &[Self::Term],
    ) -> Result<Vec<Self::Term>>;
    fn add(
        &self,
        model: &mut Model,
        sum: &[Self::Coordinate],
        term: &[Self::Coordinate],
        negative: bool,
    ) -> Result<Vec<Self::Coordinate>>;
}

impl TensorContraction for Profile {
    type Term = Word;
    type Coordinate = usize;

    fn cells(term: &Word) -> usize {
        cobar::word_cells(term)
    }

    fn project(
        &self,
        model: &mut Model,
        grade: Bidegree,
        cochain: &[Word],
    ) -> Result<(Vec<Word>, Vec<usize>)> {
        transfer::project(model, *self, grade, cochain)
    }

    fn split_closed(
        &self,
        model: &mut Model,
        grade: Bidegree,
        closed: &[Word],
        coordinates: &[usize],
    ) -> Result<(Vec<Word>, Vec<Word>)> {
        transfer::split_closed(model, *self, grade, closed, coordinates)
    }

    fn multiply(&self, model: &mut Model, left: &[Word], right: &[Word]) -> Result<Vec<Word>> {
        cobar::concat_temporary(model, left, right)
    }

    fn add(
        &self,
        model: &mut Model,
        sum: &[usize],
        term: &[usize],
        _negative: bool,
    ) -> Result<Vec<usize>> {
        linear::xor(model, sum, term)
    }
}

struct Chain<C: TensorContraction> {
    grade: Bidegree,
    terms: Rc<[C::Term]>,
    coordinates: Option<Coordinates<C::Coordinate>>,
    splitting: Option<Splitting>,
}

pub(crate) struct TensorProjection<C: TensorContraction> {
    contraction: C,
    chains: Vec<Chain<C>>,
    chain_ids: HashMap<Rc<[C::Term]>, ChainId>,
    tensors: HashMap<Vec<ChainId>, Coordinates<C::Coordinate>>,
    cells: usize,
}

impl<C: TensorContraction> TensorProjection<C> {
    pub(crate) fn new(contraction: C) -> Self {
        Self {
            contraction,
            chains: Vec::new(),
            chain_ids: HashMap::new(),
            tensors: HashMap::new(),
            cells: 0,
        }
    }

    fn retain(&mut self, added: usize) -> Result<()> {
        self.cells += added;
        if self.cells > MAX_PIVOT_CELLS {
            return Err(Error::plain("resHigher exceeds the tensor storage limit"));
        }
        Ok(())
    }

    fn intern(
        &mut self,
        model: &mut Model,
        grade: Bidegree,
        terms: Vec<C::Term>,
    ) -> Result<Option<ChainId>> {
        if terms.is_empty() {
            return Ok(None);
        }
        let cells: usize = terms.iter().map(C::cells).sum();
        model.spend(cells)?;
        if let Some(&id) = self.chain_ids.get(terms.as_slice()) {
            debug_assert_eq!(self.chains[id].grade, grade);
            return Ok(Some(id));
        }
        self.retain(1 + cells)?;
        let id = self.chains.len();
        let terms: Rc<[C::Term]> = terms.into();
        self.chain_ids.insert(terms.clone(), id);
        self.chains.push(Chain {
            grade,
            terms,
            coordinates: None,
            splitting: None,
        });
        Ok(Some(id))
    }

    pub(crate) fn project_blocks(
        &mut self,
        model: &mut Model,
        arity: usize,
        mut lift: impl FnMut(&mut Model, usize, usize) -> Result<(Bidegree, Vec<C::Term>)>,
    ) -> Result<Vec<C::Coordinate>> {
        let mut restricted = Vec::with_capacity(arity);
        for start in 0..arity {
            let mut row = Vec::with_capacity(arity - start);
            for end in start + 1..=arity {
                let (grade, terms) = lift(model, start, end)?;
                row.push(self.intern(model, grade, terms)?);
            }
            restricted.push(row);
        }
        let mut coordinates = Vec::new();
        self.partitions(model, &restricted, &mut Vec::new(), 0, &mut coordinates)?;
        Ok(coordinates)
    }

    fn partitions(
        &mut self,
        model: &mut Model,
        restricted: &[Vec<Option<ChainId>>],
        tensor: &mut Vec<ChainId>,
        start: usize,
        result: &mut Vec<C::Coordinate>,
    ) -> Result<()> {
        if start == restricted.len() {
            let projected = self.project_tensor(model, tensor)?;
            *result = self.contraction.add(model, result, &projected, false)?;
            return Ok(());
        }
        for (offset, &chain) in restricted[start].iter().enumerate() {
            model.spend(1)?;
            if let Some(chain) = chain {
                tensor.push(chain);
                self.partitions(model, restricted, tensor, start + offset + 1, result)?;
                tensor.pop();
            }
        }
        Ok(())
    }

    fn project_tensor(
        &mut self,
        model: &mut Model,
        tensor: &[ChainId],
    ) -> Result<Coordinates<C::Coordinate>> {
        model.spend(1 + tensor.len())?;
        if let [chain] = *tensor {
            return self.project_chain(model, chain);
        }
        if let Some(result) = self.tensors.get(tensor) {
            return Ok(result.clone());
        }
        let mut result = Vec::new();
        let mut transformed = tensor.to_vec();
        let mut reduced = Vec::with_capacity(tensor.len() - 1);
        let mut homotopy_prefix = false;
        // P_j = -P_(j-1) δ H_0. The adjacent products lower the tensor
        // length, so the recursion reaches the single-cochain projection.
        for (index, &chain) in tensor.iter().enumerate() {
            let (projection, homotopy) = self.split_chain(model, chain)?;
            if let Some(homotopy) = homotopy {
                // This summand of H0 has ip on the already processed prefix,
                // h in this slot, and the original entries in the suffix.
                transformed[index] = homotopy;
                let mut product_prefix = false;
                for pair in 0..tensor.len() - 1 {
                    let left_odd = odd_stem(self.chains[transformed[pair]].grade);
                    let product = self.multiply(model, transformed[pair], transformed[pair + 1])?;
                    reduced.clear();
                    reduced.extend_from_slice(&transformed[..pair]);
                    reduced.push(product);
                    reduced.extend_from_slice(&transformed[pair + 2..]);
                    let projected = self.project_tensor(model, &reduced)?;
                    // Include the global minus, the two suspended tensor
                    // prefixes, and the raw left stem from b_2 itself.
                    let negative = !(homotopy_prefix ^ product_prefix ^ left_odd);
                    result = self.contraction.add(model, &result, &projected, negative)?;
                    product_prefix ^= !left_odd;
                }
            }
            let Some(projection) = projection else {
                break;
            };
            transformed[index] = projection;
            homotopy_prefix ^= !odd_stem(self.chains[chain].grade);
        }
        self.retain(1 + tensor.len() + C::COORDINATE_CELLS * result.len())?;
        let result: Coordinates<C::Coordinate> = result.into();
        self.tensors.insert(tensor.to_vec(), result.clone());
        Ok(result)
    }

    fn project_chain(
        &mut self,
        model: &mut Model,
        id: ChainId,
    ) -> Result<Coordinates<C::Coordinate>> {
        let chain = &self.chains[id];
        if let Some(coordinates) = &chain.coordinates {
            return Ok(coordinates.clone());
        }
        let (_, coordinates) = self.contraction.project(model, chain.grade, &chain.terms)?;
        self.save_coordinates(id, coordinates)
    }

    fn save_coordinates(
        &mut self,
        id: ChainId,
        coordinates: Vec<C::Coordinate>,
    ) -> Result<Coordinates<C::Coordinate>> {
        if let Some(stored) = &self.chains[id].coordinates {
            return Ok(stored.clone());
        }
        self.retain(C::COORDINATE_CELLS * coordinates.len())?;
        let coordinates: Coordinates<C::Coordinate> = coordinates.into();
        self.chains[id].coordinates = Some(coordinates.clone());
        Ok(coordinates)
    }

    fn split_chain(&mut self, model: &mut Model, id: ChainId) -> Result<Splitting> {
        let chain = &self.chains[id];
        if let Some(splitting) = chain.splitting {
            return Ok(splitting);
        }
        let grade = chain.grade;
        let (closed, coordinates) = self.contraction.project(model, grade, &chain.terms)?;
        let coordinates = self.save_coordinates(id, coordinates)?;
        let (representative, homotopy) =
            self.contraction
                .split_closed(model, grade, &closed, &coordinates)?;
        let projection = self.intern(model, grade, representative)?;
        debug_assert!(grade.0 != 0 || homotopy.is_empty());
        let homotopy = self.intern(model, (grade.0.saturating_sub(1), grade.1), homotopy)?;
        let splitting = (projection, homotopy);
        self.chains[id].splitting = Some(splitting);
        Ok(splitting)
    }

    fn multiply(&mut self, model: &mut Model, left: ChainId, right: ChainId) -> Result<ChainId> {
        let left = &self.chains[left];
        let right = &self.chains[right];
        let overflow = || Error::plain("resHigher intermediate cobar bidegree is out of range");
        let s = left
            .grade
            .0
            .checked_add(right.grade.0)
            .ok_or_else(overflow)?;
        let t = left
            .grade
            .1
            .checked_add(right.grade.1)
            .ok_or_else(overflow)?;
        let terms = self
            .contraction
            .multiply(model, &left.terms, &right.terms)?;
        // Both factors were nonempty, and the common left length determines
        // a unique split of every concatenated word.
        Ok(self.intern(model, (s, t), terms)?.unwrap())
    }
}
