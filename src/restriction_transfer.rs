use std::collections::HashMap;
use std::rc::Rc;

use crate::adams;
use crate::cobar::{self, Word};
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::linear;
use crate::model::{Model, Value};
use crate::profile::Profile;
use crate::transfer::{self, Bidegree, Input, Intervals};

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

    // A k-entry tensor has 2^(k-1) consecutive partitions. Check this
    // before allocating even the polynomial source lift table. Partitions
    // themselves are visited one at a time and branches with zero lifts end.
    let partitions = u32::try_from(arity - 1)
        .ok()
        .and_then(|exponent| 1usize.checked_shl(exponent))
        .ok_or_else(partition_limit)?;
    model.spend(partitions).map_err(|_| partition_limit())?;

    let lifts = Intervals::build(model, source, inputs, arity)?;
    let mut projection = TensorProjection::new(target);
    let mut restricted = Vec::with_capacity(arity);
    for start in 0..arity {
        let mut row = Vec::with_capacity(arity - start);
        for end in start + 1..=arity {
            let words = target.restrict(model, lifts.lift(start, end))?;
            row.push(projection.intern(model, lifts.grade(start, end)?, words)?);
        }
        restricted.push(row);
    }
    let mut coordinates = Vec::new();
    projection.partitions(model, &restricted, &mut Vec::new(), 0, &mut coordinates)?;
    model.allocate(coordinates.len())?;
    Ok(Value::Class(coordinates.into()))
}

fn partition_limit() -> Error {
    Error::plain("resHigher exceeds the consecutive partition work limit")
}

type ChainId = usize;
type Coordinates = Rc<[usize]>;
type Splitting = (Option<ChainId>, Option<ChainId>);

struct Chain {
    grade: Bidegree,
    words: Rc<[Word]>,
    coordinates: Option<Coordinates>,
    splitting: Option<Splitting>,
}

struct TensorProjection {
    profile: Profile,
    chains: Vec<Chain>,
    chain_ids: HashMap<Rc<[Word]>, ChainId>,
    tensors: HashMap<Vec<ChainId>, Coordinates>,
    cells: usize,
}

impl TensorProjection {
    fn new(profile: Profile) -> Self {
        Self {
            profile,
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
        words: Vec<Word>,
    ) -> Result<Option<ChainId>> {
        if words.is_empty() {
            return Ok(None);
        }
        let cells: usize = words.iter().map(cobar::word_cells).sum();
        model.spend(cells)?;
        if let Some(&id) = self.chain_ids.get(words.as_slice()) {
            debug_assert_eq!(self.chains[id].grade, grade);
            return Ok(Some(id));
        }
        self.retain(1 + cells)?;
        let id = self.chains.len();
        let words: Rc<[Word]> = words.into();
        self.chain_ids.insert(words.clone(), id);
        self.chains.push(Chain {
            grade,
            words,
            coordinates: None,
            splitting: None,
        });
        Ok(Some(id))
    }

    fn partitions(
        &mut self,
        model: &mut Model,
        restricted: &[Vec<Option<ChainId>>],
        tensor: &mut Vec<ChainId>,
        start: usize,
        result: &mut Vec<usize>,
    ) -> Result<()> {
        if start == restricted.len() {
            let projected = self.project_tensor(model, tensor)?;
            *result = linear::xor(model, result, &projected)?;
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

    fn project_tensor(&mut self, model: &mut Model, tensor: &[ChainId]) -> Result<Coordinates> {
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
        // P_j = P_(j-1) b_2 H_0. The adjacent products lower the tensor
        // length, so the recursion reaches the single-cochain projection.
        for (index, &chain) in tensor.iter().enumerate() {
            let (projection, homotopy) = self.split_chain(model, chain)?;
            if let Some(homotopy) = homotopy {
                // This summand of H0 has ip on the already processed prefix,
                // h in this slot, and the original entries in the suffix.
                transformed[index] = homotopy;
                for pair in 0..tensor.len() - 1 {
                    let product = self.multiply(model, transformed[pair], transformed[pair + 1])?;
                    reduced.clear();
                    reduced.extend_from_slice(&transformed[..pair]);
                    reduced.push(product);
                    reduced.extend_from_slice(&transformed[pair + 2..]);
                    let projected = self.project_tensor(model, &reduced)?;
                    result = linear::xor(model, &result, &projected)?;
                }
            }
            let Some(projection) = projection else {
                break;
            };
            transformed[index] = projection;
        }
        self.retain(1 + tensor.len() + result.len())?;
        let result: Coordinates = result.into();
        self.tensors.insert(tensor.to_vec(), result.clone());
        Ok(result)
    }

    fn project_chain(&mut self, model: &mut Model, id: ChainId) -> Result<Coordinates> {
        let chain = &self.chains[id];
        if let Some(coordinates) = &chain.coordinates {
            return Ok(coordinates.clone());
        }
        let (_, coordinates) = transfer::project(model, self.profile, chain.grade, &chain.words)?;
        self.save_coordinates(id, coordinates)
    }

    fn save_coordinates(&mut self, id: ChainId, coordinates: Vec<usize>) -> Result<Coordinates> {
        if let Some(stored) = &self.chains[id].coordinates {
            return Ok(stored.clone());
        }
        self.retain(coordinates.len())?;
        let coordinates: Coordinates = coordinates.into();
        self.chains[id].coordinates = Some(coordinates.clone());
        Ok(coordinates)
    }

    fn split_chain(&mut self, model: &mut Model, id: ChainId) -> Result<Splitting> {
        let chain = &self.chains[id];
        if let Some(splitting) = chain.splitting {
            return Ok(splitting);
        }
        let grade = chain.grade;
        let (closed, coordinates) = transfer::project(model, self.profile, grade, &chain.words)?;
        let coordinates = self.save_coordinates(id, coordinates)?;
        let (representative, homotopy) =
            transfer::split_closed(model, self.profile, grade, &closed, &coordinates)?;
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
        let words = cobar::concat_temporary(model, &left.words, &right.words)?;
        // Both factors were nonempty, and the common left length determines
        // a unique split of every concatenated word.
        Ok(self.intern(model, (s, t), words)?.unwrap())
    }
}
