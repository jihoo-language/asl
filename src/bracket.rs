use std::collections::HashMap;
use std::rc::Rc;

use crate::cohomology::{self, Cohomology};
use crate::complex::{xor_sorted, MAX_PIVOT_CELLS};
use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};

pub(crate) fn evaluate(
    model: &mut Model,
    target: TypeId,
    inputs: &[(TypeId, Rc<[usize]>)],
) -> Result<Value> {
    let size = model.elements(target)?.len();
    model.allocate(1 + size)?;
    let mut rows = vec![Value::Point(0); size];
    if let [(left_ty, left), (right_ty, right)] = inputs {
        let product = cohomology::cup(model, *left_ty, left, *right_ty, right, target)?;
        rows[model.index(target, &product)?] = Value::Point(1);
        return Ok(Value::Map(rows.into()));
    }
    if inputs.len() == 4 && inputs.iter().all(|(_, class)| class.is_empty()) {
        mark_zero_quadruple(model, target, inputs, &mut rows)?;
        return Ok(Value::Map(rows.into()));
    }

    let Kind::Cohomology(base, _) = *model.kind(target) else {
        unreachable!()
    };
    let target_group = model.cohomology(target)?;
    let mut search = Search::new(model, base, inputs)?;
    if inputs.len() == 4 {
        if !search.fourfold_compatible(model, inputs)? {
            return Ok(Value::Map(rows.into()));
        }
        if size == 1 {
            rows[0] = Value::Point(1);
            return Ok(Value::Map(rows.into()));
        }
    }
    let first_witness = inputs.len();
    let top = search.intervals.len() - 1;
    let mut next = first_witness;
    let mut remaining = size;
    loop {
        model.spend(1)?;
        if next == top {
            let rhs = search.sum_products(model, top, &target_group)?;
            let coordinates = target_group.reduce_coordinates(model, rhs)?;
            let result = model.index(target, &Value::Class(coordinates.into()))?;
            if matches!(rows[result], Value::Point(0)) {
                rows[result] = Value::Point(1);
                remaining -= 1;
                if remaining == 0 {
                    break;
                }
            }
        } else if search.choose_primitive(model, next)? {
            next += 1;
            continue;
        }

        // A completed system or a non-exact right-hand side advances the last
        // witness with a remaining choice. Later witnesses have to be rebuilt.
        loop {
            if next == first_witness {
                return Ok(Value::Map(rows.into()));
            }
            next -= 1;
            if search.advance_choice(model, next)? {
                next += 1;
                break;
            }
        }
    }
    Ok(Value::Map(rows.into()))
}

fn mark_zero_quadruple(
    model: &mut Model,
    target: TypeId,
    inputs: &[(TypeId, Rc<[usize]>)],
    rows: &mut [Value],
) -> Result<()> {
    // With zero adjacent cochains, all gap-two choices are cocycles and the
    // gap-three equations vanish. The only surviving top term is a_02 cup a_24.
    rows[0] = Value::Point(1);
    if rows.len() == 1 {
        return Ok(());
    }
    let Kind::Cohomology(base, _) = *model.kind(target) else {
        unreachable!()
    };
    let pair_degree = |a: usize, b: usize| {
        let Kind::Cohomology(_, p) = *model.kind(inputs[a].0) else {
            unreachable!()
        };
        let Kind::Cohomology(_, q) = *model.kind(inputs[b].0) else {
            unreachable!()
        };
        // The bracket's result degree has already been checked.
        (p - 1) + q
    };
    let left_degree = pair_degree(0, 1);
    let right_degree = pair_degree(2, 3);
    let left = model.cohomology_type(base, left_degree)?;
    let right = model.cohomology_type(base, right_degree)?;
    let left_rank = model.cohomology(left)?.rank();
    let right_rank = model.cohomology(right)?.rank();
    let outer = left_rank.min(right_rank);
    let inner = left_rank.max(right_rank);
    if outer == 0 {
        return Ok(());
    }

    let entries = outer.saturating_mul(inner);
    model.allocate(entries.saturating_add(rows.len() + outer + inner))?;
    let mut basis_products = Vec::with_capacity(entries);
    for i in 0..outer {
        for j in 0..inner {
            let (a, b) = if left_rank <= right_rank {
                (i, j)
            } else {
                (j, i)
            };
            let product = cohomology::cup(model, left, &[a], right, &[b], target)?;
            // Enumeration of H uses the coordinates as a bit mask, so XOR
            // of row indices is addition of classes on this target.
            basis_products.push(model.index(target, &product)?);
        }
    }
    if basis_products.iter().all(|&product| product == 0) {
        return Ok(());
    }

    // For each choice on the smaller side, take the ordinary linear image
    // on the other side, then union those images without taking their span.
    let mut coefficients = vec![false; outer];
    let mut images = vec![0; inner];
    let mut pivots = vec![0; rows.len().ilog2() as usize];
    let mut span = Vec::with_capacity(rows.len());
    let mut remaining = rows.len() - 1;
    loop {
        if !advance_bilinear(model, &mut coefficients, &mut images, &basis_products)? {
            return Ok(());
        }
        pivots.fill(0);
        model.spend(images.len())?;
        for &image in &images {
            let mut vector = image;
            while vector != 0 {
                model.spend(1)?;
                let lead = vector.ilog2() as usize;
                if pivots[lead] == 0 {
                    pivots[lead] = vector;
                    break;
                }
                vector ^= pivots[lead];
            }
        }

        span.clear();
        span.push(0);
        for &pivot in pivots.iter().filter(|&&pivot| pivot != 0) {
            let previous = span.len();
            model.spend(previous)?;
            for index in 0..previous {
                let class = span[index] ^ pivot;
                span.push(class);
                if matches!(rows[class], Value::Point(0)) {
                    rows[class] = Value::Point(1);
                    remaining -= 1;
                    if remaining == 0 {
                        return Ok(());
                    }
                }
            }
        }
    }
}

fn advance_bilinear(
    model: &mut Model,
    coefficients: &mut [bool],
    images: &mut [usize],
    products: &[usize],
) -> Result<bool> {
    for (flag, row) in coefficients
        .iter_mut()
        .zip(products.chunks_exact(images.len()))
    {
        model.spend(1 + images.len())?;
        *flag ^= true;
        for (image, basis_image) in images.iter_mut().zip(row) {
            *image ^= basis_image;
        }
        if *flag {
            return Ok(true);
        }
    }
    Ok(false)
}

#[derive(Clone, Copy)]
struct Interval {
    from: usize,
    to: usize,
    degree: usize,
}

#[derive(Default)]
struct Witness {
    cochain: Vec<bool>,
    choice: Vec<bool>,
}

type FacePair = (usize, usize);

struct Search {
    base: TypeId,
    offsets: Vec<usize>,
    intervals: Vec<Interval>,
    witnesses: Vec<Witness>,
    products: HashMap<(usize, usize), Vec<FacePair>>,
}

impl Search {
    fn new(model: &mut Model, base: TypeId, inputs: &[(TypeId, Rc<[usize]>)]) -> Result<Self> {
        let arity = inputs.len();
        let size = arity
            .checked_add(1)
            .and_then(|endpoints| arity.checked_mul(endpoints))
            .ok_or_else(|| Error::plain("bracket interval count exceeds the storage limit"))?
            / 2;
        model.allocate(size.saturating_mul(3).saturating_add(arity))?;
        model.spend(size)?;

        let mut prefix = Vec::with_capacity(arity + 1);
        let mut degree_sum = 0;
        prefix.push(degree_sum);
        let mut witnesses = Vec::with_capacity(size - 1);
        for (ty, class) in inputs {
            let Kind::Cohomology(_, degree) = *model.kind(*ty) else {
                unreachable!()
            };
            degree_sum += degree - 1;
            prefix.push(degree_sum);
            let group = model.cohomology(*ty)?;
            model.allocate(group.cochain_len())?;
            witnesses.push(Witness {
                cochain: group.representative(model, class)?,
                choice: Vec::new(),
            });
        }
        witnesses.resize_with(size - 1, Witness::default);

        let mut offsets = Vec::with_capacity(arity);
        let mut intervals = Vec::with_capacity(size);
        for gap in 1..=arity {
            offsets.push(intervals.len());
            for from in 0..=arity - gap {
                let to = from + gap;
                let degree = 1 + (prefix[to] - prefix[from]);
                intervals.push(Interval { from, to, degree });
            }
        }

        Ok(Self {
            base,
            offsets,
            intervals,
            witnesses,
            products: HashMap::new(),
        })
    }

    fn interval_index(&self, from: usize, to: usize) -> usize {
        self.offsets[to - from - 1] + from
    }

    fn group(&self, model: &mut Model, degree: usize) -> Result<Rc<Cohomology>> {
        let ty = model.cohomology_type(self.base, degree)?;
        model.cohomology(ty)
    }

    fn fourfold_compatible(
        &mut self,
        model: &mut Model,
        inputs: &[(TypeId, Rc<[usize]>)],
    ) -> Result<bool> {
        let pairs = self.offsets[1];
        let triples = self.offsets[2];
        for index in pairs..triples {
            if !self.choose_primitive(model, index)? {
                return Ok(false);
            }
        }
        let mut target = Vec::new();
        let mut types = Vec::new();
        let mut offsets = Vec::new();
        let mut dimension = 0;
        for index in triples..self.offsets[3] {
            let ty = model.cohomology_type(self.base, self.intervals[index].degree + 1)?;
            let group = model.cohomology(ty)?;
            let rhs = self.sum_products(model, index, &group)?;
            let coordinates = group.reduce_coordinates(model, rhs)?;
            target.extend(coordinates.into_iter().map(|i| dimension + i));
            types.push(ty);
            offsets.push(dimension);
            dimension += group.rank();
        }
        if target.is_empty() {
            return Ok(true);
        }

        model.allocate(dimension)?;
        let mut pivots = vec![None::<Box<[usize]>>; dimension];
        let mut cells = 0;
        // A change h_i to a gap-two primitive changes both adjacent triple
        // obstructions. Its column must be inserted in the joint image once.
        for (i, index) in (pairs..triples).enumerate() {
            let ty = model.cohomology_type(self.base, self.intervals[index].degree)?;
            let group = model.cohomology(ty)?;
            for generator in 0..group.rank() {
                model.spend(1)?;
                let mut column = Vec::new();
                if i > 0 {
                    let (left_ty, left) = &inputs[i - 1];
                    let Value::Class(product) =
                        cohomology::cup(model, *left_ty, left, ty, &[generator], types[i - 1])?
                    else {
                        unreachable!()
                    };
                    column.extend(
                        product
                            .iter()
                            .map(|&coordinate| offsets[i - 1] + coordinate),
                    );
                }
                if i + 2 < inputs.len() {
                    let (right_ty, right) = &inputs[i + 2];
                    let Value::Class(product) =
                        cohomology::cup(model, ty, &[generator], *right_ty, right, types[i])?
                    else {
                        unreachable!()
                    };
                    column.extend(product.iter().map(|&coordinate| offsets[i] + coordinate));
                }
                reduce_columns(model, &mut column, &pivots)?;
                if let Some(&leading) = column.last() {
                    cells += column.len();
                    if cells > MAX_PIVOT_CELLS {
                        return Err(Error::plain(
                            "fourfold compatibility exceeds the pivot storage limit",
                        ));
                    }
                    model.allocate(column.len())?;
                    pivots[leading] = Some(column.into_boxed_slice());
                    reduce_columns(model, &mut target, &pivots)?;
                    if target.is_empty() {
                        return Ok(true);
                    }
                }
            }
        }
        Ok(false)
    }

    fn choose_primitive(&mut self, model: &mut Model, index: usize) -> Result<bool> {
        let degree = self.intervals[index].degree;
        let result = self.group(model, degree + 1)?;
        let rhs = self.sum_products(model, index, &result)?;
        let group = self.group(model, degree)?;
        let Some(primitive) = group.primitive(model, rhs)? else {
            return Ok(false);
        };

        let witness = &mut self.witnesses[index];
        if witness.cochain.is_empty() {
            model.allocate(group.cochain_len() + group.rank())?;
            witness.cochain.resize(group.cochain_len(), false);
            witness.choice.resize(group.rank(), false);
        }

        model.spend(witness.cochain.len() + witness.choice.len() + primitive.len())?;
        witness.cochain.fill(false);
        witness.choice.fill(false);
        for coordinate in primitive {
            witness.cochain[coordinate] = true;
        }
        Ok(true)
    }

    fn advance_choice(&mut self, model: &mut Model, index: usize) -> Result<bool> {
        let group = self.group(model, self.intervals[index].degree)?;
        let witness = &mut self.witnesses[index];
        // Adding generators in binary-counter order explores one cocycle per
        // cohomology class without requiring 2^rank to fit a machine integer.
        for (index, choice) in witness.choice.iter_mut().enumerate() {
            let generator = group.generator(index);
            model.spend(1 + generator.len())?;
            *choice ^= true;
            for &coordinate in generator {
                witness.cochain[coordinate] ^= true;
            }
            if *choice {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn sum_products(
        &mut self,
        model: &mut Model,
        index: usize,
        result: &Cohomology,
    ) -> Result<Vec<usize>> {
        if result.cochain_len() == 0 {
            return Ok(Vec::new());
        }
        let Interval { from, to, .. } = self.intervals[index];
        model.spend(result.cochain_len())?;
        let mut support = vec![false; result.cochain_len()];
        for split in from + 1..to {
            model.spend(1)?;
            let first = self.interval_index(from, split);
            let second = self.interval_index(split, to);
            let degrees = (self.intervals[first].degree, self.intervals[second].degree);
            if !self.products.contains_key(&degrees) {
                let left = self.group(model, degrees.0)?;
                let right = self.group(model, degrees.1)?;
                model.allocate(1 + 2 * result.cochain_len())?;
                let mut product = Vec::with_capacity(result.cochain_len());
                result.for_cup_faces(model, &left, &right, degrees.0, |_, front, back| {
                    product.push((front, back));
                })?;
                self.products.insert(degrees, product);
            }
            let faces = &self.products[&degrees];
            let left = &self.witnesses[first].cochain;
            let right = &self.witnesses[second].cochain;
            model.spend(faces.len())?;
            for (bit, &(a, b)) in support.iter_mut().zip(faces) {
                *bit ^= left[a] && right[b];
            }
        }
        Ok(support
            .into_iter()
            .enumerate()
            .filter_map(|(index, bit)| bit.then_some(index))
            .collect())
    }
}

fn reduce_columns(
    model: &mut Model,
    column: &mut Vec<usize>,
    pivots: &[Option<Box<[usize]>>],
) -> Result<()> {
    while let Some(pivot) = column.last().and_then(|&leading| pivots[leading].as_ref()) {
        model.spend(column.len() + pivot.len())?;
        *column = xor_sorted(column, pivot);
    }
    Ok(())
}
