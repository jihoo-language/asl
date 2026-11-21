use std::cell::OnceCell;
use std::collections::HashMap;
use std::rc::Rc;

use crate::complex::{for_each_boundary, xor_sorted, OrderComplex, Simplex, MAX_PIVOT_CELLS};
use crate::error::{Error, Result, Span};
use crate::model::{Kind, Model, TypeId, Value};
use crate::steenrod::Word;
use crate::syntax::Label;

#[derive(Default)]
pub(crate) struct Cohomology {
    simplices: Vec<Simplex>,
    indices: HashMap<Box<[usize]>, usize>,
    next: Vec<Simplex>,
    cofaces: Vec<Vec<usize>>,
    coboundaries: Vec<Option<Box<[usize]>>>,
    generator_pivots: Vec<Option<usize>>,
    representatives: Vec<Box<[usize]>>,
    primitive_images: OnceCell<ImageBasis>,
}

struct ImagePivot {
    image: Box<[usize]>,
    source: Box<[usize]>,
}

struct ImageBasis {
    pivots: Vec<Option<ImagePivot>>,
    cells: usize,
}

impl ImageBasis {
    fn new(rows: usize) -> Self {
        Self {
            pivots: (0..rows).map(|_| None).collect(),
            cells: 0,
        }
    }

    fn insert_column(
        &mut self,
        model: &mut Model,
        coordinate: usize,
        cofaces: &[usize],
    ) -> Result<Option<Vec<usize>>> {
        model.spend(1 + cofaces.len())?;
        let mut image = cofaces.to_vec();
        let mut source = vec![coordinate];
        while let Some(&leading) = image.last() {
            if let Some(pivot) = &self.pivots[leading] {
                image = symmetric_difference(model, &image, &pivot.image)?;
                source = symmetric_difference(model, &source, &pivot.source)?;
            } else {
                store_cells(&mut self.cells, image.len() + source.len())?;
                self.pivots[leading] = Some(ImagePivot {
                    image: image.into_boxed_slice(),
                    source: source.into_boxed_slice(),
                });
                return Ok(None);
            }
        }
        Ok(Some(source))
    }
}

impl Cohomology {
    pub(crate) fn build(model: &mut Model, base: TypeId, degree: usize) -> Result<Self> {
        let height = model.height(base, false)?;
        if height != usize::MAX && degree >= height {
            return Ok(Self::default());
        }
        if degree > 0 {
            let height = model.height(base, true)?;
            if height != usize::MAX && degree >= height {
                return Ok(Self::default());
            }
        }

        let mut complex = OrderComplex::new(model, base)?;
        let mut previous = Vec::new();
        let mut simplices = complex.vertices();
        for _ in 0..degree {
            let next = complex.extend(model, &simplices)?;
            previous = simplices;
            simplices = next;
        }
        let next = complex.extend(model, &simplices)?;
        let mut coboundaries = vec![None::<Box<[usize]>>; simplices.len()];
        let mut basis_cells = 0;
        // A delta column is the list of simplices containing its input face.
        let boundaries = if degree == 0 {
            Vec::new()
        } else {
            transpose(model, &previous, &simplices)?
        };
        for mut boundary in boundaries {
            while let Some(&leading) = boundary.last() {
                if let Some(pivot) = &coboundaries[leading] {
                    boundary = symmetric_difference(model, &boundary, pivot)?;
                } else {
                    store_cells(&mut basis_cells, boundary.len())?;
                    coboundaries[leading] = Some(boundary.into_boxed_slice());
                    break;
                }
            }
        }

        let cofaces = transpose(model, &simplices, &next)?;
        let mut images = ImageBasis::new(next.len());
        let mut representatives = Vec::new();
        let mut generator_pivots = vec![None; simplices.len()];
        for (coordinate, cofaces) in cofaces.iter().enumerate() {
            if coboundaries[coordinate].is_some() {
                continue;
            }

            // Nonpivot coordinates complement the coboundaries. Dependencies
            // among their delta columns give cocycles modulo coboundaries.
            if let Some(source) = images.insert_column(model, coordinate, cofaces)? {
                store_cells(&mut basis_cells, source.len())?;
                generator_pivots[coordinate] = Some(representatives.len());
                representatives.push(source.into_boxed_slice());
            }
        }

        let vertices: usize = simplices.iter().map(Vec::len).sum();
        let next_vertices: usize = next.iter().map(Vec::len).sum();
        model.allocate(
            2 * vertices + 2 * next_vertices + 4 * simplices.len() + next.len() + basis_cells,
        )?;
        let indices = simplices
            .iter()
            .enumerate()
            .map(|(index, simplex)| (simplex.clone().into_boxed_slice(), index))
            .collect();

        Ok(Self {
            simplices,
            indices,
            next,
            cofaces,
            coboundaries,
            generator_pivots,
            representatives,
            primitive_images: OnceCell::new(),
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.representatives.len()
    }

    pub(crate) fn cochain_len(&self) -> usize {
        self.simplices.len()
    }

    pub(crate) fn generator(&self, index: usize) -> &[usize] {
        &self.representatives[index]
    }

    pub(crate) fn primitive(
        &self,
        model: &mut Model,
        mut rhs: Vec<usize>,
    ) -> Result<Option<Vec<usize>>> {
        if rhs.is_empty() {
            return Ok(Some(Vec::new()));
        }
        let basis = if let Some(basis) = self.primitive_images.get() {
            basis
        } else {
            // The kernel basis from the original reduction records each dependent
            // column. Only independent columns are needed to recover primitives.
            let mut basis = ImageBasis::new(self.next.len());
            model.spend(self.cochain_len())?;
            for (coordinate, cofaces) in self.cofaces.iter().enumerate() {
                if self.coboundaries[coordinate].is_none()
                    && self.generator_pivots[coordinate].is_none()
                {
                    basis.insert_column(model, coordinate, cofaces)?;
                }
            }
            model.allocate(basis.cells + basis.pivots.len())?;
            self.primitive_images.get_or_init(|| basis)
        };
        let mut primitive = Vec::new();
        // Consecutive degrees share the same simplex order in the order complex.
        while let Some(&leading) = rhs.last() {
            let Some(pivot) = &basis.pivots[leading] else {
                return Ok(None);
            };
            rhs = symmetric_difference(model, &rhs, &pivot.image)?;
            primitive = symmetric_difference(model, &primitive, &pivot.source)?;
        }
        Ok(Some(primitive))
    }

    pub(crate) fn representative(&self, model: &mut Model, class: &[usize]) -> Result<Vec<bool>> {
        model.spend(self.simplices.len())?;
        let mut representative = vec![false; self.simplices.len()];
        for &index in class {
            let lift = &self.representatives[index];
            model.spend(lift.len())?;
            for &coordinate in lift.iter() {
                representative[coordinate] ^= true;
            }
        }
        Ok(representative)
    }

    fn coordinates(&self, model: &mut Model, support: Vec<usize>) -> Result<Value> {
        let coordinates = self.reduce_coordinates(model, support)?;
        model.allocate(coordinates.len())?;
        Ok(Value::Class(coordinates.into()))
    }

    pub(crate) fn reduce_coordinates(
        &self,
        model: &mut Model,
        mut support: Vec<usize>,
    ) -> Result<Vec<usize>> {
        // Reduce every coboundary pivot, including those below free coordinates.
        let mut reduced = Vec::new();
        while let Some(&leading) = support.last() {
            if let Some(pivot) = &self.coboundaries[leading] {
                support = symmetric_difference(model, &support, pivot)?;
            } else {
                model.spend(1)?;
                reduced.push(support.pop().unwrap());
            }
        }
        reduced.reverse();

        let mut coordinates = Vec::new();
        while let Some(&leading) = reduced.last() {
            let generator = self.generator_pivots[leading]
                .ok_or_else(|| Error::plain("cochain does not represent a closed class"))?;
            coordinates.push(generator);
            reduced = symmetric_difference(model, &reduced, &self.representatives[generator])?;
        }
        coordinates.reverse();
        Ok(coordinates)
    }

    pub(crate) fn for_cup_faces(
        &self,
        model: &mut Model,
        left: &Self,
        right: &Self,
        left_degree: usize,
        mut use_faces: impl FnMut(usize, usize, usize),
    ) -> Result<()> {
        for (coordinate, simplex) in self.simplices.iter().enumerate() {
            model.spend(simplex.len() + 1)?;
            let front = left.indices[&simplex[..=left_degree]];
            let back = right.indices[&simplex[left_degree..]];
            use_faces(coordinate, front, back);
        }
        Ok(())
    }
}

fn transpose(model: &mut Model, rows: &[Simplex], columns: &[Simplex]) -> Result<Vec<Vec<usize>>> {
    let mut transpose = vec![Vec::new(); rows.len()];
    for_each_boundary(model, rows, columns, |_, index, boundary| {
        for face in boundary {
            transpose[face].push(index);
        }
        Ok(())
    })?;
    Ok(transpose)
}

fn store_cells(cells: &mut usize, added: usize) -> Result<()> {
    *cells += added;
    if *cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "cohomology reduction exceeds the storage limit",
        ));
    }
    Ok(())
}

fn symmetric_difference(model: &mut Model, left: &[usize], right: &[usize]) -> Result<Vec<usize>> {
    model.spend(left.len() + right.len())?;
    Ok(xor_sorted(left, right))
}

pub(crate) fn zero() -> Value {
    Value::Class([].into())
}

pub(crate) fn add(model: &mut Model, left: &[usize], right: &[usize]) -> Result<Value> {
    let support = symmetric_difference(model, left, right)?;
    model.allocate(support.len())?;
    Ok(Value::Class(support.into()))
}

pub(crate) fn literal(
    model: &mut Model,
    ty: TypeId,
    chains: Vec<Vec<Label>>,
    span: Span,
) -> Result<Value> {
    let Kind::Cohomology(base, degree) = *model.kind(ty) else {
        unreachable!()
    };
    let poset = model.poset(base).ok_or_else(|| {
        Error::new(
            span,
            format!(
                "cocycle requires a named or Star space, found {}",
                model.type_name(base),
            ),
        )
    })?;
    if chains.is_empty() {
        return Ok(zero());
    }
    let mut support = Vec::with_capacity(chains.len());
    for chain in chains {
        if degree.checked_add(1) != Some(chain.len()) {
            return Err(Error::new(
                chain[0].span,
                format!(
                    "cocycle of degree {degree} cannot contain a simplex with {} points",
                    chain.len()
                ),
            ));
        }
        let mut points = Vec::with_capacity(chain.len());
        for label in chain {
            let point = model.point(base, &label)?;
            if let Some(&previous) = points.last() {
                if previous == point
                    || !model.leq(base, &Value::Point(previous), &Value::Point(point))?
                {
                    return Err(Error::new(
                        label.span,
                        format!(
                            "cocycle simplex is not a strict chain: {} < {}",
                            poset.names[previous], label.text,
                        ),
                    ));
                }
            }
            points.push(point);
        }
        support.push(points);
    }
    let data = model.cohomology(ty)?;
    model.spend(data.simplices.len() + data.next.len())?;
    let mut cochain = vec![false; data.simplices.len()];
    for simplex in support {
        model.spend(simplex.len())?;
        cochain[data.indices[simplex.as_slice()]] ^= true;
    }

    let mut coboundary = vec![false; data.next.len()];
    let mut support = Vec::new();
    for (index, coefficient) in cochain.into_iter().enumerate() {
        if coefficient {
            support.push(index);
            model.spend(data.cofaces[index].len())?;
            for &coface in &data.cofaces[index] {
                coboundary[coface] ^= true;
            }
        }
    }
    if let Some(offending) = coboundary.into_iter().position(|coefficient| coefficient) {
        let simplex: Vec<_> = data.next[offending]
            .iter()
            .map(|&point| poset.names[point].as_str())
            .collect();
        return Err(Error::new(
            span,
            format!(
                "cochain is not a cocycle: its coboundary is nonzero on {}",
                simplex.join(" < ")
            ),
        ));
    }
    data.coordinates(model, support)
}

pub(crate) fn cup(
    model: &mut Model,
    left_ty: TypeId,
    left: &[usize],
    right_ty: TypeId,
    right: &[usize],
    result_ty: TypeId,
) -> Result<Value> {
    if left.is_empty() || right.is_empty() {
        return Ok(zero());
    }
    let result = model.cohomology(result_ty)?;
    if result.rank() == 0 {
        return Ok(zero());
    }
    let Kind::Cohomology(_, p) = *model.kind(left_ty) else {
        unreachable!()
    };
    let left_group = model.cohomology(left_ty)?;
    let right_group = model.cohomology(right_ty)?;
    let front = left_group.representative(model, left)?;
    let back = right_group.representative(model, right)?;
    let mut support = Vec::new();
    result.for_cup_faces(model, &left_group, &right_group, p, |index, left, right| {
        if front[left] && back[right] {
            support.push(index);
        }
    })?;
    result.coordinates(model, support)
}

pub(crate) fn pull(
    model: &mut Model,
    function_ty: TypeId,
    table: &[Value],
    class_ty: TypeId,
    class: &[usize],
    result_ty: TypeId,
) -> Result<Value> {
    if class.is_empty() {
        return Ok(zero());
    }
    let result = model.cohomology(result_ty)?;
    if result.rank() == 0 {
        return Ok(zero());
    }
    let Kind::Arrow(_, range) = *model.kind(function_ty) else {
        unreachable!()
    };
    let group = model.cohomology(class_ty)?;
    let representative = group.representative(model, class)?;
    model.spend(table.len())?;
    let image_points = table
        .iter()
        .map(|image| model.index(range, image))
        .collect::<Result<Vec<_>>>()?;
    let mut support = Vec::new();
    let mut image = Vec::new();
    for (index, simplex) in result.simplices.iter().enumerate() {
        model.spend(simplex.len())?;
        image.clear();
        image.extend(simplex.iter().map(|&point| image_points[point]));
        // The image is weakly increasing, so equal vertices occur adjacently.
        if image.windows(2).all(|pair| pair[0] != pair[1])
            && representative[group.indices[image.as_slice()]]
        {
            support.push(index);
        }
    }
    result.coordinates(model, support)
}

pub(crate) fn act(
    model: &mut Model,
    base: TypeId,
    degree: usize,
    words: &[Word],
    class: Rc<[usize]>,
) -> Result<Value> {
    if class.is_empty() || words.is_empty() {
        return Ok(zero());
    }
    let mut support = Vec::new();
    for word in words {
        model.spend(1)?;
        let mut summand = class.clone();
        let mut current_degree = degree;
        for &index in word.iter().rev() {
            let Value::Class(next) = apply_square(model, base, current_degree, index, &summand)?
            else {
                unreachable!()
            };
            summand = next;
            if summand.is_empty() {
                break;
            }
            // Each word is homogeneous and its total action degree was checked.
            current_degree += index;
        }
        support = symmetric_difference(model, &support, &summand)?;
    }
    model.allocate(support.len())?;
    Ok(Value::Class(support.into()))
}

fn apply_square(
    model: &mut Model,
    base: TypeId,
    degree: usize,
    index: usize,
    class: &[usize],
) -> Result<Value> {
    // The action represents Sq^0 by the empty word.
    debug_assert!(index > 0);
    if index > degree {
        return Ok(zero());
    }
    let source_ty = model.cohomology_type(base, degree)?;
    let result_ty = model.cohomology_type(base, degree + index)?;
    if index == degree {
        return cup(model, source_ty, class, source_ty, class, result_ty);
    }

    let result = model.cohomology(result_ty)?;
    if result.rank() == 0 {
        return Ok(zero());
    }
    let source = model.cohomology(source_ty)?;
    let representative = source.representative(model, class)?;
    let mut support = Vec::new();
    let mut omitted = Vec::new();
    let mut front = Vec::new();
    let mut back = Vec::new();
    for (coordinate, simplex) in result.simplices.iter().enumerate() {
        omitted.clear();
        omitted.extend(0..2 * index);
        let mut value = false;
        loop {
            model.spend(simplex.len() + omitted.len())?;
            // The normalized cup-(degree-index) formula splits removed
            // vertices by the parity of each position plus its one-based rank.
            let even_count = omitted
                .iter()
                .enumerate()
                .filter(|&(rank, position)| (position + rank + 1) % 2 == 0)
                .count();
            if even_count == index {
                front.clear();
                back.clear();
                let mut rank = 0;
                for (position, &vertex) in simplex.iter().enumerate() {
                    if omitted.get(rank) == Some(&position) {
                        rank += 1;
                        if (position + rank) % 2 == 0 {
                            back.push(vertex);
                        } else {
                            front.push(vertex);
                        }
                    } else {
                        front.push(vertex);
                        back.push(vertex);
                    }
                }
                if representative[source.indices[front.as_slice()]]
                    && representative[source.indices[back.as_slice()]]
                {
                    value ^= true;
                }
            }
            if !next_subset(&mut omitted, simplex.len()) {
                break;
            }
        }
        if value {
            support.push(coordinate);
        }
    }
    result.coordinates(model, support)
}

fn next_subset(subset: &mut [usize], size: usize) -> bool {
    let offset = size - subset.len();
    for index in (0..subset.len()).rev() {
        if subset[index] < offset + index {
            subset[index] += 1;
            for following in index + 1..subset.len() {
                subset[following] = subset[following - 1] + 1;
            }
            return true;
        }
    }
    false
}

pub(crate) fn suspend(
    model: &mut Model,
    source_ty: TypeId,
    class: &[usize],
    result_ty: TypeId,
) -> Result<Value> {
    if class.is_empty() {
        return Ok(zero());
    }
    let result = model.cohomology(result_ty)?;
    if result.rank() == 0 {
        return Ok(zero());
    }
    let Kind::Cohomology(base, _) = *model.kind(source_ty) else {
        unreachable!()
    };
    // Equatorial vertices have their original indices and precede the poles.
    let north = model.elements(base)?.len();
    let source = model.cohomology(source_ty)?;
    let representative = source.representative(model, class)?;
    let mut support = Vec::new();
    for (coordinate, simplex) in result.simplices.iter().enumerate() {
        model.spend(simplex.len())?;
        let (&last, equator) = simplex.split_last().unwrap();
        if last == north && representative[source.indices[equator]] {
            support.push(coordinate);
        }
    }
    result.coordinates(model, support)
}
