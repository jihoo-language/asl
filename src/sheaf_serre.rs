use std::rc::Rc;

use crate::complex::{OrderComplex, Simplex, MAX_PIVOT_CELLS};
use crate::error::{Error, Result};
use crate::model::Model;
use crate::prime_linear::Coordinate;
use crate::sheaf::{self, Matrix, Sheaf};
use crate::sheaf_hyper::{Complex, Map, NamedSheaf};
use crate::sheaf_leray::component_name;
use crate::syntax::Label;

pub(crate) fn build(
    model: &mut Model,
    name: &Label,
    source: &Sheaf,
) -> Result<(usize, Vec<NamedSheaf>, Complex)> {
    let height = model.height(source.base(), true)?.saturating_sub(1);
    let names = (0..=height)
        .rev()
        .map(|degree| component_name(model, name, degree))
        .collect::<Result<Vec<_>>>()?;
    let mut order = if source.dims().iter().any(|&dim| dim != 0) {
        Some(OrderComplex::new(model, source.base())?)
    } else {
        None
    };
    let mut chains = order.as_ref().map_or_else(Vec::new, OrderComplex::vertices);
    chains.retain(|chain| source.dims()[chain[0]] != 0);
    let mut components = Vec::<NamedSheaf>::new();
    let mut maps = Vec::new();
    let mut previous: Option<BarLayer> = None;
    for (length, component_name) in names.into_iter().enumerate() {
        let layer = BarLayer::build(model, source, &mut chains)?;
        if let Some(previous) = &previous {
            sheaf::check_cells(layer.cells + previous.cells)?;
        }
        let dims = layer.stalks.iter().map(|stalk| stalk.size).collect();
        let sheaf =
            Sheaf::from_cover_maps(model, source.base(), source.prime, dims, |model, a, b| {
                Ok(Rc::new(layer.restriction(model, source, a, b)?))
            })?;
        let component = NamedSheaf {
            name: component_name,
            sheaf: Rc::new(sheaf),
        };
        if let Some((earlier_chain, target)) = previous.as_ref().zip(components.last()) {
            let matrix = layer.differential(model, source, earlier_chain)?;
            let d = Map::from_matrices(model, component.clone(), target.clone(), matrix)?;
            let label = Label {
                span: name.span,
                text: format!("{} bar boundary at {length}", name.text),
            };
            maps.push((label, Rc::new(d)));
        }
        components.push(component);
        previous = Some(layer);
        if length < height {
            if let Some(order) = &mut order {
                chains = order.extend(model, &chains)?;
            }
        }
    }
    // Chain length n is sheaf degree h-n.
    components.reverse();
    maps.reverse();
    let complex = Complex::build(model, components.clone(), maps)?;
    Ok((height, components, complex))
}

struct BarLayer {
    stalks: Vec<Stalk>,
    cells: usize,
}

#[derive(Default)]
struct Stalk {
    blocks: Vec<Block>,
    size: usize,
}

struct Block {
    chain: Rc<[usize]>,
    offset: usize,
}

impl BarLayer {
    fn build(model: &mut Model, source: &Sheaf, chains: &mut [Simplex]) -> Result<Self> {
        let mut cells = 0;
        let mut stalks: Vec<_> = (0..source.dims().len()).map(|_| Stalk::default()).collect();
        model.spend(
            chains.len()
                * (1 + chains.len().max(1).ilog2() as usize)
                * chains.first().map_or(1, Vec::len),
        )?;
        chains.sort_unstable();
        for chain in chains {
            let first = chain[0];
            let last = *chain.last().unwrap();
            let dimension = source.dims()[first];
            let chain: Rc<[_]> = chain.as_slice().into();
            model.spend(stalks.len() + chain.len())?;
            for (point, stalk) in stalks.iter_mut().enumerate() {
                if point == last || source.poset.strict_less(point, last) {
                    cells += 3 + chain.len();
                    sheaf::check_cells(cells)?;
                    let offset = stalk.size;
                    stalk.size = offset
                        .checked_add(dimension)
                        .filter(|&width| width <= MAX_PIVOT_CELLS)
                        .ok_or_else(|| Error::plain("Serre bar stalk exceeds the storage limit"))?;
                    stalk.blocks.push(Block {
                        chain: chain.clone(),
                        offset,
                    });
                }
            }
        }
        Ok(Self { stalks, cells })
    }

    fn restriction(&self, model: &mut Model, source: &Sheaf, a: usize, b: usize) -> Result<Matrix> {
        let bigger = &self.stalks[a];
        let smaller = &self.stalks[b];
        sheaf::check_cells(1 + 5 * smaller.size)?;
        model.spend(smaller.size)?;
        let mut columns = Vec::with_capacity(smaller.size);
        for block in &smaller.blocks {
            let coordinate = bigger.index(model, &block.chain, 0)?;
            let width = source.dims()[block.chain[0]];
            for i in 0..width {
                columns.push((coordinate + i, vec![(block.offset + i, 1)].into()));
            }
        }
        Matrix::from_columns(columns)
    }

    fn differential(
        &self,
        model: &mut Model,
        source: &Sheaf,
        lower: &Self,
    ) -> Result<Vec<Rc<Matrix>>> {
        let mut cells = 0;
        self.stalks
            .iter()
            .zip(&lower.stalks)
            .enumerate()
            .map(|(point, (stalk, target))| {
                let matrix = Matrix::from_action(model, stalk.size, |model, index| {
                    let (chain, coordinate) = stalk.generator(index);
                    boundary(model, source, point, chain, coordinate, target)
                })?;
                cells += matrix.cells();
                sheaf::check_cells(cells)?;
                Ok(Rc::new(matrix))
            })
            .collect()
    }
}

impl Stalk {
    fn generator(&self, index: usize) -> (&[usize], usize) {
        let i = self.blocks.partition_point(|block| block.offset <= index) - 1;
        let block = &self.blocks[i];
        (&block.chain, index - block.offset)
    }

    fn index(&self, model: &mut Model, chain: &[usize], coordinate: usize) -> Result<usize> {
        model.spend(1 + chain.len() * (1 + self.blocks.len().max(1).ilog2() as usize))?;
        let position = self
            .blocks
            .binary_search_by(|block| block.chain.as_ref().cmp(chain))
            .map_err(|_| Error::plain("Serre bar face is absent from its target stalk"))?;
        Ok(self.blocks[position].offset + coordinate)
    }
}

fn boundary(
    model: &mut Model,
    source: &Sheaf,
    point: usize,
    chain: &[usize],
    coordinate: usize,
    target: &Stalk,
) -> Result<Vec<Coordinate>> {
    let mut output = Vec::new();
    let last = chain.len() - 1;
    let prime = source.prime;
    let mut face = Vec::new();
    for i in 0..=last {
        model.spend(chain.len())?;
        if i == last {
            let before = chain[last - 1];
            if point != before && !source.poset.strict_less(point, before) {
                continue;
            }
        }
        face.clear();
        face.extend_from_slice(&chain[..i]);
        face.extend_from_slice(&chain[i + 1..]);
        let sign = if i % 2 == 0 { 1 } else { prime.negate(1) };
        if i == 0 {
            let rows = source.restriction(chain[0], chain[1]).column(coordinate);
            if !rows.is_empty() {
                let offset = target.index(model, &face, 0)?;
                model.spend(rows.len())?;
                output.extend(rows.iter().map(|&(j, value)| (offset + j, value)));
            }
        } else {
            let index = target.index(model, &face, coordinate)?;
            output.push((index, sign));
        }
        sheaf::check_cells(2 * output.len())?;
    }
    model.spend(output.len() * (1 + output.len().max(1).ilog2() as usize))?;
    output.sort_unstable_by_key(|&(index, _)| index);
    Ok(output)
}
