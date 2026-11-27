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
    totalize(model, name, &[source], &[])
}

pub(crate) fn build_complex(
    model: &mut Model,
    name: &Label,
    source: &Complex,
) -> Result<(usize, Vec<NamedSheaf>, Complex)> {
    let sheaves: Vec<_> = source
        .components()
        .iter()
        .map(|part| part.sheaf.as_ref())
        .collect();
    totalize(model, name, &sheaves, source.maps())
}

fn totalize(
    model: &mut Model,
    name: &Label,
    sheaves: &[&Sheaf],
    differential: &[Rc<Map>],
) -> Result<(usize, Vec<NamedSheaf>, Complex)> {
    let first = sheaves[0];
    let height = model.height(first.base(), true)?.saturating_sub(1);
    let size = height
        .checked_add(sheaves.len())
        .ok_or_else(|| Error::plain("Serre total degree is out of range"))?;
    model.spend(size)?;
    let names = (0..size)
        .rev()
        .map(|degree| component_name(model, name, degree))
        .collect::<Result<Vec<_>>>()?;
    let mut chains = Chains::new(model, sheaves)?;
    let mut components = Vec::<NamedSheaf>::new();
    let mut maps = Vec::new();
    let mut previous: Option<TotalLayer> = None;
    for (n, component_name) in (0..size).rev().zip(names) {
        let layer = TotalLayer::build(
            model,
            &mut chains,
            sheaves,
            height,
            n,
            previous.as_ref().map_or(0, |p| p.cells),
        )?;
        let sheaf = Sheaf::from_cover_maps(
            model,
            first.base(),
            first.prime,
            layer.widths.clone(),
            |model, a, b| Ok(Rc::new(layer.restriction(model, sheaves, a, b)?)),
        )?;
        let component = NamedSheaf {
            name: component_name,
            sheaf: Rc::new(sheaf),
        };
        if let Some((target, target_sheaf)) = previous.as_ref().zip(components.last()) {
            let matrices = layer.differential(model, sheaves, differential, target)?;
            let map = Map::from_matrices(model, component.clone(), target_sheaf.clone(), matrices)?;
            let label = Label {
                span: name.span,
                text: format!("{} Serre differential at {n}", name.text),
            };
            maps.push((label, Rc::new(map)));
        }
        components.push(component);
        previous = Some(layer);
    }
    components.reverse();
    maps.reverse();
    let complex = Complex::build(model, components.clone(), maps)?;
    Ok((height, components, complex))
}

struct Chains {
    order: Option<OrderComplex>,
    streams: Vec<Vec<Simplex>>,
    cells: usize,
}

impl Chains {
    fn new(model: &mut Model, source: &[&Sheaf]) -> Result<Self> {
        let order = if source
            .iter()
            .any(|sheaf| sheaf.dims().iter().any(|&dim| dim != 0))
        {
            Some(OrderComplex::new(model, source[0].base())?)
        } else {
            None
        };
        Ok(Self {
            order,
            streams: vec![Vec::new(); source.len()],
            cells: 3 * source.len(),
        })
    }

    fn take(
        &mut self,
        model: &mut Model,
        q: usize,
        length: usize,
        source: &Sheaf,
    ) -> Result<Vec<Simplex>> {
        let Some(order) = &mut self.order else {
            return Ok(Vec::new());
        };
        if length == 0 {
            model.spend(source.dims().len())?;
            let mut chains = order.vertices();
            chains.retain(|chain| source.dims()[chain[0]] != 0);
            return Ok(chains);
        }
        let previous = std::mem::take(&mut self.streams[q]);
        self.cells -= chain_cells(&previous);
        if previous.is_empty() {
            return Ok(previous);
        }
        order.extend(model, &previous)
    }

    fn retain(&mut self, q: usize, chains: Vec<Simplex>) {
        self.cells += chain_cells(&chains);
        self.streams[q] = chains;
    }
}

fn chain_cells(chains: &[Simplex]) -> usize {
    chains.len() * chains.first().map_or(0, |chain| 3 + chain.len())
}

struct TotalLayer {
    first: usize,
    sectors: Vec<Sector>,
    widths: Vec<usize>,
    cells: usize,
}

struct Sector {
    length: usize,
    start: Vec<usize>,
    layer: BarLayer,
}

impl TotalLayer {
    fn build(
        model: &mut Model,
        chains: &mut Chains,
        sheaves: &[&Sheaf],
        height: usize,
        degree: usize,
        retained: usize,
    ) -> Result<Self> {
        let first = degree.saturating_sub(height);
        let end = degree.min(sheaves.len() - 1);
        let n = sheaves[0].dims().len();
        let mut layer = Self {
            first,
            sectors: Vec::new(),
            widths: vec![0; n],
            cells: n,
        };
        for (q, &sheaf) in (first..=end).zip(&sheaves[first..=end]) {
            let length = height - (degree - q);
            let mut basis = chains.take(model, q, length, sheaf)?;
            let bar = BarLayer::build(model, sheaf, &mut basis)?;
            layer.cells += 8 + 4 * n + bar.cells;
            sheaf::check_cells(layer.cells + chains.cells + chain_cells(&basis) + retained)?;
            if length != height {
                chains.retain(q, basis);
            }
            let start = layer.widths.clone();
            model.spend(n)?;
            for (width, stalk) in layer.widths.iter_mut().zip(&bar.stalks) {
                *width += stalk.size;
                sheaf::check_cells(*width)?;
            }
            layer.sectors.push(Sector {
                length,
                start,
                layer: bar,
            });
        }
        Ok(layer)
    }

    fn sector(&self, q: usize) -> &Sector {
        &self.sectors[q - self.first]
    }

    fn restriction(
        &self,
        model: &mut Model,
        sheaves: &[&Sheaf],
        a: usize,
        b: usize,
    ) -> Result<Matrix> {
        sheaf::check_cells(1 + 5 * self.widths[b])?;
        model.spend(self.widths[b] + self.sectors.len())?;
        let mut columns = Vec::with_capacity(self.widths[b]);
        for (i, sector) in self.sectors.iter().enumerate() {
            let bigger = &sector.layer.stalks[a];
            let smaller = &sector.layer.stalks[b];
            let source = sheaves[self.first + i];
            for block in &smaller.blocks {
                let column = sector.start[a] + bigger.index(model, &block.chain, 0)?;
                let row = sector.start[b] + block.offset;
                let width = source.dims()[block.chain[0]];
                for k in 0..width {
                    columns.push((column + k, vec![(row + k, 1)].into()));
                }
            }
        }
        Matrix::from_columns(columns)
    }

    fn differential(
        &self,
        model: &mut Model,
        sheaves: &[&Sheaf],
        differential: &[Rc<Map>],
        target: &Self,
    ) -> Result<Vec<Rc<Matrix>>> {
        let mut cells = 0;
        self.widths
            .iter()
            .enumerate()
            .map(|(point, &width)| {
                let matrix = Matrix::from_action(model, width, |model, index| {
                    model.spend(1 + self.sectors.len().max(1).ilog2() as usize)?;
                    let i = self
                        .sectors
                        .partition_point(|sector| sector.start[point] <= index)
                        - 1;
                    let q = self.first + i;
                    let sector = &self.sectors[i];
                    let (chain, coordinate) =
                        sector.layer.stalks[point].generator(index - sector.start[point]);
                    let mut column = Vec::new();
                    if sector.length != 0 {
                        let bar = target.sector(q);
                        column = boundary(
                            model,
                            sheaves[q],
                            point,
                            chain,
                            coordinate,
                            &bar.layer.stalks[point],
                        )?;
                        for (row, _) in &mut column {
                            *row += bar.start[point];
                        }
                    }
                    if let Some(map) = differential.get(q) {
                        let inner = map.matrix(chain[0]).column(coordinate);
                        model.spend(1 + inner.len())?;
                        if !inner.is_empty() {
                            let next = target.sector(q + 1);
                            let offset = next.start[point]
                                + next.layer.stalks[point].index(model, chain, 0)?;
                            sheaf::check_cells(2 * (column.len() + inner.len()))?;
                            let prime = sheaves[q].prime;
                            let sign = if sector.length.is_multiple_of(2) {
                                1
                            } else {
                                prime.negate(1)
                            };
                            column.extend(
                                inner
                                    .iter()
                                    .map(|&(row, c)| (offset + row, prime.mul(sign, c))),
                            );
                        }
                    }
                    Ok(column)
                })?;
                cells += matrix.cells();
                sheaf::check_cells(cells)?;
                Ok(Rc::new(matrix))
            })
            .collect()
    }
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
