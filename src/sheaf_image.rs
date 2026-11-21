use std::collections::{btree_map::Entry, BTreeMap};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Poset, TypeId, Value};
use crate::sheaf::{self, Cohomology, Matrix, Open, Sheaf};

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Pull,
    Push(usize),
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Pull => "shpull",
            Self::Push(_) => "shpush",
        }
    }
}

pub(crate) struct Site {
    pub(crate) source: TypeId,
    pub(crate) target: TypeId,
    target_poset: Rc<Poset>,
    images: Vec<usize>,
}

impl Site {
    pub(crate) fn new(
        model: &mut Model,
        source: TypeId,
        target: TypeId,
        rows: &[Value],
    ) -> Result<Self> {
        let target_poset = model.poset(target).expect("named map target");
        model.spend(rows.len())?;
        let images = rows
            .iter()
            .map(|point| match *point {
                Value::Point(index) => index,
                _ => unreachable!("named map target values are points"),
            })
            .collect();
        Ok(Self {
            source,
            target,
            target_poset,
            images,
        })
    }

    pub(crate) fn pull(&self, model: &mut Model, sheaf: &Sheaf) -> Result<Sheaf> {
        if sheaf.base() != self.target {
            return Err(Error::plain(
                "shpull requires a sheaf on the map's named target",
            ));
        }
        model.spend(self.images.len())?;
        let dims = self
            .images
            .iter()
            .map(|&point| sheaf.dims()[point])
            .collect();
        let mut identities: BTreeMap<usize, Rc<Matrix>> = BTreeMap::new();
        Sheaf::from_cover_maps(model, self.source, sheaf.prime, dims, |model, x, z| {
            let (first, last) = (self.images[x], self.images[z]);
            if first != last {
                return Ok(sheaf.restriction_rc(first, last));
            }
            let rank = sheaf.dims()[first];
            Ok(match identities.entry(rank) {
                Entry::Occupied(entry) => entry.get().clone(),
                Entry::Vacant(entry) => entry
                    .insert(Rc::new(Matrix::identity(model, rank)?))
                    .clone(),
            })
        })
    }

    pub(crate) fn push(&self, model: &mut Model, sheaf: &Sheaf, degree: usize) -> Result<Sheaf> {
        if sheaf.base() != self.source {
            return Err(Error::plain(
                "shpush requires a sheaf on the map's named source",
            ));
        }
        if degree >= sheaf.height {
            let zero = Rc::new(Matrix::default());
            return Sheaf::from_cover_maps(
                model,
                self.target,
                sheaf.prime,
                vec![0; self.target_poset.names.len()],
                |_, _, _| Ok(zero.clone()),
            );
        }
        let opens = self.opens(model, sheaf)?;
        let mut stalks = Vec::new();
        for open in &opens.unique {
            let groups = Cohomology::default();
            let rank = if degree >= open.height {
                0
            } else {
                groups.rank(model, open, degree)?
            };
            stalks.push((groups, rank));
        }
        let dims = opens.indices.iter().map(|&open| stalks[open].1).collect();
        let mut restrictions: BTreeMap<(usize, usize), Rc<Matrix>> = BTreeMap::new();
        Sheaf::from_cover_maps(model, self.target, sheaf.prime, dims, |model, a, b| {
            let x = opens.indices[a];
            let z = opens.indices[b];
            if let Some(matrix) = restrictions.get(&(x, z)) {
                return Ok(matrix.clone());
            }
            let (from, rank) = &stalks[x];
            let (to, target_rank) = &stalks[z];
            let matrix = if x == z {
                Matrix::identity(model, *rank)?
            } else if *target_rank == 0 || *rank == 0 {
                Matrix::default()
            } else {
                let mut columns = Vec::new();
                let mut cells = 1;
                for index in 0..*rank {
                    let mut lift = from.representative_support(
                        model,
                        &opens.unique[x],
                        degree,
                        &[(index, 1)],
                    )?;
                    model.spend(lift.len())?;
                    lift.retain(|(generator, _)| opens.unique[z].contains(generator));
                    let coordinates =
                        to.class_coordinates(model, &opens.unique[z], degree, &lift)?;
                    if !coordinates.is_empty() {
                        cells += 3 + 2 * coordinates.len();
                        sheaf::check_cells(cells)?;
                        columns.push((index, coordinates.into()));
                    }
                }
                Matrix::from_columns(columns)?
            };
            let matrix = Rc::new(matrix);
            restrictions.insert((x, z), matrix.clone());
            Ok(matrix)
        })
    }

    pub(crate) fn opens<'a>(&self, model: &mut Model, sheaf: &'a Sheaf) -> Result<Opens<'a>> {
        let mut unique = Vec::new();
        let mut indices = Vec::with_capacity(self.target_poset.names.len());
        let mut masks = BTreeMap::new();
        for point in 0..self.target_poset.names.len() {
            model.spend(1 + self.images.len())?;
            let filled: Vec<_> = self
                .images
                .iter()
                .map(|&image| point == image || self.target_poset.strict_less(point, image))
                .collect();
            let index = match masks.entry(filled) {
                Entry::Occupied(entry) => *entry.get(),
                Entry::Vacant(entry) => {
                    let next = unique.len();
                    sheaf::check_cells((next + 1) * self.images.len())?;
                    unique.push(Open::new(model, sheaf, entry.key())?);
                    entry.insert(next);
                    next
                }
            };
            indices.push(index);
        }
        Ok(Opens { unique, indices })
    }
}

pub(crate) struct Opens<'a> {
    pub(crate) unique: Vec<Open<'a>>,
    pub(crate) indices: Vec<usize>,
}
