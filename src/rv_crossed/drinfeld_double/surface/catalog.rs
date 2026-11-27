use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_cyclotomic;

use super::{Double, Marking};

pub(in crate::rv_crossed::drinfeld_double) struct Catalog {
    pub(super) conductor: usize,
    pub(super) order: usize,
    pub(super) characters: Vec<Vec<usize>>,
    product: Vec<Vec<usize>>,
    inverse: Vec<usize>,
    handles: Vec<Vec<usize>>,
    counts: Vec<u128>,
}

impl Catalog {
    pub(super) fn new(model: &mut Model, double: &Double) -> Result<Self> {
        let conductor = double.conductor;
        let order = double.group.order();
        model.spend(3 * order.pow(3))?;
        let characters = (0..order)
            .map(|i| {
                let s = &double.simples[i];
                if s.flux != 0 || s.dimension != 1 {
                    return Err(Error::plain(
                        "surface electric simples are not the full character group",
                    ));
                }
                s.phases
                    .iter()
                    .map(|p| {
                        p.ok_or_else(|| Error::plain("surface electric character is undefined"))
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let index = characters
            .iter()
            .enumerate()
            .map(|(i, v)| (v.clone(), i))
            .collect::<BTreeMap<_, _>>();
        if index.len() != order || characters[0].iter().any(|&p| p != 0) {
            return Err(Error::plain(
                "surface electric catalogue lacks distinct characters or its identity",
            ));
        }
        let locate = |phase: Vec<usize>| {
            index.get(&phase).copied().ok_or_else(|| {
                Error::plain("surface residual gauge phase is not an electric character")
            })
        };
        let inverse = characters
            .iter()
            .map(|x| locate(x.iter().map(|&a| (conductor - a) % conductor).collect()))
            .collect::<Result<Vec<_>>>()?;
        let product = characters
            .iter()
            .map(|x| {
                characters
                    .iter()
                    .map(|y| locate(x.iter().zip(y).map(|(a, b)| (a + b) % conductor).collect()))
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let mut handles = vec![vec![0; order]; order];
        let mut counts = vec![0; order];
        for (a, row) in handles.iter_mut().enumerate() {
            for (b, charge) in row.iter_mut().enumerate() {
                *charge = locate(
                    (0..order)
                        .map(|t| {
                            (conductor + double.tensor_phase(t, a, b)
                                - double.tensor_phase(t, b, a))
                                % conductor
                        })
                        .collect(),
                )?;
                counts[*charge] += 1;
            }
        }
        model.allocate(3 * order.pow(2) + 10 * order + 12)?;
        Ok(Self {
            conductor,
            order,
            characters,
            product,
            inverse,
            handles,
            counts,
        })
    }

    pub(super) fn check_charge(&self, q: usize) -> Result<()> {
        if q >= self.order {
            return Err(Error::plain("surface electric charge is out of range"));
        }
        Ok(())
    }

    pub(super) fn inverse(&self, q: usize) -> usize {
        self.inverse[q]
    }

    pub(super) fn charge(&self, marks: &[usize]) -> Result<usize> {
        if !marks.len().is_multiple_of(2) || marks.iter().any(|&x| x >= self.order) {
            return Err(Error::plain(
                "surface marking requires valid pairs of handle holonomies",
            ));
        }
        Ok(marks
            .chunks_exact(2)
            .fold(0, |q, p| self.product[q][self.handles[p[0]][p[1]]]))
    }

    pub(super) fn dimensions(
        &self,
        model: &mut Model,
        double: &Double,
        mut g: usize,
    ) -> Result<Vec<u128>> {
        let genus = g;
        let mut count = vec![0; self.order];
        count[0] = 1;
        let mut matrix = self.counts.clone();
        while g != 0 {
            if g & 1 != 0 {
                count = self.convolve(model, &count, &matrix)?;
            }
            g >>= 1;
            if g != 0 {
                matrix = self.convolve(model, &matrix, &matrix)?;
            }
        }
        self.verify_verlinde(model, double, genus, &count)?;
        Ok(count)
    }

    fn convolve(&self, model: &mut Model, a: &[u128], b: &[u128]) -> Result<Vec<u128>> {
        model.spend(self.order.pow(2))?;
        let mut counts = vec![0u128; self.order];
        for (i, &x) in a.iter().enumerate().filter(|(_, x)| **x != 0) {
            for (j, &y) in b.iter().enumerate().filter(|(_, y)| **y != 0) {
                let q = self.product[i][j];
                counts[q] = counts[q]
                    .checked_add(x.checked_mul(y).ok_or_else(size)?)
                    .ok_or_else(size)?;
            }
        }
        Ok(counts)
    }

    fn verify_verlinde(
        &self,
        model: &mut Model,
        double: &Double,
        g: usize,
        counts: &[u128],
    ) -> Result<()> {
        model.spend(self.order * double.count())?;
        let power = if g == 0 {
            0
        } else {
            g.checked_sub(1)
                .and_then(|x| x.checked_mul(2))
                .and_then(|n| u32::try_from(n).ok())
                .ok_or_else(size)?
        };
        let factor = if g == 0 { self.order.pow(2) as u128 } else { 1 };
        let weights = double
            .simples
            .iter()
            .map(|s| {
                if g == 0 {
                    return Ok(s.dimension.pow(2) as i128);
                }
                if !self.order.is_multiple_of(s.dimension) {
                    return Err(Error::plain(
                        "surface Verlinde dimension has nonintegral quantum-dimension ratio",
                    ));
                }
                ((self.order / s.dimension) as i128)
                    .checked_pow(power)
                    .ok_or_else(size)
            })
            .collect::<Result<Vec<_>>>()?;
        for (q, &actual) in counts.iter().enumerate() {
            let mut terms = vec![0i128; self.conductor];
            for (a, &weight) in weights.iter().enumerate() {
                // Electric S_{qa} / S_{0a} is the unit phase on flux a.
                let cell = double.s[q][a];
                if cell.weight != double.simples[a].dimension as i128 {
                    return Err(Error::plain(
                        "surface electric character disagrees with its Verlinde coefficient",
                    ));
                }
                terms[cell.phase] = terms[cell.phase].checked_add(weight).ok_or_else(size)?;
            }
            let reduced = rv_cyclotomic::reduced(model, self.conductor, terms)?;
            let result = match reduced.as_slice() {
                [] => 0,
                &[n] if n >= 0 && (n as u128).is_multiple_of(factor) => (n as u128) / factor,
                _ => {
                    return Err(Error::plain(
                        "surface Verlinde dimension is not a nonnegative integer",
                    ))
                }
            };
            if result != actual {
                return Err(Error::plain(
                    "surface gauge-state count disagrees with the independent Verlinde dimension",
                ));
            }
        }
        Ok(())
    }
}

impl Double {
    pub(super) fn surface_catalog<'a>(&'a self, model: &mut Model) -> Result<&'a Catalog> {
        if self.surface.get().is_none() {
            let data = Catalog::new(model, self)?;
            self.surface.set(data).map_err(|_| {
                Error::plain("could not store certified surface electric characters")
            })?;
        }
        Ok(self.surface.get().expect("certified electric characters"))
    }

    pub(crate) fn surface_check_charge(&self, q: usize) -> Result<()> {
        if q >= self.group.order() {
            return Err(Error::plain("surface electric charge is out of range"));
        }
        Ok(())
    }

    pub(crate) fn surface_dimensions(&self, model: &mut Model, g: usize) -> Result<String> {
        let counts = self.surface_catalog(model)?.dimensions(model, self, g)?;
        Ok(format!(
            "rvDoubleSurfaceDimensions[{g}] = charges={counts:?}, verlinde_verified=true"
        ))
    }

    pub(crate) fn surface_charge(&self, model: &mut Model, marks: &Marking) -> Result<String> {
        let catalog = self.surface_catalog(model)?;
        let q = catalog.charge(marks)?;
        Ok(format!(
            "rvDoubleSurfaceCharge = charge={q}, conductor={}, phases={:?}",
            catalog.conductor, catalog.characters[q]
        ))
    }

    pub(crate) fn surface_sewing(&self, model: &mut Model, [g, h]: [usize; 2]) -> Result<String> {
        let catalog = self.surface_catalog(model)?;
        let a = catalog.dimensions(model, self, g)?;
        let b = catalog.dimensions(model, self, h)?;
        let dim = catalog.dimensions(model, self, g.checked_add(h).ok_or_else(size)?)?[0];
        let mut total = 0u128;
        let channels = a
            .iter()
            .enumerate()
            .map(|(q, &x)| {
                let r = catalog.inverse(q);
                let number = x.checked_mul(b[r]).ok_or_else(size)?;
                total = total.checked_add(number).ok_or_else(size)?;
                Ok((q, r, number))
            })
            .collect::<Result<Vec<_>>>()?;
        if total != dim {
            return Err(Error::plain(
                "surface sewing does not recover the closed dimension",
            ));
        }
        Ok(format!("rvDoubleSurfaceSewing[{g},{h}] = dimension={dim}, channels={channels:?}, verified=true"))
    }
}

fn size() -> Error {
    Error::plain("surface dimension exceeds the exact integer bound")
}
