use std::rc::Rc;

use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::nikulin::binary_rref;
use crate::model::Model;

use super::{Cover, Descent, Pencil, Section};

pub(crate) struct MordellImage {
    pencil: Rc<Pencil>,
    rank: usize,
    image_rows: Vec<u32>,
    pivots: Vec<usize>,
    sha_dimension: usize,
}

impl MordellImage {
    pub(crate) fn new(
        model: &mut Model,
        pencil: Rc<Pencil>,
        generators: &[Rc<Section>],
    ) -> Result<Self> {
        pencil.split_roots(model)?;
        let relative = pencil.relative(model)?;
        let rank = relative.arithmetic_mordell_weil_rank(model, &pencil)?;
        let selmer = pencil.selmer_space(model)?;
        let sha_dimension = selmer.sha_dimension(rank)?;
        let (rows, pivots) = Self::span_rows(model, &pencil, generators)?;
        if pivots.len() != rank + 2 {
            return Err(Error::plain(format!(
                "K3 pencil Mordell-Weil mod-2 image has dimension {}, expected {} from arithmetic rank {rank} and split two-torsion",
                pivots.len(), rank + 2,
            )));
        }
        model.allocate(7 + 2 * rows.len())?;
        Ok(Self {
            pencil,
            rank,
            image_rows: rows,
            pivots,
            sha_dimension,
        })
    }

    pub(crate) fn dimension_of(
        model: &mut Model,
        pencil: &Pencil,
        generators: &[Rc<Section>],
    ) -> Result<usize> {
        Ok(Self::span_rows(model, pencil, generators)?.0.len())
    }

    fn span_rows(
        model: &mut Model,
        pencil: &Pencil,
        generators: &[Rc<Section>],
    ) -> Result<(Vec<u32>, Vec<usize>)> {
        let selmer = pencil.selmer_space(model)?;
        let mut rows = Vec::new();
        let mut add_point = |model: &mut Model, point: &Section| -> Result<()> {
            let class = point.descent(model, pencil)?;
            let vector = selmer.encode_local(model, &class)?.ok_or_else(|| {
                Error::plain("K3 pencil global section failed a local Kummer test")
            })?;
            rows.push(vector);
            Ok(())
        };
        for i in 0..3 {
            let point = pencil.torsion_at(model, i)?;
            add_point(model, &point)?;
        }
        for section in generators {
            add_point(model, section)?;
        }
        model.spend(rows.len().saturating_mul(selmer.width()))?;
        let pivots = binary_rref(&mut rows, selmer.width());
        rows.truncate(pivots.len());
        Ok((rows, pivots))
    }

    pub(crate) fn report(&self, name: &str) -> String {
        format!(
            "rvK3PencilMW2[{name}] = Mordell-Weil rank={}, image dimension={}, dim Sha[2]={}, echelon masks={:?}",
            self.rank, self.pivots.len(), self.sha_dimension, self.image_rows,
        )
    }

    fn coset(&self, model: &mut Model, cover: &Cover) -> Result<Option<u32>> {
        if !Rc::ptr_eq(&self.pencil, &cover.pencil) {
            return Err(Error::plain(
                "K3 pencil covering and Mordell-Weil certificate require the same named pencil",
            ));
        }
        let selmer = self.pencil.selmer_space(model)?;
        let Some(mut vector) = selmer.encode_local(model, &cover.b)? else {
            return Ok(None);
        };
        model.spend(self.pivots.len())?;
        for (&row, &pivot) in self.image_rows.iter().zip(&self.pivots) {
            if vector & (1 << pivot) != 0 {
                vector ^= row;
            }
        }
        Ok(Some(vector))
    }

    pub(crate) fn sha_class(&self, model: &mut Model, cover: &Cover) -> Result<Descent> {
        let vector = self.coset(model, cover)?.ok_or_else(|| {
            Error::plain("K3 pencil covering has a local obstruction and defines no Sha class")
        })?;
        self.pencil.selmer_space(model)?.decode(model, vector)
    }
}

impl Cover {
    pub(crate) fn status(
        &self,
        model: &mut Model,
        name: &str,
        image: Option<&MordellImage>,
    ) -> Result<String> {
        let result = if let Some(image) = image {
            match image.coset(model, self)? {
                None => "local obstruction, rational point=no".into(),
                Some(0) => "everywhere locally soluble, in certified global image, rational point=yes, Sha coset mask=0".into(),
                Some(vector) => format!("everywhere locally soluble, outside certified global image, rational point=no, Hasse-principle failure, Sha coset mask={vector}"),
            }
        } else if self.pencil.selmer_space(model)?.contains(model, &self.b)? {
            "everywhere locally soluble, global image uncertified".into()
        } else {
            "local obstruction, rational point=no".into()
        };
        Ok(format!("rvK3PencilCoverStatus[{name}] = {result}"))
    }
}
