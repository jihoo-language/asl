use super::smith::{exterior_dimension, Spectrum};
use super::{character, Algebra, EVEN_RANK};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::nilpotent::Limit;
use crate::hyperelliptic_kummer_surface::rational_sparse::decimal;
use crate::model::Model;

struct GroupSpectrum {
    divided: Spectrum,
    cover: u128,
    factorial: u128,
}

impl GroupSpectrum {
    fn new(divided: Spectrum, m: i128) -> Self {
        let power = divided.power;
        let cover = if power == 0 { 1 } else { m.unsigned_abs() };
        // A nonzero exterior divided power has order at most 32.
        let factorial = (1..=power as u128).product();
        Self {
            divided,
            cover,
            factorial,
        }
    }

    fn torsion(&self) -> String {
        let omit_unit =
            self.divided.base.small() == Some(1) && self.cover == 1 && self.factorial == 1;
        self.divided.factor_text(omit_unit)
    }
}

impl Algebra {
    fn exterior_group_spectrum(
        &self,
        model: &mut Model,
        limit: &Limit,
        degree: usize,
        m: i128,
        order: usize,
    ) -> Result<GroupSpectrum> {
        if degree > EVEN_RANK {
            return Err(Error::plain(
                "Kuga–Satake cohomological degree must be between 0 and 64",
            ));
        }
        let spectrum = if m == 0 && order > 0 {
            Spectrum::zero(exterior_dimension(degree))
        } else {
            self.exterior_smith_spectrum(model, limit, degree, order)?
        };
        // In each singly occupied pair block the triangular integral change
        // of variables conjugates T^m-1 to mD. It applies equally to unequal d_i.
        Ok(GroupSpectrum::new(spectrum, m))
    }

    pub(crate) fn exterior_group_report(
        &self,
        model: &mut Model,
        label: &str,
        limit: &Limit,
        degree: usize,
        m: i128,
        order: usize,
    ) -> Result<String> {
        let group = self.exterior_group_spectrum(model, limit, degree, m, order)?;
        let dimension = group.divided.dimension;
        let rank = group.divided.rank();
        let free = group.divided.free();
        let base = decimal(group.divided.base);
        let power = group.divided.power;
        let cover = group.cover;
        let factorial = group.factorial;
        let entries = group.divided.factor_text(false);
        Ok(format!(
            "rvKSExteriorGroupSmith[{label},{degree},{m},{order}] = dimension={dimension}, rank={rank}, kernel={free}, cokernel_free={free}, common_base={base}, cover_base={cover}, common_power={power}, factorial={factorial}, factors=[{entries}]"
        ))
    }

    pub(crate) fn mapping_torus_report(
        &self,
        model: &mut Model,
        label: &str,
        limit: &Limit,
        m: i128,
        degree: usize,
    ) -> Result<String> {
        if degree > EVEN_RANK + 1 {
            return Err(Error::plain(
                "Kuga–Satake mapping torus degree must be between 0 and 65",
            ));
        }
        let coinvariants = if let Some(lower) = degree.checked_sub(1) {
            self.exterior_group_spectrum(model, limit, lower, m, 1)?
        } else {
            GroupSpectrum::new(Spectrum::zero(0), m)
        };
        let fiber_invariants = if m == 0 || degree == 0 || degree >= EVEN_RANK {
            exterior_dimension(degree)
        } else {
            let rank = self.action(model, limit)?.smith(model)?.len();
            let weights = character(model, rank, degree)?;
            // Every sl₂ string contributes either to relative weight zero or one.
            weights[rank] + weights.get(rank + 1).unwrap_or(&0)
        };
        let fiber_coinvariants_free = coinvariants.divided.free();
        let free = fiber_invariants + fiber_coinvariants_free;
        let base = decimal(coinvariants.divided.base);
        let cover = coinvariants.cover;
        let torsion = coinvariants.torsion();
        Ok(format!(
            "rvKSMappingTorus[{label},{m},{degree}] = free={free}, fiber_invariants={fiber_invariants}, fiber_coinvariants_free={fiber_coinvariants_free}, common_base={base}, cover_base={cover}, torsion=[{torsion}]"
        ))
    }
}
