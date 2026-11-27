use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::{Basis, Coordinate, QuotientBasis};
use crate::sheaf::{Cochains, Matrix, Sheaf};
use crate::sheaf_hyper::Complex;

pub(crate) enum StalkQuotient {
    Full(usize),
    Computed(QuotientBasis),
}

impl StalkQuotient {
    pub(crate) fn build(
        model: &mut Model,
        prime: Prime,
        rank: usize,
        incoming: Option<&Matrix>,
        outgoing: Option<&Matrix>,
    ) -> Result<Self> {
        if incoming.is_none_or(Matrix::is_zero) && outgoing.is_none_or(Matrix::is_zero) {
            return Ok(Self::Full(rank));
        }
        if rank > MAX_PIVOT_CELLS {
            return Err(Error::plain(
                "sheaf stalk quotient exceeds the storage limit",
            ));
        }
        let mut boundaries = Basis::default();
        if let Some(matrix) = incoming {
            for (_, column) in matrix.columns() {
                model.spend(1 + column.len())?;
                boundaries.insert(model, prime, column.to_vec())?;
            }
        }
        let basis = QuotientBasis::build(model, prime, boundaries, rank, |model, i| {
            let column = outgoing.map_or(&[][..], |m| m.column(i));
            model.spend(1 + column.len())?;
            Ok(column.to_vec())
        })?;
        Ok(Self::Computed(basis))
    }

    pub(crate) fn rank(&self) -> usize {
        match self {
            Self::Full(rank) => *rank,
            Self::Computed(basis) => basis.rank(),
        }
    }

    fn representative(
        &self,
        model: &mut Model,
        prime: Prime,
        index: usize,
    ) -> Result<Vec<Coordinate>> {
        match self {
            Self::Full(_) => Ok(vec![(index, 1)]),
            Self::Computed(basis) => basis.representative(model, prime, &[(index, 1)]),
        }
    }

    pub(crate) fn coordinates(
        &self,
        model: &mut Model,
        prime: Prime,
        raw: Vec<Coordinate>,
    ) -> Result<Vec<Coordinate>> {
        match self {
            Self::Full(_) => Ok(raw),
            Self::Computed(basis) => basis
                .project(model, prime, raw)?
                .ok_or_else(|| Error::plain("sheaf map does not preserve the stalk cycle space")),
        }
    }
}

pub(crate) fn induced_sheaf(
    model: &mut Model,
    original: &Sheaf,
    quotients: &[StalkQuotient],
) -> Result<Sheaf> {
    let prime = original.prime;
    let dims = quotients.iter().map(StalkQuotient::rank).collect();
    Sheaf::from_cover_maps(model, original.base(), prime, dims, |model, a, b| {
        let map = original.restriction(a, b);
        let induced = if map.is_zero() || quotients[b].rank() == 0 {
            Matrix::default()
        } else {
            Matrix::from_action(model, quotients[a].rank(), |model, index| {
                let lift = quotients[a].representative(model, prime, index)?;
                let image = map.apply(model, prime, &lift)?;
                quotients[b].coordinates(model, prime, image)
            })?
        };
        Ok(Rc::new(induced))
    })
}

pub(crate) fn cohomology(model: &mut Model, complex: &Complex, degree: usize) -> Result<Sheaf> {
    let Some(component) = complex.component_at(degree) else {
        let base = &complex.components()[0].sheaf;
        let zero = Rc::new(Matrix::default());
        return Sheaf::from_cover_maps(
            model,
            base.base(),
            base.prime,
            vec![0; base.dims().len()],
            |_, _, _| Ok(zero.clone()),
        );
    };
    let original = component.sheaf.as_ref();
    let prime = complex.prime();
    let incoming = degree.checked_sub(1).and_then(|q| complex.maps().get(q));
    let outgoing = complex.maps().get(degree);
    let quotients = original
        .dims()
        .iter()
        .enumerate()
        .map(|(point, &rank)| {
            StalkQuotient::build(
                model,
                prime,
                rank,
                incoming.map(|m| m.matrix(point)),
                outgoing.map(|m| m.matrix(point)),
            )
        })
        .collect::<Result<Vec<_>>>()?;
    induced_sheaf(model, original, &quotients)
}
