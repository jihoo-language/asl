use std::cell::OnceCell;
use std::cmp::Ordering;
use std::rc::Rc;

use crate::ec_poly::{Field, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::Matrix;
use crate::hw_closed::Curve;
use crate::hyperelliptic_jacobian::require_extension;
use crate::jacobian_level::Basis;
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

use super::{point_count, Action};
use crate::hyperelliptic_kummer_surface::{resolution::ResolvedPoint, theta, BaseOperation};

pub(crate) struct Declare(
    pub(crate) Label,
    pub(crate) Label,
    pub(crate) Label,
    pub(crate) Label,
    pub(crate) usize,
);

pub(crate) enum Expression {
    Orbit(Label, Box<ParsedExpr>),
    Read(Read, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Read {
    Split,
    Representative,
    Lift,
}

pub(crate) struct Quotient {
    pub(crate) base_name: Rc<str>,
    pub(crate) cover_name: Rc<str>,
    base: Rc<Curve>,
    cover: Rc<Curve>,
    basis: Rc<Basis>,
    mask: usize,
    count: u128,
    translation: OnceCell<Matrix<Scalar>>,
}

pub(crate) struct Orbit {
    pub(crate) representative: Rc<ResolvedPoint>,
    pub(crate) split: bool,
}

impl Orbit {
    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        self.representative.compare(model, &other.representative)
    }

    fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(3)?;
        Ok(Value::Enriques(Rc::new(self)))
    }
}

impl Quotient {
    pub(crate) fn new(
        model: &mut Model,
        base: Label,
        cover: Label,
        basis: &Label,
        mask: usize,
    ) -> Result<Self> {
        let base_curve = model
            .kummer_curve(&base.text)
            .map_err(|error| error.at(base.span))?;
        let cover_curve = model
            .kummer_curve(&cover.text)
            .map_err(|error| error.at(cover.span))?;
        if cover_curve.ring.degree() != 2 * base_curve.ring.degree() {
            return Err(Error::new(
                cover.span,
                "Enriques cover must use a quadratic extension of the base field",
            ));
        }
        require_extension(model, &base_curve, &cover_curve)
            .map_err(|error| error.at(cover.span))?;
        let basis = model.jacobian_level(basis)?;
        let action = Action::validate(model, &basis, &base_curve, mask)?;
        let count = point_count(
            base_curve.ring.size_of_residue_field()?,
            1,
            action.trace(model, 1)?,
        )?;
        model.allocate(10 + base.text.len() + cover.text.len())?;
        Ok(Self {
            base_name: base.text.into(),
            cover_name: cover.text.into(),
            base: base_curve,
            cover: cover_curve,
            basis,
            mask,
            count,
            translation: OnceCell::new(),
        })
    }

    pub(crate) fn count(&self) -> usize {
        usize::try_from(self.count).unwrap_or(usize::MAX)
    }

    pub(crate) fn orbit(
        &self,
        model: &mut Model,
        point: Rc<ResolvedPoint>,
        from_base: bool,
    ) -> Result<Value> {
        let point = if from_base {
            point
                .base_change(model, BaseOperation::Extend, &self.base, &self.cover)?
                .keep(model)?
        } else {
            point
        };
        let Some((split, ordering, mate)) = self.classify(model, &point)? else {
            return Err(Error::plain(
                "resolved point does not give a Frobenius-stable Enriques orbit",
            ));
        };
        let representative = if ordering.is_lt() {
            point
        } else {
            mate.keep(model)?
        };
        Orbit {
            representative,
            split,
        }
        .persist(model)
    }

    pub(crate) fn lift(&self, model: &mut Model, orbit: &Orbit) -> Result<Value> {
        if !orbit.split {
            return Err(Error::plain("nonsplit Enriques orbit has no rational lift"));
        }
        orbit
            .representative
            .base_change(model, BaseOperation::Descend, &self.cover, &self.base)?
            .persist(model)
    }

    pub(crate) fn elements(&self, model: &mut Model, limit: usize) -> Result<Vec<Value>> {
        if self.count > limit as u128 {
            return Err(Error::plain(
                "Enriques quotient exceeds the finite element limit",
            ));
        }
        let ty = model.kummer_k3_type(self.cover_name.clone())?;
        let cover_points = model.elements(ty)?;
        let mut values = Vec::with_capacity(self.count());
        for point in cover_points.iter() {
            let Value::KummerK3(point) = point else {
                unreachable!()
            };
            if let Some((split, Ordering::Less, _)) = self.classify(model, point)? {
                values.push(
                    Orbit {
                        representative: point.clone(),
                        split,
                    }
                    .persist(model)?,
                );
            }
        }
        if values.len() as u128 != self.count {
            return Err(Error::plain(
                "enumerated Enriques orbits disagree with the arithmetic point count",
            ));
        }
        Ok(values)
    }

    fn classify(
        &self,
        model: &mut Model,
        point: &ResolvedPoint,
    ) -> Result<Option<(bool, Ordering, ResolvedPoint)>> {
        let conjugate =
            point.base_change(model, BaseOperation::Frobenius(1), &self.cover, &self.base)?;
        let split = point.compare(model, &conjugate)?.is_eq();
        let mate = self.involution(model, point)?;
        let ordering = point.compare(model, &mate)?;
        if ordering.is_eq() {
            return Err(Error::plain(
                "Enriques involution unexpectedly fixes a resolved point",
            ));
        }
        if split || mate.compare(model, &conjugate)?.is_eq() {
            Ok(Some((split, ordering, mate)))
        } else {
            Ok(None)
        }
    }

    fn involution(&self, model: &mut Model, point: &ResolvedPoint) -> Result<ResolvedPoint> {
        let translation = self.translation(model)?;
        let switched = point.gauss_switch(model, &self.cover)?;
        let field = Field {
            ring: &self.cover.ring,
        };
        switched.projective_image(model, &self.cover, |model, coordinates| {
            theta::apply(model, field, translation, coordinates)
        })
    }

    fn translation(&self, model: &mut Model) -> Result<&Matrix<Scalar>> {
        if let Some(translation) = self.translation.get() {
            return Ok(translation);
        }
        let group = self.basis.theta(model)?;
        // A fixed theta mask gives a normalized matrix over the base.
        // Descend first, since neither named extension need contain the other.
        let preimage = self.base.ring.residue_preimage(model, &group.curve.ring)?;
        let coefficients = preimage
            .coefficients(model, group.matrix(self.mask).entries())?
            .ok_or_else(|| {
                Error::plain("Enriques translation matrix does not descend to the base field")
            })?;
        let coefficients =
            self.base
                .ring
                .embed_residue_coefficients(model, &self.cover.ring, &coefficients)?;
        model.allocate(16 * (3 + self.cover.ring.degree()))?;
        let matrix = Matrix::from_entries(4, coefficients)?;
        Ok(self.translation.get_or_init(|| matrix))
    }
}
