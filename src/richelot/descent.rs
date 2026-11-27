use std::rc::Rc;

use crate::ec_poly::{self as poly, Field};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::hyperelliptic_jacobian::{self, Divisor};
use crate::hyperelliptic_kummer_surface::{arithmetic, Point};
use crate::model::Model;
use crate::syntax::Label;

use super::correspondence::Correspondence;
use super::Splitting;

pub(crate) struct Descent {
    source: String,
    extension_target: String,
    degree: usize,
    permutation: [usize; 3],
}

impl Descent {
    pub(crate) fn define(
        model: &mut Model,
        target: Label,
        source: Label,
        extension_target: Label,
    ) -> Result<()> {
        let base = model.hww_curve(&source)?;
        if base.grade.genus != 2 || base.ring.precision() != 1 {
            return Err(Error::new(
                source.span,
                "Richelot descent source must be a smooth genus-two curve over a named precision-one field",
            ));
        }
        let isogeny = model.richelot(&extension_target)?;
        let extended_source = model.hww_curve(&Label {
            text: isogeny.source.clone(),
            span: extension_target.span,
        })?;
        let extended_equation = base
            .ring
            .embed_residue_coefficients(model, &extended_source.ring, base.polynomial())
            .map_err(|error| error.at(extension_target.span))?;
        if extended_source.polynomial() != extended_equation {
            return Err(Error::new(
                extension_target.span,
                "Richelot source equation is not the declared base extension",
            ));
        }

        let extension = model.hww_curve(&extension_target)?;
        let permutation = frobenius_permutation(
            model,
            Field {
                ring: &extension.ring,
            },
            base.ring.degree(),
            &isogeny.source_factors,
        )?
        .ok_or_else(|| {
            Error::new(
                extension_target.span,
                "Frobenius does not preserve the Richelot splitting",
            )
        })?;
        let preimage = base.ring.residue_preimage(model, &extension.ring)?;
        let polynomial = preimage
            .coefficients(model, extension.polynomial())?
            .ok_or_else(|| {
                Error::new(
                    extension_target.span,
                    "Richelot target coefficient lies outside the requested subfield",
                )
            })?;
        let descent = Self {
            source: source.text,
            extension_target: extension_target.text,
            degree: extension.ring.degree() / base.ring.degree(),
            permutation,
        };
        model.allocate(9 + descent.source.len() + descent.extension_target.len())?;
        let ring = Label {
            text: base.grade.ring.to_string(),
            span: source.span,
        };
        model.define_richelot_descent(target, &ring, polynomial, descent)
    }

    pub(crate) fn report(&self, name: &Label) -> String {
        let rational_kernel_order = 1 + self
            .permutation
            .iter()
            .enumerate()
            .filter(|&(i, &j)| i == j)
            .count();
        format!(
            "rvRichelotDescent[{}] = relative degree={}, permutation={:?}, rational kernel order={}",
            name.text, self.degree, self.permutation, rational_kernel_order
        )
    }
}

fn frobenius_permutation(
    model: &mut Model,
    field: Field<'_>,
    power: usize,
    splitting: &Splitting,
) -> Result<Option<[usize; 3]>> {
    let [a, b, c] = splitting;
    model.spend((a.len() + b.len() + c.len()) * field.cells())?;
    let monic = [
        poly::monic(model, field, a.clone())?,
        poly::monic(model, field, b.clone())?,
        poly::monic(model, field, c.clone())?,
    ];
    let mut permutation = [0; 3];
    for (index, factor) in monic.iter().enumerate() {
        let image = factor
            .iter()
            .map(|coefficient| field.ring.frobenius(model, coefficient, power))
            .collect::<Result<Vec<_>>>()?;
        let Some(target) = monic.iter().position(|other| other == &image) else {
            return Ok(None);
        };
        permutation[index] = target;
    }
    Ok(Some(permutation))
}

pub(crate) struct DescendedCorrespondence {
    registration: Rc<Descent>,
    source: Rc<Curve>,
    pub(crate) target: Rc<Curve>,
    extended: Correspondence,
}

impl DescendedCorrespondence {
    pub(crate) fn resolve(model: &Model, target: &Label) -> Result<Self> {
        let registration = model.richelot_descent(target)?;
        let source = model.hww_curve(&Label {
            text: registration.source.clone(),
            span: target.span,
        })?;
        if source.polynomial().len() != 6 || hyperelliptic_jacobian::require_curve(&source).is_err()
        {
            return Err(Error::new(
                target.span,
                "Richelot descent push requires a monic quintic source",
            ));
        }
        let curve = model
            .jacobian_curve(&target.text)
            .map_err(|error| error.at(target.span))?;
        let extended = Correspondence::resolve(
            model,
            &Label {
                text: registration.extension_target.clone(),
                span: target.span,
            },
        )?;
        Ok(Self {
            registration,
            source,
            target: curve,
            extended,
        })
    }

    pub(crate) fn source_name(&self) -> &str {
        &self.registration.source
    }

    pub(crate) fn push(&self, model: &mut Model, input: &Divisor) -> Result<Divisor> {
        let input = input.extend(model, &self.source, &self.extended.source)?;
        let output = self.extended.push(model, &input)?;
        output.descend(model, &self.extended.target, &self.target)
    }

    pub(crate) fn push_kummer(&self, model: &mut Model, input: &Point) -> Result<Point> {
        let input = input.extend(model, &self.source, &self.extended.source)?;
        let output = arithmetic::richelot(model, &self.extended, &input)?;
        output.descend(model, &self.extended.target, &self.target)
    }
}
