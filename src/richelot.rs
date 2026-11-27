use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;
use crate::syntax::Label;

pub(crate) mod correspondence;

type Literal = Vec<Vec<i128>>;
type Splitting = [Polynomial; 3];

pub(crate) enum Command {
    Define(Label, Label, [Literal; 3]),
    Dual(Label, Label),
    Data(Label),
}

pub(crate) struct Isogeny {
    pub(crate) source: String,
    pub(crate) source_factors: Splitting,
    pub(crate) target_factors: Splitting,
    pub(crate) determinant: Scalar,
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<Option<String>> {
        let (target, source, literal) = match self {
            Self::Define(target, source, literal) => (target, source, Some(literal)),
            Self::Dual(target, source) => (target, source, None),
            Self::Data(target) => {
                let isogeny = model.richelot(&target)?;
                let curve = model.hww_curve(&target)?;
                let cells = isogeny.cells(Field { ring: &curve.ring });
                model.spend(cells + curve.polynomial().len() * curve.ring.degree())?;
                return Ok(Some(format!(
                    "rvRichelotData[{}] = source={}, determinant={:?}, f={:?}, factors={:?}",
                    target.text,
                    isogeny.source,
                    isogeny.determinant,
                    curve.polynomial(),
                    isogeny.target_factors
                )));
            }
        };
        let curve = model.hww_curve(&source)?;
        if curve.ring.precision() != 1 || curve.grade.genus != 2 {
            return Err(Error::new(
                source.span,
                "Richelot source must be a smooth genus-two curve over a named precision-one field",
            ));
        }
        let field = Field { ring: &curve.ring };
        let factors = if let Some([a, b, c]) = literal {
            let make = |model: &mut Model, coefficients: Literal| -> Result<Polynomial> {
                field.check_polynomial(coefficients.len())?;
                let values = coefficients
                    .into_iter()
                    .map(|c| field.ring.literal(model, c))
                    .collect::<Result<_>>()?;
                Ok(poly::trim(values))
            };
            [make(model, a)?, make(model, b)?, make(model, c)?]
        } else {
            model.richelot(&source)?.target_factors.clone()
        };
        if factors.iter().any(|g| !(2..=3).contains(&g.len())) {
            return Err(Error::new(
                target.span,
                "Richelot splitting factors must have degree one or two",
            ));
        }
        if product(model, field, &factors)? != curve.polynomial() {
            return Err(Error::new(
                target.span,
                "Richelot factors must multiply to the source equation",
            ));
        }
        let (determinant, target_factors) =
            transform(model, field, &factors).map_err(|error| error.at(target.span))?;
        let polynomial = product(model, field, &target_factors)?;
        let data = Isogeny {
            source: source.text,
            source_factors: factors,
            target_factors,
            determinant,
        };
        model.allocate(data.cells(field))?;
        let ring = Label {
            text: curve.grade.ring.to_string(),
            span: source.span,
        };
        model.define_richelot(target, &ring, polynomial, data)?;
        Ok(None)
    }
}

impl Isogeny {
    fn cells(&self, field: Field<'_>) -> usize {
        let coefficients: usize = self
            .source_factors
            .iter()
            .chain(&self.target_factors)
            .map(Vec::len)
            .sum();
        8 + self.source.len() + (coefficients + 1) * field.cells()
    }
}

fn product(model: &mut Model, field: Field<'_>, [a, b, c]: &Splitting) -> Result<Polynomial> {
    let ab = poly::multiply(model, field, a, b)?;
    poly::multiply(model, field, &ab, c)
}

fn transform(
    model: &mut Model,
    field: Field<'_>,
    factors: &Splitting,
) -> Result<(Scalar, Splitting)> {
    let zero = field.zero();
    let [a, b, c] = factors
        .each_ref()
        .map(|g| std::array::from_fn(|i| g.get(i).unwrap_or(&zero)));
    let first = cross_product(model, field, b, c)?;
    let mut determinant = field.zero();
    for (coefficient, minor) in a.into_iter().zip(&first) {
        let term = field.multiply(model, coefficient, minor)?;
        determinant = field.add(model, &determinant, &term)?;
    }
    if field.is_zero(&determinant) {
        return Err(Error::plain(
            "Richelot determinant is zero, so the quotient is a product of elliptic curves over the algebraic closure",
        ));
    }
    let second = cross_product(model, field, c, a)?;
    let third = cross_product(model, field, a, b)?;
    let bracket = |model: &mut Model, [a, b, c]: [Scalar; 3]| -> Result<Polynomial> {
        // The ascending coefficients of G' H - G H' are (-c, 2b, -a)
        // when (a,b,c) is the cross product of the coefficient rows of G and H.
        let constant = field.ring.negate(model, &c)?;
        let linear = field.ring.scale(model, &b, 2)?;
        let quadratic = field.ring.negate(model, &a)?;
        Ok(poly::trim(vec![constant, linear, quadratic]))
    };
    let first = bracket(model, first)?;
    let inverse = field.inverse_unit(model, &determinant)?;
    let first = poly::scale(model, field, &first, &inverse)?;
    Ok((
        determinant,
        [first, bracket(model, second)?, bracket(model, third)?],
    ))
}

fn cross_product(
    model: &mut Model,
    field: Field<'_>,
    a: [&Scalar; 3],
    b: [&Scalar; 3],
) -> Result<[Scalar; 3]> {
    let mut result = std::array::from_fn(|_| field.zero());
    for (i, coordinate) in result.iter_mut().enumerate() {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        let left = field.multiply(model, a[j], b[k])?;
        let right = field.multiply(model, a[k], b[j])?;
        *coordinate = field.subtract(model, &left, &right)?;
    }
    Ok(result)
}
