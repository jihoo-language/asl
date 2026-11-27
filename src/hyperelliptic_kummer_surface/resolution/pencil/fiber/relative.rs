use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::finite_field_factor::{factor_with_multiplicity as closed_places, remove_factor};
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;

use super::super::{Parameter, Pencil};
use super::{binary_quartic, cube, field, plane_quadratics, quartic, FiberRing, Jacobian};

mod arithmetic;
pub(crate) mod sections;

const DISCRIMINANT_DEGREE: usize = 24;

pub(crate) struct RelativeJacobian(Jacobian<Polynomial>);

impl Pencil {
    pub(crate) fn relative(&self, model: &mut Model) -> Result<&RelativeJacobian> {
        if let Some(relative) = self.relative.get() {
            return Ok(relative);
        }
        let field = field(&self.curve);
        let ring = BaseRing(field);
        let [u_axis, v_axis] = self.transverse_axes(field)?;
        // Work in the symbolic chart u=1. Each term involving Y once
        // has one transverse coefficient, so the Y exponent determines
        // the original homogeneous base degree including at infinity.
        let mut transverse = std::array::from_fn(|_| ring.zero());
        transverse[u_axis] = ring.one();
        transverse[v_axis] = field.polynomial(&[0, 1]);
        let [i, j] = &self.nodes;
        let x = j.coordinates().each_ref().map(|c| constant(field, c));
        let z = i.coordinates().each_ref().map(|c| constant(field, c));
        let terms: Vec<_> = quartic(model, &self.curve)?
            .into_iter()
            .map(|(slots, c)| (slots, constant(field, &c)))
            .collect();
        let [a, b, c] = plane_quadratics(model, &ring, &terms, [&x, &transverse, &z])?;
        for (offset, coefficients) in [&a, &b, &c].into_iter().enumerate() {
            for (i, coefficient) in coefficients.iter().enumerate() {
                check_degree(coefficient, i + offset)?;
            }
        }
        let binary = binary_quartic(model, &ring, &a, &b, &c)?;
        for (i, coefficient) in binary.iter().enumerate() {
            check_degree(coefficient, i + 2)?;
        }
        let jacobian = Jacobian::new(model, &ring, &binary)?;
        let relative = RelativeJacobian(jacobian);
        relative.verify(model, &ring)?;
        let cells: usize = relative
            .forms()
            .iter()
            .map(|(f, _)| f.len() * field.cells())
            .sum();
        model.allocate(7 + cells)?;
        Ok(self.relative.get_or_init(|| relative))
    }
}

impl RelativeJacobian {
    fn forms(&self) -> [(&Polynomial, usize); 6] {
        let Jacobian {
            a2,
            a4,
            a6,
            c4,
            c6,
            discriminant,
        } = &self.0;
        [
            (a2, 4),
            (a4, 8),
            (a6, 12),
            (c4, 8),
            (c6, 12),
            (discriminant, DISCRIMINANT_DEGREE),
        ]
    }

    fn verify(&self, model: &mut Model, ring: &BaseRing<'_>) -> Result<()> {
        for (form, degree) in self.forms() {
            check_degree(form, degree)?;
        }
        if self.0.discriminant.is_empty() {
            return Err(Error::plain(
                "K3 pencil relative discriminant is identically zero",
            ));
        }
        let c4_cube = cube(model, ring, &self.0.c4)?;
        let c6_square = ring.multiply(model, &self.0.c6, &self.0.c6)?;
        let difference = ring.subtract(model, &c4_cube, &c6_square)?;
        if difference != ring.scale(model, &self.0.discriminant, 1728)? {
            return Err(Error::plain(
                "K3 pencil relative Weierstrass invariants failed their identity",
            ));
        }
        Ok(())
    }

    pub(crate) fn weierstrass(&self, pencil: &Pencil, name: &str) -> String {
        let field = field(&pencil.curve);
        let forms = self.forms();
        let [a2, a4, a6, c4, c6, _] = forms.map(|(form, degree)| homogeneous(field, form, degree));
        let infinity: Vec<_> = forms
            .iter()
            .map(|&(p, n)| coefficient(field, p, n))
            .collect();
        format!(
            "rvK3PencilWeierstrass[{name}] = a2={a2:?}, a4={a4:?}, a6={a6:?}, c4={c4:?}, c6={c6:?}, infinity (a2, a4, a6, c4, c6, discriminant)={infinity:?}"
        )
    }

    pub(crate) fn discriminant(
        &self,
        model: &mut Model,
        pencil: &Pencil,
        name: &str,
    ) -> Result<String> {
        let field = field(&pencil.curve);
        let form = &self.0.discriminant;
        let infinity = DISCRIMINANT_DEGREE + 1 - form.len();
        let places = closed_places(model, field, form)?;
        let total = infinity + places.iter().map(|(d, v, _)| d * v).sum::<usize>();
        if total != DISCRIMINANT_DEGREE {
            return Err(Error::plain(
                "K3 pencil discriminant divisor has the wrong total degree",
            ));
        }
        let coefficients = homogeneous(field, form, DISCRIMINANT_DEGREE);
        Ok(format!(
            "rvK3PencilDiscriminant[{name}] = degree={total}, coefficients={coefficients:?}, infinity order={infinity}, finite places (degree, order, monic polynomial)={places:?}"
        ))
    }

    pub(crate) fn reduction(
        &self,
        model: &mut Model,
        pencil: &Pencil,
        name: &str,
        parameter: &Parameter,
    ) -> Result<String> {
        let field = field(&pencil.curve);
        let uniformizer = local_uniformizer(model, field, parameter)?;
        let [_, v] = parameter.coordinates();
        let mut valuations = [None; 6];
        for (out, (form, degree)) in valuations.iter_mut().zip(self.forms()) {
            *out = order(
                model,
                field,
                form,
                degree,
                uniformizer.as_ref().map(|linear| linear.as_slice()),
            )?;
        }
        let [a2, a4, a6, c4, c6, discriminant] = valuations.map(display_order);
        let chart = if uniformizer.is_some() {
            format!("v/u-{v:?}")
        } else {
            "u/v".into()
        };
        let local = if field.ring.prime() == 3 {
            "Kodaira=unclassified in characteristic 3".into()
        } else {
            tame_reduction(valuations[3], valuations[4], valuations[5])?
        };
        Ok(format!(
            "rvK3FiberReduction[{name}] = uniformizer={chart}, raw orders (a2, a4, a6, c4, c6, discriminant)=[{a2}, {a4}, {a6}, {c4}, {c6}, {discriminant}], {local}"
        ))
    }
}

fn local_uniformizer(
    model: &mut Model,
    field: Field<'_>,
    parameter: &Parameter,
) -> Result<Option<[Scalar; 2]>> {
    let [u, v] = parameter.coordinates();
    if field.is_zero(u) {
        Ok(None)
    } else {
        // Parameters on this chart are normalized to [1:v].
        Ok(Some([field.ring.negate(model, v)?, field.one()]))
    }
}

fn coefficient(field: Field<'_>, p: &[Scalar], i: usize) -> Scalar {
    p.get(i).cloned().unwrap_or_else(|| field.zero())
}

fn homogeneous(field: Field<'_>, p: &[Scalar], degree: usize) -> Vec<Scalar> {
    (0..=degree).map(|i| coefficient(field, p, i)).collect()
}

fn check_degree(p: &[Scalar], degree: usize) -> Result<()> {
    if p.len() > degree + 1 {
        return Err(Error::plain(
            "K3 pencil coefficient has the wrong homogeneous degree",
        ));
    }
    Ok(())
}

fn order(
    model: &mut Model,
    field: Field<'_>,
    p: &[Scalar],
    degree: usize,
    uniformizer: Option<&[Scalar]>,
) -> Result<Option<usize>> {
    if p.is_empty() {
        return Ok(None);
    }
    if let Some(uniformizer) = uniformizer {
        let mut p = p.to_vec();
        let multiplicity = remove_factor(model, field, &mut p, uniformizer)?;
        Ok(Some(multiplicity))
    } else {
        Ok(Some(degree + 1 - p.len()))
    }
}

fn display_order(order: Option<usize>) -> String {
    order.map_or_else(|| "inf".into(), |n| n.to_string())
}

#[derive(Clone, Copy)]
enum Kodaira {
    Good,
    Multiplicative(usize),
    Star(usize),
    Other(&'static str),
}

impl Kodaira {
    fn symbol(self) -> String {
        match self {
            Self::Good => "I0".into(),
            Self::Multiplicative(n) => format!("I{n}"),
            Self::Star(n) => format!("I{n}*"),
            Self::Other(symbol) => symbol.into(),
        }
    }
}

struct MinimalReduction {
    scales: usize,
    c4: Option<usize>,
    c6: Option<usize>,
    delta: usize,
    j_order: Option<i64>,
    kodaira: Kodaira,
}

impl MinimalReduction {
    fn new(c4: Option<usize>, c6: Option<usize>, discriminant: Option<usize>) -> Result<Self> {
        let failure = || Error::plain("K3 pencil has inconsistent minimal Weierstrass valuations");
        let raw = discriminant.ok_or_else(failure)?;
        // Short model in residue characteristic at least 5. Its A and B
        // differ from c4 and c6 by units, so these are the maximal scalings.
        let scales = (c4.unwrap_or(usize::MAX) / 4).min(c6.unwrap_or(usize::MAX) / 6);
        let delta = raw
            .checked_sub(12usize.checked_mul(scales).ok_or_else(failure)?)
            .ok_or_else(failure)?;
        let c4 = c4.map(|n| n - 4 * scales);
        let c6 = c6.map(|n| n - 6 * scales);
        let j_order = c4.map(|n| (3 * n) as i64 - delta as i64);
        let kodaira = match (delta, c4) {
            (0, _) => Kodaira::Good,
            (n, Some(0)) => Kodaira::Multiplicative(n),
            (2, _) => Kodaira::Other("II"),
            (3, _) => Kodaira::Other("III"),
            (4, _) => Kodaira::Other("IV"),
            (6, _) => Kodaira::Star(0),
            (n, _) if n >= 7 && j_order.is_some_and(|v| v < 0) => Kodaira::Star(n - 6),
            (8, _) => Kodaira::Other("IV*"),
            (9, _) => Kodaira::Other("III*"),
            (10, _) => Kodaira::Other("II*"),
            _ => return Err(failure()),
        };
        Ok(Self {
            scales,
            c4,
            c6,
            delta,
            j_order,
            kodaira,
        })
    }
}

fn tame_reduction(
    c4: Option<usize>,
    c6: Option<usize>,
    discriminant: Option<usize>,
) -> Result<String> {
    let reduction = MinimalReduction::new(c4, c6, discriminant)?;
    let MinimalReduction {
        scales,
        c4,
        c6,
        delta,
        j_order,
        kodaira,
    } = reduction;
    let j = j_order.map_or_else(|| "inf".into(), |v| v.to_string());
    let kodaira = kodaira.symbol();
    Ok(format!(
        "scale={scales}, minimal orders (c4, c6, discriminant)=[{}, {}, {delta}], j order={j}, geometric Kodaira={kodaira}",
        display_order(c4), display_order(c6)
    ))
}

fn constant(field: Field<'_>, value: &Scalar) -> Polynomial {
    if field.is_zero(value) {
        Vec::new()
    } else {
        vec![value.clone()]
    }
}

struct BaseRing<'a>(Field<'a>);

impl CoefficientRing for BaseRing<'_> {
    type Element = Polynomial;

    fn cells(&self) -> usize {
        (DISCRIMINANT_DEGREE + 1) * self.0.cells()
    }
    fn zero(&self) -> Polynomial {
        Vec::new()
    }
    fn one(&self) -> Polynomial {
        vec![self.0.one()]
    }
    fn is_zero(&self, p: &Polynomial) -> bool {
        p.is_empty()
    }
    fn is_unit(&self, p: &Polynomial) -> bool {
        p.len() == 1
    }
    fn add(&self, model: &mut Model, p: &Polynomial, q: &Polynomial) -> Result<Polynomial> {
        poly::sum(model, self.0, p, q, false)
    }
    fn subtract(&self, model: &mut Model, p: &Polynomial, q: &Polynomial) -> Result<Polynomial> {
        poly::sum(model, self.0, p, q, true)
    }
    fn multiply(&self, model: &mut Model, p: &Polynomial, q: &Polynomial) -> Result<Polynomial> {
        poly::multiply(model, self.0, p, q)
    }
    fn inverse_unit(&self, model: &mut Model, p: &Polynomial) -> Result<Polynomial> {
        if let [scalar] = p.as_slice() {
            return Ok(vec![self.0.inverse_unit(model, scalar)?]);
        }
        Err(Error::plain(
            "K3 base polynomial unit must be a nonzero constant",
        ))
    }
}

impl FiberRing for BaseRing<'_> {
    fn scale(&self, model: &mut Model, p: &Polynomial, factor: u64) -> Result<Polynomial> {
        let factor = self.0.signed(i128::from(factor));
        poly::scale(model, self.0, p, &factor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rv_unram::Ring;

    #[test]
    fn closed_places_keep_a_nonsplit_factor_with_inseparable_multiplicity() {
        let mut model = Model::new();
        let ring =
            Ring::new(&mut model, 3, 1, vec![1, 0, 1]).unwrap_or_else(|error| panic!("{error}"));
        let field = Field { ring: &ring };
        let nonsquare = vec![1, 1];
        assert_eq!(
            ring.power(&mut model, &nonsquare, 4)
                .unwrap_or_else(|error| panic!("{error}")),
            field.signed(-1)
        );
        let negative = ring
            .negate(&mut model, &nonsquare)
            .unwrap_or_else(|error| panic!("{error}"));
        let linear = vec![negative.clone(), field.one()];
        let quadratic = vec![negative, field.zero(), field.one()];

        let mut polynomial = poly::multiply(&mut model, field, &linear, &linear)
            .unwrap_or_else(|error| panic!("{error}"));
        for _ in 0..6 {
            polynomial = poly::multiply(&mut model, field, &polynomial, &quadratic)
                .unwrap_or_else(|error| panic!("{error}"));
        }
        assert_eq!(
            closed_places(&mut model, field, &polynomial).unwrap_or_else(|error| panic!("{error}")),
            [(1, 2, linear), (2, 6, quadratic)]
        );
    }

    #[test]
    fn minimalization_and_j_order_select_the_geometric_type() {
        let scaled_good =
            tame_reduction(Some(4), Some(6), Some(12)).unwrap_or_else(|error| panic!("{error}"));
        assert!(scaled_good.starts_with("scale=1, minimal orders (c4, c6, discriminant)=[0, 0, 0]"));
        assert!(scaled_good.ends_with("geometric Kodaira=I0"));
        let potentially_multiplicative =
            tame_reduction(Some(2), Some(3), Some(9)).unwrap_or_else(|error| panic!("{error}"));
        assert!(potentially_multiplicative.ends_with("j order=-3, geometric Kodaira=I3*"));
        let potentially_good =
            tame_reduction(Some(3), Some(5), Some(9)).unwrap_or_else(|error| panic!("{error}"));
        assert!(potentially_good.ends_with("j order=0, geometric Kodaira=III*"));
    }
}
