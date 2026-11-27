use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::finite_field_factor::{factor_squarefree, factor_with_multiplicity};
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;

use super::{field, from_polynomial, on_curve, Pencil, Section};

impl Pencil {
    pub(crate) fn torsion_roots(&self, model: &mut Model) -> Result<&[Polynomial]> {
        if let Some(roots) = self.torsion_roots.get() {
            return Ok(roots);
        }
        let field = field(&self.curve);
        let relative = self.relative(model)?;
        let cubic = vec![
            relative.0.a6.clone(),
            relative.0.a4.clone(),
            relative.0.a2.clone(),
            vec![field.one()],
        ];
        let roots = rational_roots(model, field, cubic)?;
        let y = from_polynomial(field, Vec::new());
        for root in &roots {
            let x = from_polynomial(field, root.clone());
            if !on_curve(model, field, relative, &x, &y)? {
                return Err(Error::plain(
                    "K3 pencil two-torsion candidate failed its equation",
                ));
            }
        }
        let cells: usize = roots
            .iter()
            .map(|root| 1 + root.len() * field.cells())
            .sum();
        model.allocate(1 + cells)?;
        Ok(self.torsion_roots.get_or_init(|| roots))
    }

    pub(crate) fn torsion_at(&self, model: &mut Model, index: usize) -> Result<Section> {
        let field = field(&self.curve);
        let roots = self.torsion_roots(model)?;
        let root = roots
            .get(index)
            .ok_or_else(|| Error::plain("K3 pencil affine two-torsion index is out of range"))?;
        Ok(Section::Affine(
            from_polynomial(field, root.clone()),
            from_polynomial(field, Vec::new()),
        ))
    }

    pub(crate) fn torsion_report(&self, model: &mut Model, name: &str) -> Result<String> {
        let roots = self.torsion_roots(model)?;
        Ok(format!(
            "rvK3PencilTorsion2[{name}] = order={}, affine x-roots={roots:?}",
            1 + roots.len()
        ))
    }
}

fn rational_roots(
    model: &mut Model,
    field: Field<'_>,
    mut cubic: Vec<Polynomial>,
) -> Result<Vec<Polynomial>> {
    let mut roots = Vec::new();
    if cubic[0].is_empty() {
        roots.push(Vec::new());
        cubic.remove(0);
    }
    if cubic[0].is_empty() {
        return Err(Error::plain("generic K3 pencil cubic has a repeated root"));
    }
    for divisor in small_monic_divisors(model, field, &cubic[0])? {
        for scalar in scalar_roots(model, field, &cubic, &divisor)? {
            if scalar == field.zero() {
                return Err(Error::plain(
                    "K3 pencil monic-divisor search found a zero scale",
                ));
            }
            roots.push(poly::scale(model, field, &divisor, &scalar)?);
        }
    }
    roots.sort_unstable();
    if !matches!(roots.len(), 0 | 1 | 3) || roots.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(Error::plain(
            "generic K3 pencil two-torsion roots are inconsistent",
        ));
    }
    Ok(roots)
}

fn small_monic_divisors(
    model: &mut Model,
    field: Field<'_>,
    constant: &[Scalar],
) -> Result<Vec<Polynomial>> {
    // Integrality gives polynomial x. A degree > 4 makes the highest x
    // power strictly dominate all others, also after peeling the zero root.
    const X_DEGREE: usize = 4;
    let mut divisors = vec![vec![field.one()]];
    for (degree, multiplicity, factor) in factor_with_multiplicity(model, field, constant)? {
        for i in 0..divisors.len() {
            let mut divisor = divisors[i].clone();
            let remaining_degree = X_DEGREE - (divisor.len() - 1);
            for _ in 0..multiplicity.min(remaining_degree / degree) {
                divisor = poly::multiply(model, field, &divisor, &factor)?;
                field.check_polynomial((divisors.len() + 1) * (X_DEGREE + 1))?;
                model.spend(divisor.len() * field.cells())?;
                divisors.push(divisor.clone());
            }
        }
    }
    Ok(divisors)
}

fn scalar_roots(
    model: &mut Model,
    field: Field<'_>,
    cubic: &[Polynomial],
    divisor: &[Scalar],
) -> Result<Vec<Scalar>> {
    // Form the coefficient polynomials in u for f(u d(t)). Their gcd
    // captures precisely the scalar choices that make the whole sum zero.
    let mut coefficients = Vec::with_capacity(cubic.len());
    let mut power = vec![field.one()];
    for (exponent, term) in cubic.iter().enumerate() {
        if exponent != 0 {
            power = poly::multiply(model, field, &power, divisor)?;
        }
        coefficients.push(poly::multiply(model, field, &power, term)?);
    }
    let zero = field.zero();
    let width = coefficients.iter().map(Vec::len).max().unwrap_or(0);
    let mut gcd = Vec::new();
    for i in 0..width {
        let polynomial = poly::trim(
            coefficients
                .iter()
                .map(|row| row.get(i).unwrap_or(&zero).clone())
                .collect(),
        );
        if polynomial.is_empty() {
            continue;
        }
        gcd = poly::monic_gcd(model, field, gcd, polynomial)?;
        if gcd.len() == 1 {
            return Ok(Vec::new());
        }
    }
    if gcd.is_empty() {
        return Err(Error::plain(
            "K3 pencil scalar root search lost its monic leading term",
        ));
    }
    // A multiple root of the scalar gcd would give a multiple root of
    // the generic cubic, already excluded by its nonzero discriminant.
    factor_squarefree(model, field, gcd)?
        .into_iter()
        .filter(|factor| factor.len() == 2)
        .map(|linear| field.ring.negate(model, &linear[0]))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rv_unram::Ring;

    #[test]
    fn zero_root_peels_and_scalar_gcd_finds_two_scales_of_one_divisor() {
        let mut model = Model::new();
        let ring =
            Ring::new(&mut model, 5, 1, vec![0, 1]).unwrap_or_else(|error| panic!("{error}"));
        let field = Field { ring: &ring };
        // x (x - d) (x - 2d), d=t²+1. The two nonzero roots have
        // the same monic divisor but different leading coefficients.
        let polynomial = vec![
            Vec::new(),
            field.polynomial(&[2, 0, 4, 0, 2]),
            field.polynomial(&[-3, 0, -3]),
            vec![field.one()],
        ];
        let roots =
            rational_roots(&mut model, field, polynomial).unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(
            roots,
            [
                Vec::new(),
                field.polynomial(&[1, 0, 1]),
                field.polynomial(&[2, 0, 2])
            ]
        );

        let nonsplit = vec![
            Vec::new(),
            field.polynomial(&[-2]),
            Vec::new(),
            vec![field.one()],
        ];
        assert_eq!(
            rational_roots(&mut model, field, nonsplit).unwrap_or_else(|error| panic!("{error}")),
            [Vec::<Scalar>::new()],
        );
    }
}
