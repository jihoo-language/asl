use std::rc::Rc;

use crate::dwork_k3::fricke::{certify_endpoint_quadratic, ENDPOINT_PRODUCT, ENDPOINT_SUM};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::{cm_neighbors, cm_ring::Order};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_unram::Ring;

use super::{decimal, monic_norm, Integer, Poly, Polynomial as ClassPolynomial, LIMBS};

mod exact;
mod residue_field;

pub(crate) use residue_field::{endpoint_report, pair_coefficient};

pub(crate) struct Polynomial {
    support: Poly,
    kind: Kind,
}

enum Kind {
    Norm {
        discriminant: i128,
        strata: Vec<(usize, Poly)>,
    },
    Edge {
        sources: [Rc<ClassPolynomial>; 2],
        strata: Vec<(usize, usize, Poly)>,
        horizontal: Option<Horizontal>,
    },
}

struct Horizontal {
    distinct: Poly,
    fixed: Poly,
    paired: Poly,
}

struct Norm {
    raw: Poly,
    squarefree: Poly,
    strata: Vec<(usize, Poly)>,
    remainder: [Poly; 2],
}

impl Norm {
    fn new(model: &mut Model, source: &ClassPolynomial) -> Result<Self> {
        certify_endpoint_quadratic()?;
        let sum = exact::integers(&ENDPOINT_SUM);
        let product = exact::integers(&ENDPOINT_PRODUCT);

        // Horner in the quadratic algebra J^2 = S J - P.
        let [mut a, mut b] = [exact::zero(), exact::zero()];
        for &c in source.coefficients.0.iter().rev() {
            let constant = Poly(vec![c]);
            let next_a = constant.subtract(&product.multiply(model, &b)?)?;
            let next_b = exact::plus(&a, &sum.multiply(model, &b)?)?;
            (a, b) = (next_a, next_b);
        }
        let a2 = a.multiply(model, &a)?;
        let ab = a.multiply(model, &b)?;
        let ab = sum.multiply(model, &ab)?;
        let b2 = b.multiply(model, &b)?;
        let b2 = product.multiply(model, &b2)?;
        let raw = exact::plus(&exact::plus(&a2, &ab)?, &b2)?;
        if !raw.is_monic() || exact::degree(&raw) != 3 * source.degree() {
            return Err(certificate());
        }

        // In the other quotient algebra this is the determinant of
        // multiplication by J² - S J + P on Z[H][J]/H_D.
        let independent = monic_norm(model, &source.coefficients, 3, |i| {
            let s = ENDPOINT_SUM.get(i).copied().unwrap_or(0);
            Ok(vec![
                Integer::from_i128(ENDPOINT_PRODUCT[i]),
                Integer::from_i128(-s),
                Integer::from_i128(i128::from(i == 0)),
            ])
        })?;
        if !independent.same(&raw) {
            return Err(certificate());
        }

        let strata = exact::strata(model, &raw)?;
        let squarefree = exact::radical(model, &strata)?;
        Ok(Self {
            raw,
            squarefree,
            strata,
            remainder: [a, b],
        })
    }
}

impl Polynomial {
    pub(crate) fn norm(model: &mut Model, source: &ClassPolynomial) -> Result<Self> {
        let norm = Norm::new(model, source)?;
        Self::keep(
            model,
            norm.raw,
            Kind::Norm {
                discriminant: source.discriminant,
                strata: norm.strata,
            },
        )
    }

    pub(crate) fn edge(
        model: &mut Model,
        left: Rc<ClassPolynomial>,
        right: Rc<ClassPolynomial>,
    ) -> Result<Self> {
        let l = Norm::new(model, &left)?;
        let ld = left.discriminant;
        let rd = right.discriminant;
        let (support, strata, horizontal) = if ld == rd {
            if !left.same(&right) {
                return Err(certificate());
            }
            let horizontal = horizontal(model, &l)?;
            let support = horizontal
                .distinct
                .multiply(model, &horizontal.fixed)?
                .multiply(model, &horizontal.paired)?;
            let strata = orientations(model, &support, &l.strata, &l.strata)?;
            (support, strata, Some(horizontal))
        } else {
            let r = Norm::new(model, &right)?;
            let support = exact::gcd(model, &l.squarefree, &r.squarefree)?;
            certify_endpoints(model, &support, &l, &r, &left, &right)?;
            let strata = orientations(model, &support, &l.strata, &r.strata)?;
            (support, strata, None)
        };
        let counts = directed_degrees(&strata);
        if counts
            != [
                expected_edges(model, &left, rd)?,
                expected_edges(model, &right, ld)?,
            ]
        {
            return Err(certificate());
        }
        Self::keep(
            model,
            support,
            Kind::Edge {
                sources: [left, right],
                strata,
                horizontal,
            },
        )
    }

    fn keep(model: &mut Model, support: Poly, kind: Kind) -> Result<Self> {
        let stored = match &kind {
            Kind::Norm { strata, .. } => strata.iter().map(|(_, p)| p.0.len()).sum::<usize>(),
            Kind::Edge {
                strata, horizontal, ..
            } => {
                let count = strata.iter().map(|(_, _, p)| p.0.len()).sum::<usize>();
                count
                    + horizontal.as_ref().map_or(0, |h| {
                        h.distinct.0.len() + h.fixed.0.len() + h.paired.0.len()
                    })
            }
        };
        model.allocate((support.0.len() + stored) * (LIMBS + 1))?;
        Ok(Self { support, kind })
    }

    pub(crate) fn coefficient(&self, index: usize) -> String {
        decimal(self.support.coefficient(index))
    }

    pub(crate) fn report(&self, name: &str) -> String {
        let degree = exact::degree(&self.support);
        let coefficients = exact::show(&self.support);
        let detail = match &self.kind {
            Kind::Norm {
                discriminant,
                strata,
            } => {
                let strata_text = strata
                    .iter()
                    .map(|(m, factor)| format!("({m}, {})", exact::show(factor)))
                    .collect::<Vec<_>>()
                    .join(", ");
                format!("order={discriminant}, directed multiplicity strata=[{strata_text}], quadratic norm and Bareiss resultant agree")
            }
            Kind::Edge {
                sources,
                strata,
                horizontal,
            } => {
                let [a, b] = sources;
                let strata_text = strata
                    .iter()
                    .map(|(m, n, factor)| format!("({m}, {n}, {})", exact::show(factor)))
                    .collect::<Vec<_>>()
                    .join(", ");
                let mut result = format!(
                    "orders=({}, {}), directed degrees={:?}, directed multiplicity strata=[{strata_text}]",
                    a.discriminant, b.discriminant, directed_degrees(strata),
                );
                if let Some(h) = horizontal {
                    result.push_str(&format!(
                        ", distinct endpoints={}, Fricke-fixed loops={}, nonfixed paired loops={}",
                        exact::show(&h.distinct),
                        exact::show(&h.fixed),
                        exact::show(&h.paired),
                    ));
                }
                result.push_str(", exact CM endpoint and neighbor counts certified");
                result
            }
        };
        format!(
            "rvFrickeCMPolyData[{name}] = degree={degree}, coefficients={coefficients}, {detail}"
        )
    }

    pub(crate) fn eval(&self, model: &mut Model, ring: &Ring, h: &[u64]) -> Result<Vec<u64>> {
        residue_field::eval(model, ring, &self.support, h)
    }

    pub(crate) fn modular_coefficient(
        &self,
        model: &mut Model,
        ring: &Ring,
        i: usize,
    ) -> Result<Vec<u64>> {
        let coefficient = super::residue(model, self.support.coefficient(i), ring.modulus())?;
        Ok(ring.scalar(coefficient))
    }

    pub(crate) fn factors(&self, model: &mut Model, name: &str, ring: &Ring) -> Result<String> {
        residue_field::factors_report(model, self, name, ring)
    }

    pub(crate) fn endpoints(
        &self,
        model: &mut Model,
        name: &str,
        ring: &Ring,
        h: &[u64],
    ) -> Result<String> {
        residue_field::edge_endpoints(model, self, name, ring, h)
    }
}

fn horizontal(model: &mut Model, norm: &Norm) -> Result<Horizontal> {
    let s = exact::integers(&ENDPOINT_SUM);
    let p = exact::integers(&ENDPOINT_PRODUCT);
    let [a, b] = &norm.remainder;
    let delta = s.multiply(model, &s)?.subtract(&exact::scale(&p, 4)?)?;
    let discriminant = exact::gcd(model, &norm.squarefree, &delta)?;
    let both = exact::gcd(model, &norm.squarefree, a)?;
    let both = exact::gcd(model, &both, b)?;
    let repeated_both = exact::gcd(model, &both, &discriminant)?;
    let distinct = both.exact_divide(model, &repeated_both)?;

    // A repeated endpoint is j=S/2. 2A+S B is twice its class residual.
    let at_double = exact::plus(&exact::scale(a, 2)?, &s.multiply(model, b)?)?;
    let loops = exact::gcd(model, &discriminant, &at_double)?;
    // t=±64 are the W2-fixed points, so H=24±128.
    let fixed_parameters =
        exact::integers(&[104, 1]).multiply(model, &exact::integers(&[-152, 1]))?;
    let fixed = exact::gcd(model, &loops, &fixed_parameters)?;
    let paired = loops.exact_divide(model, &fixed)?;
    Ok(Horizontal {
        distinct,
        fixed,
        paired,
    })
}

fn certify_endpoints(
    model: &mut Model,
    support: &Poly,
    left: &Norm,
    right: &Norm,
    d: &ClassPolynomial,
    e: &ClassPolynomial,
) -> Result<()> {
    if !exact::unit(&exact::gcd(model, &d.coefficients, &e.coefficients)?) {
        return Err(certificate());
    }
    if exact::unit(support) {
        return Ok(());
    }
    let s = exact::integers(&ENDPOINT_SUM);
    let p = exact::integers(&ENDPOINT_PRODUCT);
    let delta = s.multiply(model, &s)?.subtract(&exact::scale(&p, 4)?)?;
    let [a, b] = &left.remainder;
    let [c, d] = &right.remainder;
    for denominator in [b, d, &delta] {
        if !exact::unit(&exact::gcd(model, support, denominator)?) {
            return Err(certificate());
        }
    }
    // In Q[H]/support the two roots are -A/B and -C/D.
    // Clear their denominators in the sum and product identities.
    let denominators = b.multiply(model, d)?;
    let sum = exact::plus(&a.multiply(model, d)?, &c.multiply(model, b)?)?;
    let sum = exact::plus(&sum, &s.multiply(model, &denominators)?)?;
    let product = a
        .multiply(model, c)?
        .subtract(&p.multiply(model, &denominators)?)?;
    for residual in [sum, product] {
        residual.exact_divide(model, support)?;
    }
    Ok(())
}

fn orientations(
    model: &mut Model,
    support: &Poly,
    left: &[(usize, Poly)],
    right: &[(usize, Poly)],
) -> Result<Vec<(usize, usize, Poly)>> {
    let mut result = Vec::new();
    let mut combined = Poly::one();
    for (lm, l) in left {
        for (rm, r) in right {
            let shared = exact::gcd(model, l, r)?;
            let factor = exact::gcd(model, &shared, support)?;
            if !exact::unit(&factor) {
                combined = combined.multiply(model, &factor)?;
                result.push((*lm, *rm, factor));
            }
        }
    }
    if !support.same(&combined) {
        return Err(certificate());
    }
    Ok(result)
}

fn directed_degrees(strata: &[(usize, usize, Poly)]) -> [usize; 2] {
    let mut degrees = [0; 2];
    for (a, b, factor) in strata {
        let degree = exact::degree(factor);
        degrees[0] += a * degree;
        degrees[1] += b * degree;
    }
    degrees
}

fn expected_edges(model: &mut Model, source: &ClassPolynomial, target: i128) -> Result<usize> {
    let order = Order::new(source.discriminant)?;
    let classes = order.elements(model, source.degree())?;
    if classes.len() != source.degree() {
        return Err(certificate());
    }
    let mut total = 0;
    for c in classes {
        let counts = cm_neighbors::counts(model, c, Prime::validated(2))?;
        total += if target == source.discriminant {
            counts.horizontal
        } else if source.discriminant.checked_mul(4) == Some(target) {
            counts.descending
        } else if target.checked_mul(4) == Some(source.discriminant) {
            counts.ascending
        } else {
            0
        };
    }
    Ok(total)
}

fn certificate() -> Error {
    Error::plain("Fricke CM norm or edge failed its exact certificate")
}
