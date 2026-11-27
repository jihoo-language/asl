use crate::ec_poly::{self as poly, Field, Scalar};
use crate::error::{Error, Result};
use crate::finite_field_factor::distinct_degrees;
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;
use crate::rv_qp_quad::Wide;

use super::super::super::super::super::k3;
use super::{
    closed_places, coefficient, cube, field, local_uniformizer, order, BaseRing, Kodaira,
    MinimalReduction, Parameter, Pencil, RelativeJacobian, DISCRIMINANT_DEGREE,
};

mod bsd;

type Uniformizer = [Scalar; 2];

struct Components {
    cycles: Vec<usize>,
    tamagawa: usize,
    split: Option<bool>,
}

struct BadPlace {
    uniformizer: Option<Uniformizer>,
    kind: Kodaira,
    action: Components,
}

struct LFunction {
    q: u64,
    polynomial: Vec<Integer>,
    conductor: usize,
    ns_rank: usize,
    rank: usize,
    trivial_rank: usize,
    node_orbits: Vec<usize>,
    places: Vec<BadPlace>,
}

impl Kodaira {
    fn conductor_order(self) -> usize {
        match self {
            Self::Good => 0,
            Self::Multiplicative(_) => 1,
            Self::Star(_) | Self::Other(_) => 2,
        }
    }
}

impl RelativeJacobian {
    pub(crate) fn euler(
        &self,
        model: &mut Model,
        pencil: &Pencil,
        name: &str,
        parameter: &Parameter,
    ) -> Result<String> {
        let field = field(&pencil.curve);
        require_tame(field)?;
        let q = field.ring.size_of_residue_field()?;
        let parameter = local_uniformizer(model, field, parameter)?;
        let parameter = parameter.as_ref();
        let reduction = self.minimal_at(model, field, parameter)?;
        let components = self.components(model, field, q, parameter, &reduction)?;
        let euler = match (reduction.kodaira, &components) {
            (Kodaira::Good, _) => {
                let (a, b) = self.short_coefficients(model, field, parameter, &reduction, false)?;
                let cubic = [b, a, field.zero(), field.one()];
                let trace = field.ring.residue_quadratic_trace(model, |model, x| {
                    poly::evaluate(model, field, &cubic, x)
                })?;
                vec![1, -i128::from(trace), i128::from(q)]
            }
            (Kodaira::Multiplicative(_), Some(components)) => {
                let split = components.split.ok_or_else(inconsistent_components)?;
                vec![1, if split { -1 } else { 1 }]
            }
            _ => vec![1],
        };
        let symbol = reduction.kodaira.symbol();
        let conductor = reduction.kodaira.conductor_order();
        let action = if let Some(components) = components {
            let kind = match components.split {
                Some(true) => ", split=yes",
                Some(false) => ", split=no",
                None => "",
            };
            format!(
                "{kind}, component cycles={:?}, Tamagawa={}",
                components.cycles, components.tamagawa
            )
        } else {
            ", component action=uncomputed".into()
        };
        Ok(format!(
            "rvK3FiberEuler[{name}] = q={q}, geometric Kodaira={symbol}, factor={euler:?}, conductor order={conductor}{action}"
        ))
    }

    pub(crate) fn l_function(
        &self,
        model: &mut Model,
        pencil: &Pencil,
        name: &str,
    ) -> Result<String> {
        let LFunction {
            q,
            polynomial,
            conductor,
            ns_rank,
            rank,
            trivial_rank,
            ..
        } = self.l_data(model, pencil)?;
        let degree = polynomial.len() - 1;
        let coefficients = polynomial
            .into_iter()
            .map(|coefficient| coefficient.signed().ok_or_else(polynomial_overflow))
            .collect::<Result<Vec<_>>>()?;
        Ok(format!(
            "rvK3PencilLFunction[{name}] = q={q}, coefficients={coefficients:?}, degree={degree}, conductor degree={conductor}, NS rank={ns_rank}, trivial rank={trivial_rank}, Mordell-Weil rank={rank}"
        ))
    }

    pub(super) fn arithmetic_mordell_weil_rank(
        &self,
        model: &mut Model,
        pencil: &Pencil,
    ) -> Result<usize> {
        Ok(self.l_data(model, pencil)?.rank)
    }

    fn l_data(&self, model: &mut Model, pencil: &Pencil) -> Result<LFunction> {
        let field = field(&pencil.curve);
        require_tame(field)?;
        self.require_nonconstant_j(model, field)?;
        let q = field.ring.size_of_residue_field()?;
        let finite = closed_places(model, field, &self.0.discriminant)?;
        let mut places = Vec::with_capacity(finite.len() + 1);
        if self.0.discriminant.len() <= DISCRIMINANT_DEGREE {
            places.push(None);
        }
        for (degree, _, polynomial) in finite {
            if degree != 1 {
                return Err(Error::plain(
                    "K3 pencil L function currently requires rational discriminant support",
                ));
            }
            let factor: Uniformizer = polynomial
                .try_into()
                .map_err(|_| inconsistent_components())?;
            places.push(Some(factor));
        }

        let mut local = Vec::with_capacity(places.len());
        let mut conductor = 0;
        for uniformizer in places {
            let place = uniformizer.as_ref();
            let reduction = self.minimal_at(model, field, place)?;
            let action = self
                .components(model, field, q, place, &reduction)?
                .ok_or_else(|| {
                    Error::plain(format!(
                    "K3 pencil L function has no component action for geometric Kodaira type {}",
                    reduction.kodaira.symbol()
                ))
                })?;
            local.push(BadPlace {
                uniformizer,
                kind: reduction.kodaira,
                action,
            });
            conductor += reduction.kodaira.conductor_order();
        }

        // The two abelian eigenvalues q and the branch-induced exceptional
        // orbits complete the actual Kummer H² polynomial det(1 - T F_q)
        // in the point-count convention.
        let cross = k3::cross_quartic(q, pencil.curve.weil(model)?)?;
        let mut polynomial: Vec<_> = cross.into_iter().map(Integer::from_i128).collect();
        let node_orbits = k3::node_orbits(model, &pencil.curve, q)?;
        for length in [1, 1].into_iter().chain(node_orbits.iter().copied()) {
            polynomial = multiply_cycle(model, &polynomial, q, length)?;
        }
        if polynomial.len() != 23 {
            return Err(inconsistent_components());
        }
        let ns_rank = root_multiplicity(model, &polynomial, q)?;
        let component_cycles = local
            .iter()
            .flat_map(|place| place.action.cycles.iter().copied());
        let trivial_rank = 2 + component_cycles.clone().count();
        for length in [1, 1].into_iter().chain(component_cycles) {
            polynomial = divide_cycle(model, &polynomial, q, length)?.ok_or_else(|| {
                Error::plain("K3 pencil component factor did not divide H² exactly")
            })?;
        }
        let rank = root_multiplicity(model, &polynomial, q)?;
        let degree = polynomial.len() - 1;
        if ns_rank != rank + trivial_rank || conductor.checked_sub(4) != Some(degree) {
            return Err(Error::plain(
                "K3 pencil L function failed its conductor or rank check",
            ));
        }
        Ok(LFunction {
            q,
            polynomial,
            conductor,
            ns_rank,
            rank,
            trivial_rank,
            node_orbits,
            places: local,
        })
    }

    fn require_nonconstant_j(&self, model: &mut Model, field: Field<'_>) -> Result<()> {
        let c4_cube = cube(model, &BaseRing(field), &self.0.c4)?;
        let discriminant = &self.0.discriminant;
        let constant = if c4_cube.is_empty() {
            true
        } else if c4_cube.len() != discriminant.len() {
            false
        } else {
            poly::scale(model, field, &c4_cube, discriminant.last().unwrap())?
                == poly::scale(model, field, discriminant, c4_cube.last().unwrap())?
        };
        if constant {
            return Err(Error::plain("K3 pencil L function requires nonconstant j"));
        }
        Ok(())
    }

    fn minimal_at(
        &self,
        model: &mut Model,
        field: Field<'_>,
        place: Option<&Uniformizer>,
    ) -> Result<MinimalReduction> {
        let c4 = order(model, field, &self.0.c4, 8, place)?;
        let c6 = order(model, field, &self.0.c6, 12, place)?;
        let discriminant = order(
            model,
            field,
            &self.0.discriminant,
            DISCRIMINANT_DEGREE,
            place,
        )?;
        MinimalReduction::new(c4, c6, discriminant)
    }

    fn short_coefficients(
        &self,
        model: &mut Model,
        field: Field<'_>,
        place: Option<&Uniformizer>,
        reduction: &MinimalReduction,
        residual_star: bool,
    ) -> Result<(Scalar, Scalar)> {
        let scale = reduction.scales;
        let star = usize::from(residual_star);
        let c4 = leading_residue(model, field, &self.0.c4, 8, place, 4 * scale + 2 * star)?;
        let c6 = leading_residue(model, field, &self.0.c6, 12, place, 6 * scale + 3 * star)?;
        let inverse48 = field.inverse_unit(model, &field.signed(-48))?;
        let inverse864 = field.inverse_unit(model, &field.signed(-864))?;
        let a = field.multiply(model, &c4, &inverse48)?;
        let b = field.multiply(model, &c6, &inverse864)?;
        Ok((a, b))
    }

    fn components(
        &self,
        model: &mut Model,
        field: Field<'_>,
        q: u64,
        place: Option<&Uniformizer>,
        reduction: &MinimalReduction,
    ) -> Result<Option<Components>> {
        let components = match reduction.kodaira {
            Kodaira::Good => Components {
                cycles: Vec::new(),
                tamagawa: 1,
                split: None,
            },
            Kodaira::Multiplicative(n) => {
                let c6 =
                    leading_residue(model, field, &self.0.c6, 12, place, 6 * reduction.scales)?;
                if field.is_zero(&c6) {
                    return Err(inconsistent_components());
                }
                let tangent = field.ring.negate(model, &c6)?;
                let split = field.ring.power(model, &tangent, (q - 1) / 2)? == field.one();
                let tamagawa = if split {
                    n
                } else {
                    1 + usize::from(n.is_multiple_of(2))
                };
                let mut cycles = vec![1; tamagawa - 1];
                cycles.extend(std::iter::repeat_n(2, (n - tamagawa) / 2));
                Components {
                    cycles,
                    tamagawa,
                    split: Some(split),
                }
            }
            Kodaira::Star(0) => {
                let (a, b) = self.short_coefficients(model, field, place, reduction, true)?;
                let cubic = vec![b, a, field.zero(), field.one()];
                let derivative = poly::derivative(model, field, &cubic)?;
                if !poly::coprime(model, field, cubic.clone(), derivative)? {
                    return Err(inconsistent_components());
                }
                // The central component is fixed. The three other
                // nonidentity components are the roots of the cubic.
                let mut cycles = vec![1];
                let mut tamagawa = 1;
                distinct_degrees(model, field, Some(q), cubic, |_, block, degree| {
                    let count = (block.len() - 1) / degree;
                    cycles.extend(std::iter::repeat_n(degree, count));
                    if degree == 1 {
                        tamagawa += count;
                    }
                    Ok(())
                })?;
                if cycles.iter().sum::<usize>() != 4 || !matches!(tamagawa, 1 | 2 | 4) {
                    return Err(inconsistent_components());
                }
                Components {
                    cycles,
                    tamagawa,
                    split: None,
                }
            }
            Kodaira::Star(_) | Kodaira::Other(_) => return Ok(None),
        };
        Ok(Some(components))
    }
}

fn require_tame(field: Field<'_>) -> Result<()> {
    if field.ring.prime() < 5 {
        return Err(Error::plain(
            "K3 pencil arithmetic Euler and L reports require characteristic at least five",
        ));
    }
    Ok(())
}

fn leading_residue(
    model: &mut Model,
    field: Field<'_>,
    polynomial: &[Scalar],
    degree: usize,
    place: Option<&Uniformizer>,
    power: usize,
) -> Result<Scalar> {
    if polynomial.is_empty() {
        return Ok(field.zero());
    }
    if let Some(linear) = place {
        let mut quotient = polynomial.to_vec();
        for _ in 0..power {
            quotient = poly::exact_divide(model, field, &quotient, linear)?;
        }
        let root = field.ring.negate(model, &linear[0])?;
        poly::evaluate(model, field, &quotient, &root)
    } else {
        let index = degree
            .checked_sub(power)
            .ok_or_else(inconsistent_components)?;
        Ok(coefficient(field, polynomial, index))
    }
}

// The degree-22 numerator can contain coefficients of size q^22 before
// cancellation. Keep the exact P2 and all successive quotients wide.
const LIMBS: usize = 48;
type Integer = Wide<LIMBS>;

fn cycle_power(q: u64, length: usize) -> Result<Integer> {
    if length == 0 || length > 22 {
        return Err(inconsistent_components());
    }
    (0..length).try_fold(Integer::from_i128(1), |power, _| {
        power
            .multiply(i128::from(q))
            .map_err(|_| polynomial_overflow())
    })
}

fn multiply_cycle(model: &mut Model, p: &[Integer], q: u64, length: usize) -> Result<Vec<Integer>> {
    let power = cycle_power(q, length)?;
    model.spend(p.len() * (LIMBS + LIMBS * LIMBS))?;
    let mut product = p.to_vec();
    product.resize(p.len() + length, Integer::from_i128(0));
    for (i, &coefficient) in p.iter().enumerate() {
        let term = coefficient
            .multiply_wide(power)
            .map_err(|_| polynomial_overflow())?;
        product[i + length] = product[i + length]
            .add(term.negate())
            .map_err(|_| polynomial_overflow())?;
    }
    Ok(product)
}

fn divide_cycle(
    model: &mut Model,
    p: &[Integer],
    q: u64,
    length: usize,
) -> Result<Option<Vec<Integer>>> {
    if length >= p.len() {
        return Ok(None);
    }
    let power = cycle_power(q, length)?;
    model.spend(p.len() * (LIMBS + LIMBS * LIMBS))?;
    let quotient_len = p.len() - length;
    let mut quotient: Vec<Integer> = Vec::with_capacity(quotient_len);
    for (i, &coefficient) in p.iter().enumerate() {
        let mut current = coefficient;
        if let Some(previous) = i.checked_sub(length).and_then(|j| quotient.get(j)) {
            let term = previous
                .multiply_wide(power)
                .map_err(|_| polynomial_overflow())?;
            current = current.add(term).map_err(|_| polynomial_overflow())?;
        }
        if i < quotient_len {
            quotient.push(current);
        } else if !current.is_zero() {
            return Ok(None);
        }
    }
    Ok(Some(quotient))
}

fn root_multiplicity(model: &mut Model, polynomial: &[Integer], q: u64) -> Result<usize> {
    let mut polynomial = polynomial.to_vec();
    let mut rank = 0;
    while let Some(quotient) = divide_cycle(model, &polynomial, q, 1)? {
        rank += 1;
        polynomial = quotient;
    }
    Ok(rank)
}

fn inconsistent_components() -> Error {
    Error::plain("K3 pencil has inconsistent arithmetic component data")
}

fn polynomial_overflow() -> Error {
    Error::plain("K3 pencil L coefficient exceeds the exact integer range")
}

#[cfg(test)]
mod tests {
    use super::super::super::{invariants, Jacobian};
    use super::*;
    use crate::rv_unram::Ring;

    #[test]
    fn residual_cubic_distinguishes_all_star_component_actions() {
        let mut model = Model::new();
        let ring =
            Ring::new(&mut model, 5, 1, vec![0, 1]).unwrap_or_else(|error| panic!("{error}"));
        let field = Field { ring: &ring };
        let place = [field.zero(), field.one()];
        for (a, b, expected, tamagawa) in [
            (-1, 0, vec![1, 1, 1, 1], 4),
            (2, 0, vec![1, 1, 2], 2),
            (1, 1, vec![1, 3], 1),
        ] {
            let a2 = Vec::new();
            let a4 = field.polynomial(&[0, 0, a]);
            let a6 = field.polynomial(&[0, 0, 0, b]);
            let (c4, c6, discriminant) = invariants(&mut model, &BaseRing(field), &a2, &a4, &a6)
                .unwrap_or_else(|error| panic!("{error}"));
            let relative = RelativeJacobian(Jacobian {
                a2,
                a4,
                a6,
                c4,
                c6,
                discriminant,
            });
            let reduction = relative
                .minimal_at(&mut model, field, Some(&place))
                .unwrap_or_else(|error| panic!("{error}"));
            assert!(matches!(reduction.kodaira, Kodaira::Star(0)));
            let components = relative
                .components(&mut model, field, 5, Some(&place), &reduction)
                .unwrap_or_else(|error| panic!("{error}"))
                .unwrap();
            assert_eq!(components.cycles, expected);
            assert_eq!(components.tamagawa, tamagawa);
        }
    }
}
