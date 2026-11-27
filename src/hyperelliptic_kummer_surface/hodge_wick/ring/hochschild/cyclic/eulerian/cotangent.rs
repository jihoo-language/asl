use std::fmt::Write;

use super::{storage, work_limit, BTreeMap, Error, Model, Result};

pub(crate) type Variable = (u128, u128);

pub(crate) enum Command {
    Dimension(i128, i128, i128, Vec<Variable>),
    Differential(Vec<Variable>, Vec<u128>, Vec<u128>, Vec<usize>),
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<String> {
        match self {
            Self::Dimension(degree, index, weight, variables) => {
                let dim = dimension(model, degree, index, weight, &variables)?;
                Ok(format!(
                    "rvMonomialCICotangentDim[{degree}, {index}, {weight}] = {dim}"
                ))
            }
            Self::Differential(variables, monomial, relations, wedge) => {
                let value = differential(model, &variables, &monomial, &relations, &wedge)?;
                Ok(format!("rvMonomialCICotangentDiff = {value}"))
            }
        }
    }
}

fn validate(variables: &[Variable]) -> Result<()> {
    if variables
        .iter()
        .any(|&(power, weight)| power == 0 || weight == 0)
    {
        Err(Error::plain(
            "monomial complete intersection needs positive relation powers and variable weights",
        ))
    } else {
        Ok(())
    }
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
struct Component {
    hodge: u128,
    forms: u128,
    weight: u128,
}

impl Component {
    const ZERO: Self = Self {
        hodge: 0,
        forms: 0,
        weight: 0,
    };

    fn subtract(self, other: Self) -> Option<Self> {
        Some(Self {
            hodge: self.hodge.checked_sub(other.hodge)?,
            forms: self.forms.checked_sub(other.forms)?,
            weight: self.weight.checked_sub(other.weight)?,
        })
    }
}

// A one-variable exterior slice contains gamma_j(xi) and
// gamma_(j-1)(xi) dx. The only nonzero differential coefficient
// multiplies the constant of the former by r*x^(r-1).
fn survivors(hodge: u128, forms: u128, (power, _): Variable) -> Option<(u128, u128, u128)> {
    match (hodge, forms) {
        (0, 0) => Some((0, 0, power - 1)),
        (_, 0) if power > 1 => Some((hodge.checked_mul(power)?, 1, power - 1)),
        (_, 1) if hodge > 0 && power > 1 => Some((
            (hodge - 1).checked_mul(power)?.checked_add(1)?,
            0,
            power - 2,
        )),
        _ => None,
    }
}

fn local_dimension(grade: Component, variable: Variable) -> u128 {
    let Some((base, min, max)) = survivors(grade.hodge, grade.forms, variable) else {
        return 0;
    };
    if !grade.weight.is_multiple_of(variable.1) {
        return 0;
    }
    let exponent = (grade.weight / variable.1).checked_sub(base);
    u128::from(exponent.is_some_and(|exponent| (min..=max).contains(&exponent)))
}

fn local_terms(model: &mut Model, cap: Component, variable: Variable) -> Result<Vec<Component>> {
    let bound = cap.weight / variable.1;
    let max_hodge = cap.hodge.min(bound / variable.0 + 1);
    let max_hodge = usize::try_from(max_hodge).map_err(|_| work_limit())?;
    storage(max_hodge.saturating_add(1))?;
    let mut components = Vec::new();
    for hodge in 0..=max_hodge as u128 {
        model.spend(1)?;
        for forms in 0..=cap.forms.min(1) {
            if let Some((base, min, max)) = survivors(hodge, forms, variable) {
                if base > bound {
                    continue;
                }
                for b in min..=max.min(bound - base) {
                    storage((components.len() + 1).saturating_mul(8))?;
                    model.spend(1)?;
                    components.push(Component {
                        hodge,
                        forms,
                        weight: (base + b) * variable.1,
                    });
                }
            }
        }
    }
    Ok(components)
}

fn dimension(
    model: &mut Model,
    degree: i128,
    index: i128,
    weight: i128,
    variables: &[Variable],
) -> Result<u128> {
    validate(variables)?;
    let (Ok(degree), Ok(hodge), Ok(weight)) = (
        u128::try_from(degree),
        u128::try_from(index),
        u128::try_from(weight),
    ) else {
        return Ok(0);
    };
    if degree < hodge || degree > 2 * hodge {
        return Ok(0);
    }
    let cap = Component {
        hodge,
        forms: 2 * hodge - degree,
        weight,
    };
    if cap.forms > variables.len() as u128 {
        return Ok(0);
    }
    let Some((&last, earlier)) = variables.split_last() else {
        return Ok(u128::from(cap == Component::ZERO));
    };
    let mut state = BTreeMap::from([(Component::ZERO, 1u128)]);
    for &variable in earlier {
        let terms = local_terms(model, cap, variable)?;
        let mut next = BTreeMap::<Component, u128>::new();
        for (sum, coefficient) in state {
            let remaining = cap.subtract(sum).unwrap();
            model.spend(terms.len())?;
            for &term in &terms {
                if remaining.subtract(term).is_none() {
                    continue;
                }
                let grade = Component {
                    hodge: sum.hodge + term.hodge,
                    forms: sum.forms + term.forms,
                    weight: sum.weight + term.weight,
                };
                let count = next.entry(grade).or_default();
                *count = count
                    .checked_add(coefficient)
                    .ok_or_else(dimension_overflow)?;
                storage(next.len().saturating_mul(10))?;
            }
        }
        state = next;
        if state.is_empty() {
            return Ok(0);
        }
    }
    state.into_iter().try_fold(0u128, |dim, (sum, count)| {
        model.spend(1)?;
        let term = if local_dimension(cap.subtract(sum).unwrap(), last) > 0 {
            count
        } else {
            0
        };
        dim.checked_add(term).ok_or_else(dimension_overflow)
    })
}

fn dimension_overflow() -> Error {
    Error::plain(
        "monomial complete-intersection cotangent dimension exceeds the supported integer range",
    )
}

fn differential(
    model: &mut Model,
    variables: &[Variable],
    monomial: &[u128],
    relations: &[u128],
    wedge: &[usize],
) -> Result<String> {
    validate(variables)?;
    if monomial.len() != variables.len() || relations.len() != variables.len() {
        return Err(Error::plain(
            "cotangent term exponents have the wrong number of variables",
        ));
    }
    if monomial
        .iter()
        .zip(variables)
        .any(|(&power, relation)| power >= relation.0)
    {
        return Err(Error::plain(
            "cotangent coefficient monomial exceeds a relation power",
        ));
    }
    if wedge.iter().any(|&i| i >= variables.len())
        || !wedge.windows(2).all(|pair| pair[0] < pair[1])
    {
        return Err(Error::plain(
            "cotangent wedge indices must be distinct, in range, and increasing",
        ));
    }
    let mut result = "[".to_owned();
    let mut separator = "";
    for (i, (&exponent, &divided)) in monomial.iter().zip(relations).enumerate() {
        let Err(place) = wedge.binary_search(&i) else {
            continue;
        };
        if exponent != 0 || divided == 0 {
            continue;
        }
        model.spend(1 + variables.len() + wedge.len())?;
        let mut powers = monomial.to_vec();
        powers[i] = variables[i].0 - 1;
        let mut divided_powers = relations.to_vec();
        divided_powers[i] -= 1;
        let mut forms = wedge.to_vec();
        forms.insert(place, i);
        let sign = if place.is_multiple_of(2) { "" } else { "-" };
        let coefficient = variables[i].0;
        write!(
            result,
            "{separator}({sign}{coefficient}, {powers:?}, {divided_powers:?}, {forms:?})"
        )
        .unwrap();
        if result.len() >= crate::model::MAX_SHOWN_BYTES {
            return Err(Error::plain(
                "monomial complete-intersection cotangent report exceeds the output limit",
            ));
        }
        separator = ", ";
    }
    result.push(']');
    Ok(result)
}
