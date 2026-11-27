use crate::error::Result;
use crate::model::{Model, Value};
use crate::rv_ore::{self, arithmetic, Powers, Support};
use crate::rv_qconnection::{self, PeriodScale, Polynomial, Rational};
use crate::rv_qgr::Grade;

type Blocks = Vec<(usize, Polynomial)>;

fn times_linear(model: &mut Model, polynomial: &mut Polynomial, offset: i128) -> Result<bool> {
    model.spend(polynomial.len())?;
    if let Some(&leading) = polynomial.last() {
        rv_qconnection::storage(polynomial.len() + 1)?;
        let length = polynomial.len();
        polynomial.push(leading);
        for i in (1..length).rev() {
            let Some(coefficient) = polynomial[i]
                .checked_mul(offset)
                .and_then(|value| value.checked_add(polynomial[i - 1]))
            else {
                return Ok(false);
            };
            polynomial[i] = coefficient;
        }
        let Some(constant) = polynomial[0].checked_mul(offset) else {
            return Ok(false);
        };
        polynomial[0] = constant;
    }
    Ok(true)
}

pub(crate) fn factorial(model: &mut Model, height: usize, input: &Support) -> Result<Value> {
    if height == 0 {
        return Ok(Value::RVOre(input.clone()));
    }
    let mut output = rv_ore::Sum::default();
    let mut polynomial = vec![1];
    let mut previous = 0;
    for &(powers, coefficient) in input.iter() {
        let degree = powers.q.checked_mul(height).ok_or_else(arithmetic)?;
        rv_qconnection::storage(degree.checked_add(1).ok_or_else(arithmetic)?)?;
        powers.theta.checked_add(degree).ok_or_else(arithmetic)?;
        for q in previous + 1..=powers.q {
            for _ in 0..height {
                if !times_linear(model, &mut polynomial, q as i128)? {
                    return Err(arithmetic());
                }
            }
        }
        previous = powers.q;
        for (i, &factor) in polynomial.iter().enumerate() {
            output.add(
                model,
                Powers {
                    q: powers.q,
                    theta: powers.theta + i,
                },
                coefficient.checked_mul(factor).ok_or_else(arithmetic)?,
            )?;
        }
    }
    output.finish(model)
}

fn grouped(model: &mut Model, input: &Support) -> Result<Blocks> {
    let mut blocks = Vec::new();
    let mut cells: usize = 0;
    for piece in input.chunk_by(|left, right| left.0.q == right.0.q) {
        let last = piece.last().unwrap().0;
        let length = last.theta.checked_add(1).ok_or_else(arithmetic)?;
        cells = cells.checked_add(length).ok_or_else(arithmetic)?;
        rv_qconnection::storage(cells)?;
        model.spend(length)?;
        let mut polynomial = vec![0; length];
        for &(powers, coefficient) in piece {
            polynomial[powers.theta] = coefficient;
        }
        blocks.push((last.q, polynomial));
    }
    Ok(blocks)
}

fn shifted(model: &mut Model, polynomial: &[i128], offset: i128) -> Result<Option<Polynomial>> {
    if offset == 0 {
        model.spend(polynomial.len())?;
        return Ok(Some(polynomial.to_vec()));
    }
    let mut result = vec![*polynomial.last().unwrap()];
    for &coefficient in polynomial[..polynomial.len() - 1].iter().rev() {
        if !times_linear(model, &mut result, offset)? {
            return Ok(None);
        }
        let Some(constant) = result[0].checked_add(coefficient) else {
            return Ok(None);
        };
        result[0] = constant;
    }
    Ok(Some(result))
}

fn common_factor(
    model: &mut Model,
    blocks: &[(usize, Polynomial)],
    origin: usize,
) -> Result<Polynomial> {
    if blocks.is_empty() || blocks.iter().any(|(_, p)| p.len() == 1) {
        return Ok(vec![1]);
    }
    let mut factor = Vec::new();
    for (power, polynomial) in blocks {
        model.spend(polynomial.len())?;
        let mut primitive = polynomial.clone();
        rv_qconnection::primitive(&mut primitive);
        let offset = origin as i128 - *power as i128;
        let Some(shifted) = shifted(model, &primitive, offset)? else {
            return rv_qconnection::translated_gcd(model, blocks, origin);
        };
        factor = rv_qconnection::polynomial_gcd(model, factor, shifted)?;
        if factor == [1] {
            break;
        }
    }
    Ok(factor)
}

fn divide(
    model: &mut Model,
    blocks: &mut [(usize, Polynomial)],
    factor: &[i128],
    origin: usize,
) -> Result<()> {
    if factor == [1] {
        return Ok(());
    }
    for (power, polynomial) in blocks {
        let offset = *power as i128 - origin as i128;
        *polynomial = if let Some(divisor) = shifted(model, factor, offset)? {
            rv_qconnection::exact_divide(model, std::mem::take(polynomial), &divisor)?
        } else {
            rv_qconnection::translated_quotient(model, polynomial, factor, offset)?
        };
    }
    Ok(())
}

fn ore(model: &mut Model, coefficients: Blocks) -> Result<Value> {
    let mut output = rv_ore::Sum::default();
    for (q, polynomial) in coefficients {
        for (theta, coefficient) in polynomial.into_iter().enumerate() {
            output.add(model, Powers { q, theta }, coefficient)?;
        }
    }
    output.finish(model)
}

pub(crate) fn euler_factor(model: &mut Model, input: &Support) -> Result<Value> {
    let blocks = grouped(model, input)?;
    let factor = common_factor(model, &blocks, 0)?;
    ore(model, vec![(0, factor)])
}

pub(crate) fn euler_reduce(model: &mut Model, input: &Support) -> Result<Value> {
    let mut blocks = grouped(model, input)?;
    if let [(_, polynomial)] = blocks.as_mut_slice() {
        // Translation preserves primitive content. The entire
        // normalized block is the divisor, leaving just its content.
        model.spend(polynomial.len())?;
        let content = rv_qconnection::content(polynomial);
        let negative = *polynomial.last().unwrap() < 0;
        let coefficient = if let Ok(content) = i128::try_from(content) {
            if negative {
                -content
            } else {
                content
            }
        } else {
            // Content 2^127 occurs only when all nonzero coefficients
            // equal i128::MIN, so the leading coefficient is negative.
            debug_assert!(negative);
            i128::MIN
        };
        polynomial.clear();
        polynomial.push(coefficient);
    } else {
        let origin = blocks
            .iter()
            .min_by_key(|(_, p)| p.len())
            .map_or(0, |(power, _)| *power);
        let factor = common_factor(model, &blocks, origin)?;
        divide(model, &mut blocks, &factor, origin)?;
    }
    ore(model, blocks)
}

fn root_at_one(polynomial: &[i128]) -> bool {
    // Keep the running sum in range by always taking the next opposite sign.
    let mut positive = polynomial.iter().copied().filter(|&c| c > 0);
    let mut negative = polynomial.iter().copied().filter(|&c| c < 0);
    let mut sum = 0i128;
    loop {
        let next = if sum > 0 {
            negative.next()
        } else {
            positive.next()
        };
        let Some(next) = next else {
            return sum == 0 && negative.next().is_none();
        };
        sum += next;
    }
}

fn positive_root(model: &mut Model, polynomial: &[i128], root: u128) -> Result<bool> {
    model.spend(polynomial.len() * 2)?;
    if root == 1 {
        return Ok(root_at_one(polynomial));
    }
    // Verify from the constant term using exact divisions by the candidate.
    // For root >= 2 all intermediate magnitudes are strictly below 2^128.
    let mut magnitude = polynomial[0].unsigned_abs();
    let mut negative = polynomial[0] < 0;
    for &coefficient in &polynomial[1..] {
        if !magnitude.is_multiple_of(root) {
            return Ok(false);
        }
        magnitude /= root;
        let next = coefficient.unsigned_abs();
        if negative == (coefficient < 0) {
            magnitude = magnitude.checked_add(next).ok_or_else(arithmetic)?;
        } else if magnitude < next {
            magnitude = next - magnitude;
            negative = !negative;
        } else {
            magnitude -= next;
        }
    }
    Ok(magnitude == 0)
}

fn roots(model: &mut Model, polynomial: &[i128]) -> Result<Vec<usize>> {
    model.spend(polynomial.len())?;
    let mut roots = Vec::new();
    let nonzero = polynomial.iter().position(|&c| c != 0).unwrap();
    if nonzero != 0 {
        roots.push(0);
    }
    let polynomial = &polynomial[nonzero..];
    if polynomial.len() == 1 {
        return Ok(roots);
    }
    let (leading, lower) = polynomial.split_last().unwrap();
    let bound = 1 + lower
        .iter()
        .map(|c| c.unsigned_abs() / leading.unsigned_abs())
        .max()
        .unwrap();
    let constant = polynomial[0].unsigned_abs();
    let trial = constant.isqrt();
    let mut check = |model: &mut Model, candidate| -> Result<()> {
        if positive_root(model, polynomial, candidate)? {
            let degree = usize::try_from(candidate).map_err(|_| arithmetic())?;
            rv_qconnection::check_period_degree(degree)?;
            roots.push(degree);
        }
        Ok(())
    };
    // Cauchy's bound and the rational-root theorem exhaust the candidates.
    for candidate in 1..=bound.min(trial) {
        model.spend(1)?;
        if constant.is_multiple_of(candidate) {
            check(model, candidate)?;
            let paired = constant / candidate;
            if bound > trial && paired <= bound && paired != candidate {
                check(model, paired)?;
            }
        }
    }
    roots.sort_unstable();
    Ok(roots)
}

fn evaluate(model: &mut Model, polynomial: &[i128], x: usize) -> Result<i128> {
    model.spend(polynomial.len())?;
    polynomial.iter().rev().try_fold(0i128, |previous, &c| {
        previous
            .checked_mul(x as i128)
            .and_then(|n| n.checked_add(c))
            .ok_or_else(arithmetic)
    })
}

pub(crate) fn rhs(model: &mut Model, grade: Grade) -> Result<String> {
    if grade.rank.min(grade.width()) == 1 {
        // The factorial transform is theta^n (1-q), and G=1/(1-q).
        return Ok("{0: 1}".into());
    }
    let Value::RVOre(operator) = rv_qconnection::annihilator(model, grade)? else {
        unreachable!()
    };
    let Value::RVOre(transformed) = factorial(model, grade.ambient, &operator)? else {
        unreachable!()
    };
    let mut blocks = grouped(model, &transformed)?;
    let factor = common_factor(model, &blocks, 0)?;
    let degrees = roots(model, &factor)?;
    let Some(&last) = degrees.last() else {
        return Ok("{}".into());
    };
    divide(model, &mut blocks, &factor, 0)?;
    let period = rv_qconnection::period_coefficients(model, grade, last, PeriodScale::Factorial)?;
    let mut forcing = Vec::new();
    for degree in degrees {
        let mut value = Rational::ZERO;
        for (q, polynomial) in blocks.iter().take_while(|(q, _)| *q <= degree) {
            let input = degree - q;
            let factor = evaluate(model, polynomial, input)?;
            value = value.add(period[input].scale(factor, 1)?)?;
        }
        if !value.is_zero() {
            forcing.push(format!("{degree}: {value}"));
        }
    }
    Ok(format!("{{{}}}", forcing.join(", ")))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn success<T>(value: Result<T>) -> T {
        value.unwrap_or_else(|error| panic!("{error}"))
    }

    #[test]
    fn exact_integer_roots_include_positive_and_zero_resonances() {
        let mut model = Model::new();
        // t^2 (t-3)^2 (t+2)
        assert_eq!(success(roots(&mut model, &[0, 0, 18, -3, -4, 1])), [0, 3]);
        assert_eq!(success(roots(&mut model, &[-77, 4, 1])), [7]);
        assert_eq!(success(roots(&mut model, &[1, 1])), []);
        assert!(roots(&mut model, &[-5000, 1])
            .unwrap_err()
            .message
            .contains("limit of 4096"));
        assert!(success(positive_root(
            &mut model,
            &[i128::MAX, 1, i128::MIN],
            1
        )));
        assert!(!success(positive_root(
            &mut model,
            &[i128::MAX, 1, 1, i128::MIN],
            1
        )));
        assert!(success(positive_root(
            &mut model,
            &[i128::MIN, 1],
            1u128 << 127
        )));
    }
}
