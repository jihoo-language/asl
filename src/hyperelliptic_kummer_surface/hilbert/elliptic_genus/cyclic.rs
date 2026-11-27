use std::collections::{btree_map::Entry, BTreeMap};
use std::mem::size_of;

use crate::error::{Error, Result};
use crate::model::Model;

use super::{
    add, add_theta, multiply, multiply_factor, overflow, storage, Cutoff, Integer, Product, Seed,
    INTEGER_WORDS, MAX_PIVOT_CELLS,
};

const REAL_MAX_DEGREE: usize = 3;
type Algebra = [Integer; REAL_MAX_DEGREE];
type AlgebraSeries = BTreeMap<(usize, i128), Algebra>;

const ALGEBRA_WORDS: usize = size_of::<Algebra>().div_ceil(size_of::<usize>());
const ALGEBRA_TERM_LIMIT: usize = MAX_PIVOT_CELLS / (ALGEBRA_WORDS + 8);

pub(super) fn series(model: &mut Model, cutoff: Cutoff, order: usize) -> Result<Product> {
    let mut factors: Vec<(usize, Seed)> = Vec::new();
    for d in (1..=order.min(cutoff.p)).filter(|d| order.is_multiple_of(*d)) {
        let degree = (cutoff.p / d) * (cutoff.q / d);
        // T_d uses the primitive trace of g^d. Since d divides the order,
        // this is the trace for the exact order h/d.
        let mut cycles = if d == order {
            let mut seed = Seed::new();
            add_theta(model, &mut seed, degree, false)?;
            add_theta(model, &mut seed, degree, true)?;
            seed
        } else {
            FixedPointOrbit::new(order / d).seed(model, degree)?
        };
        // T_d = Σ_{s|d} s k_s. Invert this divisor sum over the integers.
        for (s, seed) in &factors {
            if d.is_multiple_of(*s) {
                for (&key, &value) in seed.range(..=(degree, i128::MAX)) {
                    let term = multiply(model, value, Integer::from_i128(-(*s as i128)))?;
                    add(model, &mut cycles, key, term)?;
                }
            }
        }
        for value in cycles.values_mut() {
            *value = divide(
                model,
                *value,
                d as u128,
                "cyclic DMVV multiplicity is not integral",
            )?;
        }
        factors.push((d, cycles));
    }

    let mut product = Product::from([((0, 0, 0), Integer::from_i128(1))]);
    for (d, seed) in factors {
        for m in 1..=cutoff.p / d {
            for j in 0..=cutoff.q / d {
                for (&(_, r), &count) in seed.range((m * j, i128::MIN)..=(m * j, i128::MAX)) {
                    multiply_factor(
                        model,
                        &mut product,
                        cutoff,
                        (d * m, d * j, d as i128 * r),
                        count.negate(),
                        false,
                    )?;
                }
            }
        }
    }
    Ok(product)
}

struct FixedPointOrbit {
    minimal: &'static [i128],
    trace_weights: [i128; REAL_MAX_DEGREE],
    denominator: u128,
}

impl FixedPointOrbit {
    fn new(order: usize) -> Self {
        // x = ζ + ζ⁻¹. The geometric fixed points have the indicated
        // multiplicity, spread uniformly over the real conjugates.
        let (minimal, fixed): (&[i128], i128) = match order {
            2 => (&[2, 1], 8),
            3 => (&[1, 1], 6),
            4 => (&[0, 1], 4),
            5 => (&[-1, 1, 1], 4),
            6 => (&[-1, 1], 2),
            7 => (&[-1, -2, 1, 1], 3),
            8 => (&[-2, 0, 1], 2),
            _ => unreachable!("fixed-point orbit requires a symplectic order from 2 to 8"),
        };
        let degree = minimal.len() - 1;
        let mut traces = [0i128; 2 * REAL_MAX_DEGREE - 1];
        traces[0] = degree as i128;
        // Newton power sums give Tr(x^k) without adjoining any complex roots.
        for k in 1..=2 * (degree - 1) {
            traces[k] = -(1..=k.min(degree))
                .map(|j| {
                    let factor = if j == k { k as i128 } else { traces[k - j] };
                    minimal[degree - j] * factor
                })
                .sum::<i128>();
        }
        // (2-x) R(x) = f(2) in the real cyclotomic algebra.
        let mut quotient = [0i128; REAL_MAX_DEGREE];
        let mut previous = 0;
        for j in (1..=degree).rev() {
            previous = 2 * previous + minimal[j];
            quotient[j - 1] = previous;
        }
        let denominator = (2 * previous + minimal[0]) as u128;
        let trace_weights = std::array::from_fn(|i| {
            if i >= degree {
                return 0;
            }
            let trace: i128 = (0..degree).map(|j| quotient[j] * traces[i + j]).sum();
            trace * (fixed / degree as i128)
        });
        Self {
            minimal,
            trace_weights,
            denominator,
        }
    }

    fn x(&self) -> Algebra {
        if self.minimal.len() == 2 {
            scalar(-self.minimal[0])
        } else {
            [
                Integer::from_i128(0),
                Integer::from_i128(1),
                Integer::from_i128(0),
            ]
        }
    }

    fn multiply(&self, model: &mut Model, left: Algebra, right: Algebra) -> Result<Algebra> {
        let degree = self.minimal.len() - 1;
        let mut powers = [Integer::from_i128(0); 2 * REAL_MAX_DEGREE - 1];
        for (i, a) in left.into_iter().take(degree).enumerate() {
            for (j, b) in right.into_iter().take(degree).enumerate() {
                model.spend(1)?;
                if !a.is_zero() && !b.is_zero() {
                    let term = multiply(model, a, b)?;
                    powers[i + j] = powers[i + j].add(term).map_err(|_| overflow())?;
                }
            }
        }
        for power in (degree..=2 * (degree - 1)).rev() {
            let value = powers[power];
            if !value.is_zero() {
                for (i, &coefficient) in self.minimal[..degree].iter().enumerate() {
                    let term = multiply(model, value, Integer::from_i128(-coefficient))?;
                    let target = &mut powers[power - degree + i];
                    *target = target.add(term).map_err(|_| overflow())?;
                }
            }
            powers[power] = Integer::from_i128(0);
        }
        Ok(std::array::from_fn(|i| powers[i]))
    }

    fn trace(&self, model: &mut Model, value: Algebra) -> Result<Integer> {
        let mut trace = Integer::from_i128(0);
        for (coefficient, weight) in value.into_iter().zip(self.trace_weights) {
            let term = multiply(model, coefficient, Integer::from_i128(weight))?;
            trace = trace.add(term).map_err(|_| overflow())?;
        }
        divide(
            model,
            trace,
            self.denominator,
            "local K3 Lefschetz trace is not integral",
        )
    }

    fn seed(&self, model: &mut Model, degree: usize) -> Result<Seed> {
        // The oscillator slices have at most five Laurent exponents each,
        // and the denominator recurrence holds one scalar algebra element.
        if (degree + 1).saturating_mul(6 * ALGEBRA_WORDS) > MAX_PIVOT_CELLS {
            return Err(storage());
        }
        let factors = self.oscillator(model, degree)?;
        let x = self.x();
        // The origin factor is y+y⁻¹−x. Its denominator 2−x is incorporated
        // by the trace weights after the whole local product.
        let mut local = AlgebraSeries::from([((0, -1), scalar(1)), ((0, 1), scalar(1))]);
        add_algebra(model, &mut local, (0, 0), x.map(Integer::negate))?;
        for m in 1..=degree {
            model.spend(local.len())?;
            let inputs = local
                .range(..=(degree - m, i128::MAX))
                .map(|(&key, &coefficient)| (key, coefficient))
                .collect::<Vec<_>>();
            for ((d, r), coefficient) in inputs {
                let limit = (degree - d) / m;
                for (k, slice) in factors.iter().enumerate().take(limit + 1).skip(1) {
                    for (j, &factor) in slice.iter().enumerate() {
                        model.spend(1)?;
                        if nonzero(&factor) {
                            let term = self.multiply(model, coefficient, factor)?;
                            add_algebra(model, &mut local, (d + m * k, r + j as i128 - 2), term)?;
                        }
                    }
                }
            }
        }
        let mut seed = Seed::new();
        for (key, value) in local {
            let coefficient = self.trace(model, value)?;
            add(model, &mut seed, key, coefficient)?;
        }
        Ok(seed)
    }

    fn oscillator(&self, model: &mut Model, degree: usize) -> Result<Vec<[Algebra; 5]>> {
        let one = scalar(1);
        let x = self.x();
        let twice_x = sum(model, x, x)?;
        let square = self.multiply(model, x, x)?;
        let middle = sum(model, square, scalar(2))?;
        // (1-xz+z²)^(-2); obtain the reciprocal by a four-term recurrence.
        let recurrence = [twice_x, middle.map(Integer::negate), twice_x, scalar(-1)];
        let mut reciprocal = Vec::with_capacity(degree + 1);
        reciprocal.push(one);
        for k in 1..=degree {
            let mut coefficient = scalar(0);
            for (index, &factor) in recurrence.iter().take(k).enumerate() {
                let term = self.multiply(model, factor, reciprocal[k - 1 - index])?;
                coefficient = sum(model, coefficient, term)?;
            }
            reciprocal.push(coefficient);
        }
        // (1-xyz+y²z²)(1-xy⁻¹z+y⁻²z²).
        let numerator = [one, x.map(Integer::negate), one];
        let mut factors = vec![[scalar(0); 5]; degree + 1];
        for (a, &left) in numerator.iter().enumerate() {
            for (b, &right) in numerator.iter().enumerate() {
                let product = self.multiply(model, left, right)?;
                for k in a + b..=degree {
                    let term = self.multiply(model, product, reciprocal[k - a - b])?;
                    let entry = &mut factors[k][2 + a - b];
                    *entry = sum(model, *entry, term)?;
                }
            }
        }
        Ok(factors)
    }
}

fn scalar(value: i128) -> Algebra {
    [
        Integer::from_i128(value),
        Integer::from_i128(0),
        Integer::from_i128(0),
    ]
}

fn nonzero(value: &Algebra) -> bool {
    value.iter().any(|coefficient| !coefficient.is_zero())
}

fn sum(model: &mut Model, mut left: Algebra, right: Algebra) -> Result<Algebra> {
    model.spend(REAL_MAX_DEGREE * INTEGER_WORDS)?;
    for (a, b) in left.iter_mut().zip(right) {
        *a = a.add(b).map_err(|_| overflow())?;
    }
    Ok(left)
}

fn add_algebra(
    model: &mut Model,
    series: &mut AlgebraSeries,
    index: (usize, i128),
    value: Algebra,
) -> Result<()> {
    model.spend(1)?;
    if !nonzero(&value) {
        return Ok(());
    }
    let full = series.len() >= ALGEBRA_TERM_LIMIT;
    match series.entry(index) {
        Entry::Occupied(mut entry) => {
            let coefficient = sum(model, *entry.get(), value)?;
            if nonzero(&coefficient) {
                *entry.get_mut() = coefficient;
            } else {
                entry.remove();
            }
        }
        Entry::Vacant(entry) => {
            if full {
                return Err(storage());
            }
            entry.insert(value);
        }
    }
    Ok(())
}

fn divide(model: &mut Model, value: Integer, divisor: u128, diagnostic: &str) -> Result<Integer> {
    let (magnitude, remainder) = value.magnitude.divide(model, divisor)?;
    if remainder != 0 {
        return Err(Error::plain(diagnostic));
    }
    Ok(Integer { magnitude, ..value })
}
