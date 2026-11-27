use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::Model;
use crate::rv_ore::Powers;

use super::{integrity, Operator};

const MAX_DEGREE: usize = 12;
type Taylor = [Rational; 3];
type Series = Vec<Rational>;

const ZERO: Rational = Rational::ZERO;
const ONE: Rational = Rational::integer(1);

fn degree_ok(degree: usize) -> Result<()> {
    if degree > MAX_DEGREE {
        return Err(Error::plain(format!(
            "Dwork mirror jet degree must be at most {MAX_DEGREE}"
        )));
    }
    Ok(())
}

fn range(_: Error) -> Error {
    Error::plain("Dwork mirror exact rational arithmetic exceeds the signed 128-bit range")
}

fn plus(x: Rational, y: Rational) -> Result<Rational> {
    x.add(y).map_err(range)
}

fn times(x: Rational, y: Rational) -> Result<Rational> {
    x.multiply(y).map_err(range)
}

fn ratio(n: i128, d: i128) -> Result<Rational> {
    Rational::ratio(n, d).map_err(range)
}

fn scale(x: Rational, n: i128) -> Result<Rational> {
    times(x, Rational::integer(n))
}

fn minus(x: Rational, y: Rational) -> Result<Rational> {
    plus(x, scale(y, -1)?)
}

fn display(x: Rational) -> String {
    if x.denominator == 1 {
        x.numerator.to_string()
    } else {
        format!("{}/{}", x.numerator, x.denominator)
    }
}

fn display_row(x: &[Rational]) -> String {
    format!(
        "[{}]",
        x.iter()
            .copied()
            .map(display)
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn taylor_product(x: Taylor, y: Taylor, rank: usize) -> Result<Taylor> {
    let mut output = [ZERO; 3];
    for (d, coefficient) in output.iter_mut().take(rank).enumerate() {
        for i in 0..=d {
            *coefficient = plus(*coefficient, times(x[i], y[d - i])?)?;
        }
    }
    Ok(output)
}

struct Periods {
    family: Operator,
    coefficients: Vec<Taylor>,
}

impl Periods {
    fn new(model: &mut Model, family: Operator, degree: usize) -> Result<Self> {
        degree_ok(degree)?;
        let rank = family.rank();
        let (slope, offsets) = family.factors();
        let mut previous = [ONE, ZERO, ZERO];
        let mut coefficients = vec![previous];
        for n in 1..=degree {
            model.spend(rank.pow(3))?;
            let n = n as i128;
            let r = rank as u32;
            let mut numerator = [Rational::integer(4), ZERO, ZERO];
            for &offset in offsets {
                // The z-shift applies P(theta) at n - 1 + s.
                let linear = [
                    Rational::integer(slope * (n - 1) + offset),
                    Rational::integer(slope),
                    ZERO,
                ];
                numerator = taylor_product(numerator, linear, rank)?;
            }
            let r_int = i128::from(r);
            let reciprocal = [
                ratio(1, n.pow(r))?,
                ratio(-r_int, n.pow(r + 1))?,
                ratio(r_int * (r_int + 1) / 2, n.pow(r + 2))?,
            ];
            let step = taylor_product(numerator, reciprocal, rank)?;
            previous = taylor_product(previous, step, rank)?;
            coefficients.push(previous);
        }
        let periods = Self {
            family,
            coefficients,
        };
        for j in 0..rank {
            let residual = periods.jet(j)?.apply(model, family.terms())?;
            if !residual.is_zero() {
                return Err(integrity());
            }
        }
        Ok(periods)
    }

    fn row(&self, index: usize) -> Series {
        self.coefficients.iter().map(|c| c[index]).collect()
    }

    fn jet(&self, solution: usize) -> Result<LogJet> {
        if solution >= self.family.rank() {
            return Err(Error::plain(
                "Dwork mirror solution index is outside its operator rank",
            ));
        }
        Ok(LogJet(
            self.coefficients
                .iter()
                .map(|c| std::array::from_fn(|j| solution.checked_sub(j).map_or(ZERO, |i| c[i])))
                .collect(),
        ))
    }
}

struct LogJet(Vec<Taylor>);

impl LogJet {
    fn is_zero(&self) -> bool {
        self.0.iter().flatten().all(|&c| c == ZERO)
    }

    fn apply(
        &self,
        model: &mut Model,
        operator: impl IntoIterator<Item = (Powers, i128)>,
    ) -> Result<Self> {
        let mut output = vec![[ZERO; 3]; self.0.len()];
        for (Powers { q, theta }, scale_by) in operator {
            let count = self.0.len().saturating_sub(q);
            let work = count
                .checked_mul(theta.saturating_add(1))
                .ok_or_else(operator_limit)?;
            model.spend(work).map_err(|_| operator_limit())?;
            for (n, (&value, output)) in self.0.iter().zip(output.iter_mut().skip(q)).enumerate() {
                let mut c = value;
                for _ in 0..theta {
                    if c == [ZERO; 3] {
                        break;
                    }
                    // [ell^j/j!] theta(z^n c(ell)) = n c_j + c_(j+1).
                    c = differentiate_log(c, n)?;
                }
                for (target, value) in output.iter_mut().zip(c) {
                    *target = plus(*target, scale(value, scale_by)?)?;
                }
            }
        }
        Ok(Self(output))
    }

    fn shift(&self, shift: Rational, half_square: Rational) -> Result<Self> {
        let shifted = self.0.iter().map(|&[a, b, c]| {
            Ok([
                plus(plus(a, times(shift, b)?)?, times(half_square, c)?)?,
                plus(b, times(shift, c)?)?,
                c,
            ])
        });
        Ok(Self(shifted.collect::<Result<_>>()?))
    }

    fn add_scaled(&mut self, other: &Self, scalar: Rational) -> Result<()> {
        for (left, right) in self.0.iter_mut().zip(&other.0) {
            for (x, &y) in left.iter_mut().zip(right) {
                *x = plus(*x, times(y, scalar)?)?;
            }
        }
        Ok(())
    }

    fn report(&self) -> String {
        if self.is_zero() {
            return "0".into();
        }
        let rows = std::array::from_fn::<_, 3, _>(|j| {
            display_row(&self.0.iter().map(|c| c[j]).collect::<Vec<_>>())
        });
        format!(
            "coefficients of 1, log z, (log z)^2/2 = [{}, {}, {}]",
            rows[0], rows[1], rows[2]
        )
    }
}

fn differentiate_log([a, b, c]: Taylor, n: usize) -> Result<Taylor> {
    let n = n as i128;
    Ok([plus(scale(a, n)?, b)?, plus(scale(b, n)?, c)?, scale(c, n)?])
}

fn operator_limit() -> Error {
    Error::plain("Dwork logarithmic Ore application exceeds the computation limit")
}

pub(super) fn apply(
    model: &mut Model,
    family: Operator,
    solution: usize,
    degree: usize,
    terms: &[(Powers, i128)],
) -> Result<String> {
    let periods = Periods::new(model, family, degree)?;
    Ok(periods
        .jet(solution)?
        .apply(model, terms.iter().copied())?
        .report())
}

fn theta(series: &[Rational]) -> Result<Series> {
    series
        .iter()
        .enumerate()
        .map(|(n, &c)| scale(c, n as i128))
        .collect()
}

fn product(model: &mut Model, a: &[Rational], b: &[Rational]) -> Result<Series> {
    let size = a.len();
    debug_assert_eq!(size, b.len());
    model.spend(size * (size + 1) / 2)?;
    (0..size)
        .map(|d| (0..=d).try_fold(ZERO, |sum, i| plus(sum, times(a[i], b[d - i])?)))
        .collect()
}

fn quotient_by_unit(model: &mut Model, a: &[Rational], b: &[Rational]) -> Result<Series> {
    if b[0] != ONE {
        return Err(integrity());
    }
    model.spend(a.len() * a.len())?;
    let mut result = Vec::new();
    for (n, &c) in a.iter().enumerate() {
        let mut coefficient = c;
        for (i, &q) in result.iter().enumerate() {
            coefficient = minus(coefficient, times(q, b[n - i])?)?;
        }
        result.push(coefficient);
    }
    Ok(result)
}

struct Variation {
    f: Series,
    g: Series,
    h: Series,
    u: Series,
    v: Series,
    log_ratio: Series,
    wronskian: Series,
}

impl Variation {
    fn new(model: &mut Model, degree: usize) -> Result<Self> {
        let pf = Periods::new(model, Operator::PicardFuchs, degree)?;
        let cl = Periods::new(model, Operator::Clausen, degree)?;
        let (f, g, h) = (pf.row(0), pf.row(1), pf.row(2));
        let (u, v) = (cl.row(0), cl.row(1));
        let twice_h = h.iter().map(|&x| scale(x, 2)).collect::<Result<Vec<_>>>()?;
        if product(model, &u, &u)? != f
            || product(model, &u, &v)? != g
            || product(model, &v, &v)? != twice_h
            || product(model, &f, &twice_h)? != product(model, &g, &g)?
        {
            return Err(integrity());
        }
        let log_ratio = quotient_by_unit(model, &g, &f)?;
        let mut logarithmic_derivative = theta(&log_ratio)?;
        logarithmic_derivative[0] = ONE;
        let wronskian = product(model, &f, &logarithmic_derivative)?;
        let u_dv = product(model, &u, &theta(&v)?)?;
        let v_du = product(model, &v, &theta(&u)?)?;
        let mut expected = ONE;
        for (n, &coefficient) in wronskian.iter().enumerate() {
            if n != 0 {
                expected = times(expected, ratio(128 * (2 * n as i128 - 1), n as i128)?)?;
            }
            if coefficient != expected || coefficient != plus(f[n], minus(u_dv[n], v_du[n])?)? {
                return Err(integrity());
            }
        }
        Ok(Self {
            f,
            g,
            h,
            u,
            v,
            log_ratio,
            wronskian,
        })
    }
}

pub(super) fn period_report(model: &mut Model, degree: usize) -> Result<String> {
    let v = Variation::new(model, degree)?;
    Ok(format!(
        "rvDworkMirrorPeriods[{degree}] = f={}, g={}, h={}, u={}, v={}, W={}, certified Ore equations and symmetric-square identities",
        display_row(&v.f), display_row(&v.g), display_row(&v.h),
        display_row(&v.u), display_row(&v.v), display_row(&v.wronskian),
    ))
}

pub(super) fn monodromy(model: &mut Model, degree: usize, shift: i128) -> Result<String> {
    let rank = Operator::PicardFuchs.rank();
    let periods = Periods::new(model, Operator::PicardFuchs, degree)?;
    let jets = (0..rank)
        .map(|i| periods.jet(i))
        .collect::<Result<Vec<_>>>()?;
    let k = Rational::integer(shift);
    let half_square = times(k, ratio(shift, 2)?)?;
    let mut matrix = [[ZERO; 3]; 3];
    for (j, jet) in jets.iter().enumerate() {
        let shifted = jet.shift(k, half_square)?;
        let column = shifted.0[0];
        let mut recomposed = LogJet(vec![[ZERO; 3]; degree + 1]);
        for (i, (&c, basis)) in column.iter().zip(&jets).enumerate() {
            matrix[i][j] = c;
            recomposed.add_scaled(basis, c)?;
        }
        if recomposed.0 != shifted.0
            || !shifted
                .apply(model, Operator::PicardFuchs.terms())?
                .is_zero()
        {
            return Err(integrity());
        }
    }
    let rows = matrix
        .iter()
        .map(|row| display_row(row))
        .collect::<Vec<_>>();
    Ok(format!(
        "rvDworkMirrorMonodromy[{degree}, {shift}] = [{}, {}, {}], normalized log generator has one unipotent block of size 3",
        rows[0], rows[1], rows[2],
    ))
}

fn exponential(model: &mut Model, input: &[Rational]) -> Result<Series> {
    if input[0] != ZERO {
        return Err(integrity());
    }
    let derivative = theta(input)?;
    model.spend(input.len() * input.len())?;
    let mut result = vec![ONE];
    for n in 1..input.len() {
        let coefficient = (1..=n).try_fold(ZERO, |sum, i| {
            plus(sum, times(derivative[i], result[n - i])?)
        })?;
        result.push(times(coefficient, ratio(1, n as i128)?)?);
    }
    if theta(&result)? != product(model, &result, &derivative)? {
        return Err(integrity());
    }
    Ok(result)
}

fn compose(model: &mut Model, outer: &[Rational], inner: &[Rational]) -> Result<Series> {
    debug_assert_eq!(inner[0], ZERO);
    let mut result = vec![ZERO; inner.len()];
    result[0] = outer[0];
    let mut power = inner.to_vec();
    // A monomial of degree j starts in degree j after substitution.
    // Building from low powers keeps large coefficients in their true grade.
    for (j, &coefficient) in outer.iter().enumerate().skip(1) {
        for (value, &p) in result[j..].iter_mut().zip(&power[j..]) {
            *value = plus(*value, times(coefficient, p)?)?;
        }
        if j + 1 < outer.len() {
            power = product(model, &power, inner)?;
        }
    }
    Ok(result)
}

fn revert(model: &mut Model, series: &[Rational]) -> Result<Series> {
    if series[0] != ZERO || series[1] != ONE {
        return Err(integrity());
    }
    let mut inverse = vec![ZERO; series.len()];
    inverse[1] = ONE;
    for d in 2..series.len() {
        let residual = compose(model, &series[..=d], &inverse[..=d])?[d];
        inverse[d] = scale(residual, -1)?;
    }
    let mut identity = vec![ZERO; series.len()];
    identity[1] = ONE;
    if compose(model, series, &inverse)? != identity
        || compose(model, &inverse, series)? != identity
    {
        return Err(integrity());
    }
    Ok(inverse)
}

pub(super) fn mirror_map(model: &mut Model, degree: usize) -> Result<String> {
    let periods = Variation::new(model, degree)?;
    let q_over_z = exponential(model, &periods.log_ratio)?;
    let q = std::iter::once(ZERO)
        .chain(q_over_z.iter().copied())
        .collect::<Vec<_>>();
    let z = revert(model, &q)?;
    Ok(format!(
        "rvDworkMirrorMap[{degree}] = log(q/z)={}, q(z)/z={}, z(q)/q={}, certified both compositions",
        display_row(&periods.log_ratio), display_row(&q_over_z), display_row(&z[1..]),
    ))
}
