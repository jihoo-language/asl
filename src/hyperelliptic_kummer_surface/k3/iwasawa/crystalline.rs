use crate::elliptic_modular::decimal;
use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::hyperelliptic_kummer_surface::newton::formal_height;
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::wide_integer::Natural;

use super::super::artin_tate::artin_invariant;
use super::{inconsistent, modular, Breakpoint, Spectrum};

// The linear term and mu are at most 135 bits. This stores them together
// with an additional 128-bit exponentiation factor when it is available.
type Arithmetic = Wide<9>;

pub(super) struct Crystalline {
    prime: u128,
    height: Option<usize>,
    mu: Arithmetic,
    mu_per_degree: i64,
    field_valuation: i64,
    negative_sum: i64,
    negative: Vec<modular::Slope>,
    positive: Vec<modular::Slope>,
    sigma: Option<u32>,
}

impl Crystalline {
    pub(super) fn new(
        model: &mut Model,
        curve: &Curve,
        spectrum: &Spectrum<'_>,
        base: u128,
        prime: u128,
        (rank, geometric_rank): (usize, usize),
    ) -> Result<Self> {
        let height = formal_height(model, curve)?;
        let field_valuation = modular::valuation(prime, u128::from(spectrum.q))
            .filter(|&valuation| valuation > 0)
            .ok_or_else(inconsistent)?;
        let slopes = modular::unshifted(
            model,
            &spectrum.abelian.transcendental,
            prime,
            field_valuation,
        )?;
        let expected = if let Some(height) = height {
            // For a K3 of finite height h there are h slopes 1-1/h.
            let count = slopes
                .negative
                .iter()
                .map(|slope| slope.count)
                .sum::<usize>();
            if count != height
                || !slopes
                    .negative
                    .iter()
                    .all(|&slope| slope.matches_height(field_valuation, height))
                || geometric_rank == 22
            {
                return Err(height_inconsistency());
            }
            0
        } else {
            if geometric_rank != 22 || !slopes.negative.is_empty() {
                return Err(height_inconsistency());
            }
            field_valuation
        };
        let mu_per_degree = field_valuation + slopes.negative_sum;
        if mu_per_degree != expected {
            return Err(height_inconsistency());
        }
        let mu = Natural::from_small(base)
            .multiply_u128(mu_per_degree as u128)
            .ok_or_else(inconsistent)?;
        let mu = Arithmetic {
            magnitude: mu,
            negative: false,
        };
        let sigma = if rank == 22 {
            let sigma = artin_invariant(model, curve)?;
            if mu
                .magnitude
                .cmp(Natural::from_small(u128::from(2 * sigma)))
                .is_lt()
            {
                return Err(Error::plain(
                    "Kummer K3 crystalline Iwasawa base order conflicts with the Artin invariant",
                ));
            }
            Some(sigma)
        } else {
            None
        };
        Ok(Self {
            prime,
            height,
            mu,
            mu_per_degree,
            field_valuation,
            negative_sum: slopes.negative_sum,
            negative: slopes.negative,
            positive: slopes.positive,
            sigma,
        })
    }

    pub(super) fn normalization(&self) -> (i64, i64) {
        (self.field_valuation, self.negative_sum)
    }

    pub(super) fn negative_count(&self) -> usize {
        self.negative.iter().map(|slope| slope.count).sum()
    }

    pub(super) fn has_absolute(&self) -> bool {
        self.sigma.is_some()
    }

    pub(super) fn report(
        &self,
        group: &str,
        base: u128,
        level: u128,
        breakpoint: &Breakpoint,
    ) -> Result<String> {
        let scaled = if self.negative.is_empty() {
            Vec::new()
        } else {
            let exponent = self
                .prime
                .checked_pow(breakpoint.index() as u32)
                .and_then(|power| power.checked_mul(base))
                .ok_or_else(inconsistent)?;
            self.negative
                .iter()
                .map(|&slope| slope.scaled(exponent, true))
                .collect::<Result<Vec<_>>>()?
        };
        let source_slopes = self
            .negative
            .iter()
            .map(|&slope| slope.scaled(base, true))
            .chain(self.positive.iter().map(|&slope| slope.scaled(base, false)))
            .collect::<Result<Vec<_>>>()?;
        let height = self
            .height
            .map_or_else(|| "infinite".into(), |h| h.to_string());
        let relative = breakpoint.relative(level)?.resize::<9>()?;
        let from_base = breakpoint.relative_valuation(level)?.resize::<9>()?;
        let at_breakpoint = self
            .exponential(breakpoint.index())
            .ok_or_else(inconsistent)?;
        let relative = self.render(level, relative.add(at_breakpoint.negate())?)?;
        let from_base = self.render(level, from_base.add(self.mu.negate())?)?;
        let mut output = format!(
            "formal Brauer height={height}, nonzero v_p(u) at base=[{}], negative Newton valuations=[{}], positive determinant valuations through breakpoint={:?}, mu={}, lambda={}, relative v_ell({group}) from breakpoint={relative}, relative v_ell({group}) from base={from_base}",
            source_slopes.join(", "), scaled.join(", "), breakpoint.prefix, decimal(self.mu), breakpoint.lambda,
        );
        if let Some(sigma) = self.sigma {
            let discr = self.prime.checked_pow(2 * sigma).ok_or_else(inconsistent)?;
            let correction = Arithmetic::from_i128(-2 * i128::from(sigma));
            let base_order = decimal(self.mu.add(correction)?);
            let absolute = self.render(level, correction)?;
            output.push_str(&format!(", Artin invariant={sigma}, NS discriminant=-{discr}, base v_ell({group})={base_order}, v_ell({group})={absolute}"));
        }
        Ok(output)
    }

    pub(super) fn absolute_valuation(
        &self,
        model: &mut Model,
        spectrum: &Spectrum<'_>,
        base_order: u128,
        level: u128,
        breakpoint: &Breakpoint,
    ) -> Result<String> {
        let reference = modular::certify(
            model,
            spectrum.q,
            &spectrum.abelian.transcendental,
            1,
            self.prime,
            0,
        )?;
        let known = modular::valuation(self.prime, base_order).ok_or_else(inconsistent)?;
        let fixed = known + breakpoint.prefix[0] - self.mu_per_degree - reference.determinant;
        let offset = breakpoint
            .relative_valuation(level)?
            .resize::<9>()?
            .add(Arithmetic::from_i128(i128::from(fixed)))?;
        if let Some(exponential) = self.exponential(level) {
            let order = exponential.add(offset)?;
            if order.negative && !order.is_zero() {
                return Err(inconsistent());
            }
        }
        self.render(level, offset)
    }

    fn exponential(&self, level: u128) -> Option<Arithmetic> {
        if self.mu.is_zero() {
            return Some(self.mu);
        }
        let power = u32::try_from(level)
            .ok()
            .and_then(|k| self.prime.checked_pow(k))?;
        let magnitude = self.mu.magnitude.multiply_u128(power)?;
        Some(Arithmetic {
            magnitude,
            negative: false,
        })
    }

    fn render(&self, level: u128, correction: Arithmetic) -> Result<String> {
        if let Some(value) = self.exponential(level) {
            return Ok(decimal(value.add(correction)?));
        }
        let base = format!("{}*{}^{level}", decimal(self.mu), self.prime);
        if correction.is_zero() {
            Ok(base)
        } else if correction.negative {
            Ok(format!("{base} - {}", decimal(correction.negate())))
        } else {
            Ok(format!("{base} + {}", decimal(correction)))
        }
    }
}

fn height_inconsistency() -> Error {
    Error::plain("Kummer K3 crystalline Iwasawa slopes disagree with the formal Brauer height")
}
