use std::cmp::Ordering;

use crate::error::Result;
use crate::model::Model;
use crate::rv_ore::gcd;

#[derive(Clone, Copy, Eq, PartialEq)]
pub(crate) struct Natural<const N: usize>([u32; N]);

impl<const N: usize> Default for Natural<N> {
    fn default() -> Self {
        Self([0; N])
    }
}

impl<const N: usize> Natural<N> {
    pub(crate) fn from_small(mut n: u128) -> Self {
        let mut value = Self::default();
        for digit in &mut value.0[..4] {
            *digit = n as u32;
            n >>= 32;
        }
        value
    }

    pub(crate) fn small(self) -> Option<u128> {
        if self.0[4..].iter().any(|&digit| digit != 0) {
            return None;
        }
        Some(
            self.0[..4]
                .iter()
                .rev()
                .fold(0, |acc, &digit| (acc << 32) | u128::from(digit)),
        )
    }

    pub(crate) fn signed(self, negative: bool) -> Option<i128> {
        let value = self.small()?;
        if negative && value == 1u128 << 127 {
            return Some(i128::MIN);
        }
        let value = i128::try_from(value).ok()?;
        Some(if negative { -value } else { value })
    }

    pub(crate) fn product(left: i128, right: i128) -> Self {
        let left = Self::from_small(left.unsigned_abs());
        let right = Self::from_small(right.unsigned_abs());
        let mut result = Self::default();
        for (i, &a) in left.0[..4].iter().enumerate() {
            let mut carry = 0;
            for (j, &b) in right.0[..4].iter().enumerate() {
                let digit = &mut result.0[i + j];
                carry += u64::from(*digit) + u64::from(a) * u64::from(b);
                *digit = carry as u32;
                carry >>= 32;
            }
            result.0[i + 4] = carry as u32;
        }
        result
    }

    pub(crate) fn bits(self) -> usize {
        self.0
            .iter()
            .rposition(|&digit| digit != 0)
            .map_or(0, |i| 32 * i + (32 - self.0[i].leading_zeros()) as usize)
    }

    pub(crate) fn cmp(self, other: Self) -> Ordering {
        self.0.iter().rev().cmp(other.0.iter().rev())
    }

    pub(crate) fn add_multiple(&mut self, other: Self, multiplier: u32) -> Option<()> {
        let mut carry = 0;
        for (digit, other) in self.0.iter_mut().zip(other.0) {
            carry += u64::from(*digit) + u64::from(other) * u64::from(multiplier);
            *digit = carry as u32;
            carry >>= 32;
        }
        (carry == 0).then_some(())
    }

    pub(crate) fn difference(mut self, other: Self) -> Self {
        let mut borrow = 0i64;
        for (digit, other) in self.0.iter_mut().zip(other.0) {
            let remaining = i64::from(*digit) - i64::from(other) - borrow;
            *digit = remaining as u32;
            borrow = i64::from(remaining < 0);
        }
        debug_assert_eq!(borrow, 0);
        self
    }

    pub(crate) fn residue(self, p: u32) -> u32 {
        self.0.iter().rev().fold(0, |r, &digit| {
            ((u64::from(r) << 32 | u64::from(digit)) % u64::from(p)) as u32
        })
    }

    pub(crate) fn divide(self, model: &mut Model, divisor: u128) -> Result<(Self, u128)> {
        model.spend(N)?;
        if let Some(value) = self.small() {
            return Ok((Self::from_small(value / divisor), value % divisor));
        }
        if divisor == 1 {
            return Ok((self, 0));
        }
        // Keep the high bit of the shifted remainder separately. On
        // that branch subtraction brings its magnitude back into u128.
        let width = self.bits();
        model.spend(width)?;
        let mut quotient = Self::default();
        let mut remainder = 0u128;
        for bit in (0..width).rev() {
            let (i, j) = (bit / 32, bit % 32);
            let high = remainder >> 127;
            remainder = (remainder << 1) | u128::from((self.0[i] >> j) & 1);
            if high != 0 || remainder >= divisor {
                remainder = remainder.wrapping_sub(divisor);
                quotient.0[i] |= 1 << j;
            }
        }
        Ok((quotient, remainder))
    }
}

impl<const N: usize> Natural<N> {
    pub(crate) fn is_zero(self) -> bool {
        self.0.iter().all(|&digit| digit == 0)
    }

    pub(crate) fn multiply_u128(self, factor: u128) -> Option<Self> {
        let rhs = Self::from_small(factor);
        let width = rhs.bits().div_ceil(32);
        let mut result = Self::default();
        for (i, &a) in self.0.iter().enumerate().filter(|&(_, &a)| a != 0) {
            let mut carry = 0;
            for (j, &b) in rhs.0[..width].iter().enumerate() {
                let partial = u64::from(a) * u64::from(b) + carry;
                if let Some(digit) = result.0.get_mut(i + j) {
                    let value = partial + u64::from(*digit);
                    *digit = value as u32;
                    carry = value >> 32;
                } else if partial != 0 {
                    return None;
                }
            }
            if carry != 0 {
                *result.0.get_mut(i + width)? = carry as u32;
            }
        }
        Some(result)
    }

    pub(crate) fn divide_u32(self, divisor: u32) -> (Self, u32) {
        let mut quotient = self;
        let mut remainder = 0u64;
        for digit in quotient.0.iter_mut().rev() {
            let value = (remainder << 32) | u64::from(*digit);
            *digit = (value / u64::from(divisor)) as u32;
            remainder = value % u64::from(divisor);
        }
        (quotient, remainder as u32)
    }

    pub(crate) fn gcd(mut self, model: &mut Model, mut rhs: Self) -> Result<Self> {
        if self.is_zero() {
            return Ok(rhs);
        }
        if rhs.is_zero() {
            return Ok(self);
        }
        if let Some(d) = self.small().or_else(|| rhs.small()) {
            let larger = if self.small().is_some() { rhs } else { self };
            let (_, remainder) = larger.divide(model, d)?;
            return Ok(Self::from_small(gcd(d, remainder)));
        }
        let common = self.trailing_zeros().min(rhs.trailing_zeros());
        self = self.shift_right(self.trailing_zeros());
        rhs = rhs.shift_right(rhs.trailing_zeros());
        loop {
            model.spend(4 * N)?;
            if let Some(d) = rhs.small() {
                let (_, r) = self.divide(model, d)?;
                return Ok(Self::from_small(gcd(d, r)).shift_left(common));
            }
            if self.cmp(rhs).is_gt() {
                (self, rhs) = (rhs, self);
            }
            rhs = rhs.difference(self);
            if rhs.is_zero() {
                return Ok(self.shift_left(common));
            }
            rhs = rhs.shift_right(rhs.trailing_zeros());
        }
    }

    pub(crate) fn divide_exact(self, model: &mut Model, divisor: Self) -> Result<Option<Self>> {
        if let Some(divisor) = divisor.small() {
            let (quotient, remainder) = self.divide(model, divisor)?;
            return Ok((remainder == 0).then_some(quotient));
        }
        if self.is_zero() {
            return Ok(Some(self));
        }
        let Some(offset) = self.bits().checked_sub(divisor.bits()) else {
            return Ok(None);
        };
        model.spend((offset + 1) * 3 * N)?;
        let mut remainder = self;
        let mut quotient = Self::default();
        let mut shifted = divisor.shift_left(offset);
        for bit in (0..=offset).rev() {
            if !remainder.cmp(shifted).is_lt() {
                remainder = remainder.difference(shifted);
                quotient.0[bit / 32] |= 1 << (bit % 32);
            }
            shifted = shifted.shift_right(1);
        }
        Ok(remainder.is_zero().then_some(quotient))
    }

    fn trailing_zeros(self) -> usize {
        self.0
            .iter()
            .position(|&digit| digit != 0)
            .map_or(32 * N, |i| 32 * i + self.0[i].trailing_zeros() as usize)
    }

    fn shift_left(self, offset: usize) -> Self {
        debug_assert!(self.bits() + offset <= 32 * N);
        let whole = offset / 32;
        let amount = offset % 32;
        let mut result = Self::default();
        for (i, &digit) in self.0.iter().take(N - whole).enumerate() {
            result.0[i + whole] |= digit << amount;
            if amount != 0 && i + whole + 1 < N {
                result.0[i + whole + 1] |= digit >> (32 - amount);
            }
        }
        result
    }

    fn shift_right(self, offset: usize) -> Self {
        let whole = offset / 32;
        let amount = offset % 32;
        let mut result = Self::default();
        for (i, &digit) in self.0.iter().enumerate().skip(whole) {
            result.0[i - whole] |= digit >> amount;
            if amount != 0 && i > whole {
                result.0[i - whole - 1] |= digit << (32 - amount);
            }
        }
        result
    }
}
