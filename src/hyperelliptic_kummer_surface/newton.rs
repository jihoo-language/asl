use std::fmt::Write;
use std::mem::size_of;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::generalized_kummer::cohomology::{self, dual_degree, geometric_torsion_terms};
use crate::hw_closed::Curve;
use crate::model::Model;
use crate::syntax::Label;

use super::hilbert;

type SurfaceSlopes = [[usize; 5]; 5];

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Kind {
    K3Hilbert,
    JacobianHilbert,
    Generalized,
}

pub(crate) enum Command {
    Polygon(Label, Kind, u128, u128),
    Height(Label),
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<String> {
        match self {
            Self::Height(name) => {
                let curve = model
                    .kummer_curve(&name.text)
                    .map_err(|error| error.at(name.span))?;
                let h1 = h1_slopes(model, &curve)?;
                let height = match h1[0] + h1[1] {
                    0 => "1",
                    1 => "2",
                    2 => "infinite",
                    _ => return Err(integrity()),
                };
                Ok(format!("rvKummerK3Height[{}] = {height}", name.text))
            }
            Self::Polygon(name, kind, n, i) => {
                let curve = model
                    .kummer_curve(&name.text)
                    .map_err(|error| error.at(name.span))?;
                let bins = polygon(model, &curve, kind, n, i)?;
                let operation = match kind {
                    Kind::K3Hilbert => "rvKummerHilbNewton",
                    Kind::JacobianHilbert => "rvJacHilbNewton",
                    Kind::Generalized => "rvGenKummerNewton",
                };
                Ok(format!("{operation}[{}, {n}, {i}] = [{bins}]", name.text))
            }
        }
    }
}

fn polygon(model: &mut Model, curve: &Curve, kind: Kind, n: u128, i: u128) -> Result<String> {
    let length = if kind == Kind::Generalized {
        if n == 0 {
            return Err(Error::plain("generalized Kummer length must be positive"));
        }
        if n.is_multiple_of(u128::from(curve.ring.prime())) {
            return Err(Error::plain("generalized Kummer Newton polygon requires length prime to the residue characteristic"));
        }
        n - 1
    } else {
        n
    };
    let Some((grade, twist)) = dual_degree(length, i) else {
        return Ok(String::new());
    };
    if kind == Kind::K3Hilbert && !i.is_multiple_of(2) {
        return Ok(String::new());
    }
    if grade == 0 {
        return Ok(format!("{}: 1", label(twist, 0)));
    }
    let expected = match kind {
        Kind::K3Hilbert => hilbert::betti(model, n, i)?,
        Kind::JacobianHilbert => cohomology::hilbert_betti(model, n, i)?,
        Kind::Generalized => cohomology::betti(model, n, i)?,
    };
    if expected == 0 {
        return Ok(String::new());
    }
    let h1 = h1_slopes(model, curve)?;
    if h1 == [1; 4] {
        return Ok(format!("{}: {expected}", label(0, i)));
    }
    let jacobian = exterior_slopes(h1);
    let surface = if kind == Kind::K3Hilbert {
        let mut k3 = [[0; 5]; 5];
        k3[0][0] = 1;
        k3[2] = jacobian[2];
        k3[2][2] += 16;
        k3[4][4] = 1;
        k3
    } else {
        jacobian
    };
    // Except for H⁰ at length one, each generator's degree is at least
    // its length. More copies of H⁰ just pad a low-degree class.
    let length = usize::try_from(n.min(grade)).map_err(|_| storage_limit())?;
    let grade = usize::try_from(grade).map_err(|_| storage_limit())?;
    let fock = Fock::new(model, &surface, length, grade)?;
    let combined;
    let row = if kind == Kind::Generalized {
        let n = usize::try_from(n).ok().filter(|&n| n <= grade);
        combined = generalized_row(model, &fock, &jacobian, n)?;
        &combined[grade * fock.width..(grade + 1) * fock.width]
    } else {
        fock.grade_row(length, grade)
    };
    let total = row
        .iter()
        .try_fold(0u128, |sum, &term| sum.checked_add(term))
        .ok_or_else(overflow)?;
    if total != expected {
        return Err(integrity());
    }
    let mut rendered = String::new();
    for (bin, &count) in row.iter().enumerate().filter(|(_, count)| **count != 0) {
        if !rendered.is_empty() {
            rendered.push_str(", ");
        }
        write!(rendered, "{}: {count}", label(twist, bin as u128)).unwrap();
    }
    Ok(rendered)
}

fn label(twist: u128, slope_twice: u128) -> String {
    let whole = twist + slope_twice / 2;
    if slope_twice.is_multiple_of(2) {
        whole.to_string()
    } else if let Some(numerator) = whole.checked_mul(2).and_then(|n| n.checked_add(1)) {
        format!("{numerator}/2")
    } else {
        format!("{whole}+1/2")
    }
}

fn h1_slopes(model: &mut Model, curve: &Curve) -> Result<[usize; 4]> {
    let p = u128::from(curve.ring.prime());
    let e = curve.ring.degree() as i32;
    let coefficients = curve.weil(model)?;
    let mut hull: Vec<(i32, i32)> = Vec::with_capacity(5);
    for (i, &coefficient) in coefficients.iter().enumerate().filter(|(_, c)| **c != 0) {
        let mut value = coefficient.unsigned_abs();
        let mut valuation = 0;
        while value.is_multiple_of(p) {
            model.spend(1)?;
            value /= p;
            valuation += 1;
        }
        let i = i as i32;
        while let [.., (x0, y0), (x1, y1)] = hull.as_slice() {
            if (y1 - y0) * (i - x1) < (valuation - y1) * (x1 - x0) {
                break;
            }
            hull.pop();
        }
        hull.push((i, valuation));
    }
    let mut slopes = [0; 4];
    for pair in hull.windows(2) {
        let [(x0, y0), (x1, y1)] = pair else {
            unreachable!()
        };
        let numerator = 2 * (y1 - y0);
        let denominator = (x1 - x0) * e;
        if numerator < 0 || numerator > 2 * denominator || numerator % denominator != 0 {
            return Err(Error::plain(
                "Jacobian Newton hull has an invalid genus-two slope",
            ));
        }
        slopes[*x0 as usize..*x1 as usize].fill((numerator / denominator) as usize);
    }
    if slopes[0] + slopes[3] != 2 || slopes[1] + slopes[2] != 2 {
        return Err(integrity());
    }
    Ok(slopes)
}

fn exterior_slopes(h1: [usize; 4]) -> SurfaceSlopes {
    let mut powers = [[0; 5]; 5];
    for subset in 0u32..16 {
        let slope: usize = h1
            .iter()
            .enumerate()
            .filter(|(index, _)| subset & (1 << index) != 0)
            .map(|(_, slope)| slope)
            .sum();
        powers[subset.count_ones() as usize][slope] += 1;
    }
    powers
}

struct Fock {
    length: usize,
    top: usize,
    width: usize,
    stride: usize,
    values: Vec<u128>,
}

impl Fock {
    fn new(model: &mut Model, surface: &SurfaceSlopes, length: usize, top: usize) -> Result<Self> {
        let width = top.saturating_mul(2).saturating_add(1);
        let stride = width.saturating_mul(top.saturating_add(1));
        let cells = stride.saturating_mul(length.saturating_add(2));
        if cells.saturating_mul(size_of::<u128>().div_ceil(size_of::<usize>())) > MAX_PIVOT_CELLS {
            return Err(storage_limit());
        }
        model.spend(cells)?;
        let mut result = Self {
            length,
            top,
            width,
            stride,
            values: vec![0; (length + 1) * stride],
        };
        result.values[0] = 1;
        for m in 1..=length {
            let shift = 2 * (m - 1);
            for (b, slopes) in surface.iter().enumerate() {
                let grade = b + shift;
                if grade > top {
                    break;
                }
                for (slope, &multiplicity) in slopes.iter().enumerate() {
                    for _ in 0..multiplicity {
                        result.oscillator(model, m, grade, slope + shift, !b.is_multiple_of(2))?;
                    }
                }
            }
        }
        Ok(result)
    }

    fn oscillator(
        &mut self,
        model: &mut Model,
        length: usize,
        grade: usize,
        slope: usize,
        odd: bool,
    ) -> Result<()> {
        let last = self.length - length;
        // Ascending length permits another symmetric copy. Descending length
        // uses each odd generator at most once, as for exterior powers.
        for index in 0..=last {
            let l = if odd { last - index } else { index };
            for i in 0..=self.top - grade {
                let from = l * self.stride + i * self.width;
                let to = (l + length) * self.stride + (i + grade) * self.width + slope;
                model.spend(2 * i + 1)?;
                for bin in 0..=2 * i {
                    self.values[to + bin] = self.values[to + bin]
                        .checked_add(self.values[from + bin])
                        .ok_or_else(overflow)?;
                }
            }
        }
        Ok(())
    }

    fn grade_row(&self, length: usize, grade: usize) -> &[u128] {
        let start = length * self.stride + grade * self.width;
        &self.values[start..start + self.width]
    }
}

fn generalized_row(
    model: &mut Model,
    fock: &Fock,
    jacobian: &SurfaceSlopes,
    actual_length: Option<usize>,
) -> Result<Vec<u128>> {
    model.spend(fock.stride)?;
    let start = fock.length * fock.stride;
    let mut quotient = fock.values[start..start + fock.stride].to_vec();
    if let Some(n) = actual_length {
        for (k, shift, count) in geometric_torsion_terms(model, n, fock.top)? {
            for grade in 0..=fock.top - shift {
                let start = (grade + shift) * fock.width + shift;
                model.spend(2 * grade + 1)?;
                for (bin, &source) in fock.grade_row(k, grade)[..=2 * grade].iter().enumerate() {
                    let term = source.checked_mul(count).ok_or_else(overflow)?;
                    quotient[start + bin] = quotient[start + bin]
                        .checked_add(term)
                        .ok_or_else(overflow)?;
                }
            }
        }
    }
    // Divide the polynomial in cohomological degree, retaining every slope.
    // The constant coefficient of H*(A) is 1, so no division of integers is needed.
    for grade in 1..=fock.top {
        for (b, slopes) in jacobian.iter().enumerate().take(grade + 1).skip(1) {
            let from = (grade - b) * fock.width;
            for (slope, &multiplicity) in slopes.iter().enumerate().filter(|(_, m)| **m != 0) {
                let to = grade * fock.width + slope;
                model.spend(2 * (grade - b) + 1)?;
                for bin in 0..=2 * (grade - b) {
                    let term = quotient[from + bin]
                        .checked_mul(multiplicity as u128)
                        .ok_or_else(overflow)?;
                    quotient[to + bin] =
                        quotient[to + bin].checked_sub(term).ok_or_else(integrity)?;
                }
            }
        }
    }
    Ok(quotient)
}

fn storage_limit() -> Error {
    Error::plain("crystalline Newton calculation exceeds the coefficient storage limit")
}

fn overflow() -> Error {
    Error::plain("crystalline Newton multiplicity exceeds the checked unsigned 128-bit range")
}

fn integrity() -> Error {
    Error::plain("crystalline Newton calculation failed its multiplicity check")
}
