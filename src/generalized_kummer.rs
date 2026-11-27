use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::hyperelliptic_jacobian::{self as jacobian, Divisor};
use crate::hyperelliptic_kummer;
use crate::hyperelliptic_kummer_surface::k3::{node_orbits, powered_weil_coefficients};
use crate::jacobian_level::Basis;
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear;
use crate::rv_ore::gcd;
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod cohomology;
pub(crate) mod hodge;
pub(crate) mod inversion;

pub(crate) enum Command {
    Hilbert(Label, u128, u128),
    Zero(Label, u128),
    Fiber(Label, u128, Box<ParsedExpr>),
    Level(Label, Label, u128),
    Betti(u128, u128),
    Trace(Label, Label, u128, u128),
    HilbertBetti(u128, u128),
    HilbertTrace(Label, u128, u128, u128),
    TranslateChar(u128, u128, u128),
    TranslateFixed(u128, u128),
    TranslateFrob(Label, Label, u128, u128, Box<ParsedExpr>),
    TranslatePoints(Label, Label, u128, Box<ParsedExpr>),
    InversionEuler(u128, u128),
    InversionFrob(Label, Label, u128, u128, Box<ParsedExpr>),
    InversionPoints(Label, Label, u128, Box<ParsedExpr>),
    StratumFrob(Label, Label, Label, u128, u128, Box<ParsedExpr>),
}

pub(crate) fn require_curve(curve: &Curve) -> Result<()> {
    if curve.grade.genus != 2 {
        return Err(Error::plain(
            "Jacobian Hilbert schemes require a genus-two curve",
        ));
    }
    Ok(())
}

pub(crate) fn hilbert(
    model: &mut Model,
    curve: &Curve,
    length: u128,
    degree: u128,
) -> Result<u128> {
    if length == 0 {
        return Ok(1);
    }
    let (_, counts) = abelian_counts(model, curve, length, degree)?;
    Ok(*counts.last().expect("Hilbert count includes length zero"))
}

pub(crate) fn fiber(
    model: &mut Model,
    curve: &Curve,
    length: u128,
    target: Option<&Divisor>,
) -> Result<u128> {
    if length == 1 {
        return Ok(1);
    }
    let (q, hilbert) = abelian_counts(model, curve, length, 1)?;
    let order = hilbert[1];
    let relevant = gcd(length, order) as usize;
    let mut statistics: Vec<_> = (2..=relevant)
        .filter(|d| relevant.is_multiple_of(*d))
        .map(|d| DivisorStatistic {
            degree: d,
            kernel: 0,
            in_image: false,
        })
        .collect();
    let target = target.filter(|d| d.u.len() != 1);
    if let Some(two) = statistics.first_mut().filter(|entry| entry.degree == 2) {
        two.kernel = node_orbits(model, curve, q as u64)?
            .into_iter()
            .filter(|&d| d == 1)
            .count() as u128;
        two.in_image = match target {
            Some(divisor) => hyperelliptic_kummer::support(model, curve, divisor)?.is_trivial(),
            None => true,
        };
    }
    let first_unmeasured = usize::from(relevant.is_multiple_of(2));
    let to_enumerate = &mut statistics[first_unmeasured..];
    if !to_enumerate.is_empty() {
        collect_statistics(model, curve, order, target, to_enumerate)?;
    }
    let n = hilbert.len() - 1;
    let mut positive = hilbert[n];
    let mut negative = 0u128;
    let mut weights = vec![(1, 1i128)];
    for entry in statistics {
        // Möbius inversion on divisors: the sums for characters of order
        // dividing e are kernel(e) times the indicator of the image of [e].
        let mut exact = if entry.in_image {
            i128::try_from(entry.kernel).map_err(|_| overflow())?
        } else {
            0
        };
        for &(proper, weight) in &weights {
            if entry.degree.is_multiple_of(proper) {
                exact = exact.checked_sub(weight).ok_or_else(overflow)?;
            }
        }
        weights.push((entry.degree, exact));
        if exact != 0 {
            let k = n / entry.degree;
            let contribution = q
                .checked_pow((n - k) as u32)
                .and_then(|factor| factor.checked_mul(hilbert[k]))
                .and_then(|factor| factor.checked_mul(exact.unsigned_abs()))
                .ok_or_else(overflow)?;
            let sum = if exact > 0 {
                &mut positive
            } else {
                &mut negative
            };
            *sum = sum.checked_add(contribution).ok_or_else(overflow)?;
        }
    }
    let numerator = positive.checked_sub(negative).ok_or_else(integrity)?;
    exact_divide(numerator, order)
}

pub(crate) fn level(
    model: &mut Model,
    basis: &Basis,
    base: &Curve,
    degree: u128,
    translation: Option<&Divisor>,
) -> Result<(u128, u128)> {
    require_level_basis(basis)?;
    basis.require_base_extension(model, base)?;
    let coordinates = translation
        .map(|t| basis.translation_coordinates(model, t))
        .transpose()?;
    let ell = u128::from(basis.prime.get());
    let (q, hilbert) = abelian_counts(model, base, ell, degree)?;
    let n = ell as usize;
    let fixed = torsion_fixed(model, basis, base, degree as usize, coordinates.as_deref())?;
    // The correction is integral because H_1 is already the group order.
    let main = exact_divide(hilbert[n], hilbert[1])?;
    let tate = q.checked_pow((n - 1) as u32).ok_or_else(overflow)?;
    let count = if fixed == 0 {
        main.checked_sub(tate).ok_or_else(integrity)?
    } else {
        let correction = tate.checked_mul(fixed - 1).ok_or_else(overflow)?;
        main.checked_add(correction).ok_or_else(overflow)?
    };
    Ok((count, fixed))
}

fn require_level_basis(basis: &Basis) -> Result<()> {
    if basis.len() != 4 {
        return Err(Error::plain(
            "generalized Kummer level requires a genus-two torsion basis",
        ));
    }
    Ok(())
}

fn torsion_fixed(
    model: &mut Model,
    basis: &Basis,
    base: &Curve,
    degree: usize,
    coordinates: Option<&[u32]>,
) -> Result<u128> {
    let frobenius = basis.frobenius(model, base, degree)?;
    affine_fixed(
        model,
        basis.prime,
        |row, column| *frobenius.matrix.entry(row, column) as u32,
        coordinates,
    )
}

fn affine_fixed(
    model: &mut Model,
    prime: Prime,
    entry: impl Fn(usize, usize) -> u32,
    coordinates: Option<&[u32]>,
) -> Result<u128> {
    let mut image = prime_linear::Basis::default();
    // The fixed points of x ↦ Lx + t solve (I − L)x = t.
    for column in 0..4 {
        let vector = (0..4)
            .filter_map(|row| {
                let value = prime.add(u32::from(row == column), prime.negate(entry(row, column)));
                (value != 0).then_some((row, value))
            })
            .collect();
        image.insert(model, prime, vector)?;
    }
    if let Some(coordinates) = coordinates {
        let translation = coordinates
            .iter()
            .copied()
            .enumerate()
            .filter(|&(_, value)| value != 0)
            .collect();
        if !image.reduce(model, prime, translation)?.is_empty() {
            return Ok(0);
        }
    }
    Ok(u128::from(prime.get()).pow((4 - image.rank()) as u32))
}

pub(crate) fn require_translation_order(n: u128, h: u128) -> Result<()> {
    if n == 0 {
        return Err(Error::plain("generalized Kummer length must be positive"));
    }
    if h == 0 || !n.is_multiple_of(h) {
        return Err(Error::plain(
            "translation order must be a positive divisor of the generalized Kummer length",
        ));
    }
    Ok(())
}

pub(crate) fn translation_fixed_locus(
    model: &mut Model,
    n: u128,
    h: u128,
) -> Result<(u128, u128, u128)> {
    require_translation_order(n, h)?;
    let overflow =
        || Error::plain("translation fixed-locus invariant exceeds unsigned 128-bit range");
    let n_cubed = n.checked_pow(3).ok_or_else(overflow)?;
    let components = h.checked_pow(3).ok_or_else(overflow)?;
    let length = n / h;
    let root = length.isqrt();
    let work = usize::try_from(root)
        .map_err(|_| Error::plain("translation fixed-locus divisor sum exceeds the work limit"))?;
    model.spend(work)?;
    let mut sigma = 0u128;
    for d in 1..=root {
        if length.is_multiple_of(d) {
            let complement = length / d;
            sigma = sigma.checked_add(d).ok_or_else(overflow)?;
            if d != complement {
                sigma = sigma.checked_add(complement).ok_or_else(overflow)?;
            }
        }
    }
    let euler = n_cubed.checked_mul(sigma).ok_or_else(overflow)?;
    Ok((components, length - 1, euler))
}

struct DivisorStatistic {
    degree: usize,
    kernel: u128,
    in_image: bool,
}

fn collect_statistics(
    model: &mut Model,
    curve: &Curve,
    order: u128,
    target: Option<&Divisor>,
    statistics: &mut [DivisorStatistic],
) -> Result<()> {
    let mut seen = 0;
    jacobian::enumerate(model, curve, |model, divisor| {
        seen += 1;
        for entry in &mut *statistics {
            let multiple = divisor.scale(model, curve, entry.degree as i128)?;
            entry.kernel += u128::from(multiple.u.len() == 1);
            if !entry.in_image {
                entry.in_image = match target {
                    Some(target) => multiple.compare(model, target)?.is_eq(),
                    None => true,
                };
            }
        }
        Ok(())
    })?;
    if seen == order {
        Ok(())
    } else {
        Err(Error::plain(
            "enumerated Jacobian order disagrees with the Weil polynomial",
        ))
    }
}

fn abelian_counts(
    model: &mut Model,
    curve: &Curve,
    length: u128,
    degree: u128,
) -> Result<(u128, Vec<u128>)> {
    let q = u128::from(curve.ring.size_of_residue_field()?);
    let r = u32::try_from(degree).map_err(|_| overflow())?;
    let n = u32::try_from(length).map_err(|_| overflow())?;
    let horizon = r.checked_mul(n).ok_or_else(overflow)?;
    // The largest field used in the Weil recursion must support its square.
    let max_field = q.checked_pow(horizon).ok_or_else(overflow)?;
    max_field.checked_mul(max_field).ok_or_else(overflow)?;
    let field = q.checked_pow(r).ok_or_else(overflow)?;
    let weil = curve.weil(model)?;
    let mut points = Vec::with_capacity(n as usize);
    let mut order = 1u128;
    for i in 1..=n {
        order = order.checked_mul(field).ok_or_else(overflow)?;
        let (a, b) = powered_weil_coefficients(model, weil, (r * i) as usize, overflow)?;
        let correction = (order as i128 + 1)
            .checked_mul(a)
            .and_then(|term| b.checked_sub(term))
            .and_then(|term| term.checked_add(1))
            .ok_or_else(overflow)?;
        let count = (order * order)
            .checked_add_signed(correction)
            .ok_or_else(overflow)?;
        if count == 0 {
            return Err(integrity());
        }
        points.push(count);
    }
    Ok((field, hilbert_counts(model, field, &points)?))
}

pub(crate) fn hilbert_counts(model: &mut Model, q: u128, abelian: &[u128]) -> Result<Vec<u128>> {
    let n = abelian.len();
    if n.saturating_add(1).saturating_mul(4) > crate::complex::MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "Jacobian Hilbert-scheme series exceeds the storage limit",
        ));
    }
    model.spend(n + 1)?;
    let mut ghosts = vec![0u128; n + 1];
    for m in 1..=n {
        for j in (m..=n).step_by(m) {
            model.spend(1)?;
            let s = j / m;
            let term = q
                .checked_pow((j - s) as u32)
                .and_then(|factor| factor.checked_mul(abelian[s - 1]))
                .and_then(|factor| factor.checked_mul(m as u128))
                .ok_or_else(overflow)?;
            ghosts[j] = ghosts[j].checked_add(term).ok_or_else(overflow)?;
        }
    }
    let mut result = vec![1u128];
    for i in 1..=n {
        model.spend(i)?;
        let mut numerator = 0u128;
        for j in 1..=i {
            let term = ghosts[j].checked_mul(result[i - j]).ok_or_else(overflow)?;
            numerator = numerator.checked_add(term).ok_or_else(overflow)?;
        }
        result.push(exact_divide(numerator, i as u128)?);
    }
    Ok(result)
}

fn exact_divide(numerator: u128, denominator: u128) -> Result<u128> {
    if numerator.is_multiple_of(denominator) {
        Ok(numerator / denominator)
    } else {
        Err(integrity())
    }
}

fn integrity() -> Error {
    Error::plain("Jacobian Hilbert-scheme calculation failed its exact count check")
}

fn overflow() -> Error {
    Error::plain("Jacobian Hilbert-scheme count exceeds the checked 128-bit arithmetic range")
}
