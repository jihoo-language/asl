use std::rc::Rc;

use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::jacobian_level::Basis;
use crate::model::Model;
use crate::syntax::Label;

use super::{
    k3,
    nikulin::{basis_for, NodePermutation},
};

pub(crate) mod points;

const NODES: usize = 16;
const DIVISORS: usize = 1 + NODES;
const INVARIANT_RANK: usize = 10;

pub(crate) enum Query {
    Points(u128),
    Zeta,
}

pub(crate) struct Command(
    pub(crate) Label,
    pub(crate) Label,
    pub(crate) usize,
    pub(crate) Query,
);

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<String> {
        let Self(name, base, mask, query) = self;
        let (curve, action) = Action::resolve(model, &name, &base, mask)?;
        let q = curve.ring.size_of_residue_field()?;
        let name = name.text;
        let base = base.text;
        match query {
            Query::Points(degree) => {
                let trace = action.trace(model, degree)?;
                let total = point_count(q, degree, trace)?;
                let orbits = action.frobenius.orbit_lengths();
                let weil = curve.weil(model)?;
                let (_, _, cover) = k3::point_counts(model, q, weil, &orbits, degree)?;
                if !cover.is_multiple_of(2) {
                    return Err(Error::plain(
                        "a free Enriques switch must pair the rational cover points",
                    ));
                }
                let split = cover / 2;
                let nonsplit = total.checked_sub(split).ok_or_else(|| {
                    Error::plain("Enriques point count is smaller than its split fiber count")
                })?;
                let twisted = nonsplit.checked_mul(2).ok_or_else(point_overflow)?;
                Ok(format!("rvKummerEnriquesPoints[{name}, {base}, {mask}, {degree}] = points={total}, normalized trace={trace}, split={split}, nonsplit={nonsplit}, twisted cover={twisted}"))
            }
            Query::Zeta => {
                let polynomial = action.polynomial(model)?;
                let square = u128::from(q) * u128::from(q);
                Ok(format!("rvKummerEnriquesZeta[{name}, {base}, {mask}] = 1/((1-T)R({q}*T)(1-{square}*T)), R={polynomial:?}"))
            }
        }
    }
}

pub(super) struct Action {
    frobenius: NodePermutation,
    twice_switch: [[i8; DIVISORS]; DIVISORS],
}

impl Action {
    fn basis(model: &Model, name: &Label) -> Result<Rc<Basis>> {
        basis_for(
            model,
            name,
            "Enriques quotients require a certified genus-two level-2 basis",
        )
    }

    pub(super) fn resolve(
        model: &mut Model,
        name: &Label,
        base: &Label,
        mask: usize,
    ) -> Result<(Rc<Curve>, Self)> {
        let basis = Self::basis(model, name)?;
        let curve = model
            .kummer_curve(&base.text)
            .map_err(|error| error.at(base.span))?;
        let action = Self::validate(model, &basis, &curve, mask)?;
        Ok((curve, action))
    }

    fn validate(model: &mut Model, basis: &Basis, curve: &Curve, mask: usize) -> Result<Self> {
        let frobenius = NodePermutation::new(model, basis, curve, 1)?;
        let mut odd = 0u16;
        for node in 0..NODES {
            if basis.binary_combination(model, node as u8)?.u.len() <= 2 {
                odd |= 1 << node;
            }
        }
        if odd.count_ones() != 6 || frobenius.on_truth_table(odd) != odd {
            return Err(Error::plain(
                "expected six base-Frobenius-stable nodes on the odd trope",
            ));
        }
        if odd & (1 << mask) != 0 {
            return Err(Error::plain(
                "Enriques quotient requires an even theta characteristic",
            ));
        }
        if usize::from(frobenius.image(mask)) != mask {
            return Err(Error::plain(
                "Enriques switch mask must be fixed by base Frobenius",
            ));
        }
        Self::new(model, frobenius, odd, mask)
    }

    fn new(model: &mut Model, frobenius: NodePermutation, odd: u16, mask: usize) -> Result<Self> {
        // Use (H, E_0, ..., E_15) on the divisor cohomology, with the Tate
        // factor removed from Frobenius. Store 2σ so every entry is integral.
        model.spend(DIVISORS * DIVISORS)?;
        let mut twice_switch = [[0; DIVISORS]; DIVISORS];
        twice_switch[0][0] = 6;
        for node in 0..NODES {
            // σ H = 3H − Σ E_v, and σ E_v = T_{v+a}.
            twice_switch[node + 1][0] = -2;
            twice_switch[0][node + 1] = 1;
            for (other, row) in twice_switch[1..].iter_mut().enumerate() {
                // 2T_b = H − Σ_{v⊕b ∈ S0} E_v.
                row[node + 1] = -i8::from(odd & (1 << (other ^ node ^ mask)) != 0);
            }
        }
        let action = Self {
            frobenius,
            twice_switch,
        };
        // Pullback identifies the invariant cohomology with H² of the quotient.
        // Its rank is ten here, so the divisor span accounts for all of H².
        if action.trace(model, 0)? != INVARIANT_RANK as i128 {
            return Err(Error::plain(
                "Enriques invariant divisor space must have rank ten",
            ));
        }
        Ok(action)
    }

    fn trace(&self, model: &mut Model, power: u128) -> Result<i128> {
        model.spend(NODES * (2 * NODES + 1))?;
        let mut fixed = 0;
        let mut twisted = i128::from(self.twice_switch[0][0]);
        for node in 0..NODES {
            let image = self.frobenius.power_image(node, power);
            fixed += i128::from(node == image);
            twisted += i128::from(self.twice_switch[node + 1][image + 1]);
        }
        // Tr F^r = 1+f on the full span. A switched diagonal entry is −1/2
        // precisely when v⊕F^r(v)⊕a belongs to S0. Thus 2 Tr(F^r σ)=6−m,
        // and the invariant projector (1+σ)/2 gives (8+2f−m)/4.
        let numerator = 2 * (1 + fixed) + twisted;
        if numerator % 4 != 0 {
            return Err(Error::plain(
                "Enriques invariant trace division is not exact",
            ));
        }
        Ok(numerator / 4)
    }

    fn polynomial(&self, model: &mut Model) -> Result<[i128; INVARIANT_RANK + 1]> {
        let mut traces = [0; INVARIANT_RANK];
        for (i, trace) in traces.iter_mut().enumerate() {
            *trace = self.trace(model, i as u128 + 1)?;
        }
        let mut coefficients = [0; INVARIANT_RANK + 1];
        coefficients[0] = 1;
        for degree in 1..=INVARIANT_RANK {
            model.spend(degree)?;
            let numerator = -coefficients[..degree]
                .iter()
                .zip(traces[..degree].iter().rev())
                .map(|(coefficient, trace)| coefficient * trace)
                .sum::<i128>();
            let denominator = degree as i128;
            if numerator % denominator != 0 {
                return Err(Error::plain("Enriques zeta Newton division is not exact"));
            }
            coefficients[degree] = numerator / denominator;
        }
        Ok(coefficients)
    }
}

fn point_overflow() -> Error {
    Error::plain("Enriques point count exceeds unsigned 128-bit range")
}

fn point_count(q: u64, degree: u128, trace: i128) -> Result<u128> {
    let exponent = u32::try_from(degree).map_err(|_| point_overflow())?;
    let field_order = u128::from(q)
        .checked_pow(exponent)
        .ok_or_else(point_overflow)?;
    let square = field_order
        .checked_mul(field_order)
        .ok_or_else(point_overflow)?;
    let signed_order = i128::try_from(field_order).map_err(|_| point_overflow())?;
    let correction = signed_order
        .checked_mul(trace)
        .and_then(|trace| trace.checked_add(1));
    square
        .checked_add_signed(correction.ok_or_else(point_overflow)?)
        .ok_or_else(point_overflow)
}
