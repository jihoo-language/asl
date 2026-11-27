use super::{top_pair, Algebra, Form, Terms, EVEN_RANK};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::nilpotent::Limit;
use crate::hyperelliptic_kummer_surface::qz::Scalar;
use crate::hyperelliptic_kummer_surface::rational_sparse::{Integer, Magnitude, WORK_LIMBS};
use crate::model::Model;
use crate::syntax::{Label, ParsedExpr};

pub(crate) const FIBER_DEGREE: usize = EVEN_RANK;

pub(crate) struct Specification {
    pub(crate) name: Label,
    pub(crate) cover: i128,
    pub(crate) bounds: [u128; 2],
    pub(crate) inputs: Box<[ParsedExpr; 4]>,
}

pub(crate) struct Witness<'a> {
    pub(crate) class: &'a Form,
    pub(crate) preimage: &'a Form,
    pub(crate) bound: u128,
    pub(crate) degree: usize,
}

impl Algebra {
    pub(crate) fn torsion_link(
        &self,
        model: &mut Model,
        limit: &Limit,
        cover: i128,
        left: Witness<'_>,
        right: Witness<'_>,
    ) -> Result<Scalar> {
        if left.bound == 0 || right.bound == 0 {
            return Err(Error::plain(
                "Kuga–Satake linking witness bounds must be positive",
            ));
        }
        if cover == 0 && (!left.class.terms().is_empty() || !right.class.terms().is_empty()) {
            return Err(Error::plain(
                "identity monodromy at m=0 has only the zero Wang torsion class",
            ));
        }
        self.check_linking_witness(model, limit, cover, &left, "left")?;
        self.check_linking_witness(model, limit, cover, &right, "right")?;
        let terms = left
            .preimage
            .terms()
            .iter()
            .map(|&(mask, value)| Ok((mask, value.resize()?)));
        let numerator = top_pair(model, terms, right.class)?;
        Scalar::from_wide(model, numerator, Magnitude::from_small(left.bound))
    }

    fn check_linking_witness(
        &self,
        model: &mut Model,
        limit: &Limit,
        cover: i128,
        witness: &Witness<'_>,
        side: &str,
    ) -> Result<()> {
        let mut difference = if cover == 0
            || matches!(witness.degree, 0 | EVEN_RANK)
            || witness.preimage.terms().is_empty()
        {
            Terms::new()
        } else {
            self.action(model, limit)?.first_difference(
                model,
                cover,
                witness.preimage,
                witness.degree > EVEN_RANK / 2,
            )?
        };
        model.spend(witness.class.terms().len() * WORK_LIMBS)?;
        for &(mask, value) in witness.class.terms() {
            let value: Integer = value.resize()?;
            let magnitude = value
                .magnitude
                .multiply_u128(witness.bound)
                .ok_or_else(Form::arithmetic)?;
            let value = Integer {
                magnitude,
                negative: value.negative,
            };
            Form::accumulate(&mut difference, mask, value.negate())?;
        }
        if difference.is_empty() {
            Ok(())
        } else {
            Err(Error::plain(format!("{side} Kuga–Satake linking witness must satisfy (T^m−1) times preimage = bound times class")))
        }
    }
}

pub(crate) fn certificate(value: Scalar, bounds: [u128; 2]) -> String {
    let [left, right] = bounds.map(|bound| value.order_matches(bound));
    format!(
        "rvKSLinkCertificate = value={}, pairing_order={}, left_exact={left}, right_exact={right}",
        value.text(),
        value.order(),
    )
}
