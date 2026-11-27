use std::collections::BTreeMap;
use std::rc::Rc;

use super::smith::exterior_dimension;
use super::{wedge_negative, Algebra, Form, Terms, EVEN_RANK};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::nilpotent::Limit;
use crate::hyperelliptic_kummer_surface::rational_sparse::{
    Coefficients, Domain, Integer, Magnitude, WORK_LIMBS,
};
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

pub(crate) const TOP_DEGREE: usize = EVEN_RANK + 1;
const ETA: u128 = 1 << EVEN_RANK;
type CochainTerms = BTreeMap<u128, Integer>;

pub(crate) enum Expression {
    Fiber(Label, i128, Box<ParsedExpr>),
    Eta(Label, i128),
    Cup(Box<[ParsedExpr; 2]>),
    Differential(Box<ParsedExpr>),
    Scale(i128, u128, Box<ParsedExpr>),
    Massey(Box<[ParsedExpr; 5]>),
    EtaMassey(Box<ParsedExpr>),
}

pub(crate) enum Command {
    Closed(Box<ParsedExpr>),
    CohomEqual(Box<[ParsedExpr; 3]>),
    Coefficient(u64, bool, Box<ParsedExpr>),
    Certificate(Box<ParsedExpr>),
    Data(Label, i128, usize),
}

pub(crate) enum Nomizu {}

impl Domain for Nomizu {
    const NAME: &'static str = "Kuga–Satake Nomizu cochain";
}

pub(crate) type Cochain = Coefficients<Nomizu, u128>;

impl Cochain {
    pub(crate) fn eta() -> Self {
        Self::singleton(ETA, 1)
    }

    pub(crate) fn fiber(model: &mut Model, input: &Form) -> Result<Self> {
        model.spend(input.terms().len() * WORK_LIMBS)?;
        let terms = input
            .terms()
            .iter()
            .map(|&(mask, value)| Ok((u128::from(mask), value.resize()?)))
            .collect::<Result<CochainTerms>>()?;
        Self::from_terms(model, &terms, Magnitude::from_small(1))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        self.charge_storage(model)?;
        Ok(Value::KSNomizu(Rc::new(self)))
    }

    pub(crate) fn entry(&self, model: &mut Model, mask: u64, eta: bool) -> Result<String> {
        let key = u128::from(mask) | if eta { ETA } else { 0 };
        self.coefficient(model, key)
    }

    pub(crate) fn cup(&self, model: &mut Model, right: &Self) -> Result<Self> {
        if self.terms().is_empty() || right.terms().is_empty() {
            return Ok(Self::zero());
        }
        let mut output = CochainTerms::new();
        for &(a, coefficient) in self.terms() {
            let coefficient: Integer = coefficient.resize()?;
            let mask = a as u64;
            let odd = !mask.count_ones().is_multiple_of(2);
            model.spend(right.terms().len() * (WORK_LIMBS + EVEN_RANK))?;
            for &(b, value) in right.terms() {
                if a & b != 0 {
                    continue;
                }
                let mut product = coefficient.multiply_wide(value.resize()?)?;
                if wedge_negative(mask, b as u64) ^ (odd && b & ETA != 0) {
                    product = product.negate();
                }
                Self::accumulate(&mut output, a | b, product)?;
            }
        }
        Self::from_terms(model, &output, self.product_denominator(right)?)
    }

    fn from_terms(model: &mut Model, terms: &CochainTerms, denominator: Magnitude) -> Result<Self> {
        Self::normalized(
            model,
            terms.iter().map(|(&key, value)| (key, value)),
            denominator,
        )
    }
}

impl Algebra {
    pub(crate) fn nomizu_differential(
        &self,
        model: &mut Model,
        limit: &Limit,
        m: i128,
        input: &Cochain,
    ) -> Result<Cochain> {
        if m == 0 {
            return Ok(Cochain::zero());
        }
        let image = self.horizontal_derivative(model, limit, input)?;
        model.spend(image.len() * WORK_LIMBS)?;
        let output = image
            .into_iter()
            .map(|(mask, value)| Ok((ETA | u128::from(mask), value.multiply(m)?)))
            .collect::<Result<CochainTerms>>()?;
        Cochain::from_terms(model, &output, input.denominator_wide()?)
    }

    pub(crate) fn nomizu_is_closed(
        &self,
        model: &mut Model,
        limit: &Limit,
        m: i128,
        input: &Cochain,
    ) -> Result<bool> {
        Ok(m == 0 || self.horizontal_derivative(model, limit, input)?.is_empty())
    }

    pub(crate) fn nomizu_cohom_equal(
        &self,
        model: &mut Model,
        limit: &Limit,
        m: i128,
        inputs: [&Cochain; 3],
    ) -> Result<bool> {
        let [left, right, preimage] = inputs;
        for input in [left, right] {
            if !self.nomizu_is_closed(model, limit, m, input)? {
                return Err(Error::plain(
                    "rvKSNomizuCohomEqual requires closed cochains",
                ));
            }
        }
        let differential = self.nomizu_differential(model, limit, m, preimage)?;
        let corrected = right.add(model, &differential)?;
        Ok(left.compare(model, &corrected)?.is_eq())
    }

    pub(crate) fn nomizu_massey(
        &self,
        model: &mut Model,
        limit: &Limit,
        m: i128,
        degree: i128,
        inputs: [&Cochain; 5],
    ) -> Result<Cochain> {
        let [a, b, c, x, y] = inputs;
        for input in [a, b, c] {
            if !self.nomizu_is_closed(model, limit, m, input)? {
                return Err(Error::plain(
                    "rvKSNomizuMassey requires three closed inputs",
                ));
            }
        }
        for (left, right, primitive, side) in [(a, b, x, "first"), (b, c, y, "second")] {
            let product = left.cup(model, right)?;
            let differential = self.nomizu_differential(model, limit, m, primitive)?;
            if !product.compare(model, &differential)?.is_eq() {
                return Err(Error::plain(format!(
                    "rvKSNomizuMassey {side} nullhomotopy has the wrong differential"
                )));
            }
        }
        let ay = a.cup(model, y)?;
        let xc = x.cup(model, c)?;
        let xc = xc.signed_map(model, |_| degree.rem_euclid(2) == 0)?;
        ay.add(model, &xc)
    }

    pub(crate) fn nomizu_eta_massey(
        &self,
        model: &mut Model,
        limit: &Limit,
        m: i128,
        f: &Cochain,
    ) -> Result<Cochain> {
        self.check_eta_massey(model, limit, m, f)?;
        let cochains = f
            .terms()
            .iter()
            .map(|&(key, entry)| {
                let value: Integer = entry.resize()?;
                Ok((key | ETA, if m < 0 { value.negate() } else { value }))
            })
            .collect::<Result<CochainTerms>>()?;
        let denominator = f
            .denominator_wide()?
            .multiply_u128(m.unsigned_abs())
            .ok_or_else(Cochain::arithmetic)?;
        Cochain::from_terms(model, &cochains, denominator)
    }

    pub(crate) fn nomizu_certificate(
        &self,
        model: &mut Model,
        limit: &Limit,
        m: i128,
        f: &Cochain,
    ) -> Result<String> {
        self.check_eta_massey(model, limit, m, f)?;
        // On degree two, the detector D applied to the eta coefficient
        // kills boundaries and both full indeterminacy summands. It sends
        // the representative eta*f/m to the verified nonzero D(f)/m.
        Ok(
            "rvKSNomizuMasseyCertificate = triple=<eta,eta,Df>, degree=2, excludes_zero=true"
                .into(),
        )
    }

    fn check_eta_massey(
        &self,
        model: &mut Model,
        limit: &Limit,
        m: i128,
        f: &Cochain,
    ) -> Result<()> {
        if m == 0 {
            return Err(Error::plain(
                "the eta Massey certificate requires nonzero monodromy exponent",
            ));
        }
        if f.terms()
            .iter()
            .any(|&(key, _)| key & ETA != 0 || key.count_ones() != 1)
        {
            return Err(Error::plain(
                "the eta Massey certificate requires a pure fiber one-form",
            ));
        }
        if self.horizontal_derivative(model, limit, f)?.is_empty() {
            return Err(Error::plain(
                "the eta Massey certificate requires D(f) to be nonzero",
            ));
        }
        Ok(())
    }

    fn horizontal_derivative(
        &self,
        model: &mut Model,
        limit: &Limit,
        input: &Cochain,
    ) -> Result<Terms> {
        model.spend(input.terms().len() * WORK_LIMBS)?;
        let horizontal = input
            .terms()
            .iter()
            .filter(|&&(key, _)| key & ETA == 0)
            .map(|&(key, value)| Ok((key as u64, value.resize()?)))
            .collect::<Result<Terms>>()?;
        match horizontal.first_key_value() {
            None | Some((&0, _)) | Some((&u64::MAX, _)) => Ok(Terms::new()),
            Some(_) => self.action(model, limit)?.derivative(model, &horizontal),
        }
    }

    pub(crate) fn nomizu_data(
        &self,
        model: &mut Model,
        label: &str,
        limit: &Limit,
        m: i128,
        degree: usize,
    ) -> Result<String> {
        if degree > TOP_DEGREE {
            return Err(Error::plain(
                "Kuga–Satake Nomizu degree must be between 0 and 65",
            ));
        }
        let horizontal = exterior_dimension(degree);
        let kernel = self.rational_fiber_invariants(model, limit, m, degree)?;
        let (vertical, image) = if let Some(previous) = degree.checked_sub(1) {
            let vertical = exterior_dimension(previous);
            let image = vertical - self.rational_fiber_invariants(model, limit, m, previous)?;
            (vertical, image)
        } else {
            (0, 0)
        };
        let cochains = horizontal + vertical;
        let closed = kernel + vertical;
        let cohomology = closed - image;
        Ok(format!("rvKSNomizuData[{label},{m},{degree}] = cochains={cochains}, closed={closed}, exact={image}, cohomology={cohomology}"))
    }
}
