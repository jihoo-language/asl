use super::{
    check_report, differential_matrix, BaseHom, BaseKey, BaseMorphism, Cochain, Error, Hom, Key,
    Model, Rational, Result,
};

pub(crate) fn include(model: &mut Model, value: &BaseMorphism) -> Result<Cochain> {
    let mut result = Cochain::default();
    for (&key, &coefficient) in &value.terms {
        for q in [false, true] {
            result.insert(
                model,
                Key {
                    row: key.row,
                    column: key.column,
                    p_out: key.slot,
                    q_out: q,
                    q_odd: false,
                    exponents: [key.exponent, 0, 0],
                },
                coefficient,
            )?;
        }
    }
    Ok(result)
}

pub(crate) fn project(model: &mut Model, value: &Cochain) -> Result<BaseMorphism> {
    let mut result = BaseMorphism::default();
    for (&key, &coefficient) in &value.terms {
        if !key.q_odd && !key.q_out && key.exponents[1..] == [0, 0] {
            result.insert(
                model,
                BaseKey {
                    row: key.row,
                    column: key.column,
                    slot: key.p_out,
                    exponent: key.exponents[0],
                },
                coefficient,
            )?;
        }
    }
    Ok(result)
}

impl Hom<'_> {
    pub(crate) fn homotopy(&self, model: &mut Model, value: &Cochain) -> Result<Cochain> {
        let mut output = Cochain::default();
        for (&key, &coefficient) in &value.terms {
            let [x, u, v] = key.exponents;
            let mut result = key;
            let negative_q = match (key.q_odd, key.q_out) {
                (false, false) if u > 0 => {
                    result.q_odd = true;
                    result.exponents = [x, u - 1, v];
                    false
                }
                (false, false) if v > 0 => {
                    result.q_odd = true;
                    result.q_out = true;
                    result.exponents = [x, 0, v - 1];
                    false
                }
                (true, true) if u > 0 => {
                    result.q_odd = false;
                    result.exponents = [x, u - 1, v];
                    true
                }
                _ => continue,
            };
            let negative = negative_q ^ self.odd ^ key.q_odd;
            output.insert(
                model,
                result,
                if negative {
                    coefficient.negate()
                } else {
                    coefficient
                },
            )?;
        }
        Ok(output)
    }

    pub(crate) fn class(&self, model: &mut Model, value: &Cochain) -> Result<Cochain> {
        self.require_cycle(model, value)?;
        let base = project(model, value)?;
        let reduced = self.base()?.class(model, &base)?;
        include(model, &reduced)
    }

    pub(crate) fn primitive(&self, model: &mut Model, value: &Cochain) -> Result<Cochain> {
        self.require_cycle(model, value)?;
        let base = project(model, value)?;
        let base_primitive = self.base()?.primitive(model, &base)?;
        let h = self.homotopy(model, value)?;
        let lifted = include(model, &base_primitive)?;
        let result = h.add(model, &lifted)?;
        self.verify_boundary(model, value, &result)?;
        Ok(result)
    }

    fn verify_boundary(
        &self,
        model: &mut Model,
        boundary: &Cochain,
        primitive: &Cochain,
    ) -> Result<()> {
        let inverse = Self {
            odd: !self.odd,
            ..*self
        };
        let check = inverse.differential(model, primitive)?;
        if !check.compare(model, boundary)?.is_eq() {
            return Err(Error::plain(
                "Knörrer homotopy failed its boundary equation",
            ));
        }
        Ok(())
    }

    pub(crate) fn cone(&self, model: &mut Model, value: &Cochain) -> Result<String> {
        if self.odd {
            return Err(Error::plain("Knörrer cone requires an even cochain"));
        }
        self.require_cycle(model, value)?;
        let projected = project(model, value)?;
        let included = include(model, &projected)?;
        let negative = included.scale(model, (-1, 1))?;
        let difference = value.add(model, &negative)?;
        let h = self.homotopy(model, value)?;
        // The upper-right equation for the cone isomorphism with off-diagonal h
        // is D_T h + I(projected) = value - h D_S.
        self.verify_boundary(model, &difference, &h)?;
        let (_, base) =
            BaseHom::new(self.source, self.target, false)?.cone(model, &projected, true)?;
        let report = format!(
            "homotopy={} base_cone={} tensor_cone=verified",
            h.text(model, true)?,
            base.unwrap()
        );
        check_report(&report)?;
        Ok(report)
    }

    pub(crate) fn trace(&self, model: &mut Model, value: &Cochain) -> Result<String> {
        self.require_cycle(model, value)?;
        if !self.odd {
            return Ok(String::from("0"));
        }
        let matrix = differential_matrix(model, self.source)?;
        let derivatives = [
            differentiate(model, &matrix, 0)?,
            differentiate(model, &matrix, 1)?,
            differentiate(model, &matrix, 2)?,
        ];
        let mut alternating = Cochain::default();
        for (permutation, negative) in [
            ([0, 1, 2], false),
            ([1, 2, 0], false),
            ([2, 0, 1], false),
            ([0, 2, 1], true),
            ([1, 0, 2], true),
            ([2, 1, 0], true),
        ] {
            let [a, b, c] = permutation;
            let pair = derivatives[a].compose(model, true, &derivatives[b], true)?;
            let term = pair.compose(model, false, &derivatives[c], true)?;
            for (key, coefficient) in term.terms {
                alternating.insert(
                    model,
                    key,
                    if negative {
                        coefficient.negate()
                    } else {
                        coefficient
                    },
                )?;
            }
        }
        let composed = value.compose(model, true, &alternating, true)?;
        let mut trace = Rational::zero();
        for (key, coefficient) in composed.terms {
            if key.row == key.column
                && !key.q_odd
                && key.exponents == [(self.source.d - 2) as u128, 0, 0]
            {
                let coefficient = if key.p_out ^ key.q_out {
                    coefficient.negate()
                } else {
                    coefficient
                };
                trace = trace.add(model, coefficient)?;
            }
        }
        // The ordered residue (d x^(d-1), v, u) contributes a minus sign.
        let normalization = Rational::integer(-6 * self.source.d as i128).inverse()?;
        trace.multiply(model, normalization)?.text(model)
    }
}

fn differentiate(model: &mut Model, value: &Cochain, variable: usize) -> Result<Cochain> {
    let mut result = Cochain::default();
    for (&key, &coefficient) in &value.terms {
        let exponent = key.exponents[variable];
        if exponent != 0 {
            let mut derivative = key;
            derivative.exponents[variable] -= 1;
            let value = coefficient.multiply(model, Rational::integer(exponent as i128))?;
            result.insert(model, derivative, value)?;
        }
    }
    Ok(result)
}
