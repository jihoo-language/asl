use std::cell::OnceCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::binary_linear::check_storage;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::rv_as_odd_cartier::Cartier;
use crate::{rv_cyclotomic as cyclo, rv_newton};

type Character = Rc<[Vec<i128>]>;

#[derive(Clone, Copy)]
pub(crate) enum Query {
    Reduced,
    Basis,
    CartierRank(usize),
    ANumber,
    Genus,
    PRank,
    Character(usize),
    CharacterNewton(usize),
    LPolynomial,
    Newton,
    Points(usize),
}

impl Query {
    fn name(self) -> &'static str {
        match self {
            Self::Reduced => "rvASOddReduced",
            Self::Basis => "rvASOddBasis",
            Self::CartierRank(_) => "rvASOddCartierRank",
            Self::ANumber => "rvASOddaNumber",
            Self::Genus => "rvASOddGenus",
            Self::PRank => "rvASOddPrank",
            Self::Character(_) => "rvASOddCharacterL",
            Self::CharacterNewton(_) => "rvASOddCharacterNewton",
            Self::LPolynomial => "rvASOddLPolynomial",
            Self::Newton => "rvASOddNewton",
            Self::Points(_) => "rvASOddPoints",
        }
    }
}

pub(crate) struct Cover {
    prime: Prime,
    terms: Vec<(usize, u32)>,
    factor_degree: usize,
    genus: usize,
    cartier: OnceCell<Rc<Cartier>>,
    character: OnceCell<Character>,
    numerator: OnceCell<Rc<[i128]>>,
}

impl Cover {
    pub(crate) fn build(model: &mut Model, p: usize, terms: Vec<(usize, i128)>) -> Result<Self> {
        if p < 3 || p.is_multiple_of(2) {
            return Err(Error::plain("rvasodd requires an odd prime"));
        }
        cyclo::check_order(p)?;
        let prime = Prime::new(model, p)?;
        let mut normalized = BTreeMap::new();
        for (mut exponent, coefficient) in terms {
            let coefficient = coefficient.rem_euclid(i128::from(prime.get())) as u32;
            model.spend(1 + normalized.len().max(1).ilog2() as usize)?;
            // In F_p[x]/(G^p-G) as an additive quotient, x^(pe) and x^e
            // are equivalent. Constants survive an F_p-defined change of y.
            while exponent != 0 && exponent.is_multiple_of(p) {
                exponent /= p;
                model.spend(1)?;
            }
            let old = normalized.entry(exponent).or_insert(0);
            *old = prime.add(*old, coefficient);
        }
        let terms: Vec<_> = normalized.into_iter().filter(|&(_, c)| c != 0).collect();
        let d = terms.last().map_or(0, |&(exponent, _)| exponent);
        if d == 0 {
            return Err(Error::plain(
                "rvasodd requires a nonconstant reduced polynomial",
            ));
        }
        let factor_degree = d - 1;
        let genus = ((p - 1) / 2)
            .checked_mul(factor_degree)
            .ok_or_else(|| Error::plain("odd Artin–Schreier genus is out of range"))?;
        model.allocate(8 + 2 * terms.len())?;
        Ok(Self {
            prime,
            terms,
            factor_degree,
            genus,
            cartier: OnceCell::new(),
            character: OnceCell::new(),
            numerator: OnceCell::new(),
        })
    }

    pub(crate) fn prime(&self) -> Prime {
        self.prime
    }

    pub(crate) fn genus(&self) -> usize {
        self.genus
    }

    pub(crate) fn degree(&self) -> usize {
        self.factor_degree + 1
    }

    pub(crate) fn terms(&self) -> &[(usize, u32)] {
        &self.terms
    }

    pub(crate) fn cartier(&self, model: &mut Model) -> Result<Rc<Cartier>> {
        if let Some(cartier) = self.cartier.get() {
            return Ok(cartier.clone());
        }
        let cartier = Rc::new(Cartier::build(model, self)?);
        Ok(self.cartier.get_or_init(|| cartier).clone())
    }

    pub(crate) fn order(&self) -> usize {
        self.prime.get() as usize
    }

    fn checked_character(&self, character: usize) -> Result<()> {
        if !(1..self.order()).contains(&character) {
            return Err(Error::plain(
                "odd Artin–Schreier character index must lie in 1..p-1",
            ));
        }
        Ok(())
    }

    pub(crate) fn coefficient(&self, model: &mut Model, j: usize, k: usize) -> Result<Value> {
        self.checked_character(j)?;
        if k > self.factor_degree {
            return Ok(cyclo::zero());
        }
        if k == 0 {
            return cyclo::persist(model, vec![1]);
        }
        let factor = self.primitive(model)?;
        let coefficient = self.conjugate(model, j, &factor[k])?;
        cyclo::persist(model, coefficient)
    }

    fn conjugate(
        &self,
        model: &mut Model,
        character: usize,
        coefficient: &[i128],
    ) -> Result<Vec<i128>> {
        if character == 1 || coefficient.len() <= 1 {
            model.spend(coefficient.len())?;
            return Ok(coefficient.to_vec());
        }
        cyclo::substitute_root(model, self.order(), coefficient, character)
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str, query: Query) -> Result<String> {
        let operation = query.name();
        let value = match query {
            Query::Reduced => {
                model.spend(self.terms.len())?;
                format!("{:?}", self.terms)
            }
            Query::Basis => format!("{:?}", self.cartier(model)?.basis(model)?),
            Query::CartierRank(power) => {
                let rank = self.cartier(model)?.rank(model, self, power)?;
                return Ok(format!("{operation}[{name}, {power}] = {rank}"));
            }
            Query::ANumber => (self.genus - self.cartier(model)?.rank(model, self, 1)?).to_string(),
            Query::Genus => self.genus.to_string(),
            Query::PRank => "0".to_owned(),
            Query::Character(j) => {
                self.checked_character(j)?;
                let factor = self.primitive(model)?;
                let mut conjugate = Vec::with_capacity(factor.len());
                for coefficient in factor.iter() {
                    let mut image = self.conjugate(model, j, coefficient)?;
                    model.spend(self.order())?;
                    image.resize(self.order() - 1, 0);
                    conjugate.push(image);
                }
                return Ok(format!("{operation}[{name}, {j}] = {conjugate:?}"));
            }
            Query::CharacterNewton(j) => {
                self.checked_character(j)?;
                // The valuation extends uniquely to Q(zeta_p), so all
                // cyclotomic conjugates have the same Newton polygon.
                let factor = self.primitive(model)?;
                let valuations = factor
                    .iter()
                    .map(|a| uniformizer_order(model, self.prime, a))
                    .collect::<Result<Vec<_>>>()?;
                let slopes = rv_newton::scaled_slopes(model, valuations, self.order() - 1)?;
                return Ok(format!(
                    "{operation}[{name}, {j}] = [{}]",
                    slopes.join(", ")
                ));
            }
            Query::LPolynomial => format!("{:?}", self.l_polynomial(model)?),
            Query::Newton => {
                let polynomial = self.l_polynomial(model)?;
                let valuations = polynomial
                    .iter()
                    .map(|&c| p_order(model, self.prime, c))
                    .collect::<Result<Vec<_>>>()?;
                let slopes = rv_newton::scaled_slopes(model, valuations, 1)?;
                format!("[{}]", slopes.join(", "))
            }
            Query::Points(m) => {
                let points = self.points(model, m)?;
                return Ok(format!("{operation}[{name}, {m}] = {points}"));
            }
        };
        Ok(format!("{operation}[{name}] = {value}"))
    }

    fn primitive(&self, model: &mut Model) -> Result<Character> {
        if let Some(character) = self.character.get() {
            return Ok(character.clone());
        }
        let p = self.order();
        let n = self.factor_degree;
        let histograms = trace_histograms(model, self.prime, &self.terms, n)?;
        let mut sums = vec![Vec::new()];
        for histogram in histograms.into_iter().skip(1) {
            model.spend(p)?;
            let counts = histogram
                .into_iter()
                .map(|entries| entries as i128)
                .collect();
            sums.push(cyclo::reduced(model, p, counts)?);
        }
        let mut factor = vec![vec![1]];
        for k in 1..=n {
            let mut coefficient = Vec::new();
            for j in 1..=k {
                let product = cyclo::product_coefficients(model, p, &sums[j], &factor[k - j])?;
                cyclo::accumulate(model, &mut coefficient, &product, false)?;
            }
            model.spend(coefficient.len())?;
            for value in &mut coefficient {
                if *value % k as i128 != 0 {
                    return Err(Error::plain(
                        "odd Artin–Schreier character Newton division is not integral",
                    ));
                }
                *value /= k as i128;
            }
            factor.push(coefficient);
        }
        if factor[n].is_empty() {
            return Err(Error::plain(
                "odd Artin–Schreier character factor has the wrong degree",
            ));
        }
        model.allocate(factor.iter().map(|c| 3 + 2 * c.len()).sum())?;
        let factor: Character = factor.into();
        Ok(self.character.get_or_init(|| factor).clone())
    }

    fn l_polynomial(&self, model: &mut Model) -> Result<Rc<[i128]>> {
        if let Some(polynomial) = self.numerator.get() {
            return Ok(polynomial.clone());
        }
        let p = self.order();
        let leading = u32::try_from(self.genus)
            .ok()
            .and_then(|g| (p as i128).checked_pow(g))
            .ok_or_else(coefficient_error)?;
        let mut polynomial = vec![vec![1]];
        if self.genus != 0 {
            let factor = self.primitive(model)?;
            for j in 1..p {
                let conjugate = factor
                    .iter()
                    .map(|coefficient| self.conjugate(model, j, coefficient))
                    .collect::<Result<Vec<_>>>()?;
                polynomial = polynomial_product(model, p, &polynomial, &conjugate)?;
            }
        }
        let polynomial = polynomial
            .into_iter()
            .map(|coefficient| match coefficient.as_slice() {
                [] => Ok(0),
                &[c] => Ok(c),
                _ => Err(Error::plain(
                    "odd Artin–Schreier character norm did not descend to integers",
                )),
            })
            .collect::<Result<Vec<_>>>()?;
        if polynomial.len() != 2 * self.genus + 1 || polynomial.last() != Some(&leading) {
            return Err(Error::plain(
                "odd Artin–Schreier numerator disagrees with the curve genus",
            ));
        }
        model.allocate(2 * polynomial.len())?;
        Ok(self.numerator.get_or_init(|| polynomial.into()).clone())
    }

    fn points(&self, model: &mut Model, extension: usize) -> Result<u128> {
        if extension == 0 {
            return Err(Error::plain(
                "odd Artin–Schreier extension degree must be positive",
            ));
        }
        let p = self.order();
        let base = u32::try_from(extension)
            .ok()
            .and_then(|m| (p as u128).checked_pow(m))
            .and_then(|power| power.checked_add(1))
            .ok_or_else(points_error)?;
        if self.genus == 0 {
            return Ok(base);
        }
        let factor = self.primitive(model)?;
        let sum = power_trace(model, p, &factor, extension)?;
        // Q(zeta_p)/Q has trace p-1 on 1 and -1 on the other basis powers.
        let mut trace = sum
            .first()
            .copied()
            .unwrap_or(0)
            .checked_mul((p - 1) as i128)
            .ok_or_else(coefficient_error)?;
        for &c in sum.iter().skip(1) {
            trace = trace.checked_sub(c).ok_or_else(coefficient_error)?;
        }
        let result = if trace >= 0 {
            base.checked_add(trace as u128)
        } else {
            base.checked_sub(trace.unsigned_abs())
        };
        result.ok_or_else(points_error)
    }
}

fn power_trace(
    model: &mut Model,
    p: usize,
    factor: &[Vec<i128>],
    extension: usize,
) -> Result<Vec<i128>> {
    let degree = factor.len() - 1;
    check_storage((degree + 1).saturating_mul(3 + 2 * p))?;
    let mut window = vec![Vec::new(); degree];
    for power in 1..=extension {
        model.spend(1)?;
        let mut sum = Vec::new();
        if let Some(coefficient) = factor.get(power) {
            model.spend(coefficient.len())?;
            sum = coefficient
                .iter()
                .map(|&c| c.checked_mul(power as i128).ok_or_else(coefficient_error))
                .collect::<Result<Vec<_>>>()?;
        }
        // Negative Frobenius power trace S_k = k c_k - sum c_i S_(k-i).
        for (i, c) in factor
            .iter()
            .enumerate()
            .take(power.min(degree + 1))
            .skip(1)
        {
            let product = cyclo::product_coefficients(model, p, c, &window[(power - i) % degree])?;
            cyclo::accumulate(model, &mut sum, &product, true)?;
        }
        if power == extension {
            return Ok(sum);
        }
        window[power % degree] = sum;
    }
    unreachable!()
}

fn polynomial_product(
    model: &mut Model,
    p: usize,
    left: &[Vec<i128>],
    right: &[Vec<i128>],
) -> Result<Vec<Vec<i128>>> {
    let length = left.len() + right.len() - 1;
    check_storage((left.len() + right.len() + length).saturating_mul(3 + 2 * p))?;
    model.spend(left.len() * right.len())?;
    let mut product = vec![Vec::new(); length];
    for (i, x) in left.iter().enumerate() {
        for (j, y) in right.iter().enumerate() {
            let coefficient = cyclo::product_coefficients(model, p, x, y)?;
            cyclo::accumulate(model, &mut product[i + j], &coefficient, false)?;
        }
    }
    Ok(product)
}

fn p_order(model: &mut Model, prime: Prime, value: i128) -> Result<Option<i64>> {
    model.spend(1)?;
    if value == 0 {
        return Ok(None);
    }
    let p = i128::from(prime.get());
    let mut value = value;
    let mut order = 0;
    while value % p == 0 {
        model.spend(1)?;
        value /= p;
        order += 1;
    }
    Ok(Some(order))
}

fn uniformizer_order(model: &mut Model, prime: Prime, coefficient: &[i128]) -> Result<Option<i64>> {
    let mut common = None;
    for &value in coefficient {
        if let Some(order) = p_order(model, prime, value)? {
            common = Some(common.map_or(order, |old: i64| old.min(order)));
        }
    }
    let Some(common) = common else {
        return Ok(None);
    };
    let p = i128::from(prime.get());
    let divisor = p.checked_pow(common as u32).ok_or_else(coefficient_error)?;
    model.spend(coefficient.len())?;
    let mut residue: Vec<_> = coefficient
        .iter()
        .map(|&a| (a / divisor).rem_euclid(p) as u32)
        .collect();
    // Divide by x-1 over F_p until the remainder is nonzero. The degree
    // is less than p-1, so its first nonzero shifted coefficient gives
    // the residual uniformizer order without computing large binomials.
    let mut order = common * (i64::from(prime.get()) - 1);
    while !residue.is_empty() {
        model.spend(residue.len())?;
        for i in (1..residue.len()).rev() {
            residue[i - 1] = prime.add(residue[i - 1], residue[i]);
        }
        if residue[0] != 0 {
            return Ok(Some(order));
        }
        residue.rotate_left(1);
        residue.pop();
        order += 1;
    }
    Err(Error::plain(
        "odd Artin–Schreier normalized cyclotomic residue vanished",
    ))
}

fn coefficient_error() -> Error {
    Error::plain("odd Artin–Schreier cyclotomic coefficient is out of range")
}

fn points_error() -> Error {
    Error::plain("odd Artin–Schreier point count is out of range")
}

fn trace_histograms(
    model: &mut Model,
    prime: Prime,
    terms: &[(usize, u32)],
    bound: usize,
) -> Result<Vec<Vec<usize>>> {
    let p = prime.get() as usize;
    // Prepay at least one operation per monic candidate before building any
    // fields or search tables. Only degrees 1..d-1 are needed for the factor.
    check_storage((bound + 1).saturating_mul(p + 3))?;
    let mut order = 1usize;
    let mut orders = Vec::new();
    for _ in 1..=bound {
        order = order.checked_mul(p).ok_or_else(|| {
            Error::plain("odd Artin–Schreier trace extension size is out of range")
        })?;
        model.spend(order)?;
        orders.push(order);
    }
    let mut histograms = vec![vec![0; p]; bound + 1];
    let mut divisors: Vec<Vec<u32>> = Vec::new();
    let mut divisor_cells = 0;
    let trace_cells = terms.last().unwrap().0 + 1;
    let mut traces = vec![0; trace_cells];
    for (index, &count) in orders.iter().enumerate() {
        let m = index + 1;
        let mut monic = vec![1; m + 1];
        let mut remainder = vec![0; m + 1];
        'candidate: for mut code in 0..count {
            if m > 1 && code.is_multiple_of(p) {
                continue;
            }
            model.spend(m)?;
            for coefficient in &mut monic[..m] {
                *coefficient = (code % p) as u32;
                code /= p;
            }
            for divisor in &divisors {
                if divisor.len() - 1 > m / 2 {
                    break;
                }
                if divides(model, prime, &monic, divisor, &mut remainder)? {
                    continue 'candidate;
                }
            }
            let trace = trace_on_closed_point(model, prime, terms, &monic, &mut traces)? as usize;
            // A point of degree m contributes m conjugate affine points
            // over each degree k multiple of m, with trace multiplied by k/m.
            for k in (m..=bound).step_by(m) {
                model.spend(1)?;
                let image = trace * (k / m) % p;
                histograms[k][image] += m;
            }
            if m <= bound / 2 {
                divisor_cells += 3 + monic.len();
                check_storage(divisor_cells + histograms.len() * (p + 3) + traces.len())?;
                divisors.push(monic.clone());
            }
        }
    }
    for (histogram, order) in histograms.iter().skip(1).zip(orders) {
        model.spend(p)?;
        if histogram.iter().sum::<usize>() != order {
            return Err(Error::plain(
                "odd Artin–Schreier trace counts disagree with the field size",
            ));
        }
    }
    Ok(histograms)
}

fn divides(
    model: &mut Model,
    prime: Prime,
    monic: &[u32],
    divisor: &[u32],
    remainder: &mut [u32],
) -> Result<bool> {
    model.spend(1 + monic.len())?;
    remainder.copy_from_slice(monic);
    let degree = divisor.len() - 1;
    for i in (degree..monic.len()).rev() {
        let factor = remainder[i];
        if factor == 0 {
            continue;
        }
        model.spend(degree)?;
        for (c, &d) in remainder[i - degree..i].iter_mut().zip(divisor) {
            *c = prime.add(*c, prime.negate(prime.mul(factor, d)));
        }
    }
    Ok(remainder[..degree]
        .iter()
        .all(|&coefficient| coefficient == 0))
}

fn trace_on_closed_point(
    model: &mut Model,
    prime: Prime,
    terms: &[(usize, u32)],
    monic: &[u32],
    traces: &mut [u32],
) -> Result<u32> {
    let m = monic.len() - 1;
    let p = prime.get() as usize;
    traces[0] = (m % p) as u32;
    for k in 1..traces.len() {
        let mut sum = if k <= m {
            prime.mul((k % p) as u32, monic[m - k])
        } else {
            0
        };
        let limit = m.min(k - 1);
        model.spend(1 + limit)?;
        for j in 1..=limit {
            sum = prime.add(sum, prime.mul(monic[m - j], traces[k - j]));
        }
        traces[k] = prime.negate(sum);
    }
    model.spend(terms.len())?;
    Ok(terms.iter().fold(0, |trace, &(e, c)| {
        prime.add(trace, prime.mul(c, traces[e]))
    }))
}
