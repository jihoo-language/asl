use std::cell::OnceCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::binary_linear::check_storage;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::rv_as_odd_cartier::Cartier;
use crate::rv_as_odd_derham::DeRham;
use crate::rv_as_odd_tate::Tate;
use crate::{rv_cyclotomic as cyclo, rv_newton};

type Character = Rc<[Vec<i128>]>;

#[derive(Clone, Copy)]
pub(crate) enum Query {
    Reduced,
    Basis,
    CartierRank(usize),
    ANumber,
    FinalType,
    Gram,
    TateDim(i64),
    TateJordan,
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
            Self::FinalType => "rvASOddEO",
            Self::Gram => "rvASOddGram",
            Self::TateDim(_) => "rvASOddTatedim",
            Self::TateJordan => "rvASOddTateJordan",
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
    derham: OnceCell<Rc<DeRham>>,
    tate: OnceCell<Rc<Tate>>,
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
            derham: OnceCell::new(),
            tate: OnceCell::new(),
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

    pub(crate) fn derham(&self, model: &mut Model) -> Result<Rc<DeRham>> {
        if let Some(derham) = self.derham.get() {
            return Ok(derham.clone());
        }
        let derham = Rc::new(DeRham::build(model, self)?);
        Ok(self.derham.get_or_init(|| derham).clone())
    }

    pub(crate) fn tate(&self, model: &mut Model) -> Result<Rc<Tate>> {
        if let Some(tate) = self.tate.get() {
            return Ok(tate.clone());
        }
        let tate = Rc::new(Tate::build(model, self)?);
        Ok(self.tate.get_or_init(|| tate).clone())
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
            Query::FinalType => format!("{:?}", self.derham(model)?.final_type(model, self)?),
            Query::Gram => format!("{:?}", self.derham(model)?.gram(model, self)?),
            Query::TateDim(degree) => {
                let rank = self.tate(model)?.group(model, self, degree)?.rank();
                return Ok(format!("{operation}[{name}, {degree}] = {rank}"));
            }
            Query::TateJordan => format!("{:?}", self.tate(model)?.jordan(model, self)?),
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
        let n = self.factor_degree;
        let factor = stationary_factor(model, self.prime, &self.terms, self.degree())?;
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

fn stationary_factor(
    model: &mut Model,
    prime: Prime,
    terms: &[(usize, u32)],
    degree: usize,
) -> Result<Vec<Vec<i128>>> {
    let p = prime.get() as usize;
    let slots = degree / 2;
    let last = degree - 1;
    let constant = terms
        .first()
        .filter(|&&(exponent, _)| exponent == 0)
        .map_or(0, |&(_, c)| c);

    // The Euler coefficient of degree k is the sum of the character
    // over all monic polynomials of degree k, including repeated roots.
    // In its reversed polynomial A(t)=1+a_1 t+...+a_k t^k, the sum
    // of the j-th powers of the roots is -[t^j] t A'(t)/A(t).
    // Split A=L+H at floor(d/2). Terms quadratic in H vanish
    // through degree d. The high part of t A'/A is t(H/L)′
    // in this range, so each high variable enters linearly.
    let cells = degree
        .saturating_mul(p + 4)
        .saturating_add(2 * slots)
        .saturating_add(3usize.saturating_mul(terms.len()));
    check_storage(cells)?;
    let population = u32::try_from(slots)
        .ok()
        .and_then(|h| p.checked_pow(h))
        .ok_or_else(|| {
            Error::plain("odd Artin–Schreier stationary phase exceeds the computation limit")
        })?;
    model.spend(population)?;
    model.spend(terms.len())?;
    let positive: Vec<_> = terms
        .iter()
        .filter(|&&(j, _)| j != 0)
        .map(|&(j, c)| (j, c, prime.mul((j % p) as u32, c)))
        .collect();
    let mut histograms = vec![vec![0usize; p]; degree];
    let mut low = vec![0; slots];
    let mut derivative = vec![0; slots];
    let mut reciprocal = vec![0; degree];
    reciprocal[0] = 1;

    for mut index in 0..population {
        model.spend(slots)?;
        let mut low_degree = 0;
        for (i, (coefficient, derivative)) in low.iter_mut().zip(&mut derivative).enumerate() {
            *coefficient = (index % p) as u32;
            index /= p;
            *derivative = prime.mul(((i + 1) % p) as u32, *coefficient);
            if *coefficient != 0 {
                low_degree = i + 1;
            }
        }

        // Constant f contributes k*f_0, with k the number of roots
        // including multiplicities. The remaining phase only depends
        // on L, even when the low polynomial has smaller degree.
        for i in 1..degree {
            let width = i.min(low_degree);
            model.spend(1 + width)?;
            let coefficient = low[..width]
                .iter()
                .zip(reciprocal[i - width..i].iter().rev())
                .fold(0, |sum, (&a, &b)| prime.add(sum, prime.mul(a, b)));
            reciprocal[i] = prime.negate(coefficient);
        }
        let mut phase = 0;
        for &(j, c, _) in &positive {
            let width = j.min(low_degree);
            model.spend(1 + width)?;
            let sum = derivative[..width]
                .iter()
                .zip(reciprocal[j - width..j].iter().rev())
                .fold(0, |sum, (&a, &b)| prime.add(sum, prime.mul(a, b)));
            phase = prime.add(phase, prime.negate(prime.mul(c, sum)));
        }

        // Orthogonality kills every k that reaches the first nonzero
        // coefficient of an available high variable. Up to that limit,
        // each high variable contributes a factor of p.
        let mut end = last;
        let mut tail = positive.as_slice();
        for i in slots + 1..=last {
            while let Some(((j, _, _), rest)) = tail.split_first() {
                if *j >= i {
                    break;
                }
                model.spend(1)?;
                tail = rest;
            }
            model.spend(1 + tail.len())?;
            let coefficient = tail.iter().fold(0, |sum, &(j, _, weight)| {
                prime.add(sum, prime.mul(weight, reciprocal[j - i]))
            });
            if coefficient != 0 {
                end = i - 1;
                break;
            }
        }
        model.spend(end + 1 - low_degree)?;
        let mut phase = prime.add(phase, prime.mul((low_degree % p) as u32, constant));
        for histogram in &mut histograms[low_degree..=end] {
            histogram[phase as usize] += 1;
            phase = prime.add(phase, constant);
        }
    }

    let mut factor = Vec::with_capacity(degree);
    let mut weight = 1i128;
    for (k, histogram) in histograms.into_iter().enumerate() {
        if k > slots {
            weight = weight
                .checked_mul(p as i128)
                .ok_or_else(coefficient_error)?;
        }
        model.spend(histogram.len())?;
        let counts = histogram
            .into_iter()
            .map(|count| {
                (count as i128)
                    .checked_mul(weight)
                    .ok_or_else(coefficient_error)
            })
            .collect::<Result<Vec<_>>>()?;
        factor.push(cyclo::reduced(model, p, counts)?);
    }
    Ok(factor)
}
