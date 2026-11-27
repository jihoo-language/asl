use std::cell::OnceCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::binary_linear::check_storage;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::rv_as_characters::{coefficient_error, Characters, Query as CharacterQuery};
use crate::rv_as_odd_cartier::Cartier;
use crate::rv_as_odd_derham::DeRham;
use crate::rv_as_odd_tate::Tate;
use crate::rv_cyclotomic as cyclo;

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
    characters: Characters,
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
            characters: Characters::new(prime, factor_degree, genus),
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

    pub(crate) fn coefficient(&self, model: &mut Model, j: usize, k: usize) -> Result<Value> {
        self.characters.coefficient(model, j, k, |model| {
            stationary_factor(model, self.prime, &self.terms, self.degree())
        })
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
                return self.character_query(model, name, operation, CharacterQuery::Character(j))
            }
            Query::CharacterNewton(j) => {
                return self.character_query(
                    model,
                    name,
                    operation,
                    CharacterQuery::CharacterNewton(j),
                )
            }
            Query::LPolynomial => {
                return self.character_query(model, name, operation, CharacterQuery::LPolynomial)
            }
            Query::Newton => {
                return self.character_query(model, name, operation, CharacterQuery::Newton)
            }
            Query::Points(m) => {
                return self.character_query(model, name, operation, CharacterQuery::Points(m))
            }
        };
        Ok(format!("{operation}[{name}] = {value}"))
    }

    fn character_query(
        &self,
        model: &mut Model,
        name: &str,
        operation: &str,
        query: CharacterQuery,
    ) -> Result<String> {
        self.characters
            .report(model, name, operation, query, |model| {
                stationary_factor(model, self.prime, &self.terms, self.degree())
            })
    }
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
