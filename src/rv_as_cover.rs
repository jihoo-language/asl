use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::binary_field::BinaryField;
use crate::binary_linear::check_storage;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_as_cartier;
use crate::rv_as_derham;
use crate::rv_curve::{affine_form, affine_gcd, binary_remainder};
use crate::rv_newton;
use crate::rv_poly::{self, Polynomial};

#[derive(Clone, Copy)]
pub(crate) enum Query {
    Genus,
    PRank,
    Points(usize),
    LPolynomial,
    Newton,
    CartierRank(usize),
    ANumber,
    FinalType,
    PoleOrders,
}

impl Query {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Genus => "rvASGenus",
            Self::PRank => "rvASPrank",
            Self::Points(_) => "rvASPoints",
            Self::LPolynomial => "rvASLPolynomial",
            Self::Newton => "rvASNewton",
            Self::CartierRank(_) => "rvASCartierRank",
            Self::ANumber => "rvASaNumber",
            Self::FinalType => "rvASEO",
            Self::PoleOrders => "rvASPoleOrders",
        }
    }
}

pub(crate) struct Cover {
    numerator: Polynomial,
    denominator: Polynomial,
    infinity_value: Option<bool>,
    genus: usize,
    p_rank: usize,
    counts: RefCell<BTreeMap<usize, usize>>,
    l_polynomial: RefCell<Option<Rc<[i128]>>>,
    radical: Polynomial,
    pole_half: Polynomial,
    cartier_denominator: RefCell<Option<Polynomial>>,
}

impl Cover {
    pub(crate) fn build(
        model: &mut Model,
        width: usize,
        forms: Vec<(usize, Polynomial)>,
    ) -> Result<Self> {
        if width != 2 {
            return Err(Error::plain(
                "rvascover requires exactly two parameter variables",
            ));
        }
        let [(numerator_degree, numerator), (degree, denominator)] = forms
            .try_into()
            .map_err(|_| Error::plain("rvascover requires a numerator and a denominator"))?;
        if degree == 0 || numerator_degree != degree {
            return Err(Error::plain(
                "rvascover forms must have the same positive degree",
            ));
        }
        if denominator.is_empty() {
            return Err(Error::plain("rvascover denominator must be nonzero"));
        }
        model.spend(rv_poly::cells(&numerator) + rv_poly::cells(&denominator))?;
        let numerator = affine_form(&numerator);
        let denominator = affine_form(&denominator);
        let infinity_order = degree - denominator.last().unwrap()[0];
        let numerator_at_infinity = numerator.last().is_some_and(|term| term[0] == degree);
        if infinity_order > 0 && !numerator_at_infinity {
            return Err(Error::plain(
                "rvascover forms have a common projective zero at [0:1]",
            ));
        }
        let gcd = affine_gcd(model, &[numerator.clone(), denominator.clone()])?;
        if !is_one(&gcd) {
            return Err(Error::plain(
                "rvascover forms have common projective zeros on the chart s != 0",
            ));
        }
        if infinity_order > 0 && infinity_order.is_multiple_of(2) {
            return Err(Error::plain(
                "Artin–Schreier pole at [0:1] must have odd order",
            ));
        }
        let (radical, pole_half) = odd_radical(model, &denominator)?;
        let poles = radical.last().unwrap()[0] + usize::from(infinity_order > 0);
        let genus = ((degree as u128 + poles as u128) / 2 - 1) as usize;
        model.allocate(
            7 + [&numerator, &denominator, &radical, &pole_half]
                .into_iter()
                .map(|p| rv_poly::cells(p))
                .sum::<usize>(),
        )?;
        Ok(Self {
            numerator,
            denominator,
            infinity_value: (infinity_order == 0).then_some(numerator_at_infinity),
            genus,
            p_rank: poles - 1,
            counts: RefCell::new(BTreeMap::new()),
            l_polynomial: RefCell::new(None),
            radical,
            pole_half,
            cartier_denominator: RefCell::new(None),
        })
    }

    pub(crate) fn genus(&self) -> usize {
        self.genus
    }

    pub(crate) fn p_rank(&self) -> usize {
        self.p_rank
    }

    pub(crate) fn has_polynomial_chart(&self) -> bool {
        is_one(&self.denominator)
    }

    pub(crate) fn polynomial_chart(&self) -> Result<&Polynomial> {
        if self.has_polynomial_chart() {
            Ok(&self.numerator)
        } else {
            Err(Error::plain("Artin–Schreier de Rham requires a one-pole polynomial cover with affine denominator 1"))
        }
    }

    pub(crate) fn geometric_pole_profile(&self, model: &mut Model) -> Result<Vec<(usize, usize)>> {
        // A finite pole of order 2e+1 occurs once in R and to order e in H,
        // where the denominator is R H^2. Split R by the binary digits of e.
        let mut factors = BTreeMap::from([(1usize, self.radical.clone())]);
        let mut remaining = self.pole_half.clone();
        let mut bit = 2usize;
        while !is_one(&remaining) {
            let (layer, next) = squarefree_layer(model, &remaining)?;
            if !is_one(&layer) {
                let cells: usize = factors
                    .values()
                    .map(|factor| 4 + rv_poly::cells(factor))
                    .sum();
                let mut next_cells = 0;
                let mut divided = BTreeMap::new();
                for (order, factor) in factors {
                    let present = affine_gcd(model, &[factor.clone(), layer.clone()])?;
                    if is_one(&present) {
                        next_cells += 4 + rv_poly::cells(&factor);
                        divided.insert(order, factor);
                    } else {
                        let shifted = order.checked_add(bit).ok_or_else(pole_profile_error)?;
                        let absent = rv_poly::exact_divide(model, &factor, &present)?;
                        next_cells += 4 + rv_poly::cells(&present);
                        divided.insert(shifted, present);
                        if !is_one(&absent) {
                            next_cells += 4 + rv_poly::cells(&absent);
                            divided.insert(order, absent);
                        }
                    }
                    check_storage(cells + next_cells)?;
                }
                factors = divided;
            }
            remaining = next;
            if !is_one(&remaining) {
                bit = bit.checked_mul(2).ok_or_else(pole_profile_error)?;
            }
        }
        // The factors are squarefree. Degree counts roots over the algebraic
        // closure, including every member of a nonrational Frobenius orbit.
        let mut profile = BTreeMap::new();
        for (order, factor) in factors {
            let degree = factor.last().unwrap()[0];
            if degree != 0 {
                profile.insert(order, degree);
            }
        }
        if self.infinity_value.is_none() {
            let order = self.numerator.last().unwrap()[0] - self.denominator.last().unwrap()[0];
            *profile.entry(order).or_default() += 1;
        }
        let poles: usize = profile.values().copied().sum();
        let Some(ordinary) = poles.checked_sub(1) else {
            return Err(pole_profile_error());
        };
        let local_genera: u128 = profile
            .iter()
            .map(|(&m, &count)| ((m - 1) / 2) as u128 * count as u128)
            .sum();
        if ordinary != self.p_rank || local_genera + ordinary as u128 != self.genus as u128 {
            return Err(pole_profile_error());
        }
        Ok(profile.into_iter().collect())
    }

    pub(crate) fn cartier_denominator(&self, model: &mut Model) -> Result<Polynomial> {
        if let Some(denominator) = self.cartier_denominator.borrow().as_ref() {
            return Ok(denominator.clone());
        }
        // q=R*H² for the radical R and the stored square root H. Hence
        // sqrt(q*R)=R*H, avoiding a needless square and square-root pass.
        let denominator = rv_poly::multiply(model, &self.radical, &self.pole_half)?;
        model.allocate(rv_poly::cells(&denominator))?;
        *self.cartier_denominator.borrow_mut() = Some(denominator.clone());
        Ok(denominator)
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str, query: Query) -> Result<String> {
        let value = match query {
            Query::Genus => self.genus.to_string(),
            Query::PRank => self.p_rank.to_string(),
            Query::Points(extension) => {
                let points = self.points(model, extension)?;
                return Ok(format!("rvASPoints[{name}, {extension}] = {points}"));
            }
            Query::CartierRank(power) => {
                let rank = rv_as_cartier::rank(model, self, power)?;
                return Ok(format!("rvASCartierRank[{name}, {power}] = {rank}"));
            }
            Query::ANumber => (self.genus - rv_as_cartier::rank(model, self, 1)?).to_string(),
            Query::FinalType => format!("{:?}", rv_as_derham::final_type(model, self)?),
            Query::PoleOrders => format!("{:?}", self.geometric_pole_profile(model)?),
            Query::LPolynomial => format!("{:?}", self.l_polynomial(model)?),
            Query::Newton => {
                let polynomial = self.l_polynomial(model)?;
                format!(
                    "[{}]",
                    rv_newton::integral_slopes(model, &polynomial)?.join(", ")
                )
            }
        };
        Ok(format!("{}[{name}] = {value}", query.name()))
    }

    fn points(&self, model: &mut Model, extension: usize) -> Result<usize> {
        if let Some(&points) = self.counts.borrow().get(&extension) {
            return Ok(points);
        }
        let field = point_count_field(model, extension)?;
        let mut points = match self.infinity_value {
            None => 1,
            Some(value) => 2 * usize::from(!value || extension.is_multiple_of(2)),
        };
        for parameter in 0..field.order() {
            let denominator = affine_value(model, field, &self.denominator, parameter)?;
            if denominator == 0 {
                points += 1;
                continue;
            }
            let numerator = affine_value(model, field, &self.numerator, parameter)?;
            model.spend(field.degree())?;
            let value = field.multiply(numerator, field.inverse(model, denominator)?);
            if trace_zero(model, field, value)? {
                points += 2;
            }
        }
        model.allocate(3)?;
        self.counts.borrow_mut().insert(extension, points);
        Ok(points)
    }

    fn l_polynomial(&self, model: &mut Model) -> Result<Rc<[i128]>> {
        if let Some(polynomial) = self.l_polynomial.borrow().as_ref() {
            return Ok(polynomial.clone());
        }
        let coefficients = integral_l_polynomial(
            model,
            self.genus,
            "Artin–Schreier L-polynomial",
            |model, extension| self.points(model, extension),
        )?;
        model.allocate(1 + 2 * coefficients.len())?;
        let polynomial: Rc<[i128]> = coefficients.into();
        *self.l_polynomial.borrow_mut() = Some(polynomial.clone());
        Ok(polynomial)
    }
}

pub(crate) fn integral_l_polynomial(
    model: &mut Model,
    genus: usize,
    description: &str,
    mut points: impl FnMut(&mut Model, usize) -> Result<usize>,
) -> Result<Vec<i128>> {
    if genus >= usize::BITS as usize {
        return Err(Error::plain(format!(
            "{description} exceeds the point-count extension range"
        )));
    }
    let mut sums = vec![0i128; genus + 1];
    for (extension, sum) in sums.iter_mut().enumerate().skip(1) {
        *sum = 1 + extension_size(extension)? as i128 - points(model, extension)? as i128;
    }
    let mut coefficients = vec![0i128; 2 * genus + 1];
    coefficients[0] = 1;
    // Newton identities for L(T)=product_i(1-alpha_i T).
    for degree in 1..=genus {
        model.spend(degree)?;
        let mut numerator = 0i128;
        for extension in 1..=degree {
            let contribution = coefficients[degree - extension]
                .checked_mul(sums[extension])
                .ok_or_else(|| coefficient_limit(description))?;
            numerator = numerator
                .checked_sub(contribution)
                .ok_or_else(|| coefficient_limit(description))?;
        }
        if numerator % degree as i128 != 0 {
            return Err(Error::plain(format!(
                "{description} Newton identity did not give an integer coefficient"
            )));
        }
        coefficients[degree] = numerator / degree as i128;
    }
    // Complete the coefficients by the genus-g functional equation.
    for degree in 0..genus {
        coefficients[2 * genus - degree] = coefficients[degree]
            .checked_mul(1i128 << (genus - degree))
            .ok_or_else(|| coefficient_limit(description))?;
    }
    Ok(coefficients)
}

fn coefficient_limit(description: &str) -> Error {
    Error::plain(format!("{description} coefficient is out of range"))
}

pub(crate) fn point_count_field(model: &mut Model, extension: usize) -> Result<BinaryField> {
    let field_size = extension_size(extension)?;
    // Budget the full scan before field construction. Sparse evaluations and
    // arithmetic over the field charge their work separately.
    model.spend(field_size.saturating_mul(1 + extension))?;
    BinaryField::of_degree(model, extension)
}

fn extension_size(extension: usize) -> Result<usize> {
    if !(1..usize::BITS as usize).contains(&extension) {
        return Err(Error::plain("point-count extension degree is out of range"));
    }
    Ok(1 << extension)
}

fn pole_profile_error() -> Error {
    Error::plain("Artin–Schreier geometric pole orders do not agree with the curve dimensions")
}

fn is_one(polynomial: &Polynomial) -> bool {
    polynomial.len() == 1 && polynomial[0][0] == 0
}

fn odd_radical(model: &mut Model, denominator: &Polynomial) -> Result<(Polynomial, Polynomial)> {
    let (radical, half) = squarefree_layer(model, denominator)?;
    let mut remaining = half.clone();
    while !is_one(&remaining) {
        let (layer, next) = squarefree_layer(model, &remaining)?;
        if !binary_remainder(model, &radical, &layer)?.is_empty() {
            return Err(Error::plain(
                "Artin–Schreier finite poles must all have odd order",
            ));
        }
        remaining = next;
    }
    Ok((radical, half))
}

fn squarefree_layer(
    model: &mut Model,
    polynomial: &Polynomial,
) -> Result<(Polynomial, Polynomial)> {
    model.spend(rv_poly::cells(polynomial))?;
    let derivative = polynomial
        .iter()
        .filter(|term| term[0] & 1 != 0)
        .map(|term| [term[0] - 1].into())
        .collect();
    let square = affine_gcd(model, &[polynomial.clone(), derivative])?;
    let odd = rv_poly::exact_divide(model, polynomial, &square)?;
    // Every irreducible factor of the derivative gcd occurs with even
    // multiplicity. Taking its square root peels another binary digit.
    model.spend(rv_poly::cells(&square))?;
    let remaining = square.iter().map(|term| [term[0] / 2].into()).collect();
    Ok((odd, remaining))
}

pub(crate) fn affine_value(
    model: &mut Model,
    field: BinaryField,
    polynomial: &Polynomial,
    point: usize,
) -> Result<usize> {
    let Some((leading, tail)) = polynomial.split_last() else {
        return Ok(0);
    };
    let mut power = leading[0];
    let mut value = 1;
    model.spend(field.degree() * polynomial.len())?;
    for term in tail.iter().rev() {
        value = field.multiply(value, field.power(model, point, (power - term[0]) as u128)?) ^ 1;
        power = term[0];
    }
    Ok(field.multiply(value, field.power(model, point, power as u128)?))
}

fn trace_zero(model: &mut Model, field: BinaryField, mut value: usize) -> Result<bool> {
    if value <= 1 {
        return Ok(value == 0 || field.degree().is_multiple_of(2));
    }
    let mut trace = value;
    model.spend(field.degree() * (field.degree() - 1))?;
    for _ in 1..field.degree() {
        value = field.multiply(value, value);
        trace ^= value;
    }
    Ok(trace == 0)
}
