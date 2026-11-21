use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::binary_field::BinaryField;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_as_cover::{affine_value, point_count_field};
use crate::rv_newton;
use crate::rv_poly::{self, Polynomial};

const MAX_HEIGHT: usize = 8;
const MAX_CELLS: usize = 300_000;

type Character = Rc<[Vec<i128>]>;
type Histogram = Rc<[usize]>;

#[derive(Clone, Copy)]
pub(crate) enum Query {
    Breaks,
    Genus,
    Points(usize),
    CharacterL(usize),
    CharacterNewton(usize),
    CharacterBaseL(usize, usize),
    FrobPoints(usize),
    LPolynomial,
    Newton,
}

impl Query {
    fn name(self) -> &'static str {
        match self {
            Self::Breaks => "rvWittTowerBreaks",
            Self::Genus => "rvWittTowerGenus",
            Self::Points(_) => "rvWittTowerPoints",
            Self::CharacterL(_) => "rvWittTowerCharacterL",
            Self::CharacterNewton(_) => "rvWittTowerCharacterNewton",
            Self::CharacterBaseL(..) => "rvWittTowerCharacterBaseL",
            Self::FrobPoints(_) => "rvWittTowerFrobPoints",
            Self::LPolynomial => "rvWittTowerLPolynomial",
            Self::Newton => "rvWittTowerNewton",
        }
    }
}

pub(crate) struct Tower {
    coordinates: Vec<Polynomial>,
    breaks: Vec<usize>,
    genus: usize,
    traces: RefCell<BTreeMap<(usize, usize), Histogram>>,
    characters: RefCell<Vec<Option<Character>>>,
    numerator: RefCell<Option<Rc<[i128]>>>,
}

impl Tower {
    pub(crate) fn build(model: &mut Model, coordinates: Vec<Polynomial>) -> Result<Self> {
        if !(1..=MAX_HEIGHT).contains(&coordinates.len()) {
            return Err(Error::plain(format!(
                "rvwitttower requires between 1 and {MAX_HEIGHT} coordinates"
            )));
        }
        let mut breaks = Vec::with_capacity(coordinates.len());
        let mut different = 0;
        let mut previous = 0usize;
        for (index, coordinate) in coordinates.iter().enumerate() {
            let degree = coordinate.last().map_or(0, |term| term[0]);
            if degree.is_multiple_of(2) && (index == 0 || degree != 0) {
                return Err(Error::plain(format!(
                    "rvwitttower coordinate {index} must have odd degree{}",
                    if index == 0 {
                        " and be nonconstant"
                    } else {
                        " when nonconstant"
                    },
                )));
            }
            let doubled = previous.checked_mul(2).ok_or_else(|| {
                Error::plain("Witt tower upper ramification break is out of range")
            })?;
            let jump = doubled.max(degree);
            different += (1u128 << index) * (jump as u128 + 1);
            breaks.push(jump);
            previous = jump;
        }
        let genus = (different + 2 - (2u128 << coordinates.len())) / 2;
        let genus =
            usize::try_from(genus).map_err(|_| Error::plain("Witt tower genus is out of range"))?;
        model.allocate(
            4 + 2 * coordinates.len()
                + coordinates
                    .iter()
                    .map(|poly| rv_poly::cells(poly))
                    .sum::<usize>(),
        )?;
        let characters = RefCell::new(vec![None; coordinates.len()]);
        Ok(Self {
            coordinates,
            breaks,
            genus,
            traces: RefCell::new(BTreeMap::new()),
            characters,
            numerator: RefCell::new(None),
        })
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str, query: Query) -> Result<String> {
        let value = match query {
            Query::Breaks => format!("{:?}", self.breaks),
            Query::Genus => self.genus.to_string(),
            Query::Points(extension) => {
                let height = self.coordinates.len();
                let histogram = self.histogram(model, height, extension)?;
                let size = 1u128 + (histogram[0] as u128) * (1u128 << height);
                return Ok(format!("{}[{name}, {extension}] = {size}", query.name()));
            }
            Query::CharacterL(height) => {
                let character = self.character_l(model, height)?;
                return Ok(format!(
                    "{}[{name}, {height}] = {character:?}",
                    query.name()
                ));
            }
            Query::CharacterNewton(height) => {
                let character = self.character_l(model, height)?;
                let slopes = character_slopes(model, &character)?.join(", ");
                return Ok(format!("{}[{name}, {height}] = [{slopes}]", query.name()));
            }
            Query::CharacterBaseL(height, extension) => {
                check_positive_extension(extension)?;
                let character = self.character_l(model, height)?;
                let base = character_base_l(model, &character, extension)?;
                return Ok(format!(
                    "{}[{name}, {height}, {extension}] = {base:?}",
                    query.name()
                ));
            }
            Query::FrobPoints(extension) => {
                let points = self.frobenius_points(model, extension)?;
                return Ok(format!("{}[{name}, {extension}] = {points}", query.name()));
            }
            Query::LPolynomial => format!("{:?}", self.l_polynomial(model)?),
            Query::Newton => {
                let coefficients = self.l_polynomial(model)?;
                format!(
                    "[{}]",
                    rv_newton::integral_slopes(model, &coefficients)?.join(", ")
                )
            }
        };
        Ok(format!("{}[{name}] = {value}", query.name()))
    }

    fn histogram(&self, model: &mut Model, height: usize, extension: usize) -> Result<Histogram> {
        let size = 1usize << height;
        // A trace in a higher quotient determines every lower trace by reduction.
        let present = self.traces.borrow();
        if let Some(histogram) = present.get(&(height, extension)) {
            return Ok(histogram.clone());
        }
        let mut counts = vec![0; size];
        if let Some(histogram) = ((height + 1)..=self.coordinates.len())
            .find_map(|higher| present.get(&(higher, extension)))
        {
            model.spend(histogram.len())?;
            for (trace, &count) in histogram.iter().enumerate() {
                counts[trace & (size - 1)] += count;
            }
        } else {
            let field = point_count_field(model, extension)?;
            let mut trace = TeichTrace::new(field, height)?;
            for x in 0..field.order() {
                let mut sum = 0;
                // Verschiebung introduces inverse Frobenius on the
                // Teichmuller factor, which leaves its trace unchanged.
                for (index, coordinate) in self.coordinates[..height].iter().enumerate() {
                    let value = affine_value(model, field, coordinate, x)?;
                    sum += usize::from(trace.trace(model, value)?) << index;
                }
                counts[sum & (size - 1)] += 1;
            }
        }
        drop(present);
        model.allocate(2 + size)?;
        let counts: Histogram = counts.into();
        self.traces
            .borrow_mut()
            .insert((height, extension), counts.clone());
        Ok(counts)
    }

    fn character_l(&self, model: &mut Model, height: usize) -> Result<Character> {
        if !(1..=self.coordinates.len()).contains(&height) {
            return Err(Error::plain(
                "Witt character quotient height is out of range",
            ));
        }
        if let Some(polynomial) = self.characters.borrow()[height - 1].as_ref() {
            return Ok(polynomial.clone());
        }
        let degree = self.breaks[height - 1] - 1;
        if degree >= usize::BITS as usize {
            return Err(Error::plain(
                "Witt character degree exceeds the point-count extension range",
            ));
        }
        let width = 1usize << (height - 1);
        check_cells(4 * (degree + 1) * width)?;
        let mut sums = vec![vec![0; width]];
        for extension in 1..=degree {
            let histogram = self.histogram(model, height, extension)?;
            model.spend(2 * width)?;
            // In the power basis, zeta^width=-1.
            sums.push(
                histogram[..width]
                    .iter()
                    .zip(&histogram[width..])
                    .map(|(&plus, &minus)| plus as i128 - minus as i128)
                    .collect(),
            );
        }
        let polynomial = from_negative_power_sums(model, &sums)?;
        model.allocate(1 + (1 + 2 * width) * polynomial.len())?;
        let polynomial: Character = polynomial.into();
        self.characters.borrow_mut()[height - 1] = Some(polynomial.clone());
        Ok(polynomial)
    }

    fn frobenius_points(&self, model: &mut Model, extension: usize) -> Result<u128> {
        check_positive_extension(extension)?;
        if extension >= u128::BITS as usize {
            return Err(Error::plain(
                "Witt tower Frobenius point-count extension degree exceeds 127",
            ));
        }
        let mut negative_trace = 0i128;
        for height in (1..=self.coordinates.len()).rev() {
            let character = self.character_l(model, height)?;
            if character.len() == 1 {
                continue;
            }
            let sums = negative_power_sums(model, &character, extension, 1)?;
            // In the power basis, the field trace of 1 is the degree
            // and every other basis element has field trace zero.
            let trace = sums[1][0]
                .checked_mul(character[0].len() as i128)
                .ok_or_else(coefficient_error)?;
            negative_trace = negative_trace
                .checked_add(trace)
                .ok_or_else(coefficient_error)?;
        }
        let base = (1u128 << extension) + 1;
        let points = if negative_trace >= 0 {
            base.checked_add(negative_trace as u128)
        } else {
            base.checked_sub(negative_trace.unsigned_abs())
        };
        points.ok_or_else(|| Error::plain("Witt tower Frobenius point count is out of range"))
    }

    fn l_polynomial(&self, model: &mut Model) -> Result<Rc<[i128]>> {
        if let Some(polynomial) = self.numerator.borrow().as_ref() {
            return Ok(polynomial.clone());
        }
        // The top coefficient is 2^g, so this query has no representable
        // result once the sign bit would be needed.
        if self.genus >= i128::BITS as usize - 1 {
            return Err(coefficient_error());
        }
        let mut polynomial = vec![vec![1]];
        for height in (1..=self.coordinates.len()).rev() {
            let character = self.character_l(model, height)?;
            let norm = cyclotomic_norm(model, &character)?;
            polynomial = product(model, &polynomial, &norm)?;
        }
        let polynomial: Vec<_> = polynomial
            .into_iter()
            .map(|coefficient| coefficient[0])
            .collect();
        if polynomial.len() != 2 * self.genus + 1 {
            return Err(Error::plain(
                "Witt tower character degrees disagree with the genus",
            ));
        }
        if polynomial.last() != Some(&(1i128 << self.genus)) {
            return Err(Error::plain(
                "Witt tower numerator has the wrong leading coefficient for its genus",
            ));
        }
        model.allocate(1 + 2 * polynomial.len())?;
        let polynomial: Rc<[_]> = polynomial.into();
        *self.numerator.borrow_mut() = Some(polynomial.clone());
        Ok(polynomial)
    }
}

fn check_positive_extension(extension: usize) -> Result<()> {
    if extension == 0 {
        return Err(Error::plain(
            "Witt tower Frobenius extension degree must be positive",
        ));
    }
    Ok(())
}

fn character_slopes(model: &mut Model, character: &[Vec<i128>]) -> Result<Vec<String>> {
    let width = character[0].len();
    let valuations = character
        .iter()
        .map(|a| uniformizer_order(model, a))
        .collect::<Result<Vec<_>>>()?;
    rv_newton::scaled_slopes(model, valuations, width)
}

// The local uniformizer is pi = zeta - 1, with ord_pi(2) = width.
fn uniformizer_order(model: &mut Model, a: &[i128]) -> Result<Option<i64>> {
    let width = a.len();
    model.spend(width)?;
    let Some(k) = a
        .iter()
        .filter(|&&c| c != 0)
        .map(|c| c.trailing_zeros())
        .min()
    else {
        return Ok(None);
    };
    let mut parity: Vec<_> = a.iter().map(|&c| ((c >> k) & 1) != 0).collect();
    // After dividing out 2^k, the coefficient of pi^j modulo 2 is
    // the parity sum over the supersets of j, by Lucas's theorem.
    for direction in 0..width.ilog2() {
        model.spend(width)?;
        let bit = 1 << direction;
        for i in 0..width {
            if i & bit == 0 {
                parity[i] ^= parity[i | bit];
            }
        }
    }
    let first = parity
        .into_iter()
        .position(|odd| odd)
        .ok_or_else(|| Error::plain("normalized Witt character coefficient is zero modulo 2"))?;
    Ok(Some(i64::from(k) * width as i64 + first as i64))
}

fn character_base_l(
    model: &mut Model,
    character: &[Vec<i128>],
    extension: usize,
) -> Result<Vec<Vec<i128>>> {
    if extension == 1 || character.len() == 1 {
        model.spend(character[0].len() * character.len())?;
        return Ok(character.to_vec());
    }
    let degree = character.len() - 1;
    let sums = negative_power_sums(model, character, extension, degree)?;
    from_negative_power_sums(model, &sums)
}

// Recover -sum alpha_i^(k*stride) for 1 <= k <= count from
// L(T) = product (1 - alpha_i T). Keep only one recurrence window.
fn negative_power_sums(
    model: &mut Model,
    character: &[Vec<i128>],
    stride: usize,
    count: usize,
) -> Result<Vec<Vec<i128>>> {
    let degree = character.len() - 1;
    debug_assert!(degree != 0 && stride != 0);
    let width = character[0].len();
    let last = stride
        .checked_mul(count)
        .ok_or_else(|| Error::plain("Witt tower Frobenius power index is out of range"))?;
    model.spend(last)?;
    check_cells((1 + 2 * width) * (degree + count + 2))?;
    let mut window = vec![vec![0; width]; degree];
    let mut sums = Vec::with_capacity(count + 1);
    sums.push(vec![0; width]);
    for power in 1..=last {
        let mut sum = vec![0; width];
        model.spend(width)?;
        if let Some(coefficient) = character.get(power) {
            for (entry, &c) in sum.iter_mut().zip(coefficient) {
                *entry = c.checked_mul(power as i128).ok_or_else(coefficient_error)?;
            }
        }
        for (i, c) in character
            .iter()
            .enumerate()
            .take(power.min(degree + 1))
            .skip(1)
        {
            let previous = &window[(power - i) % degree];
            accumulate_product(model, &mut sum, c, previous, true)?;
        }
        if power.is_multiple_of(stride) {
            model.spend(width)?;
            sums.push(sum.clone());
        }
        window[power % degree] = sum;
    }
    Ok(sums)
}

fn from_negative_power_sums(model: &mut Model, sums: &[Vec<i128>]) -> Result<Vec<Vec<i128>>> {
    let width = sums[0].len();
    check_cells((2 + 4 * width) * sums.len())?;
    let mut one = vec![0; width];
    one[0] = 1;
    let mut polynomial = vec![one];
    for n in 1..sums.len() {
        let mut numerator = vec![0; width];
        for k in 1..=n {
            accumulate_product(model, &mut numerator, &sums[k], &polynomial[n - k], false)?;
        }
        let divisor = n as i128;
        for coefficient in &mut numerator {
            if *coefficient % divisor != 0 {
                return Err(Error::plain(
                    "Witt tower Newton identity is not integral in the power basis",
                ));
            }
            *coefficient /= divisor;
        }
        polynomial.push(numerator);
    }
    Ok(polynomial)
}

struct TeichTrace {
    field: BinaryField,
    height: usize,
    reduction_mask: u32,
    modulus_tail: Vec<usize>,
    traces: Vec<Option<u16>>,
}

impl TeichTrace {
    fn new(field: BinaryField, height: usize) -> Result<Self> {
        check_cells(field.order().saturating_add(8 * field.degree()))?;
        Ok(Self {
            field,
            height,
            reduction_mask: (1u32 << height) - 1,
            modulus_tail: (0..field.degree())
                .filter(|&bit| field.modulus() & (1 << bit) != 0)
                .collect(),
            traces: vec![None; field.order()],
        })
    }

    fn trace(&mut self, model: &mut Model, value: usize) -> Result<u16> {
        if let Some(trace) = self.traces[value] {
            return Ok(trace);
        }
        let e = self.field.degree();
        model.spend(1 + 2 * e * e)?;
        let mut lifted: Vec<_> = (0..e).map(|bit| ((value >> bit) & 1) as u32).collect();
        // Raising to Q repeatedly fixes the Teichmuller part and removes
        // one 2-adic digit of the principal unit per round.
        for _ in 0..e * (self.height - 1) {
            lifted = self.square(model, &lifted)?;
        }
        let mut sum = vec![0u32; e];
        let mut residue = value;
        for step in 0..e {
            for (entry, &coefficient) in sum.iter_mut().zip(&lifted) {
                *entry = (*entry + coefficient) & self.reduction_mask;
            }
            if step + 1 < e {
                lifted = self.square(model, &lifted)?;
            }
        }
        if sum[1..].iter().any(|&value| value != 0) {
            return Err(Error::plain(
                "Witt tower trace did not land in the scalar ring",
            ));
        }
        let trace = sum[0] as u16;
        for _ in 0..e {
            self.traces[residue] = Some(trace);
            residue = self.field.multiply(residue, residue);
            if residue == value {
                break;
            }
        }
        Ok(trace)
    }

    fn square(&self, model: &mut Model, coefficients: &[u32]) -> Result<Vec<u32>> {
        let e = self.field.degree();
        model.spend(e * (1 + e + self.modulus_tail.len()))?;
        let mask = self.reduction_mask;
        let mut square = vec![0u32; 2 * e - 1];
        for (i, &a) in coefficients.iter().enumerate() {
            square[2 * i] = (square[2 * i] + a * a) & mask;
            for (j, &b) in coefficients[..i].iter().enumerate() {
                square[i + j] = (square[i + j] + 2 * a * b) & mask;
            }
        }
        for high in (e..square.len()).rev() {
            let top = square[high];
            for &term in &self.modulus_tail {
                let coefficient = &mut square[high - e + term];
                *coefficient = coefficient.wrapping_sub(top) & mask;
            }
        }
        square.truncate(e);
        Ok(square)
    }
}

fn cyclotomic_norm(model: &mut Model, input: &[Vec<i128>]) -> Result<Vec<Vec<i128>>> {
    let mut polynomial = input.to_vec();
    while polynomial[0].len() > 1 {
        let mut conjugate = polynomial.clone();
        for coefficient in &mut conjugate {
            for value in coefficient.iter_mut().skip(1).step_by(2) {
                *value = value.checked_neg().ok_or_else(coefficient_error)?;
            }
        }
        // The relative norm takes zeta to both zeta and -zeta. Its odd
        // powers vanish, leaving the power basis in zeta^2.
        polynomial = product(model, &polynomial, &conjugate)?;
        for coefficient in &mut polynomial {
            if coefficient
                .iter()
                .skip(1)
                .step_by(2)
                .any(|&value| value != 0)
            {
                return Err(Error::plain(
                    "Witt tower relative norm is outside the smaller cyclotomic ring",
                ));
            }
            *coefficient = coefficient.iter().step_by(2).copied().collect();
        }
    }
    Ok(polynomial)
}

fn product(model: &mut Model, left: &[Vec<i128>], right: &[Vec<i128>]) -> Result<Vec<Vec<i128>>> {
    let width = left[0].len();
    let degree = left.len() + right.len() - 1;
    check_cells(2 * width * (degree + left.len() + right.len()))?;
    let mut polynomial = vec![vec![0; width]; degree];
    for (i, a) in left.iter().enumerate() {
        for (j, b) in right.iter().enumerate() {
            accumulate_product(model, &mut polynomial[i + j], a, b, false)?;
        }
    }
    Ok(polynomial)
}

fn accumulate_product(
    model: &mut Model,
    out: &mut [i128],
    left: &[i128],
    right: &[i128],
    subtract: bool,
) -> Result<()> {
    let width = out.len();
    model.spend(width)?;
    for (i, &a) in left.iter().enumerate().filter(|&(_, &c)| c != 0) {
        model.spend(width)?;
        for (j, &b) in right.iter().enumerate().filter(|&(_, &c)| c != 0) {
            let product = a.checked_mul(b).ok_or_else(coefficient_error)?;
            let target = &mut out[(i + j) & (width - 1)];
            *target = if subtract ^ (i + j >= width) {
                target.checked_sub(product)
            } else {
                target.checked_add(product)
            }
            .ok_or_else(coefficient_error)?;
        }
    }
    Ok(())
}

fn check_cells(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(Error::plain(
            "Witt tower calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

fn coefficient_error() -> Error {
    Error::plain("Witt tower L-polynomial coefficient is out of range")
}

#[cfg(test)]
mod tests {
    use super::{character_base_l, Model};

    #[test]
    fn character_base_change_composes_in_the_original_power_basis() {
        let mut model = Model::new();
        let character = vec![
            vec![1, 0, 0, 0],
            vec![1, 1, 0, 0],
            vec![1, 1, 1, -1],
            vec![2, 0, 2, 0],
        ];
        let quadratic =
            character_base_l(&mut model, &character, 2).unwrap_or_else(|error| panic!("{error}"));
        let quadratic_twice =
            character_base_l(&mut model, &quadratic, 2).unwrap_or_else(|error| panic!("{error}"));
        let quartic =
            character_base_l(&mut model, &character, 4).unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(quadratic_twice, quartic);
        assert_eq!(quartic.last(), Some(&vec![64, 0, 0, 0]));
        let twice_then_thrice =
            character_base_l(&mut model, &quadratic, 3).unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(
            twice_then_thrice,
            character_base_l(&mut model, &character, 6).unwrap_or_else(|error| panic!("{error}"))
        );
    }
}
