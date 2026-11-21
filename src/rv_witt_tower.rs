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
        let mut one = vec![0; width];
        one[0] = 1;
        let mut polynomial = vec![one];
        for n in 1..=degree {
            let mut numerator = vec![0; width];
            for k in 1..=n {
                add_product(model, &mut numerator, &sums[k], &polynomial[n - k])?;
            }
            for coefficient in &mut numerator {
                let divisor = n as i128;
                if *coefficient % divisor != 0 {
                    return Err(Error::plain(
                        "Witt tower Newton identity is not integral in the power basis",
                    ));
                }
                *coefficient /= divisor;
            }
            polynomial.push(numerator);
        }
        model.allocate(1 + (1 + 2 * width) * polynomial.len())?;
        let polynomial: Character = polynomial.into();
        self.characters.borrow_mut()[height - 1] = Some(polynomial.clone());
        Ok(polynomial)
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
            add_product(model, &mut polynomial[i + j], a, b)?;
        }
    }
    Ok(polynomial)
}

fn add_product(model: &mut Model, out: &mut [i128], left: &[i128], right: &[i128]) -> Result<()> {
    let width = out.len();
    model.spend(width)?;
    for (i, &a) in left.iter().enumerate().filter(|&(_, &c)| c != 0) {
        model.spend(width)?;
        for (j, &b) in right.iter().enumerate().filter(|&(_, &c)| c != 0) {
            let product = a.checked_mul(b).ok_or_else(coefficient_error)?;
            let target = &mut out[(i + j) & (width - 1)];
            *target = if i + j < width {
                target.checked_add(product)
            } else {
                target.checked_sub(product)
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
