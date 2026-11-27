use std::rc::Rc;

use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::{Model, Value};
use crate::syntax::Label;

use super::super::eigenpacket::{Element, Packet};
use super::super::hecke::Operator;
use super::super::INTEGER_WORDS;
use super::{binomial_table, Grade, Polynomial, Sign};

pub(crate) struct Declare(
    pub(crate) Label,
    pub(crate) Label,
    pub(crate) Sign,
    pub(crate) usize,
);

pub(crate) struct Expression(pub(crate) Label, pub(crate) Query);

pub(crate) enum Query {
    Coefficient(usize),
    Ratio(usize, usize),
}

pub(crate) struct Eigenperiod {
    pub(crate) source: Label,
    pub(super) packet: Rc<Packet>,
    pub(super) grade: Grade,
    pub(super) coefficients: Vec<Element>,
}

impl Eigenperiod {
    pub(crate) fn new(
        model: &mut Model,
        source: Label,
        packet: Rc<Packet>,
        sign: Sign,
        pivot: usize,
    ) -> Result<Self> {
        match packet.degree() {
            1 => {}
            2 => {
                packet.quadratic_order(model)?;
            }
            _ => {
                return Err(Error::plain(
                    "Hecke eigenperiod requires a rational or irreducible quadratic packet",
                ))
            }
        }
        let grade = Grade::new(packet.weight, sign)?;
        let w = grade.degree();
        if pivot > w {
            return Err(Error::plain(
                "Hecke eigenperiod normalization index is outside the homogeneous degree",
            ));
        }
        if !sign.includes(pivot) {
            return Err(Error::plain(
                "Hecke eigenperiod normalization coefficient has the wrong parity",
            ));
        }

        let frame = ManinFrame::new(model, grade)?;
        let operator = packet.generator();
        let mut eigenvalue = packet.image_element(model, operator)?;
        let mut images = frame.images(model, operator)?;
        if frame.boundary_collision(model, &packet, &images, &eigenvalue)? {
            // A separator on the cusp space can still agree with the
            // Eisenstein scalar. T₂ distinguishes the boundary in this case.
            let fallback = Operator::basis(2)?;
            eigenvalue = packet.image_element(model, &fallback)?;
            images = frame.images(model, &fallback)?;
            if frame.boundary_collision(model, &packet, &images, &eigenvalue)? {
                return Err(Error::plain(
                    "T_2 did not separate the cusp eigenperiod from the Eisenstein boundary",
                ));
            }
        }

        let n = frame.vectors.len();
        let mut matrix = vec![Vec::with_capacity(n); n];
        let negative_eigenvalue = eigenvalue.negate(model)?;
        for (j, image) in images.iter().enumerate() {
            for (i, (&index, row)) in frame.coordinates.iter().zip(&mut matrix).enumerate() {
                let mut entry = rational_element(model, &packet, image.0[index])?;
                if i == j {
                    entry = entry.add(model, &negative_eigenvalue)?;
                }
                row.push(entry);
            }
        }
        let coordinates = kernel_line(model, &packet, matrix)?;
        if coordinates[..packet.degree()].iter().all(Element::is_zero) {
            return Err(Error::plain(
                "computed eigenperiod lies in the Eisenstein boundary",
            ));
        }
        let mut coefficients = homogeneous(model, &packet, &frame.vectors, &coordinates, w)?;
        if coefficients[pivot].is_zero() {
            return Err(Error::plain(
                "Hecke eigenperiod normalization coefficient vanishes",
            ));
        }
        let scale = packet.inverse_element(model, &coefficients[pivot])?;
        for coefficient in &mut coefficients {
            *coefficient = packet.multiply_elements(model, coefficient, &scale)?;
        }
        model.allocate(1 + coefficients.len() * (1 + packet.degree()) * INTEGER_WORDS)?;
        Ok(Self {
            source,
            packet,
            grade,
            coefficients,
        })
    }

    pub(crate) fn coefficient(&self, model: &mut Model, j: usize) -> Result<Value> {
        self.coefficients
            .get(j)
            .cloned()
            .unwrap_or_else(|| Element::zero(self.packet.degree()))
            .persist(model)
    }

    pub(crate) fn critical_ratio(&self, model: &mut Model, s: usize, t: usize) -> Result<Value> {
        let w = self.grade.degree();
        let top = w + 1;
        if !(1..=top).contains(&s) || !(1..=top).contains(&t) {
            return Err(Error::plain(
                "critical L-value indices must be between 1 and k−1",
            ));
        }
        if s % 2 != t % 2 {
            return Err(Error::plain("critical L-value ratios across parities require both transcendental normalizations"));
        }
        let (j_s, j_t) = (top - s, top - t);
        if !self.grade.sign.includes(j_s) {
            return Err(Error::plain(
                "critical L-value indices do not match the eigenperiod sign",
            ));
        }
        let denominator = &self.coefficients[j_t];
        if denominator.is_zero() {
            return Err(Error::plain("the denominator critical L-value vanishes"));
        }
        let binomial = binomial_table(model, w)?;
        let factor = if (s.abs_diff(t) / 2).is_multiple_of(2) {
            binomial[w][j_t]
        } else {
            -binomial[w][j_t]
        };
        let factor = self
            .packet
            .scalar_element(model, factor, binomial[w][j_s] as u128)?;
        let inverse = self.packet.inverse_element(model, denominator)?;
        let ratio = self
            .packet
            .multiply_elements(model, &self.coefficients[j_s], &inverse)?;
        self.packet
            .multiply_elements(model, &ratio, &factor)?
            .persist(model)
    }
}

struct ManinFrame {
    vectors: Vec<Polynomial>,
    coordinates: Vec<usize>,
    grade: Grade,
}

impl ManinFrame {
    fn new(model: &mut Model, grade: Grade) -> Result<Self> {
        let basis = model.period_basis(grade)?;
        let mut vectors = basis.vectors.clone();
        let mut coordinates = basis.free.clone();
        if grade.sign == Sign::Plus {
            let mut boundary = Polynomial::zero(grade);
            boundary.0[0] = Rational::integer(-1);
            boundary.0[grade.degree()] = Rational::integer(1);
            vectors.push(boundary);
            coordinates.push(grade.degree());
        }
        Ok(Self {
            vectors,
            coordinates,
            grade,
        })
    }

    fn images(&self, model: &mut Model, operator: &Operator) -> Result<Vec<Polynomial>> {
        self.vectors
            .iter()
            .map(|p| p.hecke_full(model, self.grade, operator))
            .collect()
    }

    fn boundary_collision(
        &self,
        model: &mut Model,
        packet: &Packet,
        images: &[Polynomial],
        eigenvalue: &Element,
    ) -> Result<bool> {
        if self.grade.sign == Sign::Minus {
            return Ok(false);
        }
        let boundary = self.vectors.last().unwrap();
        let image = images.last().unwrap();
        let eigen = image.0[self.grade.degree()];
        let expected = boundary.scale(model, eigen)?;
        if !image.compare(model, &expected)?.is_eq() {
            return Err(Error::plain(
                "Hecke action failed to preserve the Manin boundary",
            ));
        }
        let element = rational_element(model, packet, eigen)?;
        Ok(element.compare(model, eigenvalue)?.is_eq())
    }
}

fn kernel_line(
    model: &mut Model,
    packet: &Packet,
    mut matrix: Vec<Vec<Element>>,
) -> Result<Vec<Element>> {
    let n = matrix.len();
    let mut pivots = Vec::new();
    for column in 0..n {
        let rank = pivots.len();
        let Some(pivot) = (rank..n).find(|&row| !matrix[row][column].is_zero()) else {
            continue;
        };
        matrix.swap(rank, pivot);
        let inverse = packet.inverse_element(model, &matrix[rank][column])?;
        let normalized = matrix[rank]
            .iter()
            .map(|entry| packet.multiply_elements(model, entry, &inverse))
            .collect::<Result<Vec<_>>>()?;
        for (i, row) in matrix.iter_mut().enumerate() {
            if i == rank || row[column].is_zero() {
                continue;
            }
            let multiple = row[column].negate(model)?;
            for (entry, value) in row.iter_mut().zip(&normalized) {
                let term = packet.multiply_elements(model, value, &multiple)?;
                *entry = entry.add(model, &term)?;
            }
        }
        matrix[rank] = normalized;
        pivots.push(column);
    }
    if pivots.len() + 1 != n {
        return Err(Error::plain(
            "Hecke action on the full Manin space did not have a unique cusp eigenline",
        ));
    }
    let free = (0..n).find(|column| !pivots.contains(column)).unwrap();
    let mut output = vec![Element::zero(packet.degree()); n];
    output[free] = packet.scalar_element(model, 1, 1)?;
    for (i, &column) in pivots.iter().enumerate() {
        output[column] = matrix[i][free].negate(model)?;
    }
    Ok(output)
}

fn homogeneous(
    model: &mut Model,
    packet: &Packet,
    basis: &[Polynomial],
    coordinates: &[Element],
    degree: usize,
) -> Result<Vec<Element>> {
    let mut output = vec![Element::zero(packet.degree()); degree + 1];
    for (polynomial, coordinate) in basis.iter().zip(coordinates) {
        for (coefficient, &entry) in output.iter_mut().zip(&polynomial.0) {
            if entry.numerator != 0 && !coordinate.is_zero() {
                let factor = rational_element(model, packet, entry)?;
                let term = packet.multiply_elements(model, &factor, coordinate)?;
                *coefficient = coefficient.add(model, &term)?;
            }
        }
    }
    Ok(output)
}

fn rational_element(model: &mut Model, packet: &Packet, q: Rational) -> Result<Element> {
    packet.scalar_element(model, q.numerator, q.denominator as u128)
}
