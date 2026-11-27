use std::rc::Rc;

use crate::error::{Error, Result};
use crate::integer_smith::{self as smith, Matrix, Track};
use crate::model::{Model, Value};
use crate::rv_qp_quad::Wide;
use crate::syntax::{Label, ParsedExpr};

use super::super::super::{
    dynamics::{Generator, IntegralBasis, Word},
    polarized::{Divisor, Generator as Class, COEFFICIENTS},
};
use super::super::ResolvedPoint;
use super::Pencil;

pub(crate) enum Command {
    Define(Label, Label, u8),
    Coord(Label, Box<ParsedExpr>),
    Height(Label, Box<[ParsedExpr; 2]>),
    Lattice(Label),
    BsdFactors(Label, u8, i128),
}

pub(crate) struct Expression(
    pub(crate) Label,
    pub(crate) Shift,
    pub(crate) i128,
    pub(crate) Box<ParsedExpr>,
);

#[derive(Clone, Copy)]
pub(crate) enum Shift {
    Divisor,
    Point,
}

// A shift has terms of degree three in signed 128-bit inputs. Keep the
// quadratic and linear contributions wide until they have been combined.
type Integer = Wide<13>;

pub(crate) struct Mordell {
    pub(crate) pencil: Rc<Pencil>,
    zero: u8,
    zero_pair: i128,
    word: Word,
}

impl Mordell {
    pub(crate) fn new(model: &mut Model, pencil: Rc<Pencil>, zero: u8) -> Result<Self> {
        pencil.check_section(zero)?;
        let [i, j] = pencil.masks;
        let zero_pair = if incidence(pencil.origin, zero) & (1 << i) != 0 {
            1
        } else {
            -1
        };
        let word = Word::new(
            pencil.basis_name.text.as_str().into(),
            pencil.basis.clone(),
            pencil.origin,
            Box::new([Generator::Project(i), Generator::Project(j)]),
        );
        model.allocate(10 + pencil.basis_name.text.len())?;
        Ok(Self {
            pencil,
            zero,
            zero_pair,
            word,
        })
    }

    fn pairings(&self, model: &mut Model, divisor: &Divisor) -> Result<(Integer, Integer)> {
        model.spend(5 * 13)?;
        let ([i, j], coefficients) = (self.pencil.masks, divisor.coefficients());
        let left = Integer::from_i128(coefficients[index(i)]);
        let right = Integer::from_i128(coefficients[index(j)]);
        let degree = Integer::product(2, coefficients[0]).add(left)?.add(right)?;
        let section_pair = right.add(left.negate())?;
        Ok((degree, section_pair))
    }

    fn normalized_pairings(
        &self,
        model: &mut Model,
        divisor: &Divisor,
    ) -> Result<(Integer, Integer)> {
        let (degree, section_pair) = self.pairings(model, divisor)?;
        let twice_free = section_pair.add(degree.multiply(-self.zero_pair)?)?;
        let (magnitude, remainder) = twice_free.magnitude.divide_u32(2);
        if remainder != 0 {
            return Err(Error::plain(
                "marked Kummer Mordell–Weil coordinate must be integral",
            ));
        }
        let free = Integer {
            magnitude,
            negative: twice_free.negative,
        };
        Ok((degree, free))
    }

    pub(crate) fn coordinates(&self, model: &mut Model, divisor: &Divisor) -> Result<(i128, i128)> {
        let (degree, free) = self.normalized_pairings(model, divisor)?;
        Ok((small(degree)?, small(free)?))
    }

    pub(crate) fn height(
        &self,
        model: &mut Model,
        left: &Divisor,
        right: &Divisor,
    ) -> Result<i128> {
        let (_, left) = self.normalized_pairings(model, left)?;
        let (_, right) = self.normalized_pairings(model, right)?;
        model.spend(13 * 13)?;
        small(left.multiply_wide(right)?)
    }

    pub(crate) fn shift_divisor(
        &self,
        model: &mut Model,
        steps: i128,
        input: &Divisor,
    ) -> Result<Value> {
        let (degree, section_pair) = self.pairings(model, input)?;
        model.spend(13 * 13)?;
        // In doubled (H,E) coordinates, S=E_i−E_j and F=H−E_i−E_j.
        // Only the three displayed coefficients change.
        let slide = degree.multiply(steps)?;
        let fiber = slide
            .multiply(steps)?
            .multiply(2)?
            .add(section_pair.multiply(steps)?.negate())?;
        let mut output = *input.coefficients();
        let [i, j] = self.pencil.masks.map(index);
        output[0] = small(Integer::from_i128(output[0]).add(fiber.multiply(2)?)?)?;
        for (index, slide) in [(i, slide), (j, slide.negate())] {
            let correction = slide.add(fiber.negate())?.multiply(2)?;
            output[index] = small(Integer::from_i128(output[index]).add(correction)?)?;
        }
        Divisor::literal(model, &self.pencil.basis, output)?.persist(model)
    }

    pub(crate) fn shift_point(
        &self,
        model: &mut Model,
        steps: i128,
        point: Rc<ResolvedPoint>,
    ) -> Result<Value> {
        self.word.iterate(model, point, steps)
    }

    pub(crate) fn certified_trope_height(&self, model: &mut Model, mask: u8) -> Result<()> {
        self.pencil.check_section(self.zero)?;
        self.pencil.check_section(mask)?;
        let section = self.class(model, Class::Trope(usize::from(mask)))?;
        let (degree, coordinate) = self.coordinates(model, &section)?;
        if degree != 1 || !matches!(coordinate, -1 | 1) {
            return Err(Error::plain(
                "K3 pencil BSD section must have positive height",
            ));
        }

        // The pencil and its marking are defined over the coefficient field
        // of the very same curve used for the arithmetic L polynomial.
        let basis = &self.pencil.basis;
        let frobenius = basis.frobenius(model, &self.pencil.curve, 1)?;
        if (0..basis.len())
            .any(|i| (0..basis.len()).any(|j| *frobenius.matrix.entry(i, j) != u64::from(i == j)))
        {
            return Err(Error::plain(
                "K3 pencil BSD section requires a rational level-2 marking",
            ));
        }

        let frame = self.trivial_frame(model)?;
        if gram_order(model, &frame.classes, true)? != 1024 {
            return Err(inconsistent_frame());
        }
        let [fiber, zero, ..] = frame.classes.as_slice() else {
            return Err(inconsistent_frame());
        };
        let [i, j] = self.pencil.masks;
        let left = self.class(model, Class::Node(usize::from(i)))?;
        let right = self
            .class(model, Class::Node(usize::from(j)))?
            .scale(model, -1)?;
        let shift = left.add(model, &right)?;
        let zero_pair = zero.pair(model, &shift)?;
        let zero_correction = fiber.scale(model, -zero_pair)?;
        let orthogonal = shift.add(model, &zero_correction)?;
        if zero_pair != self.zero_pair || orthogonal.pair(model, &orthogonal)? != -4 {
            return Err(inconsistent_frame());
        }
        for class in &frame.classes {
            if orthogonal.pair(model, class)? != 0 {
                return Err(inconsistent_frame());
            }
        }
        if section.pair(model, &orthogonal)? != 2 * coordinate {
            return Err(inconsistent_frame());
        }

        // Projection from either base node preserves every plane through the
        // secant. Their product permutes the components of each geometric
        // fiber. Its verified actions on F, S and O imply G^k(O)=O+kS+uF.
        // Intersect with any component C and use a k>0 fixing C to get S.C=0.
        // Thus S-(S.O)F is orthogonal to the complete trivial lattice, even
        // when the surface has additional vertical components. The marked
        // frame and this vector span the marked rank-17 space, so projection
        // of this trope has Shioda height (section.S')²/4 = 1.
        let zero_translation = fiber.scale(model, 2 - zero_pair)?;
        let shifted_zero = zero.add(model, &shift)?.add(model, &zero_translation)?;
        let four_fibers = fiber.scale(model, 4)?;
        let shifted_shift = shift.add(model, &four_fibers)?;
        for (input, expected) in [
            (fiber, fiber),
            (&shift, &shifted_shift),
            (zero, &shifted_zero),
        ] {
            let actual = self.word.apply(model, input)?;
            if actual.coefficients() != expected.coefficients() {
                return Err(inconsistent_frame());
            }
        }
        Ok(())
    }

    pub(crate) fn zero_mask(&self) -> u8 {
        self.zero
    }

    pub(crate) fn lattice(&self, model: &mut Model, name: &str) -> Result<String> {
        let frame = self.trivial_frame(model)?;
        let trivial_det = gram_order(model, &frame.classes, true)?;
        let marking = IntegralBasis::new(model, self.pencil.origin)?;
        let marked_classes = marking
            .vectors()
            .map(|row| Divisor::literal(model, &self.pencil.basis, row))
            .collect::<Result<Vec<_>>>()?;
        let marked_det = gram_order(model, &marked_classes, false)?;
        let mut inclusion = Matrix::new(model, COEFFICIENTS, frame.classes.len())?;
        for (column, class) in frame.classes.iter().enumerate() {
            for (row, value) in marking
                .coordinates(model, class.coefficients())?
                .into_iter()
                .enumerate()
            {
                inclusion.write(model, row, column, value)?;
            }
        }
        let quotient = smith::reduce(model, inclusion, Track::Both, Track::None)?;
        let torsion_factors = quotient
            .factors
            .iter()
            .enumerate()
            .filter(|(_, d)| **d > 1)
            .map(|(i, &d)| (i, d))
            .collect::<Vec<_>>();
        if trivial_det != 1024
            || marked_det != 64
            || quotient.factors.len() != 16
            || torsion_factors.len() != 2
            || torsion_factors.iter().any(|&(_, d)| d != 2)
        {
            return Err(inconsistent_frame());
        }

        let mut torsion = Vec::new();
        let mut represented = 0u8;
        let mut free = Vec::new();
        for mask in masks(self.pencil.section_masks) {
            let trope = self.class(model, Class::Trope(usize::from(mask)))?;
            let (degree, coordinate) = self.coordinates(model, &trope)?;
            match (degree, coordinate) {
                (1, 0) => {
                    // Check that the four zero-free sections represent the
                    // four distinct classes of the actual Smith quotient.
                    let coordinates = marking.coordinates(model, trope.coefficients())?;
                    let character =
                        torsion_factors
                            .iter()
                            .enumerate()
                            .fold(0u8, |bits, (k, &(i, _))| {
                                let residue = quotient
                                    .rows
                                    .coordinate(i)
                                    .iter()
                                    .fold(0, |residue, (&j, &c)| {
                                        residue ^ ((c & 1) & (coordinates[j] & 1))
                                    });
                                bits | ((residue as u8) << k)
                            });
                    represented |= 1 << character;
                    torsion.push(mask);
                }
                (1, -1 | 1) => free.push((mask, coordinate)),
                _ => return Err(inconsistent_frame()),
            }
        }
        if torsion.len() != 4 || represented != 0b1111 || free.len() != 4 {
            return Err(inconsistent_frame());
        }

        let [i, j] = self.pencil.masks;
        let left = self.class(model, Class::Node(usize::from(i)))?;
        let right = self
            .class(model, Class::Node(usize::from(j)))?
            .scale(model, -1)?;
        let translation = left.add(model, &right)?;
        let (degree, shift) = self.coordinates(model, &translation)?;
        let shift_height = self.height(model, &translation, &translation)?;
        if (degree, shift, shift_height) != (0, -2, 4)
            || translation.pair(model, &translation)? != -4
        {
            return Err(inconsistent_frame());
        }
        for vertical in frame
            .classes
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != 1)
            .map(|(_, class)| class)
        {
            if translation.pair(model, vertical)? != 0 {
                return Err(inconsistent_frame());
            }
        }

        Ok(format!(
            "rvK3MordellLattice[{name}] = nodes={:?}, zero={}, D4 stars={:?}, A1 nodes={:?}, trivial=U+D4(-1)^2+A1(-1)^6, trivial rank=16, trivial |det|={trivial_det}, marked rank=17, marked |det|={marked_det}, quotient=Z+(Z/2)^2, torsion sections={torsion:?}, free sections={free:?}, shift free={shift}, shift height={shift_height}",
            self.pencil.masks, self.zero, frame.stars, frame.singles,
        ))
    }

    fn class(&self, model: &mut Model, class: Class) -> Result<Divisor> {
        Divisor::generator(model, &self.pencil.basis, class)
    }

    fn trivial_frame(&self, model: &mut Model) -> Result<Frame> {
        let [i, j] = self.pencil.masks;
        let endpoints = (1 << i) | (1 << j);
        let zero = incidence(self.pencil.origin, self.zero);
        let mut stars = Vec::new();
        let mut occupied = endpoints;
        let mut classes = vec![
            self.pencil.fiber_class(model)?,
            self.class(model, Class::Trope(usize::from(self.zero)))?,
        ];
        for trope in 0..16 {
            let support = incidence(self.pencil.origin, trope);
            if support & endpoints != endpoints {
                continue;
            }
            let arms = support & !endpoints;
            if occupied & arms != 0 || arms.count_ones() != 4 || (arms & zero).count_ones() != 1 {
                return Err(inconsistent_frame());
            }
            occupied |= arms;
            stars.push((trope, masks(arms).collect()));
            classes.push(self.class(model, Class::Trope(usize::from(trope)))?);
            for arm in masks(arms & !zero) {
                classes.push(self.class(model, Class::Node(usize::from(arm)))?);
            }
        }
        let singles: Vec<_> = masks(!occupied).collect();
        if stars.len() != 2 || singles.len() != 6 {
            return Err(inconsistent_frame());
        }
        for &node in &singles {
            let class = self.class(model, Class::Node(usize::from(node)))?;
            classes.push(if zero & (1 << node) == 0 {
                class
            } else {
                let neg = class.scale(model, -1)?;
                classes[0].add(model, &neg)?
            });
        }
        Ok(Frame {
            classes,
            stars,
            singles,
        })
    }
}

struct Frame {
    classes: Vec<Divisor>,
    stars: Vec<(u8, Vec<u8>)>,
    singles: Vec<u8>,
}

fn incidence(origin: u16, trope: u8) -> u16 {
    masks(origin).fold(0, |out, node| out | (1 << (node ^ trope)))
}

fn masks(set: u16) -> impl Iterator<Item = u8> {
    (0..16).filter(move |&i| set & (1 << i) != 0)
}

fn index(node: u8) -> usize {
    1 + usize::from(node)
}

fn gram_order(model: &mut Model, classes: &[Divisor], trivial: bool) -> Result<i128> {
    let mut gram = Matrix::new(model, classes.len(), classes.len())?;
    for (i, left) in classes.iter().enumerate() {
        for (j, right) in classes[i..].iter().enumerate().map(|(j, v)| (j + i, v)) {
            let value = left.pair(model, right)?;
            if trivial {
                let expected = match (i, j) {
                    (0, 0) => 0,
                    (0, 1) | (2, 3..=5) | (6, 7..=9) => 1,
                    (i, j) if i == j => -2,
                    _ => 0,
                };
                if value != expected {
                    return Err(inconsistent_frame());
                }
            }
            gram.write(model, i, j, value)?;
            if i != j {
                gram.write(model, j, i, value)?;
            }
        }
    }
    let gram = smith::reduce(model, gram, Track::None, Track::None)?;
    if gram.factors.len() != classes.len() {
        return Err(inconsistent_frame());
    }
    gram.factors.into_iter().try_fold(1i128, |order, factor| {
        order.checked_mul(factor).ok_or_else(inconsistent_frame)
    })
}

fn small(value: Integer) -> Result<i128> {
    value.signed().ok_or_else(|| {
        Error::plain("marked Kummer Mordell–Weil result exceeds signed 128-bit range")
    })
}

fn inconsistent_frame() -> Error {
    Error::plain("K3 pencil has inconsistent marked Mordell–Weil lattice data")
}
