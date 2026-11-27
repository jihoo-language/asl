use std::collections::BTreeSet;

use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::finite_field_factor::factor_with_multiplicity;
use crate::hw_algebra::CoefficientRing;
use crate::hyperelliptic_kummer_surface::nikulin::binary_rref;
use crate::model::Model;

use super::{field, Descent, Pencil};

// The finite support is contained in the three differences of torsion
// roots, each of degree at most four. There are at most 26 global bits.
const ROOT_DEGREE: usize = 4;

pub(crate) struct Selmer {
    factors: Vec<Polynomial>,
    finite_images: Vec<(Scalar, LocalImage)>,
    infinity_image: LocalImage,
    constraints: Vec<u32>,
    basis: Vec<u32>,
}

impl Pencil {
    pub(crate) fn selmer_space(&self, model: &mut Model) -> Result<&Selmer> {
        if let Some(space) = self.selmer.get() {
            return Ok(space);
        }
        let roots = self.torsion_roots(model)?;
        let [a, b, c] = roots else {
            return Err(Error::plain(
                "K3 pencil two-Selmer space requires all three affine two-torsion points",
            ));
        };
        let field = field(&self.curve);
        let space = Selmer::new(model, field, [a, b, c])?;
        model.allocate(
            9 + space.factors.len() * (6 + 3 * field.cells())
                + space.constraints.len()
                + space.basis.len(),
        )?;
        Ok(self.selmer.get_or_init(|| space))
    }

    pub(crate) fn selmer_report(&self, model: &mut Model, name: &str) -> Result<String> {
        let space = self.selmer_space(model)?;
        let dimension = space.dimension();
        let rank_and_sha = if field(&self.curve).ring.prime() < 5 {
            "Mordell-Weil rank and dim Sha[2] unavailable in characteristic 3".into()
        } else {
            let relative = self.relative(model)?;
            let rank = relative.arithmetic_mordell_weil_rank(model, self)?;
            let sha = space.sha_dimension(rank)?;
            format!("Mordell-Weil rank={rank}, dim Sha[2]={sha}")
        };
        let finite: Vec<_> = space
            .finite_images
            .iter()
            .map(|(parameter, image)| (parameter, image.codes().collect::<Vec<_>>()))
            .collect();
        let infinity: Vec<_> = space.infinity_image.codes().collect();
        Ok(format!(
            "rvK3PencilSelmer[{name}] = dim Selmer={dimension}, {rank_and_sha}, finite local images (parameter, codes)={finite:?}, infinity codes={infinity:?}"
        ))
    }
}

impl Selmer {
    pub(crate) fn dimension(&self) -> usize {
        self.basis.len()
    }

    pub(crate) fn sha_dimension(&self, rank: usize) -> Result<usize> {
        self.dimension().checked_sub(rank + 2).ok_or_else(|| {
            Error::plain("K3 pencil Selmer dimension disagrees with the Mordell-Weil rank")
        })
    }

    fn new(model: &mut Model, field: Field<'_>, roots: [&Polynomial; 3]) -> Result<Self> {
        let factors = finite_support(model, field, roots)?;
        let width = factors.len() + 1;
        let squares = ResidueSquares::new(model, field)?;
        let mut constraints = Vec::new();
        let mut finite_images = Vec::new();
        for factor in &factors {
            let parameter = field.subtract(model, &field.zero(), &factor[0])?;
            let chart = local_roots(model, field, roots, Some(&parameter))?;
            let image = LocalImage::compute(model, field, &chart, &squares)?;
            let mut symbols = vec![2]; // The constant nonsquare.
            for f in &factors {
                let value = poly::evaluate(model, field, f, &parameter)?;
                symbols.push(if field.is_zero(&value) {
                    1 // t-a is the uniformizer and its angular component is one.
                } else {
                    squares.symbol(&value, 0)
                });
            }
            image.constrain(&mut constraints, &symbols);
            finite_images.push((parameter, image));
        }
        finite_images.sort_unstable_by(|(left, _), (right, _)| left.cmp(right));

        let chart = local_roots(model, field, roots, None)?;
        let infinity_image = LocalImage::compute(model, field, &chart, &squares)?;
        let mut symbols = vec![1; width]; // A monic linear factor has order -1 at infinity.
        symbols[0] = 2;
        infinity_image.constrain(&mut constraints, &symbols);

        let original = constraints.clone();
        let pivots = binary_rref(&mut constraints, 2 * width);
        constraints.truncate(pivots.len());
        let basis: Vec<_> = (0..2 * width)
            .filter(|column| !pivots.contains(column))
            .map(|column| {
                let free_bit = 1u32 << column;
                constraints
                    .iter()
                    .zip(&pivots)
                    .filter(|(row, _)| *row & free_bit != 0)
                    .fold(free_bit, |vector, (_, pivot)| vector | 1 << pivot)
            })
            .collect();
        if basis.iter().any(|&v| !satisfies(&original, v)) {
            return Err(Error::plain("K3 pencil Selmer kernel verification failed"));
        }
        Ok(Self {
            factors,
            finite_images,
            infinity_image,
            constraints,
            basis,
        })
    }

    pub(crate) fn basis_at(&self, model: &mut Model, index: usize) -> Result<Descent> {
        let vector = self.basis.get(index).ok_or_else(|| {
            Error::plain(format!(
                "K3 pencil Selmer basis index {index} is outside dimension {}",
                self.basis.len()
            ))
        })?;
        self.decode(model, *vector)
    }

    pub(super) fn width(&self) -> usize {
        2 * (self.factors.len() + 1)
    }

    pub(super) fn decode(&self, model: &mut Model, vector: u32) -> Result<Descent> {
        let width = self.factors.len() + 1;
        let left = vector & ((1 << width) - 1);
        let right = vector >> width;
        model.spend(3 * (1 + self.factors.len()))?;
        Ok(Descent::from_coordinates([left, right, left ^ right].map(
            |bits| {
                let support = self
                    .factors
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| bits & (1 << (i + 1)) != 0)
                    .map(|(_, f)| f.clone())
                    .collect();
                (bits & 1 != 0, support)
            },
        )))
    }

    pub(crate) fn contains(&self, model: &mut Model, descent: &Descent) -> Result<bool> {
        Ok(self.encode_local(model, descent)?.is_some())
    }

    pub(super) fn encode_local(&self, model: &mut Model, descent: &Descent) -> Result<Option<u32>> {
        let mut coordinates = [0u32; 3];
        for (bits, (nonsquare, support)) in coordinates.iter_mut().zip(descent.coordinates()) {
            *bits = u32::from(nonsquare);
            for factor in support {
                model.spend(1 + self.factors.len())?;
                let Ok(i) = self.factors.binary_search(factor) else {
                    return Ok(None);
                };
                *bits |= 1 << (i + 1);
            }
        }
        let [a, b, c] = coordinates;
        let vector = a | (b << (self.factors.len() + 1));
        model.spend(self.constraints.len())?;
        Ok((a ^ b == c && satisfies(&self.constraints, vector)).then_some(vector))
    }
}

fn satisfies(rows: &[u32], vector: u32) -> bool {
    rows.iter()
        .all(|&row| (row & vector).count_ones().is_multiple_of(2))
}

fn finite_support(
    model: &mut Model,
    field: Field<'_>,
    roots: [&Polynomial; 3],
) -> Result<Vec<Polynomial>> {
    let mut support = BTreeSet::new();
    for (i, a) in roots.iter().enumerate() {
        for b in &roots[i + 1..] {
            let difference = poly::sum(model, field, a, b, true)?;
            for (degree, _, factor) in factor_with_multiplicity(model, field, &difference)? {
                if degree != 1 {
                    return Err(Error::plain(
                        "K3 pencil two-Selmer space currently requires rational bad places",
                    ));
                }
                support.insert(factor);
            }
        }
    }
    Ok(support.into_iter().collect())
}

fn local_roots(
    model: &mut Model,
    field: Field<'_>,
    roots: [&Polynomial; 3],
    parameter: Option<&Scalar>,
) -> Result<[[Scalar; ROOT_DEGREE + 1]; 3]> {
    let mut charts = std::array::from_fn(|_| std::array::from_fn(|_| field.zero()));
    for (root, chart) in roots.into_iter().zip(&mut charts) {
        if let Some(a) = parameter {
            // Horner substitution t=a+s, preserving the uniformizer powers.
            for coefficients in root.iter().rev() {
                for i in (0..chart.len()).rev() {
                    let product = field.multiply(model, &chart[i], a)?;
                    chart[i] = if i == 0 {
                        field.add(model, &product, coefficients)?
                    } else {
                        field.add(model, &product, &chart[i - 1])?
                    };
                }
            }
        } else {
            // At infinity x'=s^4 x and y'=s^6 y. The change in each
            // coordinate is a square in the completion.
            for (out, coefficient) in chart.iter_mut().rev().zip(root) {
                out.clone_from(coefficient);
            }
        }
    }
    Ok(charts)
}

struct ResidueSquares {
    prime: usize,
    squares: Vec<bool>,
}

impl ResidueSquares {
    fn new(model: &mut Model, field: Field<'_>) -> Result<Self> {
        let count = usize::try_from(field.ring.size_of_residue_field()?)
            .map_err(|_| Error::plain("K3 pencil Selmer residue field exceeds the work limit"))?;
        model.spend(count.saturating_mul(field.ring.degree()))?;
        field.check_polynomial(count)?;
        let mut table = Self {
            prime: field.ring.prime() as usize,
            squares: vec![false; count],
        };
        for i in 1..count {
            let root = field.ring.coordinate_row(i);
            let square = field.multiply(model, &root, &root)?;
            let index = table.index(&square);
            table.squares[index] = true;
        }
        Ok(table)
    }

    fn index(&self, value: &Scalar) -> usize {
        value
            .iter()
            .fold(0, |i, &digit| i * self.prime + digit as usize)
    }

    fn symbol(&self, angular: &Scalar, order: usize) -> u8 {
        (order % 2) as u8 | (u8::from(!self.squares[self.index(angular)]) << 1)
    }
}

#[derive(Clone, Copy)]
struct LocalImage(u16);

impl LocalImage {
    fn compute(
        model: &mut Model,
        field: Field<'_>,
        roots: &[[Scalar; ROOT_DEGREE + 1]; 3],
        squares: &ResidueSquares,
    ) -> Result<Self> {
        let mut image = Self(1); // O and any point whose x has negative valuation.
        for (i, a) in roots.iter().enumerate() {
            let mut symbols = [0; 3];
            for (j, b) in roots.iter().enumerate().filter(|(j, _)| *j != i) {
                let (order, (u, v)) = a
                    .iter()
                    .zip(b)
                    .enumerate()
                    .find(|(_, (u, v))| u != v)
                    .ok_or_else(|| Error::plain("K3 pencil local roots coincide"))?;
                symbols[j] = squares.symbol(&field.subtract(model, u, v)?, order);
            }
            symbols[i] = symbols[0] ^ symbols[1] ^ symbols[2];
            image.insert(symbols);
        }

        // Expand only prefixes containing two or more roots. A prefix
        // containing just e_i forces the torsion point's two other symbols,
        // so the product condition forces the entire torsion triple.
        let mut prefixes = vec![(0b111u8, [0u8; 3])];
        for order in 0..=ROOT_DEGREE {
            if prefixes.is_empty() {
                break;
            }
            let mut next = Vec::new();
            for (active, known) in prefixes {
                for i in 0..squares.squares.len() {
                    model.spend(1 + field.ring.degree())?;
                    let residue = field.ring.coordinate_row(i);
                    let mut matches = 0u8;
                    let mut symbols = known;
                    for (j, root) in roots.iter().enumerate() {
                        if active & (1 << j) == 0 {
                            continue;
                        }
                        if residue == root[order] {
                            matches |= 1 << j;
                        } else {
                            let angular = field.subtract(model, &residue, &root[order])?;
                            symbols[j] = squares.symbol(&angular, order);
                        }
                    }
                    if matches.count_ones() >= 2 {
                        next.push((matches, symbols));
                    } else if matches == 0 && symbols[0] ^ symbols[1] ^ symbols[2] == 0 {
                        image.insert(symbols);
                    }
                }
            }
            prefixes = next;
        }
        if !prefixes.is_empty()
            || image.0.count_ones() != 4
            || image
                .codes()
                .any(|a| image.codes().any(|b| image.0 & (1 << (a ^ b)) == 0))
        {
            return Err(Error::plain(
                "K3 pencil local Kummer image failed its rank check",
            ));
        }
        Ok(image)
    }

    fn insert(&mut self, symbols: [u8; 3]) {
        self.0 |= 1 << (symbols[0] | (symbols[1] << 2));
    }

    fn codes(self) -> impl Iterator<Item = u8> {
        (0..16).filter(move |&code| self.0 & (1 << code) != 0)
    }

    fn constrain(self, constraints: &mut Vec<u32>, symbols: &[u8]) {
        for check in 1u8..16 {
            if self
                .codes()
                .all(|code| (code & check).count_ones().is_multiple_of(2))
            {
                let row = [check & 3, check >> 2]
                    .into_iter()
                    .flat_map(|coordinate| symbols.iter().map(move |symbol| coordinate & symbol))
                    .enumerate()
                    .fold(0, |row, (i, symbol)| row | ((symbol.count_ones() % 2) << i));
                constraints.push(row);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rv_unram::Ring;

    #[test]
    fn membership_rejects_failed_local_conditions_and_unsupported_places() {
        let mut model = Model::new();
        let ring =
            Ring::new(&mut model, 5, 1, vec![0, 1]).unwrap_or_else(|error| panic!("{error}"));
        let field = Field { ring: &ring };
        let zero = Vec::new();
        let a = field.polynomial(&[0, 0, 1]);
        let b = field.polynomial(&[0, 0, 1, 1]);
        let space = Selmer::new(&mut model, field, [&zero, &a, &b])
            .unwrap_or_else(|error| panic!("{error}"));
        assert!(space
            .contains(&mut model, &Descent::default())
            .unwrap_or_else(|error| panic!("{error}")));
        let factor = field.polynomial(&[1, 1]);
        let invalid_at_minus_one = Descent::from_coordinates([
            (false, vec![]),
            (false, vec![factor.clone()]),
            (false, vec![factor]),
        ]);
        assert!(!space
            .contains(&mut model, &invalid_at_minus_one)
            .unwrap_or_else(|error| panic!("{error}")));
        let factor = field.polynomial(&[-1, 1]);
        let outside_support = Descent::from_coordinates([
            (false, vec![]),
            (false, vec![factor.clone()]),
            (false, vec![factor]),
        ]);
        assert!(!space
            .contains(&mut model, &outside_support)
            .unwrap_or_else(|error| panic!("{error}")));
        let nonsquare_product =
            Descent::from_coordinates([(true, vec![]), (false, vec![]), (false, vec![])]);
        assert!(!space
            .contains(&mut model, &nonsquare_product)
            .unwrap_or_else(|error| panic!("{error}")));

        let nonrational = field.polynomial(&[2, 0, 1]);
        let one = field.polynomial(&[1]);
        let failure = Selmer::new(&mut model, field, [&zero, &nonrational, &one]);
        assert!(failure
            .err()
            .expect("support contains a nonrational place")
            .to_string()
            .contains("rational bad places"));
    }
}
