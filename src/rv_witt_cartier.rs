use std::array;

use crate::binary_linear;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_poly::{self, Polynomial};

// Polynomial coefficients of 1, y, w, and yw over F2[x].
type NormalForm = [Polynomial; 4];

#[derive(Clone, Copy)]
pub(crate) struct Semigroup {
    weights: [u128; 4],
    genus: usize,
}

#[derive(Clone, Copy)]
struct Function {
    degree: usize,
    component: usize,
}

impl Semigroup {
    pub(crate) fn new(m: usize, upper: usize, genus: usize) -> Self {
        let b = 2 * upper as u128 - m as u128;
        let a = 2 * m as u128;
        Self {
            weights: [0, a, b, a + b],
            genus,
        }
    }

    pub(crate) fn gaps(self, model: &mut Model) -> Result<Vec<u128>> {
        model.spend(self.genus.saturating_mul(1 + self.genus.ilog2() as usize))?;
        binary_linear::check_storage(self.genus.saturating_mul(2))?;
        let mut gaps: Vec<_> = self
            .weights
            .into_iter()
            .flat_map(|weight| (weight % 4..weight).step_by(4))
            .collect();
        gaps.sort_unstable();
        Ok(gaps)
    }

    pub(crate) fn basis(self, model: &mut Model) -> Result<Vec<(usize, usize, usize)>> {
        model.spend(self.genus.saturating_mul(1 + self.genus.ilog2() as usize))?;
        binary_linear::check_storage(self.genus.saturating_mul(7))?;
        let Some(bound) = self.differential_bound() else {
            return Ok(Vec::new());
        };
        let mut functions: Vec<_> = self
            .weights
            .into_iter()
            .enumerate()
            .flat_map(|(component, first)| {
                (first..=bound)
                    .step_by(4)
                    .enumerate()
                    .map(move |(degree, weight)| (weight, Function { degree, component }))
            })
            .collect();
        functions.sort_unstable_by_key(|&(weight, _)| weight);
        Ok(functions
            .into_iter()
            .map(|(_, f)| (f.degree, f.component % 2, f.component / 2))
            .collect())
    }

    fn differential_bound(self) -> Option<u128> {
        (2 * self.genus as u128).checked_sub(2)
    }

    fn count_before(self, weight: u128) -> usize {
        self.weights
            .into_iter()
            .filter(|&first| first < weight)
            .map(|first| ((weight - first - 1) / 4 + 1) as usize)
            .sum()
    }

    fn function(self, weight: u128) -> Option<Function> {
        let component = self
            .weights
            .iter()
            .position(|first| first % 4 == weight % 4)?;
        let degree = usize::try_from(weight.checked_sub(self.weights[component])? / 4).ok()?;
        Some(Function { degree, component })
    }

    fn indexed(self, index: usize) -> Result<Function> {
        let Some(mut hi) = self.differential_bound().filter(|_| index < self.genus) else {
            return Err(Error::plain("Witt differential index is out of range"));
        };
        let mut lo = 0;
        while lo < hi {
            let mid = (lo + hi) / 2;
            if self.count_before(mid + 1) <= index {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        self.function(lo)
            .ok_or_else(|| Error::plain("Witt differential basis weight is missing"))
    }

    fn weight(self, degree: usize, component: usize) -> u128 {
        4 * degree as u128 + self.weights[component]
    }
}

pub(crate) struct Cartier {
    semigroup: Semigroup,
    squares: [NormalForm; 4],
    derivatives: [NormalForm; 4],
}

impl Cartier {
    pub(crate) fn build(
        model: &mut Model,
        semigroup: Semigroup,
        f: &Polynomial,
        g: &Polynomial,
    ) -> Result<Self> {
        let m = f.last().unwrap()[0];
        let (a, b) = reduced_equation(model, f, g, m, semigroup.weights[2])?;
        let sum = rv_poly::add(model, &a, &b)?;
        let fs = rv_poly::multiply(model, f, &sum)?;
        let fb = rv_poly::multiply(model, f, &b)?;
        let ys = rv_poly::add(model, &sum, &fb)?;
        let empty: Polynomial = [].into();
        let one = rv_poly::monomial(vec![0]);
        // Use y^2=y+f and w^2=w+A+B*y.
        let squares = [
            [one.clone(), empty.clone(), empty.clone(), empty.clone()],
            [f.clone(), one.clone(), empty.clone(), empty.clone()],
            [a.clone(), b.clone(), one.clone(), empty.clone()],
            [fs, ys, f.clone(), one],
        ];
        let df = derivative(model, f)?;
        let da = derivative(model, &a)?;
        let db = derivative(model, &b)?;
        let bdf = rv_poly::multiply(model, &b, &df)?;
        let dw = rv_poly::add(model, &da, &bdf)?;
        let fdb = rv_poly::multiply(model, f, &db)?;
        let d_yw = rv_poly::add(model, &dw, &db)?;
        // D_x(y)=f' and D_x(w)=A'+B'y+B*f'.
        let derivatives = [
            array::from_fn(|_| empty.clone()),
            [df.clone(), empty.clone(), empty.clone(), empty.clone()],
            [dw, db, empty.clone(), empty.clone()],
            [fdb, d_yw, df, empty],
        ];
        let cells = squares
            .iter()
            .chain(&derivatives)
            .flatten()
            .map(|poly| rv_poly::cells(poly))
            .sum();
        binary_linear::check_storage(cells)?;
        model.allocate(cells)?;
        Ok(Self {
            semigroup,
            squares,
            derivatives,
        })
    }

    pub(crate) fn column(&self, model: &mut Model, index: usize) -> Result<Vec<usize>> {
        // The four progressions give a basis index without enumerating H^0.
        model.spend(u128::BITS as usize / 2)?;
        let Function { degree, component } = self.semigroup.indexed(index)?;
        let mut remainder: NormalForm = array::from_fn(|_| [].into());
        for (slot, term) in remainder.iter_mut().zip(&self.derivatives[component]) {
            *slot = shift(model, term, degree)?;
        }
        if degree % 2 == 1 {
            let dx = rv_poly::monomial(vec![degree - 1]);
            remainder[component] = rv_poly::add(model, &remainder[component], &dx)?;
        }

        check_normal_form(&remainder)?;
        // Cartier(H dx)/dx is the unique square root of D_x(H).
        let mut result = Vec::new();
        while let Some(weight) = self.leading_weight(&remainder) {
            model.spend(1)?;
            if weight % 2 != 0 {
                return Err(square_error());
            }
            let weight = weight / 2;
            let Function { degree, component } = self
                .semigroup
                .function(weight)
                .filter(|_| {
                    self.semigroup
                        .differential_bound()
                        .is_some_and(|bound| weight <= bound)
                })
                .ok_or_else(square_error)?;
            binary_linear::check_storage(result.len() + 1)?;
            result.push(self.semigroup.count_before(weight));

            let twice = degree.checked_mul(2).ok_or_else(exponent_error)?;
            for (slot, term) in remainder.iter_mut().zip(&self.squares[component]) {
                if !term.is_empty() {
                    let shifted = shift(model, term, twice)?;
                    *slot = rv_poly::add(model, slot, &shifted)?;
                }
            }
            check_normal_form(&remainder)?;
            // A normalized monomial square has exactly twice the weight.
            // Cancelling it must strictly decrease the pole of the remainder.
            if self
                .leading_weight(&remainder)
                .is_some_and(|new| new >= 2 * weight)
            {
                return Err(square_error());
            }
        }
        result.reverse();
        Ok(result)
    }

    fn leading_weight(&self, form: &NormalForm) -> Option<u128> {
        form.iter()
            .enumerate()
            .filter_map(|(component, p)| Some(self.semigroup.weight(p.last()?[0], component)))
            .max()
    }
}

fn reduced_equation(
    model: &mut Model,
    f: &Polynomial,
    g: &Polynomial,
    m: usize,
    expected: u128,
) -> Result<(Polynomial, Polynomial)> {
    let (mut a, mut b) = (g.clone(), f.clone());
    // The poles of A and B*y in the quadratic field have opposite parities.
    while let Some(leading) = a.last().map(|p| p[0]) {
        let odd = b.last().map_or(0, |p| 2 * p[0] as u128 + m as u128);
        if 2 * leading as u128 <= odd {
            break;
        }
        if leading % 2 == 0 {
            let correction: Polynomial = vec![[leading / 2].into(), [leading].into()].into();
            a = rv_poly::add(model, &a, &correction)?;
        } else {
            let difference = leading.checked_sub(m).ok_or_else(square_error)?;
            let constant = shift(model, f, difference)?;
            a = rv_poly::add(model, &a, &constant)?;
            let even_y = rv_poly::monomial(vec![difference]);
            let y = rv_poly::monomial(vec![difference / 2]);
            b = rv_poly::add(model, &b, &even_y)?;
            b = rv_poly::add(model, &b, &y)?;
        }
    }
    if b.last().map(|p| 2 * p[0] as u128 + m as u128) != Some(expected) {
        return Err(Error::plain(
            "Witt normal form disagrees with the ramification break",
        ));
    }
    Ok((a, b))
}

fn shift(model: &mut Model, poly: &Polynomial, amount: usize) -> Result<Polynomial> {
    if poly.is_empty() || amount == 0 {
        return Ok(poly.clone());
    }
    model.spend(rv_poly::cells(poly))?;
    poly.iter()
        .map(|term| {
            let degree = term[0].checked_add(amount).ok_or_else(exponent_error)?;
            Ok([degree].into())
        })
        .collect()
}

fn derivative(model: &mut Model, poly: &Polynomial) -> Result<Polynomial> {
    model.spend(rv_poly::cells(poly))?;
    Ok(poly
        .iter()
        .filter(|term| term[0] % 2 == 1)
        .map(|term| [term[0] - 1].into())
        .collect())
}

fn check_normal_form(form: &NormalForm) -> Result<()> {
    binary_linear::check_storage(form.iter().map(|poly| rv_poly::cells(poly)).sum())
}

fn exponent_error() -> Error {
    Error::plain("Witt differential polynomial degree is out of range")
}

fn square_error() -> Error {
    Error::plain("Witt differential could not be reduced in the basis of holomorphic squares")
}
