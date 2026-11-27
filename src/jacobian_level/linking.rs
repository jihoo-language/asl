use std::rc::Rc;

use super::brauer::{exterior_square, Exterior};
use super::Basis;
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::hyperelliptic::Residues;
use crate::model::Model;
use crate::prime_linear;
use crate::rv_unram::{add_mod, mul_mod, sub_mod};
use crate::syntax::{Label, ParsedExpr};

pub(crate) type Fixed = [u64; 6];

pub(crate) enum Statement {
    Define(Label, Label, Label, usize),
    Link(Box<[ParsedExpr; 2]>),
    Form(Label),
}

pub(crate) enum Expression {
    Fixed(Label, usize, [i128; 6]),
    Scale(i128, Box<ParsedExpr>),
    Change(Label, Box<ParsedExpr>),
    Nonzero(Box<[ParsedExpr; 2]>),
}

pub(crate) struct Pairing {
    pub(super) basis: Rc<Basis>,
    pub(super) base: Label,
    pub(super) degree: usize,
    pub(super) ring: Residues,
    pub(super) difference: Exterior,
}

pub(crate) struct Change {
    pub(super) matrix: Exterior,
}

pub(crate) struct Value {
    numerator: u64,
    order: u64,
}

impl Value {
    pub(crate) fn nonzero(&self) -> bool {
        self.numerator != 0
    }

    pub(crate) fn report(&self) -> String {
        let witness = self.nonzero();
        format!(
            "rvKummerBrauerLink = numerator={}, order={}, nonalgebraic_witness={witness}",
            self.numerator, self.order,
        )
    }
}

impl Pairing {
    pub(crate) fn new(
        model: &mut Model,
        basis: Rc<Basis>,
        base: Label,
        curve: &Curve,
        degree: usize,
    ) -> Result<Self> {
        if basis.prime.get() == 2 {
            return Err(Error::plain(
                "Kummer Brauer linking requires an odd torsion prime",
            ));
        }
        if degree == 0 {
            return Err(Error::plain(
                "Kummer Brauer linking extension degree must be positive",
            ));
        }
        if basis.len() != 4 {
            return Err(Error::plain(
                "Kummer Brauer linking requires a genus-two marking",
            ));
        }
        let ring = Residues::new(basis.prime, basis.height)?;
        let frobenius = basis.frobenius(model, curve, degree)?;
        let inverse = ring.inverse_unit(model, &frobenius.multiplier)?;
        let unit = scale_matrix(
            exterior_square(&frobenius.matrix, ring.modulus),
            inverse,
            ring.modulus,
        );
        verify_isometry(model, ring.modulus, &unit)?;
        let mut difference = unit.map(|row| row.map(|v| sub_mod(0, v, ring.modulus)));
        for (i, row) in difference.iter_mut().enumerate() {
            row[i] = add_mod(row[i], 1, ring.modulus);
        }
        model.allocate(48 + base.text.len())?;
        Ok(Self {
            basis,
            base,
            degree,
            ring,
            difference,
        })
    }

    pub(crate) fn base_name(&self) -> &Label {
        &self.base
    }

    pub(crate) fn modulus(&self, level: usize) -> Result<u64> {
        if level == 0 || level > self.ring.precision {
            return Err(Error::plain(
                "Kummer Brauer fixed lift exponent must be between one and the marking height",
            ));
        }
        Ok(self.prime().pow(level as u32))
    }

    fn prime(&self) -> u64 {
        u64::from(self.ring.prime.get())
    }

    pub(crate) fn fixed(
        &self,
        model: &mut Model,
        level: usize,
        coordinates: [i128; 6],
    ) -> Result<Fixed> {
        let modulus = self.modulus(level)?;
        let vector = coordinates.map(|entry| entry.rem_euclid(i128::from(modulus)) as u64);
        self.require_fixed(model, modulus, vector)?;
        Ok(vector)
    }

    pub(super) fn require_fixed(
        &self,
        model: &mut Model,
        modulus: u64,
        vector: Fixed,
    ) -> Result<()> {
        model.spend(36)?;
        if image(&self.difference, vector, modulus) != [0; 6] {
            return Err(Error::plain(
                "Kummer Brauer lift is not fixed at the specified exponent",
            ));
        }
        Ok(())
    }

    pub(crate) fn link(
        &self,
        model: &mut Model,
        a: usize,
        v: Fixed,
        b: usize,
        w: Fixed,
    ) -> Result<Value> {
        let (raw, denominator) = self.raw_link(model, a, v, b, w)?;
        let gcd = crate::rv_ore::gcd(u128::from(raw), u128::from(denominator)) as u64;
        Ok(Value {
            numerator: raw / gcd,
            order: denominator / gcd,
        })
    }

    fn raw_link(
        &self,
        model: &mut Model,
        a: usize,
        v: Fixed,
        b: usize,
        w: Fixed,
    ) -> Result<(u64, u64)> {
        let m = self.modulus(a)?;
        let n = self.modulus(b)?;
        let level = a.checked_add(b).filter(|&level| level <= self.ring.precision)
            .ok_or_else(|| Error::plain("Kummer Brauer linking needs marking height at least the sum of the lift exponents"))?;
        let denominator = self.modulus(level)?;
        self.require_fixed(model, m, v)?;
        self.require_fixed(model, n, w)?;
        model.spend(42)?;
        let delta = image(&self.difference, w, denominator);
        let raw = intersection(v, delta, denominator);
        if !raw.is_multiple_of(m.max(n)) {
            return Err(Error::plain(
                "Kummer Brauer linking failed its denominator check",
            ));
        }
        Ok((raw, denominator))
    }

    fn coefficient(&self, model: &mut Model, v: Fixed, w: Fixed) -> Result<u32> {
        let (numerator, _) = self.raw_link(model, 1, v, 1, w)?;
        Ok((numerator / self.prime()) as u32)
    }

    pub(crate) fn form(&self, model: &mut Model, name: &str) -> Result<String> {
        if self.ring.precision < 2 {
            return Err(Error::plain(
                "Kummer Brauer residue linking form requires marking height at least two",
            ));
        }
        let before = self.fixed_basis(model)?;
        let mut gram = vec![vec![0; before.len()]; before.len()];
        let prime = self.ring.prime;
        for (i, &v) in before.iter().enumerate() {
            if self.coefficient(model, v, v)? != 0 {
                return Err(Error::plain(
                    "Kummer Brauer residue form failed alternation",
                ));
            }
            for (j, &w) in before.iter().enumerate().skip(i + 1) {
                let value = self.coefficient(model, v, w)?;
                let reverse = self.coefficient(model, w, v)?;
                if prime.add(value, reverse) != 0 {
                    return Err(Error::plain(
                        "Kummer Brauer residue form failed skew-symmetry",
                    ));
                }
                gram[i][j] = value;
                gram[j][i] = reverse;
            }
        }
        let (symplectic, radical) = self.symplectic_split(model, before.clone())?;
        Ok(format!(
            "rvKummerBrauerLinkForm[{name}] = prime={}, fixed_dimension={}, rank={}, radical_dimension={}, fixed_basis={before:?}, gram={gram:?}, symplectic_pairs={symplectic:?}, radical_basis={radical:?}",
            self.prime(), before.len(), symplectic.len() * 2, radical.len(),
        ))
    }

    fn fixed_basis(&self, model: &mut Model) -> Result<Vec<Fixed>> {
        let mut basis = prime_linear::Basis::default();
        let mut fixed = Vec::new();
        let prime = self.ring.prime;
        let modulus = self.prime();
        for col in 0..6 {
            let image = self
                .difference
                .iter()
                .enumerate()
                .filter_map(|(i, row)| {
                    let value = (row[col] % modulus) as u32;
                    (value != 0).then_some((i, value))
                })
                .collect();
            if let Some(witness) = basis.insert_with_witness(model, prime, image, vec![(col, 1)])? {
                let mut vector = [0; 6];
                for (index, coefficient) in witness {
                    vector[index] = u64::from(coefficient);
                }
                fixed.push(vector);
            }
        }
        Ok(fixed)
    }

    fn symplectic_split(
        &self,
        model: &mut Model,
        mut residual: Vec<Fixed>,
    ) -> Result<(Vec<[Fixed; 2]>, Vec<Fixed>)> {
        let p = self.prime();
        let prime = self.ring.prime;
        let mut pairs = Vec::new();
        while let Some((i, j, value)) = self.find_pair(model, &residual)? {
            let w = residual.remove(j);
            let v = residual.remove(i);
            let inverse = u64::from(prime_linear::invert(model, prime, value)?);
            let w = scaled(w, inverse, p);
            for x in &mut residual {
                let subtract = self.coefficient(model, *x, w)?;
                let add = self.coefficient(model, *x, v)?;
                *x = add_vectors(
                    add_vectors(*x, scaled(v, u64::from(prime.negate(subtract)), p), p),
                    scaled(w, u64::from(add), p),
                    p,
                );
            }
            pairs.push([v, w]);
        }
        Ok((pairs, residual))
    }

    fn find_pair(&self, model: &mut Model, basis: &[Fixed]) -> Result<Option<(usize, usize, u32)>> {
        for (i, &v) in basis.iter().enumerate() {
            for (j, &w) in basis.iter().enumerate().skip(i + 1) {
                let value = self.coefficient(model, v, w)?;
                if value != 0 {
                    return Ok(Some((i, j, value)));
                }
            }
        }
        Ok(None)
    }

    pub(crate) fn change(&self, model: &mut Model, source: &Self) -> Result<Change> {
        if self.base.text != source.base.text || self.degree != source.degree {
            return Err(Error::plain("Kummer Brauer marking change requires the same named base curve and extension degree"));
        }
        let m = self.ring.modulus;
        let matrix = exterior_change(model, &self.basis, &source.basis, self.ring)?;
        model.spend(2 * 6 * 6 * 6)?;
        if product(&self.difference, &matrix, m) != product(&matrix, &source.difference, m) {
            return Err(Error::plain(
                "Kummer Brauer marking change does not intertwine Frobenius",
            ));
        }
        model.allocate(36)?;
        Ok(Change { matrix })
    }
}

impl Change {
    pub(crate) fn apply(
        &self,
        model: &mut Model,
        target: &Pairing,
        level: usize,
        vector: Fixed,
    ) -> Result<Fixed> {
        let m = target.modulus(level)?;
        model.spend(36)?;
        let result = image(&self.matrix, vector, m);
        target.require_fixed(model, m, result)?;
        Ok(result)
    }
}

pub(crate) fn add_vectors(a: Fixed, b: Fixed, modulus: u64) -> Fixed {
    std::array::from_fn(|i| add_mod(a[i], b[i], modulus))
}

pub(crate) fn scaled(x: Fixed, n: u64, modulus: u64) -> Fixed {
    x.map(|v| mul_mod(v, n, modulus))
}

fn scale_matrix(m: Exterior, n: u64, modulus: u64) -> Exterior {
    m.map(|row| scaled(row, n, modulus))
}

pub(super) fn image(matrix: &Exterior, vector: Fixed, modulus: u64) -> Fixed {
    matrix.map(|row| dot(row, vector, modulus))
}

fn dot(a: Fixed, b: Fixed, modulus: u64) -> u64 {
    a.into_iter().zip(b).fold(0, |sum, (x, y)| {
        add_mod(sum, mul_mod(x, y, modulus), modulus)
    })
}

pub(super) fn column(matrix: &Exterior, j: usize) -> Fixed {
    std::array::from_fn(|i| matrix[i][j])
}

pub(super) fn product(a: &Exterior, b: &Exterior, modulus: u64) -> Exterior {
    a.map(|row| std::array::from_fn(|j| dot(row, column(b, j), modulus)))
}

// In the listed order (e1,e2,f1,f2), positive volume is
// e1 ^ f1 ^ e2 ^ f2. Pullback along the degree-two quotient multiplies
// cup products by two, so use half the abelian intersection form.
pub(super) fn intersection(v: Fixed, w: Fixed, modulus: u64) -> u64 {
    let mut result = 0;
    for (i, (x, y)) in v.into_iter().zip(w.into_iter().rev()).enumerate() {
        let term = mul_mod(x, y, modulus);
        // The complementary column is 5-i in the lexicographic list.
        // Only (02,13) and its transpose have positive wedge sign.
        result = if matches!(i, 1 | 4) {
            add_mod(result, term, modulus)
        } else {
            sub_mod(result, term, modulus)
        };
    }
    mul_mod(result, modulus / 2 + 1, modulus)
}

pub(super) fn exterior_change(
    model: &mut Model,
    target: &Basis,
    source: &Basis,
    ring: Residues,
) -> Result<Exterior> {
    let (matrix, _, exponent) = target.transition(model, source)?;
    let inverse = ring.inverse_unit(model, &exponent)?;
    // Transport the inverse Tate twist along with the exterior square.
    let result = scale_matrix(
        exterior_square(&matrix, ring.modulus),
        inverse,
        ring.modulus,
    );
    verify_isometry(model, ring.modulus, &result)?;
    Ok(result)
}

pub(super) fn verify_isometry(model: &mut Model, modulus: u64, matrix: &Exterior) -> Result<()> {
    for i in 0..6 {
        for j in i..6 {
            model.spend(72)?;
            let mut u = [0; 6];
            let mut v = [0; 6];
            u[i] = 1;
            v[j] = 1;
            if intersection(column(matrix, i), column(matrix, j), modulus)
                != intersection(u, v, modulus)
            {
                return Err(Error::plain(
                    "Kummer Brauer exterior matrix failed its intersection isometry check",
                ));
            }
        }
    }
    Ok(())
}
