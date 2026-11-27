use crate::dwork_k3::fricke::hecke;
use crate::ec_poly::Field;
use crate::error::Result;
use crate::finite_field_factor;
use crate::hw_algebra::CoefficientRing;
use crate::lee_linear::Rational;
use crate::model::Model;
use crate::rv_lens::gcd;
use crate::rv_qp_quad::Wide;
use crate::rv_unram::Ring;

use super::{arithmetic, certificate, Basis};

#[derive(Clone, Eq, PartialEq)]
pub(super) struct Matrix(Vec<Vec<i128>>);

impl Matrix {
    pub(super) fn hecke(
        model: &mut Model,
        basis: &Basis,
        ring: &Ring,
        level: usize,
    ) -> Result<Self> {
        let degree = check_level(basis.p, level)?;
        let equation = model.dwork_hecke_equation(level)?;
        let field = Field { ring };
        let factors = basis
            .vertices
            .iter()
            .map(|h| Ok([field.subtract(model, &field.zero(), h)?, field.one()]))
            .collect::<Result<Vec<_>>>()?;
        let mut rows = Vec::with_capacity(factors.len());
        for x in &basis.vertices {
            let mut specialization = (0..=degree)
                .map(|j| equation.specialization(model, ring, x, j))
                .collect::<Result<Vec<_>>>()?;
            if specialization.last() != Some(&field.one()) {
                return Err(certificate("monic modular specialization"));
            }
            let row = factors
                .iter()
                .map(|f| {
                    Ok(
                        finite_field_factor::remove_factor(model, field, &mut specialization, f)?
                            as i128,
                    )
                })
                .collect::<Result<Vec<_>>>()?;
            if specialization != [field.one()] || row.iter().sum::<i128>() != degree as i128 {
                return Err(certificate(
                    "complete Hecke targets with multiplicity and correct row sums",
                ));
            }
            rows.push(row);
        }
        let matrix = Self(rows);
        for (i, row) in matrix.0.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                if value != matrix.0[basis.frobenius[i]][basis.frobenius[j]] {
                    return Err(certificate("Hecke commutes with Frobenius"));
                }
            }
        }
        // Also check every individual operator, even before the relations query.
        balance(model, &[&matrix], &basis.frobenius)?;
        let n = matrix.0.len();
        model.allocate(4 + n * (3 + 2 * n))?;
        Ok(matrix)
    }

    fn size(&self) -> usize {
        self.0.len()
    }

    fn multiply(&self, model: &mut Model, other: &Self) -> Result<Self> {
        let n = self.size();
        model.spend(n * n * n)?;
        let mut rows = vec![vec![0i128; n]; n];
        for (row, a) in rows.iter_mut().zip(&self.0) {
            for (&x, b) in a.iter().zip(&other.0) {
                for (entry, &y) in row.iter_mut().zip(b) {
                    *entry = add(*entry, multiply(x, y)?)?;
                }
            }
        }
        Ok(Self(rows))
    }

    pub(super) fn push(&self, model: &mut Model, value: &[i128]) -> Result<Vec<i128>> {
        let n = self.size();
        model.spend(5 * n * n)?;
        // At most 3h_ss source vertices, each with row degree at most 12.
        // A signed 160-bit sum therefore covers every i128 input term.
        let mut output = vec![Wide::<5>::from_i128(0); n];
        for (row, &coefficient) in self.0.iter().zip(value) {
            let coefficient = Wide::from_i128(coefficient);
            for (entry, &multiplicity) in output.iter_mut().zip(row) {
                if multiplicity != 0 {
                    let term = coefficient
                        .multiply(multiplicity)
                        .map_err(|_| arithmetic())?;
                    *entry = entry.add(term).map_err(|_| arithmetic())?;
                }
            }
        }
        output
            .into_iter()
            .map(|c| c.signed().ok_or_else(arithmetic))
            .collect()
    }

    pub(super) fn spectra(
        &self,
        model: &mut Model,
        frobenius: &[usize],
        level: usize,
    ) -> Result<String> {
        let (weights, components) = balance(model, &[self], frobenius)?;
        let n = self.size();
        let pairs = frobenius
            .iter()
            .enumerate()
            .filter(|(i, j)| *i < **j)
            .count();
        let plus_rank = n - pairs;
        let mut traces = vec![0i128];
        let mut plus_traces = vec![0i128];
        let mut minus_traces = vec![0i128];
        let mut power = self.clone();
        for k in 1..=n {
            let (mut trace, mut twisted) = (0i128, 0i128);
            for (i, &j) in frobenius.iter().enumerate() {
                trace = add(trace, power.0[i][i])?;
                twisted = add(twisted, power.0[i][j])?;
            }
            let plus = add(trace, twisted)?;
            let minus = trace.checked_sub(twisted).ok_or_else(arithmetic)?;
            if plus % 2 != 0 || minus % 2 != 0 {
                return Err(certificate("integral traces on Frobenius sectors"));
            }
            traces.push(trace);
            plus_traces.push(plus / 2);
            minus_traces.push(minus / 2);
            if k < n {
                power = power.multiply(model, self)?;
            }
        }
        let full = characteristic(&traces)?;
        let plus = characteristic(&plus_traces[..=plus_rank])?;
        let minus = characteristic(&minus_traces[..=pairs])?;
        let mut recombined = vec![0; n + 1];
        for (i, &a) in plus.iter().enumerate() {
            for (j, &b) in minus.iter().enumerate() {
                recombined[i + j] = add(recombined[i + j], multiply(a, b)?)?;
            }
        }
        if full != recombined {
            return Err(certificate(
                "product of Frobenius-sector characteristic polynomials",
            ));
        }
        let row_sum = hecke::degree(level)? as i128;
        let mut cuspidal = full.clone();
        for _ in &components {
            let mut quotient = Vec::with_capacity(cuspidal.len() - 1);
            let mut previous = 0;
            for &coefficient in &cuspidal[..cuspidal.len() - 1] {
                previous = add(coefficient, multiply(row_sum, previous)?)?;
                quotient.push(previous);
            }
            if add(*cuspidal.last().unwrap(), multiply(row_sum, previous)?)? != 0 {
                return Err(certificate(
                    "one constant eigenvalue per connected component",
                ));
            }
            cuspidal = quotient;
        }
        Ok(format!(
            "B{level} = {:?}; weights = {weights:?}; components = {components:?}\n\
             characteristic (descending) = {full:?}; cuspidal = {cuspidal:?}\n\
             Frobenius +1 (rank {plus_rank}) = {plus:?}; -1 (rank {pairs}) = {minus:?}",
            self.0,
        ))
    }
}

pub(super) fn check_level(p: u64, level: usize) -> Result<usize> {
    if level as u64 == p {
        return Err(crate::error::Error::plain(
            "supersingular Fricke Brandt action requires level prime to the characteristic",
        ));
    }
    hecke::degree(level)
}

pub(super) fn relations(
    model: &mut Model,
    frobenius: &[usize],
    three: &Matrix,
    five: Option<&Matrix>,
    nine: &Matrix,
) -> Result<String> {
    let mut square = three.multiply(model, three)?;
    for (i, row) in square.0.iter_mut().enumerate() {
        row[i] -= 4;
    }
    if square != *nine {
        return Err(certificate("B3² = B9 + 4I"));
    }
    let mut operators = vec![three, nine];
    if let Some(five) = five {
        if three.multiply(model, five)? != five.multiply(model, three)? {
            return Err(certificate("B3 B5 = B5 B3"));
        }
        operators.push(five);
    }
    let (weights, components) = balance(model, &operators, frobenius)?;
    let commutation = if five.is_some() {
        "B3 B5 = B5 B3"
    } else {
        "good levels 3 and 9 in characteristic 5"
    };
    Ok(format!(
        "B3^2 = B9 + 4I; {commutation}; Frobenius commutes with all available operators\n\
         common positive detailed-balance weights = {weights:?}; components = {components:?}"
    ))
}

fn balance(
    model: &mut Model,
    matrices: &[&Matrix],
    frobenius: &[usize],
) -> Result<(Vec<i128>, Vec<Vec<usize>>)> {
    let n = frobenius.len();
    model.spend(n * n * matrices.len())?;
    let mut ratios = vec![None; n];
    let mut weights = vec![0i128; n];
    let mut components = Vec::new();
    for start in 0..n {
        if ratios[start].is_some() {
            continue;
        }
        ratios[start] = Some(Rational::integer(1));
        let mut component = vec![start];
        let mut next = 0;
        while let Some(&i) = component.get(next) {
            let weight = ratios[i].unwrap();
            for operator in matrices {
                for (j, &outward) in operator.0[i].iter().enumerate() {
                    let inward = operator.0[j][i];
                    if (outward == 0) != (inward == 0) {
                        return Err(certificate("symmetric support of Hecke neighbors"));
                    }
                    if outward == 0 {
                        continue;
                    }
                    let next_weight = weight.multiply(Rational::new(outward, inward)?)?;
                    if let Some(known) = ratios[j] {
                        if known != next_weight {
                            return Err(certificate("simultaneous positive detailed balance"));
                        }
                    } else {
                        ratios[j] = Some(next_weight);
                        component.push(j);
                    }
                }
            }
            next += 1;
        }
        component.sort_unstable();
        let mut denominator = 1i128;
        for &i in &component {
            let d = ratios[i].unwrap().denominator;
            denominator = multiply(denominator / gcd(denominator as u128, d as u128) as i128, d)?;
        }
        let mut divisor = 0u128;
        for &i in &component {
            let weight = ratios[i].unwrap();
            weights[i] = multiply(weight.numerator, denominator / weight.denominator)?;
            divisor = gcd(divisor, weights[i] as u128);
        }
        for &i in &component {
            weights[i] /= divisor as i128;
        }
        components.push(component);
    }
    for (i, &j) in frobenius.iter().enumerate() {
        if weights[i] <= 0 || weights[i] != weights[j] {
            return Err(certificate("Frobenius preserves the positive weights"));
        }
    }
    Ok((weights, components))
}

fn characteristic(traces: &[i128]) -> Result<Vec<i128>> {
    let mut coefficients = vec![1];
    for k in 1..traces.len() {
        let mut sum = 0i128;
        for j in 1..=k {
            sum = add(sum, multiply(coefficients[k - j], traces[j])?)?;
        }
        let k = k as i128;
        if sum % k != 0 {
            return Err(certificate("integral Newton characteristic coefficients"));
        }
        coefficients.push((sum / k).checked_neg().ok_or_else(arithmetic)?);
    }
    Ok(coefficients)
}

fn add(a: i128, b: i128) -> Result<i128> {
    a.checked_add(b).ok_or_else(arithmetic)
}

fn multiply(a: i128, b: i128) -> Result<i128> {
    a.checked_mul(b).ok_or_else(arithmetic)
}
