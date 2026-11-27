use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;

use super::super::super::selmer::{local_roots, ResidueSquares, ROOT_DEGREE};
use super::Pairing;

mod series;
use series::{certified_lead, Laurent, Root};

const MAX_PRECISION: usize = 128;
type RootCharts = [[Scalar; ROOT_DEGREE + 1]; 3];

pub(super) struct Witness {
    x: Laurent,
    tangents: [(i64, Scalar); 3],
}

pub(super) struct Evidence<'a> {
    witness: &'a Witness,
    right: [(i64, Scalar); 3],
    symbols: [i32; 3],
}

impl Evidence<'_> {
    pub(super) fn sign(&self) -> i32 {
        self.symbols.into_iter().product()
    }

    pub(super) fn report(&self) -> String {
        format!(
            "x-chart s^{}*{:?}, tangent (valuation, leading unit)={:?}, right (valuation, leading unit)={:?}, symbols={:?}, product={}",
            self.witness.x.order, self.witness.x.coefficients, self.witness.tangents, self.right, self.symbols, self.sign(),
        )
    }
}

pub(super) fn certify(
    model: &mut Model,
    field: Field<'_>,
    squares: &ResidueSquares,
    eta: &Scalar,
    pairing: &Pairing,
    parameter: Option<&Scalar>,
) -> Result<Witness> {
    let cover = &pairing.cover;
    let weight = if parameter.is_some() {
        0
    } else {
        ROOT_DEGREE as i64
    };
    let split_roots = cover.pencil.split_roots(model)?;
    let roots = local_roots(model, field, split_roots, parameter)?;
    let [b0, b1, b2] = &cover.representatives;
    let b = [
        chart(model, field, b0, parameter)?,
        chart(model, field, b1, parameter)?,
        chart(model, field, b2, parameter)?,
    ];
    let required = b.each_ref().map(|b| b.symbol(squares));
    let mut x = if required == [0; 3] {
        // The negative x chart yields the three trivial local classes.
        Laurent::new(-2, vec![field.one()])
    } else {
        let prefix = prefix_x(model, field, squares, eta, &roots, required)?.ok_or_else(|| {
            Error::plain("Cassels local seed search did not realize the certified Selmer class")
        })?;
        Laurent::new(0, prefix)
    };
    let e = roots.map(|root| Laurent::new(0, root.into()));
    let mut delta = differences(model, field, &x, &e)?;
    let high_order = 1 + delta.iter().map(|v| v.order).max().unwrap();
    let mut forbidden = Vec::new();
    for tangent in &pairing.tangents {
        if let Some(tangent_x) = &tangent.tangent_x {
            let numerator = chart(model, field, &tangent_x.numerator, parameter)?.shift(weight);
            let denominator = chart(model, field, &tangent_x.denominator, parameter)?;
            forbidden.push((numerator, denominator));
        }
    }
    // Each conic tangent vanishes at one projective point, which has at
    // most one affine x under the covering map. The four exact x seeds
    // have the same local square classes and cannot all be forbidden.
    let base = x.clone();
    for step in 0..=forbidden.len() {
        let mut admissible = true;
        for (numerator, denominator) in &forbidden {
            if x.product(model, field, denominator)? == *numerator {
                admissible = false;
                break;
            }
        }
        if admissible {
            break;
        }
        if step == forbidden.len() {
            return Err(Error::plain(
                "Cassels local seed selection failed to avoid tangencies",
            ));
        }
        let perturb = Laurent::new(high_order + step as i64, vec![field.one()]);
        x = base.sum(model, field, &perturb, false)?;
    }
    if x != base {
        delta = differences(model, field, &x, &e)?;
    }
    if delta.each_ref().map(|d| d.symbol(squares)) != required {
        return Err(Error::plain(
            "Cassels local seed has incorrect valuation or angular square class",
        ));
    }
    let mut roots = [
        Root::new(model, field, &delta[0], &b[0], weight)?,
        Root::new(model, field, &delta[1], &b[1], weight)?,
        Root::new(model, field, &delta[2], &b[2], weight)?,
    ];
    let tangent_data = pairing
        .tangents
        .iter()
        .map(|tangent| {
            let [a, b, c] = &tangent.coefficients;
            Ok([
                chart(model, field, a, parameter)?,
                chart(model, field, b, parameter)?,
                chart(model, field, c, parameter)?,
            ])
        })
        .collect::<Result<Vec<_>>>()?;
    let mut leading: [Option<(i64, Scalar)>; 3] = [None, None, None];
    let mut precision = 1;
    loop {
        for root in &mut roots {
            root.extend(model, field, precision)?;
        }
        for (i, value) in leading.iter_mut().enumerate() {
            if value.is_none() {
                let [j, k] = pairing.tangents[i].indices;
                *value = certified_lead(
                    model,
                    field,
                    &tangent_data[i],
                    [&roots[j], &roots[k]],
                    precision,
                )?;
            }
        }
        if leading.iter().all(Option::is_some) {
            break;
        }
        precision *= 2;
        if precision > MAX_PRECISION {
            return Err(Error::plain(
                "Cassels local tangent certification exceeded its series precision limit",
            ));
        }
    }
    let tangents = leading.map(|v| v.expect("leading terms are certified"));
    Ok(Witness { x, tangents })
}

impl Witness {
    pub(super) fn evaluate(
        &self,
        model: &mut Model,
        field: Field<'_>,
        squares: &ResidueSquares,
        right: &[Polynomial; 3],
        parameter: Option<&Scalar>,
    ) -> Result<Evidence<'_>> {
        let [a, b, c] = right;
        let right = [
            chart(model, field, a, parameter)?,
            chart(model, field, b, parameter)?,
            chart(model, field, c, parameter)?,
        ];
        let minus_one = squares.symbol(&field.signed(-1), 0) >> 1;
        let symbols = std::array::from_fn(|i| {
            let (order, ref angular) = self.tangents[i];
            let left = squares.symbol(angular, order.rem_euclid(2) as usize);
            let right = right[i].symbol(squares);
            hilbert(left, right, minus_one)
        });
        let right = right.map(|r| (r.order, r.coefficients[0].clone()));
        Ok(Evidence {
            witness: self,
            right,
            symbols,
        })
    }
}

fn differences(
    model: &mut Model,
    field: Field<'_>,
    x: &Laurent,
    roots: &[Laurent; 3],
) -> Result<[Laurent; 3]> {
    let [a, b, c] = roots;
    let differences = [
        x.sum(model, field, a, true)?,
        x.sum(model, field, b, true)?,
        x.sum(model, field, c, true)?,
    ];
    if differences.iter().any(Laurent::is_zero) {
        return Err(Error::plain("Cassels local x seed meets a branch root"));
    }
    Ok(differences)
}

fn chart(
    model: &mut Model,
    field: Field<'_>,
    polynomial: &[Scalar],
    parameter: Option<&Scalar>,
) -> Result<Laurent> {
    let Some(parameter) = parameter else {
        return Ok(Laurent::new(
            -(polynomial.len().saturating_sub(1) as i64),
            polynomial.iter().rev().cloned().collect(),
        ));
    };
    // Substitute t=a+s exactly, since multiple vanishing coefficients
    // can determine the initial term of a tangent or squareclass.
    let linear = [parameter.clone(), field.one()];
    let mut result = Vec::new();
    for coefficient in polynomial.iter().rev() {
        result = poly::multiply(model, field, &result, &linear)?;
        result = poly::sum(
            model,
            field,
            &result,
            std::slice::from_ref(coefficient),
            false,
        )?;
    }
    Ok(Laurent::new(0, result))
}

fn hilbert(left: u8, right: u8, minus_one: u8) -> i32 {
    let (la, ra) = (left & 1, right & 1);
    let power = (la * ra * minus_one) ^ ((left >> 1) * ra) ^ ((right >> 1) * la);
    if power == 0 {
        1
    } else {
        -1
    }
}

fn prefix_x(
    model: &mut Model,
    field: Field<'_>,
    squares: &ResidueSquares,
    eta: &Scalar,
    charts: &RootCharts,
    required: [u8; 3],
) -> Result<Option<Polynomial>> {
    let selector = SeedSelector {
        field,
        squares,
        eta,
        charts,
        required,
    };
    selector.extend(model, 0b111, &mut Vec::new())
}

struct SeedSelector<'a> {
    field: Field<'a>,
    squares: &'a ResidueSquares,
    eta: &'a Scalar,
    charts: &'a RootCharts,
    required: [u8; 3],
}

impl SeedSelector<'_> {
    fn extend(
        &self,
        model: &mut Model,
        active: u8,
        prefix: &mut Polynomial,
    ) -> Result<Option<Polynomial>> {
        let order = prefix.len();
        // Two distinct roots of degree at most four cannot share their
        // full length-five prefix. Singleton branches are solved directly.
        if order > ROOT_DEGREE {
            return Err(Error::plain(
                "Cassels root-prefix search exceeded the root degree",
            ));
        }
        'residues: for i in 0..self.squares.count() {
            model.spend(1 + self.field.ring.degree())?;
            let residue = self.field.ring.coordinate_row(i);
            let mut matches = 0u8;
            for (j, root) in self.charts.iter().enumerate() {
                if active & (1 << j) == 0 {
                    continue;
                }
                if residue == root[order] {
                    matches |= 1 << j;
                } else {
                    let angular = self.field.subtract(model, &residue, &root[order])?;
                    if self.squares.symbol(&angular, order) != self.required[j] {
                        continue 'residues;
                    }
                }
            }
            prefix.push(residue);
            let result = match matches.count_ones() {
                0 => Some(prefix.clone()),
                1 => {
                    let j = matches.trailing_zeros() as usize;
                    let mut jump = order + 1;
                    if (jump & 1) as u8 != self.required[j] & 1 {
                        jump += 1;
                    }
                    let mut x = self.charts[j].to_vec();
                    if jump >= x.len() {
                        x.resize(jump + 1, self.field.zero());
                    }
                    let unit = if self.required[j] >> 1 == 0 {
                        self.field.one()
                    } else {
                        self.eta.clone()
                    };
                    x[jump] = self.field.add(model, &x[jump], &unit)?;
                    Some(x)
                }
                _ => self.extend(model, matches, prefix)?,
            };
            prefix.pop();
            if result.is_some() {
                return Ok(result);
            }
        }
        Ok(None)
    }
}
