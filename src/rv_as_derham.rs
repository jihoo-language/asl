use std::collections::{BTreeMap, VecDeque};

use crate::binary_field::BinaryField;
use crate::binary_linear::{self, check_storage, Entry};
use crate::error::{Error, Result};
use crate::linear::{self, Basis};
use crate::model::{Model, Value};
use crate::rv_as_cover::Cover;
use crate::rv_poly::Polynomial;

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Include,
    Scale(usize),
    Frobenius,
    Verschiebung,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Include => "rvASDeRhamIncl",
            Self::Scale(_) => "rvASDeRhamScale",
            Self::Frobenius => "rvASF",
            Self::Verschiebung => "rvASV",
        }
    }
}

pub(crate) fn generator(model: &mut Model, genus: usize, index: usize) -> Result<Value> {
    let dimension = 2 * genus;
    if index >= dimension {
        return Err(Error::plain(format!(
            "rvASDeRhamGen index {index} is out of range for dimension {dimension}"
        )));
    }
    persist(model, vec![(index, 1)])
}

fn persist(model: &mut Model, vector: Vec<Entry>) -> Result<Value> {
    Ok(Value::ASDeRham(binary_linear::persist(model, vector)?))
}

pub(crate) fn apply(
    model: &mut Model,
    cover: &Cover,
    field: BinaryField,
    operation: Operation,
    input: Value,
) -> Result<Value> {
    if let Operation::Include = operation {
        let Value::CurveDiff(coordinates) = input else {
            unreachable!()
        };
        return Ok(Value::ASDeRham(coordinates));
    }
    let Value::ASDeRham(coordinates) = input else {
        unreachable!()
    };
    if let Operation::Scale(scalar) = operation {
        let result = binary_linear::scaled(model, field, &coordinates, scalar)?;
        return persist(model, result);
    }
    let chart = Chart::new(cover)?;
    let frobenius = matches!(operation, Operation::Frobenius);
    let twist = if frobenius { 1 } else { field.degree() - 1 };
    let mut result = Vec::new();
    for &(index, coefficient) in coordinates.iter() {
        let column = if frobenius {
            chart.frobenius(model, index)?
        } else {
            chart.verschiebung(model, index)?
        };
        if !column.is_empty() {
            let scale = field.frobenius(model, coefficient, twist)?;
            let vector: Vec<_> = column.into_iter().map(|i| (i, 1)).collect();
            result = binary_linear::add_multiple(model, field, &result, &vector, scale)?;
        }
    }
    persist(model, result)
}

pub(crate) fn pair(
    model: &mut Model,
    cover: &Cover,
    field: BinaryField,
    left: &[Entry],
    right: &[Entry],
) -> Result<Value> {
    let chart = Chart::new(cover)?;
    let mut value = 0;
    for &(i, x) in left {
        model.spend(1 + right.len())?;
        for &(j, y) in right {
            if chart.pairs(model, i, j)? {
                model.spend(field.degree())?;
                value ^= field.multiply(x, y);
            }
        }
    }
    Ok(Value::BinaryScalar(value))
}

struct Chart<'a> {
    polynomial: &'a Polynomial,
    genus: usize,
}

impl<'a> Chart<'a> {
    fn new(cover: &'a Cover) -> Result<Self> {
        Ok(Self {
            polynomial: cover.polynomial_chart()?,
            genus: cover.genus(),
        })
    }

    fn after(&self, model: &mut Model, degree: usize) -> Result<impl Iterator<Item = usize> + '_> {
        model.spend(1 + self.polynomial.len().max(1).ilog2() as usize)?;
        let offset = self.polynomial.partition_point(|term| term[0] <= degree);
        let tail = &self.polynomial[offset..];
        model.spend(tail.len())?;
        Ok(tail.iter().map(|term| term[0]).filter(|d| d % 2 != 0))
    }

    fn frobenius(&self, model: &mut Model, index: usize) -> Result<Vec<usize>> {
        model.spend(1)?;
        let g = self.genus;
        if index < g {
            return Ok(Vec::new());
        }
        let twice_j = 2 * (index - g + 1);
        // Squaring y/x^j gives y/x^(2j) + f/x^(2j). If the former
        // crosses the gap basis, the exact parts leave a holomorphic form.
        if twice_j <= g {
            return Ok(vec![g + twice_j - 1]);
        }
        Ok(self
            .after(model, twice_j)?
            .map(|d| d - 1 - twice_j)
            .collect())
    }

    fn verschiebung(&self, model: &mut Model, index: usize) -> Result<Vec<usize>> {
        model.spend(1)?;
        let g = self.genus;
        if index < g {
            return Ok(if !index.is_multiple_of(2) {
                vec![index / 2]
            } else {
                Vec::new()
            });
        }
        let j = index - g + 1;
        // The U-form in the lift of y/x^j is the polynomial part of
        // f' x^(-j) dx. Cartier selects the odd powers of x there.
        if j.is_multiple_of(2) {
            return Ok(Vec::new());
        }
        Ok(self.after(model, j + 1)?.map(|d| (d - 2 - j) / 2).collect())
    }

    fn pairs(&self, model: &mut Model, a: usize, b: usize) -> Result<bool> {
        let g = self.genus;
        match (a < g, b < g) {
            (true, true) => Ok(false),
            (true, false) => Ok(a + g == b),
            (false, true) => Ok(b + g == a),
            (false, false) => {
                // The Čech cup product takes the residue of
                // h_a*omega_V(b) + omega_U(a)*h_b at the omitted fibre.
                let index = (a - g + 1) + (b - g + 1);
                if index.is_multiple_of(2) {
                    return Ok(false);
                }
                model.spend(1 + self.polynomial.len().max(1).ilog2() as usize)?;
                Ok(self
                    .polynomial
                    .binary_search_by_key(&index, |term| term[0])
                    .is_ok())
            }
        }
    }
}

struct Matrices {
    frobenius: Vec<Vec<usize>>,
    verschiebung: Vec<Vec<usize>>,
    cells: usize,
}

impl Matrices {
    fn new(dimension: usize) -> Result<Self> {
        let cells = 6usize.saturating_mul(dimension);
        check_storage(cells)?;
        Ok(Self {
            frobenius: Vec::with_capacity(dimension),
            verschiebung: Vec::with_capacity(dimension),
            cells,
        })
    }

    fn push(&mut self, model: &mut Model, f: Vec<usize>, v: Vec<usize>) -> Result<()> {
        model.spend(1 + f.len() + v.len())?;
        self.cells += f.len() + v.len();
        check_storage(self.cells)?;
        self.frobenius.push(f);
        self.verschiebung.push(v);
        Ok(())
    }

    fn append_chart(&mut self, model: &mut Model, chart: &Chart<'_>) -> Result<()> {
        let offset = self.frobenius.len();
        for i in 0..2 * chart.genus {
            let mut f = chart.frobenius(model, i)?;
            let mut v = chart.verschiebung(model, i)?;
            for index in f.iter_mut().chain(&mut v) {
                *index += offset;
            }
            self.push(model, f, v)?;
        }
        Ok(())
    }

    fn geometric(model: &mut Model, cover: &Cover, dimension: usize) -> Result<Self> {
        let mut matrices = Self::new(dimension)?;
        let profile = cover.geometric_pole_profile(model)?;
        // Elkin–Pries Theorems 1.2–1.3 over the algebraic closure. Each
        // G_c can use the one-pole model with monomial x^(2c+1).
        for _ in 0..cover.p_rank() {
            let offset = matrices.frobenius.len();
            matrices.push(model, vec![offset], Vec::new())?;
            matrices.push(model, Vec::new(), vec![offset + 1])?;
        }
        for (order, count) in profile {
            let genus = (order - 1) / 2;
            if genus == 0 {
                continue;
            }
            let polynomial: Polynomial = [[order].into()].into();
            let chart = Chart {
                polynomial: &polynomial,
                genus,
            };
            for _ in 0..count {
                matrices.append_chart(model, &chart)?;
            }
        }
        if matrices.frobenius.len() != dimension {
            return Err(dimension_error());
        }
        Ok(matrices)
    }

    // The F2 spans extend to the canonical flag of the semilinear
    // operators over any binary extension.
    fn image(&self, model: &mut Model, input: &Basis, retained: usize) -> Result<Basis> {
        let mut image = Basis::default();
        for vector in input.vectors() {
            let mut column = Vec::new();
            for &i in vector {
                column = linear::xor(model, &column, &self.verschiebung[i])?;
                check_storage(retained + image.cells() + column.len())?;
            }
            image.insert(model, column)?;
            check_storage(retained + image.cells())?;
        }
        Ok(image)
    }

    fn preimage(&self, model: &mut Model, target: &Basis, retained: usize) -> Result<Basis> {
        let mut preimage = Basis::default();
        // Dependencies on the F-columns modulo target record a basis of
        // F^{-1}(target); only source columns receive witness coordinates.
        let mut quotient = Basis::default();
        for vector in target.vectors() {
            quotient.insert(model, vector.to_vec())?;
            check_storage(retained + quotient.cells())?;
        }
        for (i, column) in self.frobenius.iter().enumerate() {
            if let Some(vector) = quotient.insert_with_witness(model, column.clone(), vec![i])? {
                preimage.insert(model, vector)?;
            }
            check_storage(retained + quotient.cells() + preimage.cells())?;
        }
        Ok(preimage)
    }
}

struct FlagStep {
    space: Basis,
    image_rank: usize,
}

struct Flag {
    steps: BTreeMap<usize, FlagStep>,
    pending: VecDeque<usize>,
    cells: usize,
}

impl Flag {
    fn include(&mut self, model: &mut Model, space: Basis) -> Result<()> {
        let rank = space.rank();
        if let Some(existing) = self.steps.get(&rank) {
            return contained(model, &space, &existing.space);
        }
        let cells = self.cells + space.cells() + 6;
        check_storage(cells)?;
        self.cells = cells;
        self.steps.insert(
            rank,
            FlagStep {
                space,
                image_rank: 0,
            },
        );
        self.pending.push_back(rank);
        Ok(())
    }
}

fn contained(model: &mut Model, lower: &Basis, upper: &Basis) -> Result<()> {
    for vector in lower.vectors() {
        if !upper.reduce(model, vector.to_vec())?.is_empty() {
            return Err(Error::plain(
                "Artin–Schreier canonical EO subspaces are not nested",
            ));
        }
    }
    Ok(())
}

fn dimension_error() -> Error {
    Error::plain("Artin–Schreier geometric Dieudonné dimension is out of range")
}

pub(crate) fn final_type(model: &mut Model, cover: &Cover) -> Result<Vec<usize>> {
    let genus = cover.genus();
    let dimension = genus.checked_mul(2).ok_or_else(dimension_error)?;
    let matrices = if cover.has_polynomial_chart() {
        let chart = Chart::new(cover)?;
        let mut matrices = Matrices::new(dimension)?;
        matrices.append_chart(model, &chart)?;
        matrices
    } else {
        Matrices::geometric(model, cover, dimension)?
    };
    let mut flag = Flag {
        steps: BTreeMap::new(),
        pending: VecDeque::new(),
        cells: matrices.cells,
    };
    flag.include(model, Basis::default())?;
    let mut ambient = Basis::default();
    for i in 0..dimension {
        ambient.insert(model, vec![i])?;
        check_storage(flag.cells + ambient.cells())?;
    }
    flag.include(model, ambient)?;
    while let Some(rank) = flag.pending.pop_front() {
        let step = flag.steps.get_mut(&rank).unwrap();
        let image = matrices.image(model, &step.space, flag.cells)?;
        let preimage = matrices.preimage(model, &step.space, flag.cells + image.cells())?;
        step.image_rank = image.rank();
        flag.include(model, image)?;
        flag.include(model, preimage)?;
    }
    let mut previous_rank = 0;
    let mut previous = &flag.steps[&0];
    let mut sequence = Vec::with_capacity(genus);
    for (&rank, step) in flag.steps.iter().skip(1) {
        contained(model, &previous.space, &step.space)?;
        let width = rank - previous_rank;
        let rise = step.image_rank.checked_sub(previous.image_rank);
        if !matches!(rise, Some(0)) && rise != Some(width) {
            return Err(Error::plain(
                "Artin–Schreier canonical EO flag has invalid Verschiebung ranks",
            ));
        }
        for dimension in previous_rank + 1..=rank.min(genus) {
            let increment = if rise == Some(width) {
                dimension - previous_rank
            } else {
                0
            };
            sequence.push(previous.image_rank + increment);
        }
        previous_rank = rank;
        previous = step;
    }
    Ok(sequence)
}
