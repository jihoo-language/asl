use std::collections::{btree_map::Entry, BTreeMap};
use std::rc::Rc;

use crate::comodule::Generator;
use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};
use crate::q_koszul::{self, Complex, Term};
use crate::v_koszul::{self, Grade};

type Column = BTreeMap<usize, usize>;
type Coefficient = (Generator, Generator);
const ENTRY_CELLS: usize = 2;

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Family {
    pub(crate) source: TypeId,
    pub(crate) target: TypeId,
    pub(crate) indices: [usize; 2],
}

impl Family {
    pub(crate) fn new(source: TypeId, target: TypeId, selected: &[usize]) -> Result<Self> {
        let indices: [usize; 2] = selected
            .try_into()
            .map_err(|_| Error::plain("localized bars require exactly two Milnor primitives"))?;
        if indices[0] >= indices[1] {
            return Err(Error::plain(
                "Milnor primitive indices must be strictly increasing",
            ));
        }
        if indices[1] >= usize::BITS as usize {
            return Err(Error::plain("Milnor primitive degree is out of range"));
        }
        Ok(Self {
            source,
            target,
            indices,
        })
    }

    fn top(self) -> usize {
        q_koszul::weight(self.indices[1])
    }

    pub(crate) fn gap(self) -> usize {
        self.top() - q_koszul::weight(self.indices[0])
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Certificate {
    Representative,
    Nullhomotopy,
}

impl Certificate {
    pub(crate) fn ty(
        self,
        model: &mut Model,
        family: Family,
        degree: i64,
        exponent: usize,
    ) -> Result<TypeId> {
        if exponent == 0 {
            return Err(Error::plain("a torsion bar exponent must be positive"));
        }
        let (s, t) = match self {
            Self::Representative => (0, degree),
            Self::Nullhomotopy => (
                -1,
                ratio_degree(degree, exponent, family.gap(), "vbarNull")?,
            ),
        };
        model.v_koszul_type(Grade {
            source: family.source,
            target: family.target,
            selected: Rc::from(family.indices),
            s,
            t,
        })
    }
}

struct Bar {
    degree: i128,
    exponent: usize,
    representative: Column,
    nullhomotopy: Column,
}

pub(crate) struct Decomposition {
    coefficients: Vec<Coefficient>,
    free: Vec<i128>,
    torsion: Vec<Bar>,
}

impl Decomposition {
    pub(crate) fn build(model: &mut Model, family: Family) -> Result<Self> {
        let mut complex = Complex::new(model, family.source, family.target, &family.indices)?;
        let basis = complex.coefficient_basis(model)?;
        let coefficients: Vec<_> = basis
            .iter()
            .map(|term| (term.source, term.target))
            .collect();
        let grades: Vec<_> = coefficients
            .iter()
            .map(|(source, target)| target.degree as i128 - source.degree as i128)
            .collect();
        let mut matrix = Vec::with_capacity(basis.len());
        let column_cells = std::mem::size_of::<Column>().div_ceil(std::mem::size_of::<usize>());
        let mut cells = basis.len() * column_cells;
        q_koszul::check_storage(cells)?;
        for generator in &basis {
            let mut column = Column::new();
            for term in complex.differential(model, std::slice::from_ref(generator))? {
                model.spend(1 + basis.len().max(1).ilog2() as usize)?;
                let row = coefficients
                    .binary_search(&(term.source, term.target))
                    .map_err(|_| Error::plain("Koszul differential left its coefficient basis"))?;
                // Dividing d by the inverted top variable gives Q_n + w Q_i.
                debug_assert_eq!(term.powers.iter().sum::<i128>(), 1);
                let exponent = term.powers[0] as usize;
                toggle(model, &mut column, row, exponent, &mut cells)?;
            }
            matrix.push(column);
        }
        drop(basis);
        drop(complex);

        let (free, mut torsion) = smith(model, matrix, &grades, cells, column_cells)?;
        model.spend(sort_work(torsion.len()))?;
        torsion.sort_by_key(|bar| (bar.degree, bar.exponent));
        let coefficients = if torsion.is_empty() {
            Vec::new()
        } else {
            coefficients
        };
        let certificate_cells: usize = torsion
            .iter()
            .map(|bar| 6 + ENTRY_CELLS * (bar.representative.len() + bar.nullhomotopy.len()))
            .sum();
        model.allocate(1 + 4 * coefficients.len() + 2 * free.len() + certificate_cells)?;
        Ok(Self {
            coefficients,
            free,
            torsion,
        })
    }

    pub(crate) fn report(&self, model: &mut Model) -> Result<String> {
        model.spend(1 + self.free.len() + self.torsion.len())?;
        let torsion: Vec<_> = self
            .torsion
            .iter()
            .map(|bar| (bar.degree, bar.exponent))
            .collect();
        Ok(format!("free {:?}, torsion {torsion:?}", self.free))
    }

    pub(crate) fn rank(&self, model: &mut Model, family: Family, s: i64, t: i64) -> Result<usize> {
        model.spend(1 + self.free.len() + self.torsion.len())?;
        let normalized = i128::from(t) - i128::from(s) * family.top() as i128;
        let gap = family.gap() as u128;
        let multiple = |degree: i128| {
            if degree < normalized {
                return None;
            }
            let difference = degree.abs_diff(normalized);
            difference.is_multiple_of(gap).then_some(difference / gap)
        };
        let free = self
            .free
            .iter()
            .filter(|&&degree| multiple(degree).is_some())
            .count();
        let torsion = self
            .torsion
            .iter()
            .filter(|bar| multiple(bar.degree).is_some_and(|k| k < bar.exponent as u128))
            .count();
        Ok(free + torsion)
    }

    pub(crate) fn certificate(
        &self,
        model: &mut Model,
        kind: Certificate,
        degree: i64,
        exponent: usize,
        index: usize,
    ) -> Result<Value> {
        model.spend(1 + self.torsion.len())?;
        let bar = self.torsion.iter()
            .filter(|bar| bar.degree == i128::from(degree) && bar.exponent == exponent)
            .nth(index)
            .ok_or_else(|| Error::plain(format!(
                "torsion bar index {index} is out of range at degree {degree} and exponent {exponent}"
            )))?;
        let (column, s) = match kind {
            Certificate::Representative => (&bar.representative, 0i128),
            Certificate::Nullhomotopy => (&bar.nullhomotopy, -1i128),
        };
        let terms = column
            .iter()
            .map(|(&coordinate, &power)| {
                let (source, target) = self.coefficients[coordinate];
                let lower = power as i128;
                Term {
                    source,
                    target,
                    powers: [lower, s - lower].into(),
                }
            })
            .collect();
        Ok(Value::QKoszul(q_koszul::persist(model, terms)?))
    }
}

fn smith(
    model: &mut Model,
    mut matrix: Vec<Column>,
    grades: &[i128],
    mut cells: usize,
    column_cells: usize,
) -> Result<(Vec<i128>, Vec<Bar>)> {
    if matrix.iter().all(Column::is_empty) {
        model.spend(sort_work(grades.len()))?;
        let mut free = grades.to_vec();
        free.sort_unstable();
        return Ok((free, Vec::new()));
    }
    cells += 2 * grades.len() * (column_cells + ENTRY_CELLS);
    q_koszul::check_storage(cells)?;
    let unit = |index| Column::from([(index, 0)]);
    let mut target: Vec<_> = (0..grades.len()).map(unit).collect();
    let mut source = target.clone();
    let mut free = BTreeMap::<i128, i32>::new();
    let mut bars = Vec::new();
    for &degree in grades {
        *free.entry(degree).or_default() += 1;
    }

    while let Some((exponent, row, col)) = pivot(model, &matrix)? {
        let mut column = std::mem::take(&mut matrix[col]);
        column.remove(&row);
        cells -= ENTRY_CELLS;
        let mut target_lift = std::mem::take(&mut target[row]);
        let source_lift = std::mem::take(&mut source[col]);

        // Entries in the pivot row are needed before the row is cleared.
        let across = remove_row(model, &mut matrix, row)?;
        let across_cells = ENTRY_CELLS * across.len();
        let column_cells = ENTRY_CELLS * column.len();
        for (other_row, power) in column {
            let shift = power - exponent;
            for &(other_col, row_power) in &across {
                toggle(
                    model,
                    &mut matrix[other_col],
                    other_row,
                    add_powers(shift, row_power)?,
                    &mut cells,
                )?;
            }
            // The inverse row operation adds the other target basis vector
            // to the pivot target vector.
            add_scaled(
                model,
                &mut target_lift,
                &target[other_row],
                shift,
                &mut cells,
            )?;
        }
        cells -= column_cells;
        for (other_col, power) in across {
            add_scaled(
                model,
                &mut source[other_col],
                &source_lift,
                power - exponent,
                &mut cells,
            )?;
        }
        cells -= across_cells;

        // The free cokernel grades lose the degrees of the pivot source
        // vectors through the graded map onto the free image.
        *free.get_mut(&grades[row]).unwrap() -= 1;
        *free.get_mut(&grades[col]).unwrap() -= 1;
        if exponent == 0 {
            cells -= ENTRY_CELLS * (target_lift.len() + source_lift.len());
        } else {
            bars.push(Bar {
                degree: grades[row],
                exponent,
                representative: target_lift,
                nullhomotopy: source_lift,
            });
        }
    }

    let mut degrees = Vec::new();
    for (degree, count) in free {
        let count = usize::try_from(count)
            .map_err(|_| Error::plain("graded Koszul rank subtraction is inconsistent"))?;
        degrees.extend(std::iter::repeat_n(degree, count));
    }
    Ok((degrees, bars))
}

fn pivot(model: &mut Model, matrix: &[Column]) -> Result<Option<(usize, usize, usize)>> {
    let mut best = None;
    model.spend(matrix.len())?;
    for (col, entries) in matrix.iter().enumerate() {
        model.spend(entries.len())?;
        for (&row, &exponent) in entries {
            let candidate = (exponent, row, col);
            if best.is_none_or(|best| candidate < best) {
                best = Some(candidate);
            }
        }
    }
    Ok(best)
}

fn remove_row(model: &mut Model, matrix: &mut [Column], row: usize) -> Result<Vec<(usize, usize)>> {
    let mut removed = Vec::new();
    for (index, col) in matrix.iter_mut().enumerate() {
        model.spend(1 + col.len().max(1).ilog2() as usize)?;
        if let Some(power) = col.remove(&row) {
            removed.push((index, power));
        }
    }
    Ok(removed)
}

fn sort_work(size: usize) -> usize {
    size * (1 + size.max(1).ilog2() as usize)
}

fn add_scaled(
    model: &mut Model,
    target: &mut Column,
    source: &Column,
    shift: usize,
    cells: &mut usize,
) -> Result<()> {
    for (&coordinate, &power) in source {
        toggle(model, target, coordinate, add_powers(power, shift)?, cells)?;
    }
    Ok(())
}

fn add_powers(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b)
        .ok_or_else(|| Error::plain("Smith reduction exponent is out of range"))
}

fn toggle(
    model: &mut Model,
    column: &mut Column,
    coordinate: usize,
    exponent: usize,
    cells: &mut usize,
) -> Result<()> {
    model.spend(1 + column.len().max(1).ilog2() as usize)?;
    match column.entry(coordinate) {
        Entry::Occupied(entry) => {
            if *entry.get() != exponent {
                return Err(Error::plain("Smith reduction lost homogeneity"));
            }
            entry.remove();
            *cells -= ENTRY_CELLS;
        }
        Entry::Vacant(entry) => {
            *cells += ENTRY_CELLS;
            q_koszul::check_storage(*cells)?;
            entry.insert(exponent);
        }
    }
    Ok(())
}

pub(crate) fn ratio_type(model: &mut Model, ty: TypeId, k: usize) -> Result<TypeId> {
    let (grade, class) = match model.kind(ty) {
        Kind::VKoszul(grade) => (grade, false),
        Kind::VExt(grade) => (grade, true),
        _ => {
            return Err(Error::plain(
                "vRatio requires a rank-two VKoszul cochain or VExt class",
            ))
        }
    };
    let family = Family::new(grade.source, grade.target, &grade.selected)?;
    let t = ratio_degree(grade.t, k, family.gap(), "vRatio")?;
    let grade = Grade { t, ..grade.clone() };
    if class {
        model.v_ext_type(grade)
    } else {
        model.v_koszul_type(grade)
    }
}

fn ratio_degree(t: i64, k: usize, gap: usize, operation: &str) -> Result<i64> {
    let delta = (k as i128).checked_mul(gap as i128);
    delta
        .and_then(|delta| i128::from(t).checked_sub(delta))
        .and_then(|degree| i64::try_from(degree).ok())
        .ok_or_else(|| {
            Error::plain(format!(
                "{operation} signed internal degree is out of range"
            ))
        })
}

pub(crate) fn ratio(
    model: &mut Model,
    source: TypeId,
    target: TypeId,
    k: usize,
    value: Value,
) -> Result<Value> {
    if k == 0 {
        return Ok(value);
    }
    v_koszul::multiply(model, source, target, value, |model, terms| {
        model.spend(terms.iter().map(Term::cells).sum())?;
        terms
            .iter()
            .map(|term| {
                let lower = term.powers[0].checked_add(k as i128);
                let upper = term.powers[1].checked_sub(k as i128);
                let (Some(lower), Some(upper)) = (lower, upper) else {
                    return Err(Error::plain("vRatio exponent is out of range"));
                };
                Ok(Term {
                    powers: [lower, upper].into(),
                    ..term.clone()
                })
            })
            .collect()
    })
}
