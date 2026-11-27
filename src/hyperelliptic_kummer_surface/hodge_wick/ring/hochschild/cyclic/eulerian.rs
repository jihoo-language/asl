use super::complex::Axis;
use super::maps::boundary;
use super::{
    storage, work_limit, Algebra, BTreeMap, Chain, Cochain, Complex, Error, Flavor, Form, Grade,
    Key, Model, ParsedExpr, Rational, Result, SyntaxGrade, LIMBS,
};

pub(crate) mod cotangent;
mod words;

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Adams(i128),
    Project(i128),
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Adams(_) => "rvVerbitskyAdams",
            Self::Project(_) => "rvVerbitskyHodgeProject",
        }
    }

    pub(crate) fn check(self, flavor: Flavor) -> Result<()> {
        match self {
            Self::Adams(m) if m < 0 => {
                Err(Error::plain("Verbitsky Adams index must be nonnegative"))
            }
            Self::Adams(_) if flavor != Flavor::Hochschild => Err(Error::plain(
                "Verbitsky Adams operation acts on Hochschild chains or cochains",
            )),
            _ => Ok(()),
        }
    }
}

pub(crate) struct Expression(pub(crate) Operation, pub(crate) Box<ParsedExpr>);

pub(crate) struct Command(pub(crate) SyntaxGrade, pub(crate) i128);

fn transform(
    model: &mut Model,
    operation: Operation,
    power: i128,
    word: &[usize],
    transpose: bool,
) -> Result<words::Words> {
    match operation {
        Operation::Adams(m) => words::adams(model, word, m, transpose),
        Operation::Project(index) => {
            // u has Hodge index -1; a u^power column is of bar index index+power.
            let Some(bar_index) = index.checked_add(power) else {
                return Ok(words::Words::new());
            };
            words::project(model, word, bar_index, transpose)
        }
    }
}

fn chain_operator(model: &mut Model, operation: Operation, value: &Chain) -> Result<Chain> {
    let mut result = Chain::default();
    for (key, &coefficient) in &value.entries {
        for (bars, scalar) in transform(model, operation, key.power, &key.word[1..], false)? {
            let mut word = Vec::with_capacity(bars.len() + 1);
            word.push(key.word[0]);
            word.extend(bars);
            let value = coefficient.multiply(model, scalar)?;
            result.insert(
                model,
                Key {
                    power: key.power,
                    word,
                },
                value,
            )?;
        }
    }
    Ok(result)
}

impl Form {
    pub(crate) fn cyclic_hodge_dimension(
        &self,
        model: &mut Model,
        grade: &Grade,
        index: i128,
    ) -> Result<usize> {
        self.cyclic_validate(grade)?;
        if !grade.possible() {
            return Ok(0);
        }
        let algebra = self.cyclic_algebra(model, grade)?;
        let current = hodge_basis(model, &algebra, grade, index)?;
        if current.rows.is_empty() {
            return Ok(0);
        }
        let previous_grade = grade.with_offset(1)?;
        let previous = hodge_basis(model, &algebra, &previous_grade, index)?;
        let outgoing = boundary_rank(model, &algebra, grade, &current)?;
        let incoming = boundary_rank(model, &algebra, &previous_grade, &previous)?;
        current
            .rows
            .len()
            .checked_sub(outgoing + incoming)
            .ok_or_else(|| {
                Error::plain("Verbitsky Hodge ranks violate the mixed complex identities")
            })
    }

    pub(crate) fn cyclic_hodge_operator(
        &self,
        model: &mut Model,
        grade: &Grade,
        operation: Operation,
        class: bool,
        value: &Chain,
    ) -> Result<Chain> {
        let result = chain_operator(model, operation, value)?;
        if class && !result.empty() {
            let algebra = self.cyclic_algebra(model, grade)?;
            Complex::new(model, &algebra, grade)?.project(model, &result)
        } else {
            Ok(result)
        }
    }

    pub(crate) fn hochschild_hodge_operator(
        &self,
        model: &mut Model,
        grade: &super::super::Grade,
        operation: Operation,
        class: bool,
        value: &Cochain,
    ) -> Result<Cochain> {
        let mut result = Cochain::default();
        for (key, &coefficient) in &value.entries {
            for (word, scalar) in transform(model, operation, 0, &key.word, true)? {
                let value = coefficient.multiply(model, scalar)?;
                result.insert(
                    model,
                    super::super::Key {
                        word,
                        output: key.output,
                    },
                    value,
                )?;
            }
        }
        if class {
            self.hochschild_class(model, grade, &result)
        } else {
            Ok(result)
        }
    }
}

#[derive(Default)]
struct SparseBasis {
    rows: BTreeMap<Key, Chain>,
    cells: usize,
}

impl SparseBasis {
    fn insert(&mut self, model: &mut Model, mut row: Chain) -> Result<()> {
        while let Some((key, &leading)) = row.entries.first_key_value() {
            if let Some(pivot) = self.rows.get(key) {
                let scale = leading.negate();
                for (key, &coefficient) in &pivot.entries {
                    let correction = coefficient.multiply(model, scale)?;
                    row.insert(model, key.clone(), correction)?;
                }
            } else {
                let key = key.clone();
                let inverse = leading.inverse()?;
                for coefficient in row.entries.values_mut() {
                    *coefficient = coefficient.multiply(model, inverse)?;
                }
                self.cells += row.cells + key.word.len() + 12;
                storage(self.cells)?;
                self.rows.insert(key, row);
                break;
            }
        }
        Ok(())
    }
}

fn hodge_basis(
    model: &mut Model,
    algebra: &Algebra,
    grade: &Grade,
    index: i128,
) -> Result<SparseBasis> {
    let mut basis = SparseBasis::default();
    if grade.flavor == Flavor::Hochschild
        && (grade.degree == 0 && index != 0
            || grade.degree > 0 && !(1..=grade.degree).contains(&index))
    {
        return Ok(basis);
    }
    let axis = Axis::new(model, algebra, grade)?;
    for key in axis.keys() {
        let mut term = Chain::default();
        term.insert(model, key.clone(), Rational::integer(1))?;
        let projected = chain_operator(model, Operation::Project(index), &term)?;
        basis.insert(model, projected)?;
    }
    Ok(basis)
}

fn boundary_rank(
    model: &mut Model,
    algebra: &Algebra,
    grade: &Grade,
    space: &SparseBasis,
) -> Result<usize> {
    let mut image = SparseBasis::default();
    for row in space.rows.values() {
        let differential = boundary(model, algebra, grade, row)?;
        image.insert(model, differential)?;
    }
    Ok(image.rows.len())
}
