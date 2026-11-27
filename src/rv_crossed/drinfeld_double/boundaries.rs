use std::rc::Rc;

use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::qz::Scalar;
use crate::model::Model;
use crate::rv_cyclotomic;
use crate::syntax::Label;

use super::{negate_phase, Double, Sector};

mod bar;
use bar::{Solution, Subgroup, Witness};

pub(crate) type Fractions = Vec<(i128, i128)>;
pub(crate) type Chain = Vec<([usize; 3], i128)>;

pub(crate) enum Command {
    Declare(Label, Label, Vec<usize>, Option<Fractions>),
    Shift(Label, Label, Fractions),
    Data(Label),
    Witness(Label),
    Cycle(Label, Vec<usize>, Chain),
    Annulus(Label, Label),
}

pub(crate) enum Expression {
    Cochain(Label, [usize; 2]),
    Condensate(Label),
}

pub(crate) struct Boundary {
    pub(crate) double_name: Rc<str>,
    subgroup: Rc<Subgroup>,
    status: Status,
}

enum Status {
    Admissible {
        cochain: Vec<Scalar>,
        condensate: Rc<Sector>,
        dimension: usize,
    },
    Obstructed(Witness),
}

impl Boundary {
    pub(crate) fn new(
        model: &mut Model,
        double_name: Rc<str>,
        double: &Double,
        members: Vec<usize>,
        supplied: Option<Fractions>,
    ) -> Result<Self> {
        let subgroup = Rc::new(Subgroup::new(model, double, members)?);
        let status = match supplied {
            Some(dense) => {
                let cochain = subgroup.cochain(model, dense)?;
                Self::admissible(model, double, &subgroup, cochain)?
            }
            None => match subgroup.solve(model, double)? {
                Solution::Cochain(cochain) => Self::admissible(model, double, &subgroup, cochain)?,
                Solution::Obstruction(witness) => {
                    model.allocate(18 + 5 * witness.terms.len())?;
                    Status::Obstructed(witness)
                }
            },
        };
        model.allocate(6 + subgroup.cells())?;
        Ok(Self {
            double_name,
            subgroup,
            status,
        })
    }

    fn admissible(
        model: &mut Model,
        double: &Double,
        subgroup: &Subgroup,
        cochain: Vec<Scalar>,
    ) -> Result<Status> {
        subgroup.verify(model, double, &cochain)?;
        let condensate = double.condensate(model, subgroup, &cochain)?;
        model.allocate(16 * cochain.len() + 3 * condensate.0.len() + 6)?;
        Ok(Status::Admissible {
            cochain,
            condensate: Rc::new(condensate),
            dimension: double.group.order(),
        })
    }

    pub(crate) fn shifted(
        &self,
        model: &mut Model,
        double: &Double,
        shift: Fractions,
    ) -> Result<Self> {
        let (base, _) = self.parts()?;
        let values = self.subgroup.cochain(model, shift)?;
        let cochain = base
            .iter()
            .zip(values)
            .map(|(&a, b)| a.add(model, b))
            .collect::<Result<Vec<_>>>()?;
        let status = Self::admissible(model, double, &self.subgroup, cochain)?;
        model.allocate(5)?;
        Ok(Self {
            status,
            double_name: self.double_name.clone(),
            subgroup: self.subgroup.clone(),
        })
    }

    fn parts(&self) -> Result<(&[Scalar], Rc<Sector>)> {
        match &self.status {
            Status::Admissible {
                cochain,
                condensate,
                ..
            } => Ok((cochain, condensate.clone())),
            Status::Obstructed(_) => {
                Err(Error::plain("double boundary has an obstructed subgroup"))
            }
        }
    }

    pub(crate) fn cochain(&self, [g, h]: [usize; 2]) -> Result<Scalar> {
        self.subgroup.mu(self.parts()?.0, g, h)
    }

    pub(crate) fn sector(&self) -> Result<Rc<Sector>> {
        Ok(self.parts()?.1)
    }

    pub(crate) fn data(&self, name: &str) -> String {
        let group = &self.subgroup.members;
        match &self.status {
            Status::Admissible {condensate,dimension,..} => format!(
                "rvDoubleBoundaryData[{name}] = subgroup={group:?}, obstructed=false, quantum_dimension={dimension}, sectors={}, vacuum=1, modular_fixed=true",
                condensate.0.len(),
            ),
            Status::Obstructed(proof) => format!(
                "rvDoubleBoundaryData[{name}] = subgroup={group:?}, obstructed=true, witness_terms={}, pairing={}",
                proof.terms.len(), proof.value.text(),
            ),
        }
    }

    pub(crate) fn witness(&self, name: &str) -> Result<String> {
        let Status::Obstructed(proof) = &self.status else {
            return Err(Error::plain(
                "this double boundary has a trivializing cochain",
            ));
        };
        Ok(format!(
            "rvDoubleBoundaryWitness[{name}] = cycle={:?}, pairing={}",
            proof.terms,
            proof.value.text()
        ))
    }

    pub(crate) fn cycle(
        model: &mut Model,
        double: &Double,
        name: &str,
        members: Vec<usize>,
        terms: Chain,
    ) -> Result<String> {
        let group = Subgroup::new(model, double, members)?;
        let proof = group.explicit_cycle(model, double, terms)?;
        Ok(format!(
            "rvDoubleBoundaryCycle[{name}] = terms={}, pairing={}, obstructs={}",
            proof.terms.len(),
            proof.value.text(),
            proof.value != Scalar::zero()
        ))
    }

    pub(crate) fn annulus(&self, model: &mut Model, other: &Self) -> Result<i128> {
        if self.double_name != other.double_name {
            return Err(Error::plain(
                "double annulus requires the same nominal bulk double",
            ));
        }
        let (a, b) = (self.sector()?, other.sector()?);
        model.spend(a.0.len() + b.0.len())?;
        let mut value = 0;
        for &(i, m) in &a.0 {
            value += m * b.coefficient(i);
        }
        Ok(value)
    }
}

impl Double {
    fn condensate(
        &self,
        model: &mut Model,
        group: &Subgroup,
        cochain: &[Scalar],
    ) -> Result<Sector> {
        let phases = self.boundary_characters(model, group, cochain)?;
        let mut coeff = vec![0; self.count()];
        let h = group.order() as i128;
        for (i, simple) in self.simples.iter().enumerate() {
            if !group.contains(simple.flux) {
                continue;
            }
            model.spend(self.conductor + group.order())?;
            let character = &phases[group.index(simple.flux)?];
            let mut terms = vec![0; self.conductor];
            for (local, &x) in group.members.iter().enumerate() {
                if let Some(a) = simple.phases[x] {
                    let power =
                        (a + negate_phase(self.conductor, character[local])) % self.conductor;
                    terms[power] += simple.dimension as i128;
                }
            }
            let multiplicity = rv_cyclotomic::reduced(model, self.conductor, terms)?;
            match multiplicity.as_slice() {
                [] => {}
                &[x] if x >= 0 && x % h == 0 => coeff[i] = x / h,
                _ => {
                    return Err(Error::plain(
                        "double boundary projective multiplicity is not a nonnegative integer",
                    ))
                }
            }
        }
        self.verify_condensate(model, &coeff)?;
        Ok(Sector(
            coeff
                .into_iter()
                .enumerate()
                .filter(|&(_, x)| x != 0)
                .collect(),
        ))
    }

    fn boundary_characters(
        &self,
        model: &mut Model,
        group: &Subgroup,
        cochain: &[Scalar],
    ) -> Result<Vec<Vec<usize>>> {
        let n = self.conductor;
        let mut phases = Vec::with_capacity(group.order());
        for &g in &group.members {
            let mut character = Vec::with_capacity(group.order());
            for &x in &group.members {
                let value = group
                    .mu(cochain, x, g)?
                    .add(model, group.mu(cochain, g, x)?.negate())?;
                character.push(value.root_index(n).ok_or_else(|| {
                    Error::plain("boundary character phase is outside the double conductor")
                })?);
            }
            model.spend(group.order() * group.order())?;
            for (i, &x) in group.members.iter().enumerate() {
                for (j, &y) in group.members.iter().enumerate() {
                    let xy = group.index(self.group.mul(x, y))?;
                    if (character[i] + character[j]) % n
                        != (character[xy] + self.tensor_phase(g, x, y)) % n
                    {
                        return Err(Error::plain(
                            "double boundary phase fails its projective character equation",
                        ));
                    }
                }
            }
            phases.push(character);
        }
        Ok(phases)
    }

    fn verify_condensate(&self, model: &mut Model, coeff: &[i128]) -> Result<()> {
        let n = self.group.order() as i128;
        let dimension = self
            .simples
            .iter()
            .zip(coeff)
            .map(|(s, &m)| s.dimension as i128 * m)
            .sum::<i128>();
        if coeff[0] != 1 || dimension != n {
            return Err(Error::plain(
                "double condensate fails its vacuum or quantum dimension check",
            ));
        }
        for (i, row) in self.s.iter().enumerate() {
            if coeff[i] != 0 && self.twists[i] != 0 {
                return Err(Error::plain(
                    "double condensate is not fixed by the modular T matrix",
                ));
            }
            model.spend(row.len() + self.conductor)?;
            let mut terms = vec![0; self.conductor];
            for (cell, &m) in row.iter().zip(coeff) {
                terms[cell.phase] += cell.weight * m;
            }
            let reduced = rv_cyclotomic::reduced(model, self.conductor, terms)?;
            if (coeff[i] == 0 && !reduced.is_empty())
                || (coeff[i] != 0 && reduced != [n * coeff[i]])
            {
                return Err(Error::plain(
                    "double condensate is not fixed by the modular S matrix",
                ));
            }
        }
        Ok(())
    }
}
