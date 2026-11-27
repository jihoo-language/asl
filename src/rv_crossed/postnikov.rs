use std::cell::OnceCell;
use std::rc::Rc;

use super::{Choices, Crossed};
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::{self, Basis, Coordinate};
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Command {
    Declare(Label, Label, usize),
    Data(Label),
    Class(Box<ParsedExpr>),
    Coefficient(Vec<usize>, Box<ParsedExpr>),
}

pub(crate) enum Expression {
    Literal(Label, usize, Vec<usize>),
    Canonical(Label),
    FromChoices(Label, Vec<usize>, Vec<Vec<usize>>),
    Coboundary(Box<ParsedExpr>),
    Gauge(Box<[ParsedExpr; 2]>),
    GaugeWitness(Box<[ParsedExpr; 2]>),
    Cell(Label, usize, usize),
    Vertical(Box<[ParsedExpr; 2]>),
    Horizontal(Box<[ParsedExpr; 2]>),
    Associator([usize; 3], Box<ParsedExpr>),
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Cochain(pub(crate) Vec<Coordinate>);

impl Cochain {
    pub(crate) fn zero() -> Self {
        Self(Vec::new())
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(2 * self.0.len() + 1)?;
        Ok(Value::PostnikovCochain(Rc::new(self)))
    }

    pub(super) fn at(&self, i: usize) -> u32 {
        self.0
            .binary_search_by_key(&i, |&(key, _)| key)
            .map_or(0, |j| self.0[j].1)
    }
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Cell {
    pub(crate) object: usize,
    pub(crate) coordinate: u32,
}

pub(crate) struct Postnikov {
    pub(super) crossed: Rc<Crossed>,
    pub(crate) prime: Prime,
    pub(super) action: Vec<u32>,
    kernel_coordinates: Vec<Option<u32>>,
    canonical: Rc<Cochain>,
    boundaries: OnceCell<Basis>,
    rank_d3: OnceCell<usize>,
}

impl Postnikov {
    pub(crate) fn new(model: &mut Model, crossed: Rc<Crossed>, generator: usize) -> Result<Self> {
        if (crossed.order() - 1).pow(4) > MAX_PIVOT_CELLS {
            return Err(Error::plain(
                "normalized Postnikov bar complex exceeds its degree-four size limit",
            ));
        }
        let post = Self::carrier(model, crossed, generator)?;
        post.boundary_basis(model)?;
        post.rank_d3(model)?;
        Ok(post)
    }

    pub(super) fn carrier(
        model: &mut Model,
        crossed: Rc<Crossed>,
        generator: usize,
    ) -> Result<Self> {
        // Checking a supplied cocycle streams over the bar rows. Computing
        // cohomology ranks also retains a sparse elimination basis.
        if (crossed.order() - 1).pow(4) > 1_000_000 {
            return Err(Error::plain(
                "normalized Postnikov cocycle exceeds its degree-four size limit",
            ));
        }
        let size = crossed.kernel.len();
        let prime = if size == 2 {
            Prime::validated(2)
        } else {
            Prime::new(model, size).map_err(|_| {
                Error::plain("Postnikov cochains require a cyclic kernel of prime order")
            })?
        };
        if generator == 0 || crossed.boundary.get(generator) != Some(&0) {
            return Err(Error::plain(
                "Postnikov generator must be a nonidentity kernel element",
            ));
        }
        let mut kernel_coordinates = vec![None; crossed.h.order()];
        let mut power = 0;
        for i in 0..prime.get() {
            if kernel_coordinates[power].replace(i).is_some() {
                return Err(Error::plain(
                    "Postnikov generator does not generate the kernel",
                ));
            }
            power = crossed.h.mul(power, generator);
        }
        let action = (0..crossed.order())
            .map(|g| {
                let image = crossed.kernel_act(g, generator)?;
                kernel_coordinates[image]
                    .ok_or_else(|| Error::plain("skeletal action left the kernel"))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut post = Self {
            crossed,
            prime,
            action,
            kernel_coordinates,
            canonical: Rc::new(Cochain::zero()),
            boundaries: OnceCell::new(),
            rank_d3: OnceCell::new(),
        };
        let canonical = post.compute(model, &post.crossed.standard)?;
        post.require_cocycle(model, &canonical)?;
        post.canonical = Rc::new(canonical);
        Ok(post)
    }

    fn boundary_basis(&self, model: &mut Model) -> Result<&Basis> {
        if let Some(basis) = self.boundaries.get() {
            return Ok(basis);
        }
        let mut delta2 = vec![Vec::new(); self.size(2)];
        for i in 0..self.size(3) {
            for (j, coefficient) in self.row(model, 2, i)? {
                delta2[j].push((i, coefficient));
            }
        }
        let mut basis = Basis::default();
        for (i, column) in delta2.into_iter().enumerate() {
            basis.insert_with_witness(model, self.prime, column, vec![(i, 1)])?;
        }
        model.allocate(basis.cells())?;
        Ok(self.boundaries.get_or_init(|| basis))
    }

    fn rank_d3(&self, model: &mut Model) -> Result<usize> {
        if let Some(&rank) = self.rank_d3.get() {
            return Ok(rank);
        }
        if self.size(4) > MAX_PIVOT_CELLS {
            return Err(Error::plain(
                "normalized Postnikov cohomology rank exceeds its degree-four size limit",
            ));
        }
        let mut delta3 = Basis::default();
        for i in 0..self.size(4) {
            let row = self.row(model, 3, i)?;
            delta3.insert(model, self.prime, row)?;
        }
        Ok(*self.rank_d3.get_or_init(|| delta3.rank()))
    }

    pub(crate) fn cells(&self) -> usize {
        12 + self.action.len() + self.kernel_coordinates.len() + 2 * self.canonical.0.len()
    }

    pub(crate) fn order(&self) -> usize {
        self.crossed.order()
    }

    pub(crate) fn size(&self, degree: usize) -> usize {
        (self.order() - 1).pow(degree as u32)
    }

    pub(crate) fn literal(
        &self,
        model: &mut Model,
        degree: usize,
        dense: Vec<usize>,
    ) -> Result<Cochain> {
        let width = self.size(degree);
        let p = self.prime.get() as usize;
        if dense.len() != width || dense.iter().any(|&x| x >= p) {
            return Err(Error::plain(format!(
                "normalized degree-{degree} cochain requires {width} residues from 0 through {}",
                p - 1
            )));
        }
        model.spend(width)?;
        Ok(Cochain(
            dense
                .into_iter()
                .enumerate()
                .filter(|&(_, x)| x != 0)
                .map(|(i, x)| (i, x as u32))
                .collect(),
        ))
    }

    pub(crate) fn canonical(&self) -> Rc<Cochain> {
        self.canonical.clone()
    }

    pub(crate) fn with_choices(
        &self,
        model: &mut Model,
        section: Vec<usize>,
        lifts: Vec<Vec<usize>>,
    ) -> Result<Cochain> {
        let choices = self.crossed.check_choices(model, section, lifts)?;
        let cochain = self.compute(model, &choices)?;
        self.require_cocycle(model, &cochain)?;
        Ok(cochain)
    }

    fn compute(&self, model: &mut Model, choices: &Choices) -> Result<Cochain> {
        model.spend(8 * self.size(3))?;
        let mut output = Vec::new();
        for i in 0..self.size(3) {
            let [g, h, l, _] = self.tuple(3, i);
            let assoc = self.crossed.associator(choices, [g, h, l])?;
            let value = self.kernel_coordinates[assoc].ok_or_else(|| {
                Error::plain("Postnikov associator left the chosen cyclic kernel")
            })?;
            if value != 0 {
                output.push((i, value));
            }
        }
        Ok(Cochain(output))
    }

    fn index(&self, tuple: &[usize]) -> Option<usize> {
        let mut index = 0;
        for &g in tuple {
            let digit = g.checked_sub(1)?;
            index = index * (self.order() - 1) + digit;
        }
        Some(index)
    }

    fn tuple(&self, degree: usize, mut index: usize) -> [usize; 4] {
        let mut tuple = [0; 4];
        let base = self.order() - 1;
        for digit in tuple[..degree].iter_mut().rev() {
            *digit = 1 + index % base;
            index /= base;
        }
        tuple
    }

    fn row(&self, model: &mut Model, degree: usize, out: usize) -> Result<Vec<Coordinate>> {
        model.spend(degree + 2)?;
        let tuple = self.tuple(degree + 1, out);
        let mut entries = Vec::with_capacity(degree + 2);
        let mut put = |t: &[usize], coefficient| {
            if let Some(i) = self.index(t) {
                entries.push((i, coefficient));
            }
        };
        let minus = self.prime.negate(1);
        put(&tuple[1..=degree], self.action[tuple[0]]);
        for i in 0..degree {
            let mut merged = tuple;
            merged[i] = self.crossed.quotient_mul(tuple[i], tuple[i + 1]);
            merged.copy_within(i + 2..degree + 1, i + 1);
            put(
                &merged[..degree],
                if i.is_multiple_of(2) { minus } else { 1 },
            );
        }
        put(
            &tuple[..degree],
            if degree.is_multiple_of(2) { minus } else { 1 },
        );
        entries.sort_unstable_by_key(|&(i, _)| i);
        let mut row = Vec::<Coordinate>::new();
        for (i, value) in entries {
            if let Some(last) = row.last_mut().filter(|last| last.0 == i) {
                last.1 = self.prime.add(last.1, value);
            } else {
                row.push((i, value));
            }
        }
        row.retain(|&(_, value)| value != 0);
        Ok(row)
    }

    fn evaluate_row(&self, row: &[Coordinate], cochain: &Cochain) -> u32 {
        row.iter().fold(0, |sum, &(key, value)| {
            self.prime.add(sum, self.prime.mul(value, cochain.at(key)))
        })
    }

    pub(super) fn require_cocycle(&self, model: &mut Model, cochain: &Cochain) -> Result<()> {
        if cochain.0.is_empty() {
            return Ok(());
        }
        for i in 0..self.size(4) {
            if self.evaluate_row(&self.row(model, 3, i)?, cochain) != 0 {
                return Err(Error::plain(
                    "Postnikov 3-cochain fails the twisted cocycle (pentagon) equation",
                ));
            }
        }
        Ok(())
    }

    pub(crate) fn coefficient(
        &self,
        model: &mut Model,
        cochain: &Cochain,
        tuple: &[usize],
    ) -> Result<u32> {
        for &g in tuple {
            self.crossed.object(g)?;
        }
        model.spend(1 + cochain.0.len().max(1).ilog2() as usize)?;
        Ok(self.index(tuple).map_or(0, |i| cochain.at(i)))
    }

    pub(crate) fn coboundary(&self, model: &mut Model, cochain: &Cochain) -> Result<Cochain> {
        let mut result = Vec::new();
        if !cochain.0.is_empty() {
            for i in 0..self.size(3) {
                let value = self.evaluate_row(&self.row(model, 2, i)?, cochain);
                if value != 0 {
                    result.push((i, value));
                }
            }
        }
        Ok(Cochain(result))
    }

    pub(crate) fn add(
        &self,
        model: &mut Model,
        left: &Cochain,
        right: &Cochain,
    ) -> Result<Cochain> {
        Ok(Cochain(prime_linear::add_scaled(
            model, self.prime, &left.0, &right.0, 1,
        )?))
    }

    pub(crate) fn neg(&self, model: &mut Model, input: &Cochain) -> Result<Cochain> {
        model.spend(input.0.len())?;
        Ok(Cochain(
            input
                .0
                .iter()
                .map(|&(i, value)| (i, self.prime.negate(value)))
                .collect(),
        ))
    }

    pub(crate) fn gauge(
        &self,
        model: &mut Model,
        cocycle: &Cochain,
        witness: &Cochain,
    ) -> Result<Cochain> {
        self.require_cocycle(model, cocycle)?;
        let delta = self.coboundary(model, witness)?;
        self.add(model, cocycle, &delta)
    }

    pub(crate) fn gauge_witness(
        &self,
        model: &mut Model,
        first: &Cochain,
        second: &Cochain,
    ) -> Result<Cochain> {
        self.require_cocycle(model, first)?;
        self.require_cocycle(model, second)?;
        let difference =
            prime_linear::add_scaled(model, self.prime, &second.0, &first.0, self.prime.negate(1))?;
        let (rest, witness) = self
            .boundary_basis(model)?
            .reduce_with_witness(model, self.prime, difference)?;
        if !rest.is_empty() {
            return Err(Error::plain("Postnikov cocycles are not gauge equivalent"));
        }
        Ok(Cochain(witness))
    }

    pub(crate) fn class_report(&self, model: &mut Model, cochain: &Cochain) -> Result<String> {
        self.require_cocycle(model, cochain)?;
        let remainder = self
            .boundary_basis(model)?
            .reduce(model, self.prime, cochain.0.clone())?;
        Ok(format!(
            "rvPostnikovClass = nonzero={}, reduced={remainder:?}",
            !remainder.is_empty()
        ))
    }

    pub(crate) fn data(&self, model: &mut Model, name: &str) -> Result<String> {
        let rank_d3 = self.rank_d3(model)?;
        let basis = self.boundary_basis(model)?;
        let rank_d2 = basis.rank();
        let dimension = self
            .size(3)
            .checked_sub(rank_d2 + rank_d3)
            .ok_or_else(|| Error::plain("normalized Postnikov bar ranks are inconsistent"))?;
        let nonzero = !basis
            .reduce(model, self.prime, self.canonical.0.clone())?
            .is_empty();
        Ok(format!("rvPostnikovData[{name}] = quotient={}, prime={}, action={:?}, rank_d2={rank_d2}, rank_d3={rank_d3}, dim_h3={dimension}, nonzero={nonzero}", self.order(), self.prime.get(), self.action))
    }

    pub(crate) fn cell(&self, object: usize, value: usize) -> Result<Cell> {
        self.crossed.object(object)?;
        if value >= self.prime.get() as usize {
            return Err(Error::plain(
                "Postnikov cell coordinate is outside its prime field",
            ));
        }
        Ok(Cell {
            object,
            coordinate: value as u32,
        })
    }

    pub(crate) fn vertical(&self, a: Cell, b: Cell) -> Result<Cell> {
        if a.object != b.object {
            return Err(Error::plain(
                "vertical Postnikov composition requires the same object",
            ));
        }
        Ok(Cell {
            object: a.object,
            coordinate: self.prime.add(a.coordinate, b.coordinate),
        })
    }

    pub(crate) fn horizontal(&self, a: Cell, b: Cell) -> Cell {
        Cell {
            object: self.crossed.quotient_mul(a.object, b.object),
            coordinate: self.prime.add(
                a.coordinate,
                self.prime.mul(self.action[a.object], b.coordinate),
            ),
        }
    }

    pub(crate) fn associator_cell(
        &self,
        model: &mut Model,
        cochain: &Cochain,
        [g, h, l]: [usize; 3],
    ) -> Result<Cell> {
        self.require_cocycle(model, cochain)?;
        let coordinate = self.coefficient(model, cochain, &[g, h, l])?;
        let object = self
            .crossed
            .quotient_mul(self.crossed.quotient_mul(g, h), l);
        Ok(Cell { object, coordinate })
    }
}
