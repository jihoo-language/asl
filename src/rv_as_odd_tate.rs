use std::cell::{Cell, OnceCell, RefCell};
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::binary_linear::check_storage;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::odd_ext;
use crate::prime_linear::{add_scaled, Basis, Coordinate, QuotientBasis};
use crate::rv_as_odd::Cover;
use crate::rv_as_odd_cartier::Column;
use crate::rv_as_odd_derham::{DeRham, Map};

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Differential(i64),
    Class(i64),
    Representative,
    Frobenius,
    Verschiebung,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Differential(_) => "rvASOddTateD",
            Self::Class(_) => "rvASOddTateClass",
            Self::Representative => "rvASOddTateRep",
            Self::Frobenius => "rvASOddTateF",
            Self::Verschiebung => "rvASOddTateV",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Differential {
    Difference,
    Norm,
}

impl Differential {
    fn at(degree: i64) -> Self {
        if degree % 2 == 0 {
            Self::Difference
        } else {
            Self::Norm
        }
    }

    fn incoming(self) -> Self {
        match self {
            Self::Difference => Self::Norm,
            Self::Norm => Self::Difference,
        }
    }

    fn index(self) -> usize {
        match self {
            Self::Difference => 0,
            Self::Norm => 1,
        }
    }
}

pub(crate) struct Tate {
    prime: Prime,
    dimension: usize,
    derham: Rc<DeRham>,
    columns: RefCell<BTreeMap<(Differential, usize), Column>>,
    groups: [OnceCell<QuotientBasis>; 2],
    cells: Cell<usize>,
}

impl Tate {
    pub(crate) fn build(model: &mut Model, cover: &Cover) -> Result<Self> {
        let dimension = cover
            .genus()
            .checked_mul(2)
            .ok_or_else(|| Error::plain("odd Artin–Schreier de Rham dimension is out of range"))?;
        Ok(Self {
            prime: cover.prime(),
            dimension,
            derham: cover.derham(model)?,
            columns: RefCell::new(BTreeMap::new()),
            groups: [OnceCell::new(), OnceCell::new()],
            cells: Cell::new(0),
        })
    }

    fn column(
        &self,
        model: &mut Model,
        cover: &Cover,
        differential: Differential,
        index: usize,
    ) -> Result<Column> {
        // A nilpotent operator on a space of this dimension has zero
        // (p-1)-st power if the dimension is below the length of a free block.
        if differential == Differential::Norm && cover.order() > self.dimension {
            return Ok(Rc::from([]));
        }
        if let Some(column) = self.columns.borrow().get(&(differential, index)) {
            return Ok(column.clone());
        }
        let unit = [(index, 1)];
        let image = match differential {
            Differential::Difference => {
                let translated = self.derham.deck(model, cover, 1, &unit)?;
                add_scaled(model, self.prime, &translated, &unit, self.prime.negate(1))?
            }
            Differential::Norm => {
                let mut image = unit.to_vec();
                for _ in 1..cover.order() {
                    image = self.image(model, cover, Differential::Difference, &image)?;
                    if image.is_empty() {
                        break;
                    }
                }
                image
            }
        };
        let cells = 6 + 2 * image.len();
        let total = self.cells.get() + cells;
        check_storage(total)?;
        model.allocate(cells)?;
        self.cells.set(total);
        let column: Column = image.into();
        self.columns
            .borrow_mut()
            .insert((differential, index), column.clone());
        Ok(column)
    }

    fn image(
        &self,
        model: &mut Model,
        cover: &Cover,
        differential: Differential,
        input: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        let mut image = Vec::new();
        for &(index, scalar) in input {
            let column = self.column(model, cover, differential, index)?;
            image = add_scaled(model, self.prime, &image, &column, scalar)?;
        }
        Ok(image)
    }

    pub(crate) fn group(
        &self,
        model: &mut Model,
        cover: &Cover,
        degree: i64,
    ) -> Result<&QuotientBasis> {
        let outgoing = Differential::at(degree);
        let group = &self.groups[outgoing.index()];
        if let Some(quotient) = group.get() {
            return Ok(quotient);
        }
        check_storage(self.dimension)?;
        model.spend(self.dimension)?;
        let mut boundaries = Basis::default();
        for index in 0..self.dimension {
            let column = self.column(model, cover, outgoing.incoming(), index)?;
            boundaries.insert(model, self.prime, column.to_vec())?;
        }
        let quotient = QuotientBasis::build(
            model,
            self.prime,
            boundaries,
            self.dimension,
            |model, index| {
                let column = self.column(model, cover, outgoing, index)?;
                Ok(column.to_vec())
            },
        )?;
        Ok(group.get_or_init(|| quotient))
    }

    fn project(
        &self,
        model: &mut Model,
        cover: &Cover,
        degree: i64,
        raw: Vec<Coordinate>,
    ) -> Result<Vec<Coordinate>> {
        self.group(model, cover, degree)?
            .project(model, self.prime, raw)?
            .ok_or_else(|| Error::plain(format!("input is not a Tate cycle in degree {degree}")))
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        cover: &Cover,
        degree: i64,
        operation: Operation,
        input: Column,
    ) -> Result<Value> {
        if input.is_empty() {
            return Ok(Value::OddClass(input));
        }
        let result = match operation {
            Operation::Differential(_) => {
                self.image(model, cover, Differential::at(degree), &input)?
            }
            Operation::Class(_) => self.project(model, cover, degree, input.to_vec())?,
            Operation::Representative => self
                .group(model, cover, degree)?
                .representative(model, self.prime, &input)?,
            Operation::Frobenius | Operation::Verschiebung => {
                let operation = if matches!(operation, Operation::Frobenius) {
                    Map::Frobenius
                } else {
                    Map::Verschiebung
                };
                let representative = self
                    .group(model, cover, degree)?
                    .representative(model, self.prime, &input)?;
                let image = self
                    .derham
                    .image(model, cover, operation, &representative)?;
                self.project(model, cover, degree, image)?
            }
        };
        odd_ext::persist(model, result)
    }

    pub(crate) fn pair(
        &self,
        model: &mut Model,
        cover: &Cover,
        left: (i64, &[Coordinate]),
        right: (i64, &[Coordinate]),
    ) -> Result<Value> {
        let first = self
            .group(model, cover, left.0)?
            .representative(model, self.prime, left.1)?;
        let second = self
            .group(model, cover, right.0)?
            .representative(model, self.prime, right.1)?;
        self.derham.pair(model, cover, &first, &second)
    }

    pub(crate) fn jordan(&self, model: &mut Model, cover: &Cover) -> Result<Vec<usize>> {
        check_storage(self.dimension)?;
        model.spend(self.dimension)?;
        let mut images = Basis::default();
        for index in 0..self.dimension {
            let column = self.column(model, cover, Differential::Difference, index)?;
            images.insert(model, self.prime, column.to_vec())?;
        }
        let mut ranks = vec![self.dimension, images.rank()];
        while images.rank() != 0 {
            let mut next = Basis::default();
            for column in images.vectors(self.prime) {
                let image = self.image(model, cover, Differential::Difference, &column)?;
                next.insert(model, self.prime, image)?;
            }
            if next.rank() >= images.rank() {
                return Err(Error::plain(
                    "odd Artin–Schreier deck difference is not nilpotent",
                ));
            }
            images = next;
            ranks.push(images.rank());
        }
        ranks.push(0);
        let mut blocks = Vec::new();
        for (i, ranks) in ranks.windows(3).enumerate().rev() {
            let count = (ranks[0] - ranks[1])
                .checked_sub(ranks[1] - ranks[2])
                .ok_or_else(|| Error::plain("invalid ranks for the deck difference"))?;
            blocks.extend(std::iter::repeat_n(i + 1, count));
        }
        Ok(blocks)
    }
}
