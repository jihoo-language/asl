use std::rc::Rc;

use super::{
    add_scaled, column_integer, dot, dual_plane, pairing_row, polarized_basis, raw, realize,
    render_integer, Class, Column, ColumnInteger, Cusp, Integer, Lattice, COLUMN_CELLS, RAW_RANK,
};
use crate::error::{Error, Result};
use crate::model::Model;
use crate::syntax::{Label, ParsedExpr};
use crate::wide_integer::Natural;

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Power(usize),
    Exp(i128),
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Power(_) => "rvOG6NilpotentPower",
            Self::Exp(_) => "rvOG6NilpotentExp",
        }
    }
}

pub(crate) struct Expression(
    pub(crate) Label,
    pub(crate) Operation,
    pub(crate) Box<ParsedExpr>,
);

pub(crate) enum Command {
    TypeII(Label, Label, i128),
    TypeIII(Label, Label, usize, Box<ParsedExpr>),
    Data(Label),
    InWeight(Label, i128, Box<ParsedExpr>),
}

pub(crate) struct Limit {
    pub(crate) cusp: Rc<Cusp>,
    direction: Direction,
    smith: [Vec<ColumnInteger>; 2],
}

enum Direction {
    Central {
        k: i128,
        dual: Box<[Column; 2]>,
    },
    Eichler {
        index: usize,
        u: Rc<Class>,
        square: Integer,
    },
}

impl Limit {
    pub(crate) fn type_ii(
        model: &mut Model,
        lattice: &Lattice,
        cusp: Rc<Cusp>,
        k: i128,
    ) -> Result<Self> {
        if k == 0 {
            return Err(Error::plain(
                "rvog6nilpotentii requires a nonzero central coefficient",
            ));
        }
        let [e, f] = &cusp.oriented_plane;
        let dual = Box::new(dual_plane(model, &raw(e), &raw(f))?);
        Self::new(model, lattice, cusp, Direction::Central { k, dual })
    }

    pub(crate) fn type_iii(
        model: &mut Model,
        lattice: &Lattice,
        cusp: Rc<Cusp>,
        index: usize,
        u: Rc<Class>,
    ) -> Result<Self> {
        let a = cusp
            .oriented_plane
            .get(index)
            .ok_or_else(|| Error::plain("rvog6nilpotentiii axis must be 0 or 1"))?;
        if !u.wide_pair(model, &cusp.polarization)?.is_zero() || !u.wide_pair(model, a)?.is_zero() {
            return Err(Error::plain(
                "rvog6nilpotentiii requires a direction perpendicular to the polarization and the selected axis",
            ));
        }
        let square = u.wide_pair(model, &u)?;
        if square.negative || square.is_zero() {
            return Err(Error::plain(
                "rvog6nilpotentiii requires a direction of positive square",
            ));
        }
        Self::new(
            model,
            lattice,
            cusp,
            Direction::Eichler { index, u, square },
        )
    }

    fn new(
        model: &mut Model,
        lattice: &Lattice,
        cusp: Rc<Cusp>,
        direction: Direction,
    ) -> Result<Self> {
        let mut limit = Self {
            cusp,
            direction,
            smith: [Vec::new(), Vec::new()],
        };
        let h = raw(&limit.cusp.polarization);
        let basis = polarized_basis(model, lattice, &h)?;
        let mut first = Vec::with_capacity(basis.len());
        let mut second = Vec::with_capacity(basis.len());
        for column in &basis {
            first.push(limit.image(model, 1, column)?);
            second.push(limit.image(model, 2, column)?);
        }
        // H is saturated in the raw integral coordinates, so the
        // determinantal divisors of columns originating in H are also
        // their determinantal divisors in any integral basis of H.
        limit.smith[0] = smith_image(model, &first, 2)?;
        let rank2 = usize::from(matches!(limit.direction, Direction::Eichler { .. }));
        limit.smith[1] = smith_image(model, &second, rank2)?;
        Ok(limit)
    }

    pub(super) fn spin_factors(&self) -> (&Class, &Class, i128) {
        match &self.direction {
            Direction::Central { k, .. } => {
                let [e, f] = &self.cusp.oriented_plane;
                (f, e, *k)
            }
            Direction::Eichler { index, u, .. } => (u, &self.cusp.oriented_plane[*index], 1),
        }
    }

    pub(crate) fn cells(&self) -> usize {
        let direction = match self.direction {
            Direction::Central { .. } => 2 * RAW_RANK,
            Direction::Eichler { .. } => 1,
        };
        let smith_entries = self.smith.iter().map(Vec::len).sum::<usize>();
        5 + (direction + smith_entries) * COLUMN_CELLS
    }

    fn image(&self, model: &mut Model, order: usize, x: &Column) -> Result<Column> {
        let mut result = [column_integer(0); RAW_RANK];
        match (&self.direction, order) {
            (_, 0) => return Ok(*x),
            (Direction::Central { k, .. }, 1) => {
                let [e, f] = self.cusp.oriented_plane.each_ref().map(|a| raw(a));
                let ex = dot(model, &pairing_row(&e)?, x)?.multiply(*k)?;
                let fx = dot(model, &pairing_row(&f)?, x)?.multiply(*k)?.negate();
                add_scaled(model, &mut result, &f, ex)?;
                add_scaled(model, &mut result, &e, fx)?;
            }
            (Direction::Eichler { index, u, square }, 1 | 2) => {
                let a = raw(&self.cusp.oriented_plane[*index]);
                let ax = dot(model, &pairing_row(&a)?, x)?;
                if order == 1 {
                    let u = raw(u);
                    let ux = dot(model, &pairing_row(&u)?, x)?;
                    add_scaled(model, &mut result, &u, ax)?;
                    add_scaled(model, &mut result, &a, ux.negate())?;
                } else {
                    let square: ColumnInteger = square.resize()?;
                    add_scaled(model, &mut result, &a, square.multiply_wide(ax)?.negate())?;
                }
            }
            _ => {}
        }
        Ok(result)
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        operation: Operation,
        class: &Class,
    ) -> Result<Class> {
        let x = raw(class);
        match operation {
            Operation::Power(order) => realize(&self.image(model, order, &x)?),
            Operation::Exp(m) => {
                let mut result = x;
                if m != 0 {
                    let first = self.image(model, 1, &x)?;
                    add_scaled(model, &mut result, &first, column_integer(m))?;
                    if matches!(self.direction, Direction::Eichler { .. }) {
                        let second = self.image(model, 2, &x)?;
                        // The even lattice makes each coordinate of N²(x) even.
                        let half = column_integer(2);
                        let mut half_second = [column_integer(0); RAW_RANK];
                        for (out, value) in half_second.iter_mut().zip(second) {
                            *out = value
                                .divide_exact(model, half)?
                                .expect("the even square of the direction divides the second term");
                        }
                        add_scaled(
                            model,
                            &mut result,
                            &half_second,
                            ColumnInteger::product(m, m),
                        )?;
                    }
                }
                realize(&result)
            }
        }
    }

    pub(crate) fn in_weight(&self, model: &mut Model, weight: i128, class: &Class) -> Result<bool> {
        if !class.wide_pair(model, &self.cusp.polarization)?.is_zero() {
            return Err(Error::plain(
                "rvOG6NilpotentInWeight requires a class perpendicular to the polarization",
            ));
        }
        let x = raw(class);
        match (&self.direction, weight) {
            (Direction::Central { .. }, i128::MIN..=0)
            | (Direction::Eichler { .. }, i128::MIN..=-1) => {
                Ok(x.iter().all(|entry| entry.is_zero()))
            }
            (Direction::Central { dual, .. }, 1) => {
                let mut residual = x;
                for (axis, covector) in self.cusp.oriented_plane.iter().zip(dual.iter()) {
                    let coefficient = dot(model, covector, &x)?.negate();
                    add_scaled(model, &mut residual, &raw(axis), coefficient)?;
                }
                Ok(residual.iter().all(|entry| entry.is_zero()))
            }
            (Direction::Central { .. }, 2) => {
                let [e, f] = &self.cusp.oriented_plane;
                Ok(class.wide_pair(model, e)?.is_zero() && class.wide_pair(model, f)?.is_zero())
            }
            (Direction::Eichler { index, .. }, 0 | 1) => {
                on_axis(model, &x, &raw(&self.cusp.oriented_plane[*index]))
            }
            (Direction::Eichler { index, .. }, 2 | 3) => Ok(class
                .wide_pair(model, &self.cusp.oriented_plane[*index])?
                .is_zero()),
            _ => Ok(true),
        }
    }

    pub(crate) fn report(&self, name: &str) -> Result<String> {
        let (label, graded) = match self.direction {
            Direction::Central { .. } => ("II", [0, 2, 3, 2, 0]),
            Direction::Eichler { .. } => ("III", [1, 0, 5, 0, 1]),
        };
        let smith_entries = |index: usize| {
            self.smith[index]
                .iter()
                .copied()
                .map(render_integer)
                .collect::<Vec<_>>()
                .join(", ")
        };
        let index = |order: usize| -> Result<String> {
            let product = self.smith[order]
                .iter()
                .try_fold(column_integer(1), |product, &entry| {
                    product.multiply_wide(entry)
                })?;
            Ok(render_integer(product))
        };
        Ok(format!(
            "rvOG6NilpotentData[{name}] = type={label}, H_rank=7, graded={graded:?}, imN_smith=[{}], imN2_smith=[{}], imN_index={}, imN2_index={}",
            smith_entries(0), smith_entries(1), index(0)?, index(1)?
        ))
    }
}

fn on_axis(model: &mut Model, x: &Column, axis: &Column) -> Result<bool> {
    let pivot = axis
        .iter()
        .position(|entry| !entry.is_zero())
        .expect("a cusp axis is nonzero");
    let Some(coefficient) = x[pivot].divide_exact(model, axis[pivot])? else {
        return Ok(false);
    };
    let mut residual = *x;
    add_scaled(model, &mut residual, axis, coefficient.negate())?;
    Ok(residual.iter().all(|entry| entry.is_zero()))
}

fn smith_image(model: &mut Model, columns: &[Column], rank: usize) -> Result<Vec<ColumnInteger>> {
    if rank == 0 {
        return Ok(Vec::new());
    }
    let one = Natural::from_small(1);
    let mut first = Natural::from_small(0);
    for entry in columns.iter().flatten() {
        first = first.gcd(model, entry.magnitude)?;
        if first == one {
            break;
        }
    }
    if first.is_zero() {
        return Err(Error::plain(
            "OG6 nilpotent image has smaller rank than expected",
        ));
    }
    let as_integer = |magnitude| ColumnInteger {
        magnitude,
        negative: false,
    };
    let mut invariants = vec![as_integer(first)];
    if rank == 1 {
        return Ok(invariants);
    }
    let mut second = Natural::from_small(0);
    'columns: for (j, right) in columns.iter().enumerate() {
        for left in &columns[..j] {
            for (i, (&a, &b)) in left.iter().zip(right).enumerate() {
                for (&c, &d) in left[..i].iter().zip(&right[..i]) {
                    model.spend(COLUMN_CELLS)?;
                    let determinant = a.multiply_wide(d)?.add(b.multiply_wide(c)?.negate())?;
                    second = second.gcd(model, determinant.magnitude)?;
                    if second == one {
                        break 'columns;
                    }
                }
            }
        }
    }
    if second.is_zero() {
        return Err(Error::plain(
            "OG6 nilpotent image has smaller rank than expected",
        ));
    }
    let second = second
        .divide_exact(model, first)?
        .expect("Smith determinantal divisors divide");
    invariants.push(as_integer(second));
    Ok(invariants)
}
