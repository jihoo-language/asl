use std::rc::Rc;

use super::super::WEDGE_SIGNS;
use super::{abelian_range, coordinate_range, Class, Integer, Lattice, Vector, RANK};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::mukai::moduli::{
    column_integer, eliminate, ColumnInteger, COLUMN_CELLS,
};
use crate::hyperelliptic_kummer_surface::rational_sparse::decimal;
use crate::model::Model;
use crate::syntax::{Label, ParsedExpr};
use crate::wide_integer::Natural;

const RAW_RANK: usize = RANK + 1;
const CUSP_RANK: usize = 3;
type Column = [ColumnInteger; RAW_RANK];

pub(crate) enum Expression {
    Basis(Label, usize),
    Transvection(Label, usize, Box<[ParsedExpr; 2]>),
    Central(Label, i128, Box<ParsedExpr>),
}

pub(crate) enum Command {
    Define(Label, Box<[ParsedExpr; 3]>),
    Gram(Label),
}

pub(crate) struct Cusp {
    pub(crate) lattice_name: Rc<str>,
    polarization: Rc<Class>,
    oriented_plane: [Rc<Class>; 2],
    lifts: Vec<Column>,
}

impl Cusp {
    pub(crate) fn new(
        model: &mut Model,
        lattice_name: Rc<str>,
        lattice: &Lattice,
        [polarization, e, f]: [Rc<Class>; 3],
    ) -> Result<Self> {
        let h_square = polarization.wide_pair(model, &polarization)?;
        let h_raw = raw(&polarization);
        let h_content = h_raw
            .iter()
            .try_fold(Natural::from_small(0), |gcd, entry| {
                gcd.gcd(model, entry.magnitude)
            })?;
        if h_square.negative || h_square.is_zero() || h_content.small() != Some(1) {
            return Err(Error::plain(
                "rvog6cusp requires a primitive positive polarization",
            ));
        }
        for (left, right) in [
            (&e, &e),
            (&e, &f),
            (&f, &f),
            (&polarization, &e),
            (&polarization, &f),
        ] {
            if !left.wide_pair(model, right)?.is_zero() {
                return Err(Error::plain(
                    "rvog6cusp requires mutually orthogonal isotropic generators perpendicular to the polarization",
                ));
            }
        }
        let e_raw = raw(&e);
        let f_raw = raw(&f);
        require_primitive_plane(model, &e_raw, &f_raw)?;

        let mut lifts = lattice
            .columns
            .iter()
            .map(|column| {
                let mut lift = [column_integer(0); RAW_RANK];
                lift[..RANK].copy_from_slice(column);
                lift
            })
            .collect::<Vec<_>>();
        let mut exceptional = [column_integer(0); RAW_RANK];
        exceptional[RANK] = column_integer(1);
        lifts.push(exceptional);

        for vector in [&h_raw, &e_raw, &f_raw] {
            let row = pairing_row(vector)?;
            cut_kernel(model, &mut lifts, &row, false)?;
        }

        // Integral coordinate functionals dual to the primitive plane split
        // I off from the rank-five kernel. Their common kernel gives lifts
        // of a basis of the quotient, rather than an index sublattice of it.
        let [e_dual, f_dual] = dual_plane(model, &e_raw, &f_raw)?;
        for dual in [&e_dual, &f_dual] {
            cut_kernel(model, &mut lifts, dual, true)?;
        }
        debug_assert_eq!(lifts.len(), CUSP_RANK);

        Ok(Self {
            lattice_name,
            polarization,
            oriented_plane: [e, f],
            lifts,
        })
    }

    pub(crate) fn cells(&self) -> usize {
        5 + CUSP_RANK * RAW_RANK * COLUMN_CELLS
    }

    pub(crate) fn basis(&self, index: usize) -> Result<Class> {
        let lift = self
            .lifts
            .get(index)
            .ok_or_else(|| Error::plain("rvOG6CuspBasis index must be between 0 and 2"))?;
        realize(lift)
    }

    pub(crate) fn gram(&self, model: &mut Model, name: &str) -> Result<String> {
        let mut gram = [[column_integer(0); CUSP_RANK]; CUSP_RANK];
        for (i, lift) in self.lifts.iter().enumerate() {
            let row = pairing_row(lift)?;
            for (j, other) in self.lifts[..=i].iter().enumerate() {
                let entry = dot(model, &row, other)?;
                gram[i][j] = entry;
                gram[j][i] = entry;
            }
        }
        let [[a, b, c], [_, d, e], [_, _, f]] = gram;
        let minor0 = d.multiply_wide(f)?.add(e.multiply_wide(e)?.negate())?;
        let minor1 = b.multiply_wide(f)?.add(c.multiply_wide(e)?.negate())?;
        let minor2 = b.multiply_wide(e)?.add(c.multiply_wide(d)?.negate())?;
        let determinant = a
            .multiply_wide(minor0)?
            .add(b.multiply_wide(minor1)?.negate())?
            .add(c.multiply_wide(minor2)?)?;
        let rows: Vec<_> = gram
            .into_iter()
            .map(|row| {
                let entries = row.map(render_integer);
                format!("[{}]", entries.join(", "))
            })
            .collect();
        Ok(format!(
            "rvOG6CuspGram[{name}] = rank=3, signature=(0,3), determinant={}, gram=[{}]",
            render_integer(determinant),
            rows.join(", ")
        ))
    }

    pub(crate) fn transvection(
        &self,
        model: &mut Model,
        index: usize,
        u: &Class,
        x: &Class,
    ) -> Result<Class> {
        let a = self
            .oriented_plane
            .get(index)
            .ok_or_else(|| Error::plain("rvOG6CuspTransvection axis must be 0 or 1"))?;
        for axis in [
            &self.polarization,
            &self.oriented_plane[0],
            &self.oriented_plane[1],
        ] {
            if !u.wide_pair(model, axis)?.is_zero() {
                return Err(Error::plain(
                    "rvOG6CuspTransvection requires a direction perpendicular to the polarization and the cusp plane",
                ));
            }
        }
        let ax = a.wide_pair(model, x)?.resize()?;
        let ux: ColumnInteger = u.wide_pair(model, x)?.resize()?;
        let uu = u.wide_pair(model, u)?;
        let (half, remainder) = uu.magnitude.divide_u32(2);
        debug_assert_eq!(remainder, 0, "OG6 is even");
        let half: ColumnInteger = Integer {
            magnitude: half,
            negative: uu.negative,
        }
        .resize()?;
        let shift = ux.add(half.multiply_wide(ax)?)?.negate();
        let mut result = raw(x);
        add_scaled(model, &mut result, &raw(u), ax)?;
        add_scaled(model, &mut result, &raw(a), shift)?;
        realize(&result)
    }

    pub(crate) fn central(&self, model: &mut Model, k: i128, x: &Class) -> Result<Class> {
        let [e, f] = &self.oriented_plane;
        let ex = e.wide_pair(model, x)?.resize()?.multiply(k)?;
        let fx = f.wide_pair(model, x)?.resize()?.multiply(k)?.negate();
        let mut result = raw(x);
        add_scaled(model, &mut result, &raw(f), ex)?;
        add_scaled(model, &mut result, &raw(e), fx)?;
        realize(&result)
    }
}

fn require_primitive_plane(model: &mut Model, e: &Column, f: &Column) -> Result<()> {
    // w^perp ⊕ Ze is saturated in the nine raw integral coordinates.
    // The gcd of Plücker minors therefore measures this plane's index.
    let mut divisor = Natural::from_small(0);
    for (j, (&ej, &fj)) in e.iter().zip(f).enumerate() {
        for (&ei, &fi) in e[..j].iter().zip(&f[..j]) {
            model.spend(COLUMN_CELLS)?;
            let minor = ei.multiply_wide(fj)?.add(ej.multiply_wide(fi)?.negate())?;
            divisor = divisor.gcd(model, minor.magnitude)?;
        }
    }
    match divisor.small() {
        Some(1) => Ok(()),
        Some(0) => Err(Error::plain(
            "rvog6cusp requires two independent isotropic generators",
        )),
        _ => Err(Error::plain(
            "rvog6cusp requires a primitive isotropic plane",
        )),
    }
}

fn dual_plane(model: &mut Model, e: &Column, f: &Column) -> Result<[Column; 2]> {
    let mut forms = vec![[column_integer(0); RAW_RANK]; RAW_RANK];
    for (i, column) in forms.iter_mut().enumerate() {
        column[i] = column_integer(1);
    }
    let message = "rvog6cusp cannot split the primitive isotropic plane";
    let mut e_dual = eliminate(model, &mut forms, e.to_vec(), message)?;
    let row = forms
        .iter()
        .map(|form| dot(model, form, f))
        .collect::<Result<Vec<_>>>()?;
    let f_dual = eliminate(model, &mut forms, row, message)?;
    let correction = dot(model, &e_dual, f)?.negate();
    add_scaled(model, &mut e_dual, &f_dual, correction)?;
    Ok([e_dual, f_dual])
}

fn cut_kernel(
    model: &mut Model,
    columns: &mut Vec<Column>,
    covector: &Column,
    primitive: bool,
) -> Result<()> {
    let mut row = columns
        .iter()
        .map(|column| dot(model, covector, column))
        .collect::<Result<Vec<_>>>()?;
    if !primitive {
        // Dividing the covector by its full content preserves its
        // integral kernel even when its image has index greater than one.
        let divisor = row.iter().try_fold(Natural::from_small(0), |gcd, entry| {
            gcd.gcd(model, entry.magnitude)
        })?;
        if divisor.is_zero() {
            return Err(Error::plain("rvog6cusp has dependent kernel conditions"));
        }
        if divisor.small() != Some(1) {
            let divisor = ColumnInteger {
                magnitude: divisor,
                negative: false,
            };
            for entry in &mut row {
                *entry = entry
                    .divide_exact(model, divisor)?
                    .expect("gcd divides the covector");
            }
        }
    }
    eliminate(
        model,
        columns,
        row,
        "rvog6cusp cannot split the integral kernel",
    )?;
    Ok(())
}

fn raw(class: &Class) -> Column {
    let mut result = [column_integer(0); RAW_RANK];
    for (out, &coordinate) in result.iter_mut().zip(&class.perpendicular.coordinates) {
        *out = column_integer(coordinate);
    }
    result[RANK] = column_integer(class.exceptional);
    result
}

fn pairing_row(vector: &Column) -> Result<Column> {
    let mut row = [column_integer(0); RAW_RANK];
    row[0] = vector[RANK - 1].negate();
    row[RANK - 1] = vector[0].negate();
    for (i, &sign) in WEDGE_SIGNS.iter().enumerate() {
        row[i + 1] = vector[RANK - 2 - i].multiply(sign)?;
    }
    row[RANK] = vector[RANK].multiply(-2)?;
    Ok(row)
}

fn dot(model: &mut Model, row: &Column, column: &Column) -> Result<ColumnInteger> {
    model.spend(RAW_RANK * COLUMN_CELLS)?;
    let mut sum = column_integer(0);
    for (&a, &b) in row.iter().zip(column) {
        if !a.is_zero() && !b.is_zero() {
            sum = sum.add(a.multiply_wide(b)?)?;
        }
    }
    Ok(sum)
}

fn add_scaled(
    model: &mut Model,
    out: &mut Column,
    axis: &Column,
    coefficient: ColumnInteger,
) -> Result<()> {
    model.spend(RAW_RANK * COLUMN_CELLS)?;
    for (entry, &axis) in out.iter_mut().zip(axis) {
        if !axis.is_zero() && !coefficient.is_zero() {
            *entry = entry.add(axis.multiply_wide(coefficient)?)?;
        }
    }
    Ok(())
}

fn realize(column: &Column) -> Result<Class> {
    let mut vector = Vector::zero();
    for (out, &coordinate) in vector.coordinates.iter_mut().zip(column) {
        *out = coordinate.signed().ok_or_else(abelian_range)?;
    }
    let exceptional = column[RANK].signed().ok_or_else(coordinate_range)?;
    Ok(Class {
        perpendicular: Rc::new(vector),
        exceptional,
    })
}

fn render_integer(value: ColumnInteger) -> String {
    let sign = if value.negative && !value.is_zero() {
        "-"
    } else {
        ""
    };
    format!("{sign}{}", decimal(value.magnitude))
}
