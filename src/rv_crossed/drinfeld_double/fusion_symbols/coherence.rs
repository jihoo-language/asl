use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::Model;

use super::matrix::{coordinates, Matrix};
use super::{Channel, Double, Tree, Work, A, AB, B, BC, C, LEFT};

const D: Tree = Tree::Leaf(3);
const CD: Tree = Tree::Join(&C, &D);
const A_BC: Tree = Tree::Join(&A, &BC);
const BC_D: Tree = Tree::Join(&BC, &D);
const B_CD: Tree = Tree::Join(&B, &CD);
const LEFT_FOUR: Tree = Tree::Join(&LEFT, &D);
const ROTATE_LEFT: Tree = Tree::Join(&A_BC, &D);
const ROTATE_CENTER: Tree = Tree::Join(&A, &BC_D);
const RIGHT_FOUR: Tree = Tree::Join(&A, &B_CD);
const MIDDLE: Tree = Tree::Join(&AB, &CD);

impl Double {
    pub(crate) fn pentagon(
        &self,
        model: &mut Model,
        name: &str,
        [a, b, c, d, total]: [usize; 5],
    ) -> Result<String> {
        let mut work = Work::new(model, self);
        let colors = [a, b, c, d];
        work.check(&[a, b, c, d, total])?;
        let [g, h, j, k] = colors.map(|i| self.simples[i].flux);
        let gh = self.group.mul(g, h);
        let hj = self.group.mul(h, j);
        let jk = self.group.mul(j, k);
        let l = work.basis(&LEFT_FOUR, &colors, total)?;
        let n = work.basis(&ROTATE_LEFT, &colors, total)?;
        let m = work.basis(&ROTATE_CENTER, &colors, total)?;
        let r = work.basis(&RIGHT_FOUR, &colors, total)?;
        let s = work.basis(&MIDDLE, &colors, total)?;
        let ln = work.move_by(&l, &n, self.treeshift([g, h, j], false))?;
        let nm = work.move_by(&n, &m, self.treeshift([g, hj, k], false))?;
        let mr = work.move_by(&m, &r, self.treeshift([h, j, k], false))?;
        let ls = work.move_by(&l, &s, self.treeshift([gh, j, k], false))?;
        let sr = work.move_by(&s, &r, self.treeshift([g, h, jk], false))?;
        let on_long_path = mr
            .multiply(&mut work.field, &nm)?
            .multiply(&mut work.field, &ln)?;
        let on_short_path = sr.multiply(&mut work.field, &ls)?;
        if on_long_path != on_short_path {
            return Err(Error::plain("double fusion associators fail the pentagon"));
        }
        Ok(format!(
            "rvDoublePentagon[{name},{a},{b},{c},{d},{total}] = dimension={}, verified=true",
            l.len(),
        ))
    }

    pub(crate) fn hexagon(
        &self,
        model: &mut Model,
        name: &str,
        [a, b, c, total]: [usize; 4],
    ) -> Result<String> {
        let mut work = Work::new(model, self);
        let colors = [a, b, c];
        work.check(&[a, b, c, total])?;
        let ab = work.crossing(colors, total, 0)?;
        let ac_after_ab = work.crossing([b, a, c], total, 1)?;
        let bc = work.crossing(colors, total, 1)?;
        let ac_after_bc = work.crossing([a, c, b], total, 0)?;
        let first = ac_after_ab.multiply(&mut work.field, &ab)?;
        let second = ac_after_bc.multiply(&mut work.field, &bc)?;
        if first != work.block_crossing(colors, total, false)?
            || second != work.block_crossing(colors, total, true)?
        {
            return Err(Error::plain("double fusion braidings fail the hexagon"));
        }
        Ok(format!(
            "rvDoubleHexagon[{name},{a},{b},{c},{total}] = dimension={}, verified=true",
            first.cols,
        ))
    }
}

impl Work<'_> {
    fn crossing(&mut self, colors: [usize; 3], total: usize, slot: usize) -> Result<Matrix> {
        let mut output = colors;
        output.swap(slot, slot + 1);
        let (left, right) = self.triple(colors, total)?;
        let (next_left, next_right) = self.triple(output, total)?;
        let actual = if slot == 0 {
            self.local_r(&left, &next_left, [colors[0], colors[1]])?
        } else {
            let f = self.f_move(colors, &left, &right, false)?;
            let r = self.local_r(&right, &next_right, [colors[1], colors[2]])?;
            let finv = self.f_move(output, &next_left, &next_right, true)?;
            finv.multiply(&mut self.field, &r)?
                .multiply(&mut self.field, &f)?
        };
        let (_, braid) = self
            .double
            .braid(self.field.model, colors, &[(slot + 1) as i128])?;
        let expected = self.transformed_coordinates(&left, &next_left, |i| {
            let image = braid.0[i];
            (image.target, image.phase)
        })?;
        if actual != expected {
            return Err(Error::plain(
                "local F and R matrices disagree with the induced braid",
            ));
        }
        Ok(actual)
    }

    fn local_r(
        &mut self,
        source: &[Channel],
        target: &[Channel],
        colors: [usize; 2],
    ) -> Result<Matrix> {
        let mut matrix = Matrix::zero(target.len(), source.len())?;
        let mut symbols = BTreeMap::new();
        for (i, a) in source.iter().enumerate() {
            let (e, mu, nu) = a.pair();
            let r = match symbols.entry(e) {
                std::collections::btree_map::Entry::Occupied(entry) => entry.into_mut(),
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(self.r_move([colors[0], colors[1], e], false)?)
                }
            };
            for (j, b) in target.iter().enumerate() {
                let (f, lam, gam) = b.pair();
                if e == f && nu == gam {
                    matrix.set(j, i, r.entry([lam, mu])?);
                }
            }
        }
        Ok(matrix)
    }

    fn transformed_coordinates(
        &mut self,
        source: &[Channel],
        target: &[Channel],
        image: impl Fn(usize) -> (usize, usize),
    ) -> Result<Matrix> {
        let basis = target.iter().map(|m| &m.map).collect::<Vec<_>>();
        let changed = source
            .iter()
            .map(|m| m.map.transformed(&mut self.field, &image))
            .collect::<Result<Vec<_>>>()?;
        coordinates(&mut self.field, &basis, &changed.iter().collect::<Vec<_>>())
    }

    fn block_crossing(
        &mut self,
        colors: [usize; 3],
        total: usize,
        left_block: bool,
    ) -> Result<Matrix> {
        let [a, b, c] = colors;
        let d = self.double;
        let [g, h, j] = colors.map(|i| d.simples[i].flux);
        let [da, db, dc] = colors.map(|i| d.simples[i].dimension);
        let source = self.basis(&LEFT, &colors, total)?;
        if left_block {
            let target = self.basis(&LEFT, &[c, a, b], total)?;
            let combined = d.group.mul(g, h);
            let phase = d.treeshift([j, g, h], true);
            let action = &d.representations[c].action[combined];
            self.transformed_coordinates(&source, &target, |row| {
                let (ab, c) = (row / dc, row % dc);
                let p = action[c];
                (p.target * da * db + ab, (phase + p.phase) % d.conductor)
            })
        } else {
            let target = self.basis(&LEFT, &[b, c, a], total)?;
            let phase = d.treeshift([g, h, j], false) + d.tensor_phase(g, h, j);
            let (act_b, act_c) = (
                &d.representations[b].action[g],
                &d.representations[c].action[g],
            );
            self.transformed_coordinates(&source, &target, |row| {
                let (a, b, c) = (row / (db * dc), (row / dc) % db, row % dc);
                let (q, p) = (act_b[b], act_c[c]);
                (
                    (q.target * dc + p.target) * da + a,
                    (phase + q.phase + p.phase) % d.conductor,
                )
            })
        }
    }
}
