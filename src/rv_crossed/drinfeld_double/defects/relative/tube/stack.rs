use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::rv_lens::Rational;

use super::super::linear::{self, GradedMap};
use super::super::tree::{self, Tree};
use super::super::Module;
use super::build::Builder;
use super::{Element, Form};

impl Builder<'_> {
    pub(super) fn check_stack(
        &mut self,
        f: Form,
        g: Form,
        output: &Element,
        forms: &[Form],
    ) -> Result<()> {
        if f.i != g.j {
            return if output.0.is_empty() {
                Ok(())
            } else {
                Err(Error::plain(
                    "tubes with distinct glued boundaries have a nonzero product",
                ))
            };
        }
        let actual = self.stack(f, g)?;
        let mut coefficients = BTreeMap::<usize, Vec<(Form, &Rational)>>::new();
        for &(index, ref value) in &output.0 {
            let form = forms[index];
            if form.j != f.j || form.i != g.i {
                return Err(Error::plain(
                    "tube stacking changed its outer annular boundaries",
                ));
            }
            coefficients.entry(form.k).or_default().push((form, value));
        }
        for (k, built) in actual {
            let left = self.product([f.j, k])?;
            let right = self.product([k, g.i])?;
            let mut given = GradedMap::zero(&left.balanced.module, &right.balanced.module)?;
            for (form, coeff) in coefficients.remove(&k).into_iter().flatten() {
                let mut tensor = self.form_map(form)?;
                tensor.scale(&mut self.field, coeff)?;
                given.add_assign(&mut self.field, &tensor)?;
            }
            if given.grades != built.grades {
                return Err(Error::plain(
                    "exact tube product disagrees with direct annular stacking",
                ));
            }
        }
        if !coefficients.is_empty() {
            return Err(Error::plain(
                "tube product contains an unsupported fused-loop sector",
            ));
        }
        Ok(())
    }

    fn stack(&mut self, f: Form, g: Form) -> Result<BTreeMap<usize, GradedMap>> {
        let j = self.objects[f.j].clone();
        let k = self.objects[f.k].clone();
        let i = self.objects[f.i].clone();
        let l = self.objects[g.k].clone();
        let h = self.objects[g.i].clone();
        let a = [j.as_ref(), k.as_ref(), l.as_ref()];
        let b = [k.as_ref(), i.as_ref(), l.as_ref()];
        let c = [k.as_ref(), l.as_ref(), h.as_ref()];
        let (a_left, a_right, alpha_a) = self.assoc(&a)?;
        let (b_left, b_right, alpha_b) = self.assoc(&b)?;
        let (c_left, c_right, alpha_c) = self.assoc(&c)?;
        let fm = self.form_map(f)?;
        let gm = self.form_map(g)?;
        let f_ext = linear::tensor(
            &mut self.field,
            &self.ambient,
            a_left.layout(),
            b_left.layout(),
            [&fm, &self.identities[g.k]],
        )?;
        let g_ext = linear::tensor(
            &mut self.field,
            &self.ambient,
            b_right.layout(),
            c_right.layout(),
            [&self.identities[f.k], &gm],
        )?;
        // j(kl) -> (jk)l -> (ki)l -> k(il) -> k(lh) -> (kl)h.
        let alpha_a_inv = alpha_a.adjoint(&mut self.field)?;
        let alpha_c_inv = alpha_c.adjoint(&mut self.field)?;
        let first = f_ext.compose(&mut self.field, &alpha_a_inv)?;
        let first = alpha_b.compose(&mut self.field, &first)?;
        let first = g_ext.compose(&mut self.field, &first)?;
        let stacked = alpha_c_inv.compose(&mut self.field, &first)?;
        let mut result = BTreeMap::new();
        for m in self.loops(f.k, g.k)? {
            let source = self.product([f.j, m])?;
            let target = self.product([m, g.i])?;
            let channels = self.channels([f.k, g.k, m])?;
            let mut term = GradedMap::zero(&source.balanced.module, &target.balanced.module)?;
            for channel in channels.iter() {
                let pre = linear::tensor(
                    &mut self.field,
                    &self.ambient,
                    source.layout(),
                    a_right.layout(),
                    [&self.identities[f.j], &channel.inclusion],
                )?;
                let post = linear::tensor(
                    &mut self.field,
                    &self.ambient,
                    c_left.layout(),
                    target.layout(),
                    [&channel.projection, &self.identities[g.i]],
                )?;
                let mid = stacked.compose(&mut self.field, &pre)?;
                let part = post.compose(&mut self.field, &mid)?;
                term.add_assign(&mut self.field, &part)?;
            }
            result.insert(m, term);
        }
        Ok(result)
    }

    fn assoc<'a>(&mut self, inputs: &'a [&Module; 3]) -> Result<(Tree<'a>, Tree<'a>, GradedMap)> {
        let [left, right] = tree::TRIPLE
            .each_ref()
            .map(|s| Tree::build(self.field.model, &self.ambient, inputs, s));
        let (left, right) = (left?, right?);
        let map = tree::rotation(
            self.field.model,
            &self.ambient,
            inputs,
            &left,
            &right,
            [0, 1, 2, 3],
        )?;
        let matrix = GradedMap::monomial(&mut self.field, &map)?;
        Ok((left, right, matrix))
    }
}
