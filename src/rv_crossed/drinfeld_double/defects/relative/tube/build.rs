use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_lens::Rational;

use super::super::central::{Bulk, Central};
use super::super::linear::{self, ChannelMap, GradedMap};
use super::super::tensor::Layout;
use super::super::{decompose, Ambient, Balanced, Category, Double, Hom, Image, Module};
use super::matrix::{Field, Matrix};
use super::{Block, Form, MAX_TUBE};

pub(super) struct Builder<'a> {
    pub(super) ambient: Ambient<'a>,
    pub(super) field: Field<'a>,
    pub(super) category: &'a Category,
    pub(super) objects: Vec<Rc<Module>>,
    pub(super) identities: Vec<GradedMap>,
    products: BTreeMap<[usize; 2], Rc<Product>>,
    channels: BTreeMap<[usize; 3], Rc<[ChannelMap]>>,
}

pub(super) struct Product {
    pub(super) balanced: Balanced,
    factors: [Rc<Module>; 2],
    homs: Vec<Hom>,
}

impl Product {
    pub(super) fn layout(&self) -> Layout<'_> {
        Layout {
            product: &self.balanced,
            factors: self.factors.each_ref().map(Rc::as_ref),
        }
    }
}

impl<'a> Builder<'a> {
    pub(super) fn new(
        model: &'a mut Model,
        double: &'a Double,
        category: &'a Category,
    ) -> Result<Self> {
        super::super::central::self_boundary(category)?;
        if category.count() > MAX_TUBE {
            return Err(super::bound());
        }
        let ambient = Ambient::new(model, double, &[category], [])?;
        let objects = (0..category.count())
            .map(|i| ambient.induce(model, category, i).map(Rc::new))
            .collect::<Result<Vec<_>>>()?;
        let identities = objects
            .iter()
            .map(|m| {
                model.spend(m.fibers.iter().map(|&n| n * n).sum())?;
                GradedMap::identity(m)
            })
            .collect::<Result<Vec<_>>>()?;
        let field = Field::new(model, ambient.conductor);
        Ok(Self {
            ambient,
            field,
            category,
            objects,
            identities,
            products: BTreeMap::new(),
            channels: BTreeMap::new(),
        })
    }

    pub(super) fn product(&mut self, pair @ [a, b]: [usize; 2]) -> Result<Rc<Product>> {
        if let Some(product) = self.products.get(&pair) {
            return Ok(product.clone());
        }
        let factors = [self.objects[a].clone(), self.objects[b].clone()];
        let balanced = Balanced::new(self.field.model, &self.ambient, &factors[0], &factors[1])?;
        let homs = decompose(
            self.field.model,
            &self.ambient,
            self.category,
            &balanced.module,
        )?;
        let product = Rc::new(Product {
            balanced,
            factors,
            homs,
        });
        self.products.insert(pair, product.clone());
        Ok(product)
    }

    pub(super) fn channels(&mut self, key @ [a, b, c]: [usize; 3]) -> Result<Rc<[ChannelMap]>> {
        if let Some(result) = self.channels.get(&key) {
            return Ok(result.clone());
        }
        let product = self.product([a, b])?;
        let channels = linear::hom_basis(
            &mut self.field,
            &self.ambient,
            &self.objects[c],
            &product.balanced.module,
            self.category.simple(c)?.simple.flux,
            &product.homs[c],
        )?;
        let result: Rc<[ChannelMap]> = channels.into();
        self.channels.insert(key, result.clone());
        Ok(result)
    }

    pub(super) fn loops(&mut self, a: usize, b: usize) -> Result<Vec<usize>> {
        Ok(self
            .product([a, b])?
            .homs
            .iter()
            .enumerate()
            .filter_map(|(i, x)| (!x.orbits.is_empty()).then_some(i))
            .collect())
    }

    pub(super) fn forms(&mut self) -> Result<Vec<Form>> {
        let n = self.objects.len();
        let mut forms = Vec::new();
        for j in 0..n {
            for i in 0..n {
                for k in 0..n {
                    let a = self.product([j, k])?;
                    let b = self.product([k, i])?;
                    for (c, (x, y)) in a.homs.iter().zip(&b.homs).enumerate() {
                        for u in 0..x.orbits.len() {
                            for v in 0..y.orbits.len() {
                                if forms.len() >= MAX_TUBE {
                                    return Err(super::bound());
                                }
                                forms.push(Form { j, k, i, c, u, v });
                            }
                        }
                    }
                }
            }
        }
        Ok(forms)
    }

    pub(super) fn form_map(&mut self, f: Form) -> Result<GradedMap> {
        let a = self.channels([f.j, f.k, f.c])?;
        let b = self.channels([f.k, f.i, f.c])?;
        b[f.v]
            .inclusion
            .compose(&mut self.field, &a[f.u].projection)
    }

    pub(super) fn fourier(&mut self, forms: &[Form]) -> Result<(Vec<Block>, Matrix, usize)> {
        let mut fourier = Matrix::zero(forms.len(), forms.len())?;
        let mut blocks = Vec::new();
        let mut start = 0;
        let mut unit = None;
        let maps = forms
            .iter()
            .map(|&f| self.form_map(f))
            .collect::<Result<Vec<_>>>()?;
        for a in 0..self.ambient.double.count() {
            let bulk = Bulk::new(self.field.model, &self.ambient, self.category, a)?;
            let homs = decompose(self.field.model, &self.ambient, self.category, &bulk.module)?;
            let mut branches = Vec::new();
            let mut slots = Vec::new();
            for (i, h) in homs.iter().enumerate() {
                let f = linear::hom_basis(
                    &mut self.field,
                    &self.ambient,
                    &self.objects[i],
                    &bulk.module,
                    self.category.simple(i)?.simple.flux,
                    h,
                )?;
                slots.extend((0..f.len()).map(|s| [i, s]));
                branches.extend(f);
            }
            if a == 0 {
                let &[[identity, 0]] = slots.as_slice() else {
                    return Err(Error::plain(
                        "bulk vacuum does not branch to one simple annular unit",
                    ));
                };
                unit = Some(identity);
            }
            let block = Block { start, slots };
            let d = block.size();
            start += d * d;
            if d == 0 || start > forms.len() {
                return Err(Error::plain(
                    "tube annular basis and bulk block dimensions disagree",
                ));
            }
            self.representation(&mut fourier, forms, &maps, &bulk, &branches, &block)?;
            blocks.push(block);
        }
        if start != forms.len() {
            return Err(Error::plain(
                "tube annular basis and squared bulk block dimensions disagree",
            ));
        }
        Ok((blocks, fourier, unit.expect("bulk vacuum exists")))
    }

    fn representation(
        &mut self,
        fourier: &mut Matrix,
        forms: &[Form],
        maps: &[GradedMap],
        bulk: &Bulk,
        branches: &[ChannelMap],
        block: &Block,
    ) -> Result<()> {
        for (k, object) in self.objects.clone().iter().enumerate() {
            let central = Central::new(self.field.model, &self.ambient, [bulk, bulk], object)?;
            let [source, target] = &central.products;
            let bulk_first = Layout {
                product: source,
                factors: [&bulk.module, object],
            };
            let bulk_second = Layout {
                product: target,
                factors: [object, &bulk.module],
            };
            let mut outgoing = Vec::new();
            let mut incoming = Vec::new();
            for (idx, &[i, _]) in block.slots.iter().enumerate() {
                let left = self.product([i, k])?;
                let right = self.product([k, i])?;
                let h = &branches[idx].inclusion;
                let id = &self.identities[k];
                outgoing.push(linear::tensor(
                    &mut self.field,
                    &self.ambient,
                    left.layout(),
                    bulk_first,
                    [h, id],
                )?);
                incoming.push(
                    linear::tensor(
                        &mut self.field,
                        &self.ambient,
                        right.layout(),
                        bulk_second,
                        [id, h],
                    )?
                    .after_monomial(&mut self.field, &central.maps[1].0)?,
                );
            }
            for (f, (descriptor, map)) in forms
                .iter()
                .zip(maps)
                .enumerate()
                .filter(|(_, (f, _))| f.k == k)
            {
                for (c, &[i, _]) in block
                    .slots
                    .iter()
                    .enumerate()
                    .filter(|(_, slot)| slot[0] == descriptor.i)
                {
                    let _ = i;
                    let composed = incoming[c].compose(&mut self.field, map)?;
                    for (r, &[j, _]) in block
                        .slots
                        .iter()
                        .enumerate()
                        .filter(|(_, slot)| slot[0] == descriptor.j)
                    {
                        let _ = j;
                        // Categorical trace on a self-bimodule is its raw
                        // graded trace divided by |H|. The common factor
                        // cancels against the output channel's Gram norm.
                        let term = outgoing[r].inner(&mut self.field, &composed)?;
                        let value = self.field.multiply(&term, &branches[r].inverse_norm)?;
                        fourier.set(block.index(r, c), f, value);
                    }
                }
            }
        }
        Ok(())
    }

    pub(super) fn geometric_unit(&mut self, forms: &[Form], unit: usize) -> Result<Vec<Rational>> {
        let mut answer = vec![Rational::zero(); forms.len()];
        for i in 0..self.objects.len() {
            let identity = Form {
                i,
                j: i,
                k: unit,
                c: i,
                u: 0,
                v: 0,
            };
            let column = forms.iter().position(|f| *f == identity).ok_or_else(|| {
                Error::plain("annular basis lacks the unitor through the neutral loop")
            })?;
            let a = self.product([i, unit])?;
            let b = self.product([unit, i])?;
            let object = &self.objects[i];
            let left = unitor(&self.ambient, b.layout(), object, 0);
            let right = unitor(&self.ambient, a.layout(), object, 1);
            let inverse = |row: &[Image], y: Image| {
                let (target, x) = row
                    .iter()
                    .enumerate()
                    .find(|(_, x)| x.target == y.target)
                    .expect("left unitor acts bijectively");
                Image {
                    target,
                    phase: self.ambient.sum([y.phase, self.ambient.neg(x.phase)]),
                }
            };
            let rows = left
                .iter()
                .zip(&right)
                .map(|(l, r)| r.iter().map(|&x| inverse(l, x)).collect())
                .collect::<Vec<_>>();
            let raw = GradedMap::monomial(&mut self.field, &rows)?;
            let basis = self.form_map(identity)?;
            let norm = basis.inner(&mut self.field, &basis)?;
            let inverse = self.field.inverse(&norm)?;
            let coefficient = basis.inner(&mut self.field, &raw)?;
            let coefficient = self.field.multiply(&coefficient, &inverse)?;
            let mut reconstructed = basis;
            reconstructed.scale(&mut self.field, &coefficient)?;
            if reconstructed.grades != raw.grades {
                return Err(Error::plain(
                    "annular neutral loop fails the geometric unitor comparison",
                ));
            }
            answer[column] = coefficient;
        }
        Ok(answer)
    }
}

fn unitor(
    ambient: &Ambient<'_>,
    layout: Layout<'_>,
    object: &Module,
    side: usize,
) -> Vec<Vec<Image>> {
    layout
        .product
        .bases
        .iter()
        .map(|basis| {
            basis
                .iter()
                .map(|p| {
                    object.grade_action(
                        ambient,
                        side,
                        p.grades[side],
                        p.grades[1 - side],
                        p.columns[1 - side],
                    )
                })
                .collect()
        })
        .collect()
}
