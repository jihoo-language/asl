use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::sheaf::{
    self, Cochain, Cochains, Cohomology, Generator, Layer, Matrix, Sheaf, Sink, Term,
};
use crate::sheaf_page::{Grade, Page};
use crate::syntax::Label;

#[derive(Clone)]
pub(crate) struct NamedSheaf {
    pub(crate) name: Rc<str>,
    pub(crate) sheaf: Rc<Sheaf>,
}

impl NamedSheaf {
    pub(crate) fn resolve(model: &Model, name: &Label) -> Result<Self> {
        Ok(Self {
            name: name.text.as_str().into(),
            sheaf: model.sheaf(name)?,
        })
    }
}

pub(crate) struct Map {
    pub(crate) source: NamedSheaf,
    pub(crate) target: NamedSheaf,
    matrices: Vec<Rc<Matrix>>,
}

impl Map {
    pub(crate) fn build(
        model: &mut Model,
        source: NamedSheaf,
        target: NamedSheaf,
        entries: Vec<(Label, Vec<i64>)>,
    ) -> Result<Self> {
        if source.sheaf.base() != target.sheaf.base() || source.sheaf.prime != target.sheaf.prime {
            return Err(Error::plain(
                "shmap endpoints must have the same named base and prime field",
            ));
        }
        let x = &source.sheaf;
        let y = &target.sheaf;
        let mut matrices: Vec<Option<Rc<Matrix>>> = (0..x.dims().len()).map(|_| None).collect();
        let mut cells = matrices.len();
        for (point, values) in entries {
            let index = model.point(x.base(), &point)?;
            if matrices[index].is_some() {
                return Err(Error::new(
                    point.span,
                    format!("duplicate sheaf map component at '{}'", point.text),
                ));
            }
            let rows = y.dims()[index];
            let width = x.dims()[index];
            let expected = rows
                .checked_mul(width)
                .ok_or_else(|| Error::new(point.span, "sheaf map matrix size is out of range"))?;
            if values.len() != expected {
                return Err(Error::new(
                    point.span,
                    format!(
                        "sheaf map component at '{}' requires {rows} rows and {width} columns ({expected} entries), found {}",
                        point.text, values.len()
                    ),
                ));
            }
            let matrix = Matrix::literal(model, x.prime, width, values)?;
            cells += matrix.cells();
            sheaf::check_cells(cells)?;
            matrices[index] = Some(Rc::new(matrix));
        }
        let matrices: Vec<_> = matrices
            .into_iter()
            .enumerate()
            .map(|(point, matrix)| {
                matrix.ok_or_else(|| {
                    Error::plain(format!(
                        "missing sheaf map component at '{}'",
                        x.poset.names[point]
                    ))
                })
            })
            .collect::<Result<_>>()?;
        Self::from_matrices(model, source, target, matrices)
    }

    pub(crate) fn from_matrices(
        model: &mut Model,
        source: NamedSheaf,
        target: NamedSheaf,
        matrices: Vec<Rc<Matrix>>,
    ) -> Result<Self> {
        let x = &source.sheaf;
        let y = &target.sheaf;
        debug_assert!(x.base() == y.base());
        debug_assert_eq!(x.prime, y.prime);
        debug_assert_eq!(x.dims().len(), matrices.len());
        let cells = matrices.len() + matrices.iter().map(|matrix| matrix.cells()).sum::<usize>();
        sheaf::check_cells(cells)?;
        for (a, targets) in x.covers().iter().enumerate() {
            for &b in targets {
                model.spend(1)?;
                let first = Matrix::compose(model, x.prime, y.restriction(a, b), &matrices[a])?;
                let second = Matrix::compose(model, x.prime, &matrices[b], x.restriction(a, b))?;
                model.spend(first.cells() + second.cells())?;
                if first != second {
                    return Err(Error::plain(format!(
                        "sheaf map is not natural on cover '{}<{}'",
                        x.poset.names[a], x.poset.names[b]
                    )));
                }
            }
        }
        model.allocate(cells)?;
        Ok(Self {
            source,
            target,
            matrices,
        })
    }

    pub(crate) fn matrix(&self, point: usize) -> &Matrix {
        &self.matrices[point]
    }

    pub(crate) fn apply_cochain(&self, model: &mut Model, raw: &[Term]) -> Result<Value> {
        let output = self.apply_support(model, raw)?;
        sheaf::persist(model, output)
    }

    pub(crate) fn apply_class(
        &self,
        model: &mut Model,
        degree: usize,
        class: &[crate::prime_linear::Coordinate],
    ) -> Result<Value> {
        let source = &self.source.sheaf;
        let target = &self.target.sheaf;
        let raw =
            source
                .cohomology
                .representative_support(model, source.as_ref(), degree, class)?;
        let mapped = self.apply_support(model, &raw)?;
        target
            .cohomology
            .class(model, target.as_ref(), degree, &mapped)
    }

    fn apply_support(&self, model: &mut Model, raw: &[Term]) -> Result<Vec<Term>> {
        let mut output = Cochain::default();
        let prime = self.source.sheaf.prime;
        for (generator, coefficient) in raw {
            self.apply_term(
                model,
                generator,
                *coefficient,
                generator.vertical,
                &mut |model, target, c| output.add(model, prime, target, c),
            )?;
        }
        Ok(output.into_terms())
    }

    fn apply_term(
        &self,
        model: &mut Model,
        generator: &Generator,
        coefficient: u32,
        vertical: usize,
        emit: &mut Sink<'_>,
    ) -> Result<()> {
        let stalk = *generator.chain.last().unwrap();
        let column = self.matrices[stalk].column(generator.index);
        model.spend(1 + column.len())?;
        let prime = self.source.sheaf.prime;
        for &(index, scalar) in column {
            emit(
                model,
                Generator {
                    vertical,
                    chain: generator.chain.clone(),
                    index,
                },
                prime.mul(coefficient, scalar),
            )?;
        }
        Ok(())
    }
}

pub(crate) struct Complex {
    components: Vec<NamedSheaf>,
    maps: Vec<Rc<Map>>,
    pub(crate) height: usize,
    pub(crate) cohomology: Cohomology,
    pages: RefCell<BTreeMap<Grade, Rc<Page>>>,
}

impl Complex {
    pub(crate) fn build(
        model: &mut Model,
        components: Vec<NamedSheaf>,
        maps: Vec<(Label, Rc<Map>)>,
    ) -> Result<Self> {
        if components.is_empty() || maps.len() + 1 != components.len() {
            return Err(Error::plain(
                "shcomplex requires a sheaf followed by alternating maps and sheaves",
            ));
        }
        model.spend(components.len())?;
        for ((label, map), endpoints) in maps.iter().zip(components.windows(2)) {
            if map.source.name != endpoints[0].name || map.target.name != endpoints[1].name {
                return Err(Error::new(
                    label.span,
                    format!(
                        "sheaf map '{}' has endpoints {} -> {}, expected {} -> {}",
                        label.text,
                        map.source.name,
                        map.target.name,
                        endpoints[0].name,
                        endpoints[1].name
                    ),
                ));
            }
        }
        let prime = components[0].sheaf.prime;
        for consecutive in maps.windows(2) {
            let (label, before) = &consecutive[0];
            let (next, after) = &consecutive[1];
            for (point, (a, b)) in before.matrices.iter().zip(&after.matrices).enumerate() {
                model.spend(1)?;
                let composite = Matrix::compose(model, prime, b, a)?;
                if !composite.is_zero() {
                    return Err(Error::new(
                        next.span,
                        format!(
                            "sheaf complex composite '{} then {}' is nonzero at '{}'",
                            label.text, next.text, before.source.sheaf.poset.names[point]
                        ),
                    ));
                }
            }
        }
        let height = components
            .iter()
            .enumerate()
            .filter(|(_, entry)| entry.sheaf.height != 0)
            .try_fold(0, |height: usize, (q, entry)| {
                q.checked_add(entry.sheaf.height)
                    .map(|edge| height.max(edge))
                    .ok_or_else(|| Error::plain("sheaf total degree is out of range"))
            })?;
        model.allocate(2 * components.len() + maps.len())?;
        Ok(Self {
            components,
            maps: maps.into_iter().map(|(_, map)| map).collect(),
            height,
            cohomology: Cohomology::default(),
            pages: RefCell::new(BTreeMap::new()),
        })
    }

    pub(crate) fn components(&self) -> &[NamedSheaf] {
        &self.components
    }

    pub(crate) fn maps(&self) -> &[Rc<Map>] {
        &self.maps
    }

    pub(crate) fn component_at(&self, degree: usize) -> Option<&NamedSheaf> {
        self.components.get(degree)
    }

    pub(crate) fn page(&self, model: &mut Model, grade: Grade) -> Result<Rc<Page>> {
        if let Some(page) = self.pages.borrow().get(&grade) {
            return Ok(page.clone());
        }
        let page = Rc::new(Page::build(model, self, grade)?);
        model.allocate(1)?;
        self.pages.borrow_mut().insert(grade, page.clone());
        Ok(page)
    }

    pub(crate) fn layer_window(
        &self,
        model: &mut Model,
        degree: usize,
        lo: usize,
        hi: usize,
        limit: usize,
    ) -> Result<Layer> {
        let mut layer = Layer::default();
        if degree >= self.height || lo > hi {
            return Ok(layer);
        }
        let first = degree.saturating_sub(hi);
        let last = degree
            .checked_sub(lo)
            .map(|q| q.min(self.components.len() - 1));
        if let Some(last) = last.filter(|&last| first <= last) {
            for (q, component) in (first..=last).zip(&self.components[first..=last]) {
                model.spend(1)?;
                let part = component.sheaf.layer(model, degree - q, limit)?;
                layer.append_layer(model, part, q, limit)?;
            }
        }
        Ok(layer)
    }

    pub(crate) fn component(&self, degree: usize) -> Result<&NamedSheaf> {
        self.components.get(degree).ok_or_else(|| {
            Error::plain(format!(
                "sheaf complex stalk degree {degree} is out of range"
            ))
        })
    }

    pub(crate) fn poset(&self) -> &Rc<crate::model::Poset> {
        &self.components[0].sheaf.poset
    }

    pub(crate) fn rank(&self, model: &mut Model, degree: usize) -> Result<usize> {
        if degree >= self.height {
            return Ok(0);
        }
        self.cohomology.rank(model, self, degree)
    }

    pub(crate) fn embed(&self, model: &mut Model, vertical: usize, raw: &[Term]) -> Result<Value> {
        model.spend(raw.len())?;
        sheaf::persist(
            model,
            raw.iter()
                .cloned()
                .map(|(mut term, scalar)| {
                    term.vertical = vertical;
                    (term, scalar)
                })
                .collect(),
        )
    }

    pub(crate) fn part(&self, model: &mut Model, vertical: usize, total: &[Term]) -> Result<Value> {
        let start = total.partition_point(|(term, _)| term.vertical < vertical);
        let end = total.partition_point(|(term, _)| term.vertical <= vertical);
        model.spend(1 + 2 * total.len().max(1).ilog2() as usize + end - start)?;
        sheaf::persist(
            model,
            total[start..end]
                .iter()
                .cloned()
                .map(|(mut term, scalar)| {
                    term.vertical = 0;
                    (term, scalar)
                })
                .collect(),
        )
    }
}

impl Cochains for Complex {
    fn prime(&self) -> Prime {
        self.components[0].sheaf.prime
    }

    fn cofaces(
        &self,
        model: &mut Model,
        generator: &Generator,
        coefficient: u32,
        emit: &mut Sink<'_>,
    ) -> Result<()> {
        let q = generator.vertical;
        self.components[q]
            .sheaf
            .cofaces(model, generator, coefficient, emit)?;
        if let Some(map) = self.maps.get(q) {
            // Vertical maps raise q. The sign (-1)^p makes their
            // composition with the horizontal Roos differential cancel.
            let coefficient = if generator.chain.len().is_multiple_of(2) {
                self.prime().negate(coefficient)
            } else {
                coefficient
            };
            map.apply_term(model, generator, coefficient, q + 1, emit)?;
        }
        Ok(())
    }

    fn layer(&self, model: &mut Model, degree: usize, limit: usize) -> Result<Layer> {
        self.layer_window(model, degree, 0, degree, limit)
    }
}
