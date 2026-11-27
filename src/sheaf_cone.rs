use std::collections::BTreeMap;
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};
use crate::prime_linear::Coordinate;
use crate::sheaf::{self, Cochain, Cochains, Generator, Layer, Matrix, Place, Sheaf, Term};
use crate::sheaf_hyper::{Complex, Map, NamedSheaf};
use crate::sheaf_image::{Opens, Site};
use crate::sheaf_leray::{self, Input, Local};
use crate::syntax::{Label, MAX_IDENTIFIER_BYTES};

type Columns = Vec<Rc<[Term]>>;

struct Attachment {
    apex: usize,
    degree: usize,
    columns: Columns,
}

pub(crate) trait Apex {
    fn height(&self) -> Result<usize>;
    fn size(&self, degree: usize) -> usize;
    fn column(&self, degree: usize, index: usize) -> &[Term];
    fn differential(&self, _degree: usize) -> Option<Rc<Matrix>> {
        None
    }
}

impl Apex for Attachment {
    fn column(&self, degree: usize, index: usize) -> &[Term] {
        debug_assert_eq!(degree, self.degree);
        &self.columns[index]
    }

    fn size(&self, degree: usize) -> usize {
        if degree == self.degree {
            self.columns.len()
        } else {
            0
        }
    }

    fn height(&self) -> Result<usize> {
        if self.columns.is_empty() {
            return Ok(0);
        }
        self.degree
            .checked_add(1)
            .ok_or_else(|| Error::plain("cone gluing degree is out of range"))
    }
}

pub(crate) struct Cone {
    pub(crate) source: Input,
    site: Site,
    attachment: Attachment,
    pub(crate) fiber_name: Rc<str>,
    pub(crate) fiber: Rc<Complex>,
}

impl Cone {
    pub(crate) fn build(
        model: &mut Model,
        name: &Label,
        source: Input,
        site: Site,
        apex: usize,
        degree: usize,
        columns: Columns,
    ) -> Result<(Self, Complex)> {
        check_cone(model, &source, &site, apex)?;
        for (j, column) in columns.iter().enumerate() {
            if !source
                .cochains()
                .differential_support(model, column)?
                .is_empty()
            {
                return Err(Error::plain(format!(
                    "cone attachment column {j} is not closed"
                )));
            }
        }
        let attachment = Attachment {
            apex,
            degree,
            columns,
        };
        let fiber_label = derived_label(name, "fiber")?;
        model.check_new_sheaf_complex(&fiber_label)?;
        let opens = source.local_families_at(model, &site, 0..site.target_poset().names.len())?;
        let link_height = opens.unique[opens.indices[apex]].height();
        let ambient = ambient_complex(
            model,
            name,
            site.target,
            &opens,
            apex,
            &attachment,
            link_height,
        )?;
        let fiber = fiber_complex(
            model,
            &fiber_label,
            source.cochains(),
            &attachment,
            link_height,
        )?;
        model.allocate(8 + site.cells() + attachment.columns.len())?;
        Ok((
            Self {
                source,
                site,
                attachment,
                fiber_name: fiber_label.text.into(),
                fiber: Rc::new(fiber),
            },
            ambient,
        ))
    }

    pub(crate) fn degree(&self) -> usize {
        self.attachment.degree
    }

    pub(crate) fn apex_vector(&self, model: &mut Model, index: usize) -> Result<Value> {
        if index >= self.attachment.columns.len() {
            return Err(Error::plain("cone gluing apex index is outside its basis"));
        }
        sheaf::persist(model, vec![(fiber_generator(self.degree(), index), 1)])
    }

    pub(crate) fn lift_tail(
        &self,
        model: &mut Model,
        fiber_degree: usize,
        raw: &[Term],
    ) -> Result<Value> {
        if raw.is_empty() {
            return Ok(Value::SheafC([].into()));
        }
        let layer = self
            .source
            .cochains()
            .layer(model, fiber_degree - 1, MAX_PIVOT_CELLS)?;
        let mut output = Cochain::default();
        self.embed_tail(model, &layer, fiber_degree, raw, 1, &mut output)?;
        sheaf::persist(model, output.into_terms())
    }

    fn embed_tail(
        &self,
        model: &mut Model,
        layer: &Layer,
        fiber_degree: usize,
        raw: &[Term],
        scale: u32,
        output: &mut Cochain,
    ) -> Result<()> {
        let offset = self.attachment.size(fiber_degree);
        let prime = self.fiber.prime();
        for (term, value) in raw {
            let index = layer.index(model, term)? + offset;
            output.add(
                model,
                prime,
                fiber_generator(fiber_degree, index),
                prime.mul(*value, scale),
            )?;
        }
        Ok(())
    }

    pub(crate) fn tail(
        &self,
        model: &mut Model,
        fiber_degree: usize,
        raw: &[Term],
    ) -> Result<Value> {
        if raw.is_empty() {
            return Ok(Value::SheafC([].into()));
        }
        let layer = self
            .source
            .cochains()
            .layer(model, fiber_degree - 1, MAX_PIVOT_CELLS)?;
        let base = self.attachment.size(fiber_degree);
        let prime = self.fiber.prime();
        let mut output = Cochain::default();
        for (term, value) in raw {
            if let Some(index) = term.index.checked_sub(base) {
                output.add(model, prime, layer.generator(index), *value)?;
            }
        }
        sheaf::persist(model, output.into_terms())
    }
}

fn check_cone(model: &mut Model, source: &Input, site: &Site, apex: usize) -> Result<()> {
    if source.base() != site.source {
        return Err(Error::plain(
            "cone gluing coefficients must live on the named open source",
        ));
    }
    site.require_open_embedding(model, "cone gluing")?;
    let count = site.target_poset().names.len();
    let source_count = source.components()[0].sheaf.dims().len();
    model.spend(count + source_count)?;
    if source_count + 1 != count
        || (0..source_count).any(|u| !site.target_poset().strict_less(apex, site.image(u)))
    {
        return Err(Error::plain("cone gluing requires the named apex strictly below exactly the complementing open image"));
    }
    Ok(())
}

pub(crate) fn strict_ambient(
    model: &mut Model,
    name: &Label,
    source: &Input,
    site: &Site,
    apex: usize,
    attachment: &dyn Apex,
) -> Result<Complex> {
    check_cone(model, source, site, apex)?;
    let opens = source.local_families_at(model, site, 0..site.target_poset().names.len())?;
    let height = opens.unique[opens.indices[apex]].height();
    ambient_complex(model, name, site.target, &opens, apex, attachment, height)
}

fn ambient_complex(
    model: &mut Model,
    name: &Label,
    base: TypeId,
    opens: &Opens<Local<'_>>,
    apex: usize,
    attachment: &dyn Apex,
    link_height: usize,
) -> Result<Complex> {
    let height = link_height.max(attachment.height()?).max(1);
    let prime = opens.unique[0].prime();
    let count = opens.indices.len();
    bound_height(model, height, count)?;
    let mut components = Vec::<NamedSheaf>::new();
    let mut maps = Vec::new();
    let mut previous = Vec::new();
    for degree in 0..height {
        let layers = sheaf_leray::local_layers(model, opens, degree, &previous)?;
        let mut dims: Vec<_> = opens.indices.iter().map(|&i| layers[i].size).collect();
        let rank = attachment.size(degree);
        dims[apex] = rank;
        let mut restrictions = BTreeMap::<_, Rc<Matrix>>::new();
        let zero = Rc::new(Matrix::default());
        let sheaf = Sheaf::from_cover_maps(model, base, prime, dims, |model, x, y| {
            let j = opens.indices[y];
            if x == apex {
                if rank == 0 {
                    return Ok(zero.clone());
                }
                let transition = Matrix::from_action(model, rank, |model, column| {
                    index_support(
                        model,
                        attachment.column(degree, column),
                        &layers[j],
                        |term| opens.unique[j].contains(term),
                    )
                })?;
                return Ok(Rc::new(transition));
            }
            let i = opens.indices[x];
            if let Some(matrix) = restrictions.get(&(i, j)) {
                return Ok(matrix.clone());
            }
            let matrix = Rc::new(sheaf_leray::restriction(model, &layers[i], &layers[j])?);
            restrictions.insert((i, j), matrix.clone());
            Ok(matrix)
        })?;
        let component = NamedSheaf {
            name: sheaf_leray::component_name(model, name, degree)?,
            sheaf: Rc::new(sheaf),
        };
        if let Some(from) = components.last() {
            let mut cache = BTreeMap::<usize, Rc<Matrix>>::new();
            let mut blocks = Vec::with_capacity(count);
            let mut cells = count;
            for (x, &i) in opens.indices.iter().enumerate() {
                let block = if x == apex {
                    attachment
                        .differential(degree - 1)
                        .unwrap_or_else(|| zero.clone())
                } else if let Some(block) = cache.get(&i) {
                    block.clone()
                } else {
                    let block = Rc::new(Matrix::from_action(
                        model,
                        previous[i].size,
                        |model, index| {
                            opens.unique[i].indexed_column(
                                model,
                                &previous[i].generator(index),
                                &layers[i],
                            )
                        },
                    )?);
                    cache.insert(i, block.clone());
                    block
                };
                cells += block.cells();
                sheaf::check_cells(cells)?;
                blocks.push(block);
            }
            let edge = Map::from_matrices(model, from.clone(), component.clone(), blocks)?;
            maps.push((edge_label(name, degree - 1), Rc::new(edge)));
        }
        components.push(component);
        previous = layers;
    }
    Complex::build(model, components, maps)
}

fn fiber_complex(
    model: &mut Model,
    name: &Label,
    cochains: &dyn Cochains,
    attachment: &Attachment,
    link_height: usize,
) -> Result<Complex> {
    let height = (link_height + 1).max(attachment.height()?);
    let base = model.named("One").expect("builtin singleton");
    let prime = cochains.prime();
    bound_height(model, height, 1)?;
    let mut components = Vec::<NamedSheaf>::new();
    let mut maps = Vec::new();
    let mut previous = Layer::default();
    for degree in 0..height {
        let layer = if degree == 0 {
            Layer::default()
        } else {
            cochains.layer(model, degree - 1, MAX_PIVOT_CELLS)?
        };
        sheaf::check_cells(previous.cells() + layer.cells())?;
        let rank = attachment.size(degree);
        let sheaf =
            Sheaf::from_cover_maps(model, base, prime, vec![rank + layer.size], |_, _, _| {
                unreachable!("singleton has no covers")
            })?;
        let component = NamedSheaf {
            name: sheaf_leray::component_name(model, name, degree)?,
            sheaf: Rc::new(sheaf),
        };
        if let Some(from) = components.last() {
            let m = attachment.size(degree - 1);
            let matrix = Matrix::from_action(model, m + previous.size, |model, index| {
                let coordinates = if index < m {
                    index_support(model, &attachment.columns[index], &layer, |_| true)?
                } else {
                    let generator = previous.generator(index - m);
                    cochains
                        .indexed_column(model, &generator, &layer)?
                        .into_iter()
                        .map(|(i, value)| (i, prime.negate(value)))
                        .collect()
                };
                Ok(coordinates
                    .into_iter()
                    .map(|(i, value)| (i + rank, value))
                    .collect())
            })?;
            let edge = Map::from_matrices(
                model,
                from.clone(),
                component.clone(),
                vec![Rc::new(matrix)],
            )?;
            maps.push((edge_label(name, degree - 1), Rc::new(edge)));
        }
        components.push(component);
        previous = layer;
    }
    Complex::build(model, components, maps)
}

fn bound_height(model: &mut Model, height: usize, stalks: usize) -> Result<()> {
    let cells = height.saturating_mul(stalks + 1);
    sheaf::check_cells(cells)?;
    model.spend(cells)
}

fn index_support(
    model: &mut Model,
    terms: &[Term],
    layer: &Layer,
    retain: impl Fn(&Generator) -> bool,
) -> Result<Vec<Coordinate>> {
    model.spend(terms.len())?;
    let mut values = Vec::new();
    for (term, value) in terms {
        if retain(term) {
            values.push((layer.index(model, term)?, *value));
        }
    }
    model.spend(values.len() * (values.len().max(1).ilog2() as usize + 1))?;
    values.sort_unstable_by_key(|&(index, _)| index);
    Ok(values)
}

fn fiber_generator(degree: usize, index: usize) -> Generator {
    Generator {
        vertical: degree,
        place: Place::Chain([0].into()),
        index,
    }
}

fn edge_label(name: &Label, degree: usize) -> Label {
    Label {
        text: format!("{} differential at {degree}", name.text),
        span: name.span,
    }
}

fn derived_label(name: &Label, suffix: &str) -> Result<Label> {
    let text = format!("{}__{suffix}", name.text);
    if text.len() > MAX_IDENTIFIER_BYTES {
        return Err(Error::new(
            name.span,
            "cone gluing generated name exceeds the identifier length limit",
        ));
    }
    Ok(Label {
        text,
        span: name.span,
    })
}

pub(crate) struct Homotopy {
    from: Rc<Cone>,
    to: Rc<Cone>,
    corrections: Columns,
}

impl Homotopy {
    pub(crate) fn build(
        model: &mut Model,
        from: Rc<Cone>,
        to: Rc<Cone>,
        corrections: Columns,
    ) -> Result<Self> {
        let a = &from.attachment;
        let b = &to.attachment;
        let source_count = from.source.components()[0].sheaf.dims().len();
        model.spend(source_count + 1)?;
        if a.degree != b.degree
            || a.columns.len() != b.columns.len()
            || from.source.name() != to.source.name()
            || from.source.base() != to.source.base()
            || from.source.form() != to.source.form()
            || from.site.target != to.site.target
            || a.apex != b.apex
            || (0..source_count).any(|u| from.site.image(u) != to.site.image(u))
        {
            return Err(Error::plain("cone homotopy requires two attachments on the same cone with the same coefficient complex, degree and apex rank"));
        }
        if a.degree == 0 {
            model.spend(
                a.columns
                    .iter()
                    .chain(&b.columns)
                    .map(|column| column.len())
                    .sum(),
            )?;
            if !corrections.is_empty() || a.columns != b.columns {
                return Err(Error::plain("a degree-zero cone attachment admits only the zero homotopy to identical columns"));
            }
        } else if corrections.len() != a.columns.len() {
            return Err(Error::plain(format!(
                "cone homotopy requires {} correcting cochains, found {}",
                a.columns.len(),
                corrections.len()
            )));
        }
        let cochains = from.source.cochains();
        let prime = cochains.prime();
        for (j, (h, (first, last))) in corrections
            .iter()
            .zip(a.columns.iter().zip(&b.columns))
            .enumerate()
        {
            let mut result = Cochain::default();
            for (term, value) in cochains.differential_support(model, h)? {
                result.add(model, prime, term, value)?;
            }
            for (term, value) in first.iter() {
                result.add(model, prime, term.clone(), *value)?;
            }
            if result.into_terms().as_slice() != last.as_ref() {
                return Err(Error::plain(format!("cone homotopy column {j} does not satisfy D h = new attachment minus old attachment")));
            }
        }
        model.allocate(2 + corrections.len())?;
        Ok(Self {
            from,
            to,
            corrections,
        })
    }

    pub(crate) fn endpoints(&self, reverse: bool) -> (&Cone, &Cone) {
        if reverse {
            (&self.to, &self.from)
        } else {
            (&self.from, &self.to)
        }
    }

    pub(crate) fn apply_raw(
        &self,
        model: &mut Model,
        degree: usize,
        reverse: bool,
        raw: &[Term],
    ) -> Result<Value> {
        let terms = self.support(model, degree, reverse, raw)?;
        sheaf::persist(model, terms)
    }

    fn support(
        &self,
        model: &mut Model,
        degree: usize,
        reverse: bool,
        raw: &[Term],
    ) -> Result<Vec<Term>> {
        if degree != self.from.degree() || self.corrections.is_empty() || raw.is_empty() {
            model.spend(raw.len())?;
            return Ok(raw.to_vec());
        }
        let prime = self.from.fiber.prime();
        let layer = self
            .from
            .source
            .cochains()
            .layer(model, degree - 1, MAX_PIVOT_CELLS)?;
        let mut result = Cochain::default();
        for (term, value) in raw {
            result.add(model, prime, term.clone(), *value)?;
            if let Some(correction) = self.corrections.get(term.index) {
                let scalar = if reverse {
                    prime.negate(*value)
                } else {
                    *value
                };
                self.to
                    .embed_tail(model, &layer, degree, correction, scalar, &mut result)?;
            }
        }
        Ok(result.into_terms())
    }

    pub(crate) fn apply_class(
        &self,
        model: &mut Model,
        degree: usize,
        reverse: bool,
        raw: &[Coordinate],
    ) -> Result<Value> {
        let (from, to) = self.endpoints(reverse);
        let representative = from.fiber.cohomology.representative_support(
            model,
            from.fiber.as_ref(),
            degree,
            raw,
        )?;
        let image = self.support(model, degree, reverse, &representative)?;
        to.fiber
            .cohomology
            .class(model, to.fiber.as_ref(), degree, &image)
    }
}
