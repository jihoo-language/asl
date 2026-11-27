use std::collections::BTreeMap;
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::sheaf::{self, Cochain, Cochains, Generator, Layer, Matrix, Sheaf, Term};
use crate::sheaf_hyper::{Complex, Map, NamedSheaf};
use crate::sheaf_image::{Opens, Site};
use crate::syntax::{Label, MAX_IDENTIFIER_BYTES};

pub(crate) struct Leray {
    pub(crate) source: NamedSheaf,
    site: Site,
}

impl Leray {
    pub(crate) fn build(
        model: &mut Model,
        name: &Label,
        source: NamedSheaf,
        site: Site,
    ) -> Result<(Self, Vec<NamedSheaf>, Complex)> {
        let result = Self::build_through(model, name, source, site, usize::MAX)?;
        for component in &result.1 {
            model.check_new_sheaf(&Label {
                text: component.name.to_string(),
                span: name.span,
            })?;
        }
        Ok(result)
    }

    pub(crate) fn build_through(
        model: &mut Model,
        name: &Label,
        source: NamedSheaf,
        site: Site,
        maximum_degree: usize,
    ) -> Result<(Self, Vec<NamedSheaf>, Complex)> {
        if source.sheaf.base() != site.source {
            return Err(Error::plain(
                "shleray requires a sheaf on the map's named source",
            ));
        }
        let opens = site.opens(model, &source.sheaf)?;
        let height = opens
            .unique
            .iter()
            .map(|open| open.height)
            .max()
            .unwrap_or(0)
            .max(1)
            .min(maximum_degree.saturating_add(1));
        let names = (0..height)
            .map(|q| component_label(name, q).map(|label| label.text.into()))
            .collect::<Result<Vec<_>>>()?;
        let prime = source.sheaf.prime;
        let mut components = Vec::<NamedSheaf>::with_capacity(height);
        let mut maps = Vec::new();
        let mut previous = Vec::new();
        for (q, component_name) in names.into_iter().enumerate() {
            let layers = local_layers(model, &opens, q, &previous)?;
            let dims = opens.indices.iter().map(|&i| layers[i].size).collect();
            let mut restrictions: BTreeMap<(usize, usize), Rc<Matrix>> = BTreeMap::new();
            let sheaf = Sheaf::from_cover_maps(model, site.target, prime, dims, |model, a, b| {
                let (x, z) = (opens.indices[a], opens.indices[b]);
                if let Some(matrix) = restrictions.get(&(x, z)) {
                    return Ok(matrix.clone());
                }
                let matrix = Rc::new(restriction(model, &layers[x], &layers[z])?);
                restrictions.insert((x, z), matrix.clone());
                Ok(matrix)
            })?;
            let component = NamedSheaf {
                name: component_name,
                sheaf: Rc::new(sheaf),
            };
            if let Some(from) = components.last() {
                let map = differential(model, &opens, &previous, &layers, from, &component)?;
                let label = Label {
                    span: name.span,
                    text: format!("{} differential at {}", name.text, q - 1),
                };
                maps.push((label, Rc::new(map)));
            }
            components.push(component);
            previous = layers;
        }
        let complex = Complex::build(model, components.clone(), maps)?;
        model.allocate(site.cells())?;
        Ok((Self { source, site }, components, complex))
    }

    pub(crate) fn include(&self, model: &mut Model, degree: usize, raw: &[Term]) -> Result<Value> {
        if raw.is_empty() {
            return Ok(Value::SheafC([].into()));
        }
        let opens = self.site.opens(model, &self.source.sheaf)?;
        let layers = local_layers(model, &opens, degree, &[])?;
        let prime = self.source.sheaf.prime;
        let mut output = Cochain::default();
        for (point, &local) in opens.indices.iter().enumerate() {
            let chain = [point].into();
            model.spend(raw.len())?;
            for (term, coefficient) in raw {
                if opens.unique[local].contains(term) {
                    let index = layers[local].index(model, term)?;
                    let generator = Generator {
                        vertical: degree,
                        chain: Rc::clone(&chain),
                        index,
                    };
                    output.add(model, prime, generator, *coefficient)?;
                }
            }
        }
        sheaf::persist(model, output.into_terms())
    }
}

pub(crate) fn component_name(model: &Model, name: &Label, degree: usize) -> Result<Rc<str>> {
    let label = component_label(name, degree)?;
    model.check_new_sheaf(&label)?;
    Ok(label.text.into())
}

fn component_label(name: &Label, degree: usize) -> Result<Label> {
    let text = format!("{}__{degree}", name.text);
    if text.len() > MAX_IDENTIFIER_BYTES {
        return Err(Error::new(
            name.span,
            "generated sheaf component name exceeds the identifier length limit",
        ));
    }
    Ok(Label {
        text,
        span: name.span,
    })
}

fn local_layers(
    model: &mut Model,
    opens: &Opens<'_>,
    degree: usize,
    previous: &[Layer],
) -> Result<Vec<Layer>> {
    let mut cells: usize = previous.iter().map(Layer::cells).sum();
    let mut layers = Vec::new();
    for open in &opens.unique {
        let layer = open.layer(model, degree, MAX_PIVOT_CELLS)?;
        cells += layer.cells();
        sheaf::check_cells(cells)?;
        layers.push(layer);
    }
    Ok(layers)
}

fn restriction(model: &mut Model, larger: &Layer, smaller: &Layer) -> Result<Matrix> {
    // Both layers order the same original X-chain coordinates. Iterating
    // the smaller one gives increasing source coordinates in the larger.
    sheaf::check_cells(1 + 5 * smaller.size)?;
    let mut columns = Vec::with_capacity(smaller.size);
    for row in 0..smaller.size {
        let generator = smaller.generator(row);
        let column = larger.index(model, &generator)?;
        columns.push((column, vec![(row, 1)].into()));
    }
    Matrix::from_columns(columns)
}

fn differential(
    model: &mut Model,
    opens: &Opens<'_>,
    previous: &[Layer],
    current: &[Layer],
    source: &NamedSheaf,
    target: &NamedSheaf,
) -> Result<Map> {
    let mut local = Vec::new();
    let mut cells = 0;
    for (i, open) in opens.unique.iter().enumerate() {
        let matrix = if current[i].size == 0 {
            Matrix::default()
        } else {
            Matrix::from_action(model, previous[i].size, |model, index| {
                let generator = previous[i].generator(index);
                open.indexed_column(model, &generator, &current[i])
            })?
        };
        cells += matrix.cells();
        sheaf::check_cells(cells)?;
        local.push(Rc::new(matrix));
    }
    let matrices = opens.indices.iter().map(|&i| local[i].clone()).collect();
    Map::from_matrices(model, source.clone(), target.clone(), matrices)
}
