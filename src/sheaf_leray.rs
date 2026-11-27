use std::collections::BTreeMap;
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};
use crate::odd_cobar::Prime;
use crate::sheaf::{self, Cochain, Cochains, Generator, Layer, Matrix, Open, Sheaf, Sink, Term};
use crate::sheaf_hyper::{Complex, Map, NamedSheaf};
use crate::sheaf_image::{Opens, Site};
use crate::syntax::{Label, MAX_IDENTIFIER_BYTES};

#[derive(Clone, Copy)]
pub(crate) enum InputForm {
    Sheaf,
    Complex,
}

impl InputForm {
    pub(crate) fn leray_name(self) -> &'static str {
        match self {
            Self::Sheaf => "shleray",
            Self::Complex => "shleraycomplex",
        }
    }

    pub(crate) fn extra_name(self) -> &'static str {
        match self {
            Self::Sheaf => "shextra",
            Self::Complex => "shextracomplex",
        }
    }

    pub(crate) fn deligne_name(self) -> &'static str {
        match self {
            Self::Sheaf => "shdeligne",
            Self::Complex => "shdelignecomplex",
        }
    }
}

pub(crate) enum Input {
    Sheaf(NamedSheaf),
    Complex(Rc<str>, Rc<Complex>),
}

impl Input {
    pub(crate) fn resolve(model: &Model, label: &Label, form: InputForm) -> Result<Self> {
        Ok(match form {
            InputForm::Sheaf => Self::Sheaf(NamedSheaf::resolve(model, label)?),
            InputForm::Complex => {
                Self::Complex(label.text.as_str().into(), model.sheaf_complex(label)?)
            }
        })
    }

    pub(crate) fn form(&self) -> InputForm {
        match self {
            Self::Sheaf(_) => InputForm::Sheaf,
            Self::Complex(..) => InputForm::Complex,
        }
    }

    pub(crate) fn name(&self) -> &str {
        match self {
            Self::Sheaf(entry) => &entry.name,
            Self::Complex(name, _) => name,
        }
    }

    pub(crate) fn base(&self) -> TypeId {
        self.components()[0].sheaf.base()
    }

    pub(crate) fn components(&self) -> &[NamedSheaf] {
        match self {
            Self::Sheaf(entry) => std::slice::from_ref(entry),
            Self::Complex(_, complex) => complex.components(),
        }
    }

    pub(crate) fn cochains(&self) -> &dyn Cochains {
        match self {
            Self::Sheaf(entry) => entry.sheaf.as_ref(),
            Self::Complex(_, complex) => complex.as_ref(),
        }
    }

    fn local_families<'a>(&'a self, model: &mut Model, site: &Site) -> Result<Opens<Local<'a>>> {
        let mut cells = 0;
        site.open_families(model, |model, mask| {
            cells += mask.len() * self.components().len();
            sheaf::check_cells(cells)?;
            Local::new(model, self, mask)
        })
    }
}

pub(crate) struct Local<'a> {
    input: &'a Input,
    components: Vec<Open<'a>>,
    height: usize,
}

impl<'a> Local<'a> {
    pub(crate) fn new(model: &mut Model, input: &'a Input, mask: &[bool]) -> Result<Self> {
        Self::from_opens(model, input, |model, sheaf| Open::new(model, sheaf, mask))
    }

    pub(crate) fn supported(
        model: &mut Model,
        input: &'a Input,
        starts: Rc<[bool]>,
    ) -> Result<Self> {
        Self::from_opens(model, input, |model, sheaf| {
            Open::supported(model, sheaf, starts.clone())
        })
    }

    fn from_opens(
        model: &mut Model,
        input: &'a Input,
        mut make: impl FnMut(&mut Model, &'a Sheaf) -> Result<Open<'a>>,
    ) -> Result<Self> {
        let mut height = 0;
        let mut components = Vec::new();
        for (q, sheaf) in input.components().iter().enumerate() {
            let open = make(model, &sheaf.sheaf)?;
            if open.height != 0 {
                height = height.max(q + open.height);
            }
            components.push(open);
        }
        Ok(Self {
            input,
            components,
            height,
        })
    }

    pub(crate) fn contains(&self, term: &Generator) -> bool {
        // Membership in the open is independent of the input stalk degree.
        self.components[0].contains(term)
    }
}

impl Cochains for Local<'_> {
    fn prime(&self) -> Prime {
        self.input.cochains().prime()
    }

    fn cofaces(
        &self,
        model: &mut Model,
        generator: &Generator,
        coefficient: u32,
        emit: &mut Sink<'_>,
    ) -> Result<()> {
        if generator.vertical + generator.place.chain().len() >= self.height {
            return Ok(());
        }
        self.input.cochains().cofaces(
            model,
            generator,
            coefficient,
            &mut |model, target, value| {
                if self.contains(&target) {
                    emit(model, target, value)?;
                }
                Ok(())
            },
        )
    }

    fn layer(&self, model: &mut Model, degree: usize, limit: usize) -> Result<Layer> {
        let mut layer = Layer::default();
        if degree < self.height {
            for (q, open) in self.components.iter().enumerate().take(degree + 1) {
                model.spend(1)?;
                let part = open.layer(model, degree - q, limit)?;
                layer.append_layer(model, part, q, limit)?;
            }
        }
        Ok(layer)
    }
}

pub(crate) struct Leray {
    pub(crate) source: Input,
    site: Site,
}

impl Leray {
    pub(crate) fn build(
        model: &mut Model,
        name: &Label,
        source: Input,
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
        source: Input,
        site: Site,
        maximum_degree: usize,
    ) -> Result<(Self, Vec<NamedSheaf>, Complex)> {
        if source.base() != site.source {
            let form = source.form();
            let message = match form {
                InputForm::Sheaf => "shleray requires a sheaf on the map's named source",
                InputForm::Complex => "shleraycomplex requires a complex on the map's named source",
            };
            return Err(Error::plain(message));
        }
        let opens = source.local_families(model, &site)?;
        let (components, complex) = build_local_complex(
            model,
            name,
            site.target,
            source.cochains().prime(),
            &opens,
            maximum_degree,
        )?;
        model.allocate(site.cells())?;
        Ok((Self { source, site }, components, complex))
    }

    pub(crate) fn include(&self, model: &mut Model, degree: usize, raw: &[Term]) -> Result<Value> {
        if raw.is_empty() {
            return Ok(Value::SheafC([].into()));
        }
        let opens = self.source.local_families(model, &self.site)?;
        let layers = local_layers(model, &opens, degree, &[])?;
        let prime = self.source.cochains().prime();
        let mut output = Cochain::default();
        for (point, &local) in opens.indices.iter().enumerate() {
            let chain = [point].into();
            model.spend(raw.len())?;
            for (term, coefficient) in raw {
                if opens.unique[local].contains(term) {
                    let index = layers[local].index(model, term)?;
                    let generator = Generator {
                        vertical: degree,
                        place: sheaf::Place::Chain(Rc::clone(&chain)),
                        index,
                    };
                    output.add(model, prime, generator, *coefficient)?;
                }
            }
        }
        sheaf::persist(model, output.into_terms())
    }
}

pub(crate) fn build_local_complex(
    model: &mut Model,
    name: &Label,
    base: TypeId,
    prime: Prime,
    opens: &Opens<Local<'_>>,
    maximum_degree: usize,
) -> Result<(Vec<NamedSheaf>, Complex)> {
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
    let mut components = Vec::<NamedSheaf>::with_capacity(height);
    let mut maps = Vec::new();
    let mut previous = Vec::new();
    for (q, component_name) in names.into_iter().enumerate() {
        let layers = local_layers(model, opens, q, &previous)?;
        let dims = opens.indices.iter().map(|&i| layers[i].size).collect();
        let mut restrictions: BTreeMap<(usize, usize), Rc<Matrix>> = BTreeMap::new();
        let sheaf = Sheaf::from_cover_maps(model, base, prime, dims, |model, a, b| {
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
            let map = differential(model, opens, &previous, &layers, from, &component)?;
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
    Ok((components, complex))
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
    opens: &Opens<Local<'_>>,
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
    opens: &Opens<Local<'_>>,
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
