use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Poset, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::{Basis, Coordinate};
use crate::sheaf::{self, Cochain, Cochains, Layer, Matrix, Place, Sheaf, Term};
use crate::sheaf_hyper::{Complex, Map, NamedSheaf};
use crate::sheaf_image::Site;
use crate::sheaf_leray::{self, Input};
use crate::sheaf_stalk;
use crate::syntax::{Label, MAX_IDENTIFIER_BYTES};

pub(crate) struct Square {
    f: Site,
    u: Site,
    v: Site,
    fp: Site,
    pub(crate) x: Rc<Poset>,
    pub(crate) xp: Rc<Poset>,
    pub(crate) yp: Rc<Poset>,
}

impl Square {
    pub(crate) fn build(model: &mut Model, f: Site, u: Site, v: Site, fp: Site) -> Result<Self> {
        if f.source != u.target
            || fp.source != u.source
            || fp.target != v.source
            || f.target != v.target
        {
            return Err(Error::plain(
                "base change requires compatible endpoints for f, u, v and f'",
            ));
        }
        let x = model.poset(f.source).expect("base change site");
        let xp = model.poset(fp.source).expect("base change site");
        let yp = model.poset(fp.target).expect("base change site");
        model.spend(xp.names.len())?;
        for (point, label) in xp.names.iter().enumerate() {
            if f.image(u.image(point)) != v.image(fp.image(point)) {
                return Err(Error::plain(format!(
                    "base-change square does not commute at '{label}'"
                )));
            }
        }
        model.allocate(4 + f.cells() + u.cells() + v.cells() + fp.cells())?;
        Ok(Self {
            f,
            u,
            v,
            fp,
            x,
            xp,
            yp,
        })
    }

    pub(crate) fn x_index(&self, model: &Model, label: &Label) -> Result<usize> {
        model.point(self.f.source, label)
    }

    pub(crate) fn y_index(&self, model: &Model, label: &Label) -> Result<usize> {
        model.point(self.fp.target, label)
    }

    pub(crate) fn in_future(&self, y: usize, x: usize) -> bool {
        self.f.contains_in_future(self.v.image(y), x)
    }

    pub(crate) fn comma(&self, model: &mut Model, y: usize, x: usize) -> Result<Vec<usize>> {
        if !self.in_future(y, x) {
            return Err(Error::plain(
                "comma target is outside the pulled-back principal open",
            ));
        }
        model.spend(1 + self.xp.names.len())?;
        Ok((0..self.xp.names.len())
            .filter(|&xp| {
                let image = self.u.image(xp);
                self.fp.contains_in_future(y, xp) && (image == x || self.x.strict_less(image, x))
            })
            .collect())
    }
}

pub(crate) struct Comparison {
    pub(crate) from_name: Rc<str>,
    pub(crate) to_name: Rc<str>,
    pub(crate) from: Rc<Complex>,
    pub(crate) to: Rc<Complex>,
    maps: Vec<Rc<Map>>,
    homology: RefCell<BTreeMap<usize, Rc<Induced>>>,
    base: crate::model::TypeId,
}

impl Comparison {
    pub(crate) fn build(
        model: &mut Model,
        name: &Label,
        square: &Square,
        input: Input,
    ) -> Result<Self> {
        if input.base() != square.f.source {
            return Err(Error::plain(
                "base-change coefficients must live on the source of f",
            ));
        }
        let from_name = derived_name(name, "from")?;
        let to_name = derived_name(name, "to")?;
        model.check_new_sheaf_complex(&from_name)?;
        model.check_new_sheaf_complex(&to_name)?;
        let pulled = pull_input(model, &input, &square.u, square.xp.names.len(), name)?;
        let base = square.fp.target;
        let count = square.yp.names.len();
        let source =
            input.local_families_at(model, &square.f, (0..count).map(|y| square.v.image(y)))?;
        let target = pulled.local_families_at(model, &square.fp, 0..count)?;
        let prime = input.cochains().prime();
        let (_, from) =
            sheaf_leray::build_local_complex(model, &from_name, base, prime, &source, usize::MAX)?;
        let (_, to) =
            sheaf_leray::build_local_complex(model, &to_name, base, prime, &target, usize::MAX)?;
        for component in from.components().iter().chain(to.components()) {
            model.check_new_sheaf(&Label {
                text: component.name.to_string(),
                span: name.span,
            })?;
        }
        let mut maps = Vec::new();
        for (q, (a, b)) in from.components().iter().zip(to.components()).enumerate() {
            let (cols, rows) = (
                sheaf_leray::local_layers(model, &source, q, &[])?,
                sheaf_leray::local_layers(model, &target, q, &[])?,
            );
            sheaf::check_cells(cols.iter().chain(&rows).map(Layer::cells).sum())?;
            let mut unique = BTreeMap::<(usize, usize), Rc<Matrix>>::new();
            let mut matrices = Vec::with_capacity(count);
            let mut cells = count;
            for (&i, &j) in source.indices.iter().zip(&target.indices) {
                model.spend(1)?;
                let map = if let Some(map) = unique.get(&(i, j)) {
                    map.clone()
                } else {
                    let map = Rc::new(pull_matrix(model, &square.u, &cols[i], &rows[j])?);
                    cells += map.cells();
                    sheaf::check_cells(cells)?;
                    unique.insert((i, j), map.clone());
                    map
                };
                matrices.push(map);
            }
            maps.push(Rc::new(Map::from_matrices(
                model,
                a.clone(),
                b.clone(),
                matrices,
            )?));
        }
        verify_chain_map(model, prime, &from, &to, &maps)?;
        model.allocate(8 + maps.len())?;
        Ok(Self {
            from_name: from_name.text.into(),
            to_name: to_name.text.into(),
            from: Rc::new(from),
            to: Rc::new(to),
            maps,
            homology: RefCell::default(),
            base,
        })
    }

    pub(crate) fn y(&self, model: &Model, label: &Label) -> Result<usize> {
        model.point(self.base, label)
    }

    pub(crate) fn raw_map(&self, degree: usize) -> Result<Rc<Map>> {
        self.maps.get(degree).cloned().ok_or_else(|| {
            Error::plain(format!(
                "no common base-change cochain term in degree {degree}"
            ))
        })
    }

    pub(crate) fn apply_cochains(&self, model: &mut Model, raw: &[Term]) -> Result<Value> {
        let output = self.apply_support(model, raw)?;
        sheaf::persist(model, output)
    }

    fn apply_support(&self, model: &mut Model, raw: &[Term]) -> Result<Vec<Term>> {
        let prime = self.from.prime();
        let mut output = Cochain::default();
        for (term, coefficient) in raw {
            model.spend(1)?;
            let Some(map) = self.maps.get(term.vertical) else {
                continue;
            };
            let point = *term.place.chain().last().expect("nonempty Roos flag");
            for &(index, value) in map.matrix(point).column(term.index) {
                let mut target = term.clone();
                target.index = index;
                output.add(model, prime, target, prime.mul(*coefficient, value))?;
            }
        }
        Ok(output.into_terms())
    }

    pub(crate) fn apply_class(
        &self,
        model: &mut Model,
        degree: usize,
        raw: &[Coordinate],
    ) -> Result<Value> {
        if raw.is_empty() {
            return Ok(crate::odd_ext::zero());
        }
        let before =
            self.from
                .cohomology
                .representative_support(model, self.from.as_ref(), degree, raw)?;
        let after = self.apply_support(model, &before)?;
        self.to
            .cohomology
            .class(model, self.to.as_ref(), degree, &after)
    }

    pub(crate) fn induced(&self, model: &mut Model, degree: usize) -> Result<Rc<Induced>> {
        if let Some(data) = self.homology.borrow().get(&degree) {
            return Ok(data.clone());
        }
        let from = sheaf_stalk::quotients(model, &self.from, degree)?;
        let to = sheaf_stalk::quotients(model, &self.to, degree)?;
        let a = self
            .from
            .component_at(degree)
            .unwrap_or(&self.from.components()[0]);
        let b = self
            .to
            .component_at(degree)
            .unwrap_or(&self.to.components()[0]);
        let prime = self.from.prime();
        let source = sheaf_stalk::induced_sheaf(model, &a.sheaf, &from)?;
        let target = sheaf_stalk::induced_sheaf(model, &b.sheaf, &to)?;
        let zero = Rc::new(Matrix::default());
        let mut matrices = Vec::with_capacity(from.len());
        let mut retained = from.len();
        for point in 0..from.len() {
            model.spend(1)?;
            if from[point].rank() == 0 || to[point].rank() == 0 {
                matrices.push(zero.clone());
                continue;
            }
            let map = self.maps[degree].matrix(point);
            let matrix = Matrix::from_action(model, from[point].rank(), |model, index| {
                let cochain = from[point].representative(model, prime, index)?;
                let image = map.apply(model, prime, &cochain)?;
                to[point].coordinates(model, prime, image)
            })?;
            retained += matrix.cells();
            sheaf::check_cells(retained)?;
            matrices.push(Rc::new(matrix));
        }
        model.allocate(retained + 6)?;
        let data = Rc::new(Induced {
            source: Rc::new(source),
            target: Rc::new(target),
            matrices,
        });
        // Avoid retaining infinitely many identical zero groups for large degrees.
        if degree < self.from.components().len().max(self.to.components().len()) {
            self.homology.borrow_mut().insert(degree, data.clone());
        }
        Ok(data)
    }
}

fn pull_input(
    model: &mut Model,
    input: &Input,
    u: &Site,
    count: usize,
    label: &Label,
) -> Result<Input> {
    let pull_sheaf = |model: &mut Model, entry: &NamedSheaf| -> Result<NamedSheaf> {
        Ok(NamedSheaf {
            name: entry.name.clone(),
            sheaf: Rc::new(u.pull(model, &entry.sheaf)?),
        })
    };
    match input {
        Input::Sheaf(sheaf) => Ok(Input::Sheaf(pull_sheaf(model, sheaf)?)),
        Input::Complex(name, complex) => {
            let components = complex
                .components()
                .iter()
                .map(|s| pull_sheaf(model, s))
                .collect::<Result<Vec<_>>>()?;
            let mut maps = Vec::new();
            for (q, original) in complex.maps().iter().enumerate() {
                model.spend(count)?;
                let matrices = (0..count)
                    .map(|point| original.matrix_rc(u.image(point)))
                    .collect();
                let map = Map::from_matrices(
                    model,
                    components[q].clone(),
                    components[q + 1].clone(),
                    matrices,
                )?;
                maps.push((label.clone(), Rc::new(map)));
            }
            let complex = Complex::build(model, components, maps)?;
            Ok(Input::Complex(name.clone(), Rc::new(complex)))
        }
    }
}

fn pull_matrix(model: &mut Model, u: &Site, source: &Layer, target: &Layer) -> Result<Matrix> {
    let mut columns = BTreeMap::<usize, Vec<Coordinate>>::new();
    let mut cells = 1;
    for row in 0..target.size {
        let mut term = target.generator(row);
        let chain = term.place.chain();
        model.spend(1 + chain.len())?;
        let image: Vec<_> = chain.iter().map(|&point| u.image(point)).collect();
        if image.windows(2).any(|pair| pair[0] == pair[1]) {
            continue;
        }
        term.place = Place::Chain(image.into());
        let column = source.index(model, &term)?;
        model.spend(1 + columns.len().max(1).ilog2() as usize)?;
        let entries = columns.entry(column).or_default();
        cells += if entries.is_empty() { 5 } else { 2 };
        sheaf::check_cells(cells)?;
        entries.push((row, 1));
    }
    Matrix::from_columns(
        columns
            .into_iter()
            .map(|(col, rows)| (col, rows.into()))
            .collect(),
    )
}

fn verify_chain_map(
    model: &mut Model,
    prime: Prime,
    source: &Complex,
    target: &Complex,
    maps: &[Rc<Map>],
) -> Result<()> {
    for degree in 0..maps.len() {
        for point in 0..source.poset().names.len() {
            model.spend(1)?;
            let left = match target.maps().get(degree) {
                Some(d) => {
                    Matrix::compose(model, prime, d.matrix(point), maps[degree].matrix(point))?
                }
                None => Matrix::default(),
            };
            let right = match source.maps().get(degree).zip(maps.get(degree + 1)) {
                Some((d, next)) => {
                    Matrix::compose(model, prime, next.matrix(point), d.matrix(point))?
                }
                None => Matrix::default(),
            };
            if left != right {
                return Err(Error::plain(
                    "base-change map does not commute with the Roos differential",
                ));
            }
        }
    }
    Ok(())
}

pub(crate) struct Induced {
    pub(crate) source: Rc<Sheaf>,
    pub(crate) target: Rc<Sheaf>,
    matrices: Vec<Rc<Matrix>>,
}

impl Induced {
    pub(crate) fn named(
        &self,
        model: &mut Model,
        source_name: Rc<str>,
        target_name: Rc<str>,
    ) -> Result<Map> {
        let source = NamedSheaf {
            name: source_name,
            sheaf: self.source.clone(),
        };
        let target = NamedSheaf {
            name: target_name,
            sheaf: self.target.clone(),
        };
        Map::from_matrices(model, source, target, self.matrices.clone())
    }

    pub(crate) fn report(
        &self,
        model: &mut Model,
        point: usize,
        query: LinearQuery,
    ) -> Result<String> {
        let from = self.source.dims()[point];
        let to = self.target.dims()[point];
        let matrix = &self.matrices[point];
        let prime = self.source.prime;
        let mut basis = Basis::default();
        if matches!(query, LinearQuery::Kernel) {
            for col in 0..from {
                model.spend(1 + matrix.column(col).len())?;
                if let Some(witness) = basis.insert_with_witness(
                    model,
                    prime,
                    matrix.column(col).to_vec(),
                    vec![(col, 1)],
                )? {
                    return Ok(format!("shbckernel = {witness:?}"));
                }
            }
            return Ok("shbckernel = []".to_owned());
        }
        for (_, column) in matrix.columns() {
            model.spend(1 + column.len())?;
            basis.insert(model, prime, column.to_vec())?;
        }
        let rank = basis.rank();
        if matches!(query, LinearQuery::Cokernel) {
            for row in 0..to {
                model.spend(1)?;
                if !basis.has_pivot(row) {
                    return Ok(format!("shbccokernel = [({row}, 1)]"));
                }
            }
            return Ok("shbccokernel = []".to_owned());
        }
        Ok(format!(
            "shbcrank = [{from}, {to}, {rank}, {}, {}]",
            from - rank,
            to - rank
        ))
    }
}

#[derive(Clone, Copy)]
pub(crate) enum LinearQuery {
    Rank,
    Kernel,
    Cokernel,
}

pub(crate) fn derived_name(name: &Label, part: &str) -> Result<Label> {
    let text = format!("{}__{part}", name.text);
    if text.len() > MAX_IDENTIFIER_BYTES {
        return Err(Error::new(
            name.span,
            "generated base-change name exceeds the identifier length limit",
        ));
    }
    Ok(Label {
        text,
        span: name.span,
    })
}
