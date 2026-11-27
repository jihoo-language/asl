use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Poset, TypeId, Value};
use crate::odd_cobar::Prime;
use crate::odd_ext;
use crate::prime_linear::Coordinate;
use crate::sheaf::{Cochains, Cohomology, Term};
use crate::sheaf_hyper::{Complex, NamedSheaf};
use crate::sheaf_image::{Opens, Site};
use crate::sheaf_leray::{build_local_complex, Input, InputForm, Local};
use crate::syntax::Label;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) enum Mode {
    Absolute,
    Punctured,
    Relative,
}

impl Mode {
    pub(crate) fn from_label(label: &Label) -> Result<Self> {
        match label.text.as_str() {
            "absolute" => Ok(Self::Absolute),
            "punctured" => Ok(Self::Punctured),
            "relative" => Ok(Self::Relative),
            _ => Err(Error::new(
                label.span,
                "expected absolute, punctured or relative",
            )),
        }
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Absolute => "absolute",
            Self::Punctured => "punctured",
            Self::Relative => "relative",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) point: usize,
    pub(crate) mode: Mode,
    pub(crate) degree: usize,
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Rep,
    Forget,
    Restrict,
    Connect,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Rep => "shLocalRep",
            Self::Forget => "shLocalForget",
            Self::Restrict => "shLocalRestrict",
            Self::Connect => "shLocalConnect",
        }
    }

    pub(crate) fn output(self, grade: Grade) -> Result<Option<Grade>> {
        let (from, to) = match self {
            Self::Rep => return Ok(None),
            Self::Forget => (Mode::Relative, Mode::Absolute),
            Self::Restrict => (Mode::Absolute, Mode::Punctured),
            Self::Connect => (Mode::Punctured, Mode::Relative),
        };
        if grade.mode != from {
            let article = if from == Mode::Absolute { "an" } else { "a" };
            return Err(Error::plain(format!(
                "{} requires {article} {} local class",
                self.name(),
                from.name()
            )));
        }
        let degree = if matches!(self, Self::Connect) {
            grade
                .degree
                .checked_add(1)
                .ok_or_else(|| Error::plain("local connecting degree is out of range"))?
        } else {
            grade.degree
        };
        Ok(Some(Grade {
            point: grade.point,
            mode: to,
            degree,
        }))
    }
}

pub(crate) struct Supports {
    pub(crate) input: Input,
    site: Site,
    closed: Vec<bool>,
    poset: Rc<Poset>,
    groups: RefCell<BTreeMap<(usize, Mode), Rc<Cohomology>>>,
}

impl Supports {
    pub(crate) fn build(
        model: &mut Model,
        name: &Label,
        input: Input,
        site: Site,
    ) -> Result<(Self, Vec<NamedSheaf>, Complex)> {
        let operation = input.form().extra_name();
        if input.base() != site.target {
            return Err(Error::plain(format!(
                "{operation} requires coefficients on the map's named target"
            )));
        }
        let closed = site.closed_image(model, operation)?;
        let poset = model.poset(site.source).expect("named closed source");
        crate::sheaf::check_cells(
            poset.names.len() * closed.len() * (1 + input.components().len()),
        )?;
        let supports = Self {
            input,
            site,
            closed,
            poset,
            groups: RefCell::default(),
        };
        let count = supports.poset.names.len();
        let locals = (0..count)
            .map(|z| supports.local(model, z, Mode::Relative))
            .collect::<Result<_>>()?;
        let relative = Opens {
            unique: locals,
            indices: (0..count).collect(),
        };
        let (parts, complex) = build_local_complex(
            model,
            name,
            supports.site.source,
            supports.prime(),
            &relative,
            usize::MAX,
        )?;
        for part in &parts {
            model.check_new_sheaf(&Label {
                text: part.name.to_string(),
                span: name.span,
            })?;
        }
        model.allocate(supports.site.cells() + supports.closed.len() + 3)?;
        Ok((supports, parts, complex))
    }

    pub(crate) fn point(&self, model: &Model, label: &Label) -> Result<usize> {
        model.point(self.site.source, label)
    }

    pub(crate) fn point_name(&self, point: usize) -> &str {
        &self.poset.names[point]
    }

    pub(crate) fn prime(&self) -> Prime {
        self.input.cochains().prime()
    }

    pub(crate) fn cochain_type(&self, model: &mut Model, degree: usize) -> Result<TypeId> {
        let form = crate::sheaf::Form::Cochain;
        let name = self.input.name().into();
        match self.input.form() {
            InputForm::Sheaf => model.known_sheaf_type(form, name, degree),
            InputForm::Complex => model.known_sheaf_complex_type(form, name, degree),
        }
    }

    fn in_open(&self, anchor: usize, start: usize) -> bool {
        start == anchor || self.site.target_poset().strict_less(anchor, start)
    }

    fn local(&self, model: &mut Model, point: usize, mode: Mode) -> Result<Local<'_>> {
        let anchor = self.site.image(point);
        model.spend(self.closed.len())?;
        let mask: Vec<_> = self
            .closed
            .iter()
            .enumerate()
            .map(|(start, &closed)| {
                self.in_open(anchor, start)
                    && match mode {
                        Mode::Absolute => true,
                        Mode::Punctured => !closed,
                        Mode::Relative => closed,
                    }
            })
            .collect();
        match mode {
            Mode::Relative => Local::supported(model, &self.input, mask.into()),
            _ => Local::new(model, &self.input, &mask),
        }
    }

    fn groups(&self, model: &mut Model, grade: Grade) -> Result<Rc<Cohomology>> {
        let key = (grade.point, grade.mode);
        if let Some(groups) = self.groups.borrow().get(&key) {
            return Ok(groups.clone());
        }
        model.allocate(5)?;
        let groups = Rc::new(Cohomology::default());
        self.groups.borrow_mut().insert(key, groups.clone());
        Ok(groups)
    }

    pub(crate) fn rank(&self, model: &mut Model, grade: Grade) -> Result<usize> {
        let local = self.local(model, grade.point, grade.mode)?;
        self.groups(model, grade)?.rank(model, &local, grade.degree)
    }

    pub(crate) fn class(&self, model: &mut Model, grade: Grade, raw: &[Term]) -> Result<Value> {
        if raw.is_empty() {
            return Ok(odd_ext::zero());
        }
        let local = self.local(model, grade.point, grade.mode)?;
        let anchor = self.site.image(grade.point);
        let mut terms = Vec::new();
        for (term, coefficient) in raw {
            model.spend(1)?;
            if local.contains(term) {
                terms.push((term.clone(), *coefficient));
            } else if grade.mode == Mode::Relative && self.in_open(anchor, term.place.chain()[0]) {
                return Err(Error::plain(
                    "relative local cochain must vanish on the punctured open",
                ));
            }
        }
        self.groups(model, grade)?
            .class(model, &local, grade.degree, &terms)
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        grade: Grade,
        operation: Operation,
        coordinates: &[Coordinate],
    ) -> Result<Value> {
        let output = operation.output(grade)?;
        if coordinates.is_empty() {
            return Ok(if output.is_some() {
                odd_ext::zero()
            } else {
                Value::SheafC([].into())
            });
        }
        let source = self.local(model, grade.point, grade.mode)?;
        let raw = self.groups(model, grade)?.representative_support(
            model,
            &source,
            grade.degree,
            coordinates,
        )?;
        let Some(output) = output else {
            return crate::sheaf::persist(model, raw);
        };
        let terms = if matches!(operation, Operation::Connect) {
            // A punctured representative is lifted by zero on all relative chains.
            let absolute = self.local(model, grade.point, Mode::Absolute)?;
            absolute.differential_support(model, &raw)?
        } else {
            raw
        };
        self.class(model, output, &terms)
    }
}
