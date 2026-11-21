use crate::error::{Error, Result};
use crate::linear::Basis;
use crate::map_ext;
use crate::model::{Kind, Model, TypeId, Value};
use crate::profile::Profile;

#[derive(Clone, Copy)]
pub(crate) enum Direction {
    Next,
    Previous,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum PageFamily {
    May,
    Ce,
}

impl PageFamily {
    fn name(self, direction: Direction) -> &'static str {
        match (self, direction) {
            (Self::May, Direction::Next) => "mayNext",
            (Self::May, Direction::Previous) => "mayPrev",
            (Self::Ce, Direction::Next) => "ceNext",
            (Self::Ce, Direction::Previous) => "cePrev",
        }
    }

    fn expected(self) -> &'static str {
        match self {
            Self::May => "a May or MapMay class",
            Self::Ce => "a CEAn class",
        }
    }
}

impl Direction {
    pub(crate) fn result_type(
        self,
        model: &mut Model,
        ty: TypeId,
        family: PageFamily,
    ) -> Result<TypeId> {
        if family == PageFamily::May {
            if let Kind::PMay(p, r, s, t, f) = *model.kind(ty) {
                let page = self.target_page(r, family)?;
                return model.pmay_type(p, page, s, t, f);
            }
        }
        let grade = Grade::from_type(model, ty)
            .filter(|grade| grade.family.page_family() == family)
            .ok_or_else(|| {
                Error::plain(format!(
                    "{} requires {}, found {}",
                    family.name(self),
                    family.expected(),
                    model.type_name(ty),
                ))
            })?;
        self.changed_type(model, grade)
    }

    fn changed_type(self, model: &mut Model, mut grade: Grade) -> Result<TypeId> {
        grade.page = self.target_page(grade.page, grade.family.page_family())?;
        grade.ty(model)
    }

    fn target_page(self, page: usize, family: PageFamily) -> Result<usize> {
        let name = family.name(self);
        match self {
            Self::Next => page
                .checked_add(1)
                .ok_or_else(|| Error::plain(format!("{name} page is out of range"))),
            Self::Previous => page
                .checked_sub(1)
                .ok_or_else(|| Error::plain(format!("{name} requires a page greater than zero"))),
        }
    }
}

#[derive(Clone, Copy)]
enum Family {
    Scalar(Profile, usize),
    Mapping(TypeId, TypeId, i64),
    CentralExtension(usize, usize),
}

#[derive(Clone, Copy)]
struct Grade {
    family: Family,
    page: usize,
    length: usize,
    filtration: usize,
}

impl Grade {
    fn from_type(model: &Model, ty: TypeId) -> Option<Self> {
        let (family, page, length, filtration) = match *model.kind(ty) {
            Kind::MapMay(source, target, page, length, degree, filtration) => (
                Family::Mapping(source, target, degree),
                page,
                length,
                filtration,
            ),
            Kind::CeAn(n, page, length, degree, filtration) => (
                Family::CentralExtension(n, degree),
                page,
                length,
                filtration,
            ),
            ref kind => {
                let (profile, page, length, degree, filtration) = kind.may_grade()?;
                (Family::Scalar(profile, degree), page, length, filtration)
            }
        };
        Some(Self {
            family,
            page,
            length,
            filtration,
        })
    }

    fn ty(self, model: &mut Model) -> Result<TypeId> {
        match self.family {
            Family::Scalar(profile, degree) => {
                model.profile_may_type(profile, self.page, self.length, degree, self.filtration)
            }
            Family::Mapping(source, target, degree) => model.map_may_type(
                source,
                target,
                self.page,
                self.length,
                degree,
                self.filtration,
            ),
            Family::CentralExtension(n, degree) => {
                model.ce_type(n, self.page, self.length, degree, self.filtration)
            }
        }
    }

    fn incoming(self) -> Option<Self> {
        Some(Self {
            length: self.length.checked_sub(1)?,
            filtration: self.filtration.checked_sub(self.page)?,
            ..self
        })
    }
}

impl Family {
    fn page_family(self) -> PageFamily {
        match self {
            Self::Scalar(..) | Self::Mapping(..) => PageFamily::May,
            Self::CentralExtension(..) => PageFamily::Ce,
        }
    }

    fn rank(self, model: &mut Model, ty: TypeId) -> Result<usize> {
        match self {
            Self::Scalar(..) => Ok(model.may_page(ty)?.rank()),
            Self::Mapping(..) => Ok(model.map_may(ty)?.rank()),
            Self::CentralExtension(..) => Ok(model.ce_page(ty)?.rank()),
        }
    }

    fn lift_coordinates(
        self,
        model: &mut Model,
        from: TypeId,
        to: TypeId,
        class: &[usize],
        differentiate: bool,
    ) -> Result<Vec<usize>> {
        match self {
            Self::Scalar(profile, _) => {
                let words = model.may_page(from)?.representative(model, class)?;
                let words = if differentiate {
                    profile.differential_support(model, &words)?
                } else {
                    words
                };
                model.may_page(to)?.project(model, &words)
            }
            Self::Mapping(source, target, _) => {
                let terms = model.map_may(from)?.representative(model, class)?;
                let terms = if differentiate {
                    map_ext::differential_support(model, source, target, &terms)?
                } else {
                    terms
                };
                model.map_may(to)?.project(model, &terms)
            }
            Self::CentralExtension(n, _) => {
                let words = model.ce_page(from)?.representative(model, class)?;
                let words = if differentiate {
                    Profile::An(n).differential_support(model, &words)?
                } else {
                    words
                };
                model.ce_page(to)?.project(model, &words)
            }
        }
    }
}

pub(crate) struct Transition {
    page: usize,
    cycles: Basis,
    name: &'static str,
}

impl Transition {
    pub(crate) fn build(model: &mut Model, current: TypeId) -> Result<Self> {
        let grade = Grade::from_type(model, current).unwrap();
        let name = grade.family.page_family().name(Direction::Next);
        let mut cycles = Basis::default();
        if let Some(incoming) = grade.incoming() {
            let incoming = incoming.ty(model)?;
            for index in 0..grade.family.rank(model, incoming)? {
                let boundary =
                    grade
                        .family
                        .lift_coordinates(model, incoming, current, &[index], true)?;
                cycles.insert(model, boundary)?;
            }
        }

        let next = Direction::Next.changed_type(model, grade)?;
        for index in 0..grade.family.rank(model, next)? {
            // The raw next-page lift is already corrected through filtration
            // f+r. Its image on the current page is a cycle. Its witness
            // records the next-page coordinate modulo incoming boundaries.
            let earlier = grade
                .family
                .lift_coordinates(model, next, current, &[index], false)?;
            if cycles
                .insert_with_witness(model, earlier, vec![index])?
                .is_some()
            {
                return Err(Error::plain(format!("{name} has dependent cycle lifts")));
            }
        }
        model.allocate(2 + cycles.cells())?;
        Ok(Self {
            page: grade.page,
            cycles,
            name,
        })
    }

    fn project(&self, model: &mut Model, class: &[usize]) -> Result<Vec<usize>> {
        model.spend(class.len())?;
        let (remainder, next) =
            self.cycles
                .reduce_with_witness(model, class.to_vec(), Vec::new())?;
        if !remainder.is_empty() {
            return Err(Error::plain(format!(
                "{} requires a cycle for d_{}",
                self.name, self.page,
            )));
        }
        Ok(next)
    }
}

pub(crate) fn next(model: &mut Model, current: TypeId, class: &[usize]) -> Result<Value> {
    if class.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let coordinates = model.page_transition(current)?.project(model, class)?;
    model.allocate(coordinates.len())?;
    Ok(Value::Class(coordinates.into()))
}

pub(crate) fn previous(
    model: &mut Model,
    source: TypeId,
    target: TypeId,
    class: &[usize],
) -> Result<Value> {
    if class.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let grade = Grade::from_type(model, source).unwrap();
    let coordinates = grade
        .family
        .lift_coordinates(model, source, target, class, false)?;
    model.allocate(coordinates.len())?;
    Ok(Value::Class(coordinates.into()))
}
