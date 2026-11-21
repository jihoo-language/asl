use std::cmp::Ordering;
use std::collections::HashMap;
use std::rc::Rc;

use crate::adams::Adams;
use crate::cobar;
use crate::cohomology::Cohomology;
use crate::comodule;
use crate::error::{Error, Result};
use crate::map_cochain;
use crate::map_ext::{self, MapAdams};
use crate::map_may::MapMayPage;
use crate::may::MayPage;
use crate::milnor::MilnorBasis;
use crate::odd_cobar::{self, Prime};
use crate::odd_ext::{self, OddExt};
use crate::page_transition::Transition;
use crate::prime_linear::Coordinate;
use crate::profile::Profile;
use crate::steenrod::{self, Word};
use crate::syntax::Label;

const MAX_SPACE_POINTS: usize = 512;
const MAX_TYPES: usize = 4096;
const MAX_TYPE_DEPTH: usize = 64;
const MAX_ELEMENTS: usize = 4096;
const MAX_WORK: usize = 25_000_000;
const MAX_STORAGE: usize = 1_000_000;
const MAX_SHOWN_BYTES: usize = 65_536;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) struct TypeId(usize);

#[derive(Clone, PartialEq, Eq, Hash)]
pub(crate) enum Kind {
    Space(String),
    Star(TypeId, usize),
    Product(TypeId, TypeId),
    Arrow(TypeId, TypeId),
    Cohomology(TypeId, usize),
    Operations(usize),
    Cobar(usize, usize),
    PCobar(Prime, usize, usize),
    PExt(Prime, usize, usize),
    CobarAn(usize, usize, usize),
    Ext(usize, usize),
    ExtAn(usize, usize, usize),
    XExt(TypeId, usize, usize),
    MapCobar(TypeId, TypeId, usize, i64),
    MapExt(TypeId, TypeId, usize, i64),
    MapMay(TypeId, TypeId, usize, usize, i64, usize),
    May(usize, usize, usize, usize),
    MayAn(usize, usize, usize, usize, usize),
    Suspension(TypeId),
    Stable(TypeId, usize, usize),
    Bar(usize),
    CyclicH2(usize),
}

impl Kind {
    pub(crate) fn may_grade(&self) -> Option<(Profile, usize, usize, usize, usize)> {
        match *self {
            Self::May(r, s, t, f) => Some((Profile::Full, r, s, t, f)),
            Self::MayAn(n, r, s, t, f) => Some((Profile::An(n), r, s, t, f)),
            _ => None,
        }
    }

    pub(crate) fn cobar_grade(&self) -> Option<(Profile, usize, usize)> {
        match *self {
            Self::Cobar(s, t) => Some((Profile::Full, s, t)),
            Self::CobarAn(n, s, t) => Some((Profile::An(n), s, t)),
            _ => None,
        }
    }

    pub(crate) fn ext_grade(&self) -> Option<(Profile, usize, usize)> {
        match *self {
            Self::Ext(s, t) => Some((Profile::Full, s, t)),
            Self::ExtAn(n, s, t) => Some((Profile::An(n), s, t)),
            _ => None,
        }
    }
}

#[derive(Clone)]
pub(crate) enum Value {
    Point(usize),
    Pair(Rc<(Value, Value)>),
    Map(Rc<[Value]>),
    Class(Rc<[usize]>),
    Operation(Rc<[Word]>),
    Cobar(Rc<[cobar::Word]>),
    PCobar(odd_cobar::Support),
    OddClass(Rc<[Coordinate]>),
    MapCobar(map_cochain::Support),
    Suspension(SuspensionPoint),
    Bar(Rc<[usize]>),
    CyclicClass(bool),
}

#[derive(Clone)]
pub(crate) enum SuspensionPoint {
    Equator(Rc<Value>),
    North,
    South,
}

pub(crate) struct Poset {
    pub(crate) names: Vec<String>,
    lookup: HashMap<String, usize>,
    order: Vec<Vec<bool>>,
}

impl Poset {
    fn new(names: Vec<String>, order: Vec<Vec<bool>>) -> Self {
        let lookup = names
            .iter()
            .enumerate()
            .map(|(index, name)| (name.clone(), index))
            .collect();
        Self {
            names,
            lookup,
            order,
        }
    }
}

struct TypeInfo {
    kind: Kind,
    poset: Option<Rc<Poset>>,
    parent_points: Option<Rc<[usize]>>,
    elements: Option<Rc<[Value]>>,
    cohomology: Option<Rc<Cohomology>>,
    adams: Option<Rc<Adams>>,
    odd_ext: Option<Rc<OddExt>>,
    xext: Option<Rc<comodule::Adams>>,
    map_ext: Option<Rc<MapAdams>>,
    map_may: Option<Rc<MapMayPage>>,
    may_page: Option<Rc<MayPage>>,
    page_transition: Option<Rc<Transition>>,
    cardinality_bound: usize,
    heights: [Option<usize>; 2],
    inhabited: bool,
    depth: usize,
}

pub(crate) struct Model {
    types: Vec<TypeInfo>,
    type_ids: HashMap<Kind, TypeId>,
    milnor_bases: HashMap<usize, Rc<MilnorBasis>>,
    work: usize,
    storage: usize,
}

impl Model {
    pub(crate) fn new() -> Self {
        let mut model = Self {
            types: Vec::new(),
            type_ids: HashMap::new(),
            milnor_bases: HashMap::new(),
            work: MAX_WORK,
            storage: MAX_STORAGE,
        };
        model.insert(
            Kind::Space("One".into()),
            Some(Poset::new(vec!["pt".into()], vec![vec![true]])),
            None,
            true,
            0,
        );
        model.insert(
            Kind::Space("S".into()),
            Some(Poset::new(
                vec!["lo".into(), "hi".into()],
                vec![vec![true, true], vec![false, true]],
            )),
            None,
            true,
            0,
        );
        model
    }

    pub(crate) fn reset_work(&mut self) {
        self.work = MAX_WORK;
    }

    pub(crate) fn spend(&mut self, amount: usize) -> Result<()> {
        self.work = self
            .work
            .checked_sub(amount)
            .ok_or_else(|| Error::plain("statement exceeds the computation limit"))?;
        Ok(())
    }

    pub(crate) fn allocate(&mut self, amount: usize) -> Result<()> {
        self.storage = self
            .storage
            .checked_sub(amount)
            .ok_or_else(|| Error::plain("program exceeds the finite-value storage limit"))?;
        Ok(())
    }

    fn insert(
        &mut self,
        kind: Kind,
        poset: Option<Poset>,
        parent_points: Option<Rc<[usize]>>,
        inhabited: bool,
        depth: usize,
    ) -> TypeId {
        let cardinality_bound = match kind {
            Kind::Space(_) | Kind::Star(..) => poset.as_ref().unwrap().names.len(),
            Kind::Product(left, right) => self.types[left.0]
                .cardinality_bound
                .saturating_mul(self.types[right.0].cardinality_bound),
            Kind::Arrow(domain, range) => {
                let base = self.types[range.0].cardinality_bound;
                let exponent = self.types[domain.0].cardinality_bound;
                match (base, exponent) {
                    _ if !self.types[domain.0].inhabited => 1,
                    (base @ 0..=1, _) => base,
                    (_, exponent) if exponent >= usize::BITS as usize => usize::MAX,
                    _ => base.saturating_pow(exponent as u32),
                }
            }
            Kind::Suspension(base) => self.types[base.0].cardinality_bound.saturating_add(2),
            Kind::Cohomology(..) | Kind::Operations(..) | Kind::Stable(..) => usize::MAX,
            Kind::PCobar(_, s, t) | Kind::PExt(_, s, t) if s > t || (s == 0 && t != 0) => 1,
            Kind::PCobar(prime, s, t) | Kind::PExt(prime, s, t) if s == t => prime.get() as usize,
            Kind::PCobar(..) | Kind::PExt(..) => usize::MAX,
            Kind::Cobar(s, t) | Kind::Ext(s, t) | Kind::CobarAn(_, s, t) | Kind::ExtAn(_, s, t)
                if s > t || (s == 0 && t != 0) =>
            {
                1
            }
            Kind::Cobar(s, t) | Kind::Ext(s, t) | Kind::CobarAn(_, s, t) | Kind::ExtAn(_, s, t)
                if s == t =>
            {
                2
            }
            Kind::MapMay(source, target, _, s, _, f)
                if !self.types[source.0].inhabited
                    || !self.types[target.0].inhabited
                    || (s == 0 && f != 0) =>
            {
                1
            }
            Kind::MapMay(..) => usize::MAX,
            Kind::May(_, s, t, f) | Kind::MayAn(_, _, s, t, f)
                if s > t || f > t - s || (s == 0 && t != 0) =>
            {
                1
            }
            Kind::May(_, s, t, _) | Kind::MayAn(_, _, s, t, _) if s == t => 2,
            Kind::XExt(base, s, t) if s > t || !self.types[base.0].inhabited => 1,
            Kind::MapExt(source, target, ..) | Kind::MapCobar(source, target, ..)
                if !self.types[source.0].inhabited || !self.types[target.0].inhabited =>
            {
                1
            }
            Kind::Cobar(..)
            | Kind::CobarAn(..)
            | Kind::Ext(..)
            | Kind::ExtAn(..)
            | Kind::XExt(..)
            | Kind::MapCobar(..)
            | Kind::MapExt(..)
            | Kind::May(..)
            | Kind::MayAn(..) => usize::MAX,
            Kind::Bar(modulus) => support_space_size(modulus - 1),
            Kind::CyclicH2(modulus) => 2 - (modulus & 1),
        };
        let id = TypeId(self.types.len());
        self.type_ids.insert(kind.clone(), id);
        self.types.push(TypeInfo {
            kind,
            poset: poset.map(Rc::new),
            parent_points,
            elements: None,
            cohomology: None,
            adams: None,
            odd_ext: None,
            xext: None,
            map_ext: None,
            map_may: None,
            may_page: None,
            page_transition: None,
            cardinality_bound,
            heights: [None; 2],
            inhabited,
            depth,
        });
        id
    }

    fn check_type_limit(&self, depth: usize) -> Result<()> {
        if self.types.len() == MAX_TYPES {
            return Err(Error::plain("program exceeds the number of allowed types"));
        }
        if depth > MAX_TYPE_DEPTH {
            return Err(Error::plain("type nesting exceeds 64 levels"));
        }
        Ok(())
    }

    pub(crate) fn named(&self, name: &str) -> Option<TypeId> {
        self.type_ids.get(&Kind::Space(name.into())).copied()
    }

    pub(crate) fn kind(&self, id: TypeId) -> &Kind {
        &self.types[id.0].kind
    }

    pub(crate) fn intrinsically_discrete(&self, ty: TypeId) -> bool {
        matches!(
            self.kind(ty),
            Kind::Cohomology(..)
                | Kind::Operations(..)
                | Kind::Cobar(..)
                | Kind::PCobar(..)
                | Kind::PExt(..)
                | Kind::CobarAn(..)
                | Kind::Ext(..)
                | Kind::ExtAn(..)
                | Kind::XExt(..)
                | Kind::MapCobar(..)
                | Kind::MapExt(..)
                | Kind::MapMay(..)
                | Kind::May(..)
                | Kind::MayAn(..)
                | Kind::Stable(..)
                | Kind::Bar(..)
                | Kind::CyclicH2(..)
        )
    }

    pub(crate) fn define_space(&mut self, name: &Label, chains: &[Vec<Label>]) -> Result<()> {
        if self.named(&name.text).is_some() {
            return Err(Error::new(
                name.span,
                format!("space '{}' is already defined", name.text),
            ));
        }
        self.check_type_limit(0)?;

        let mut labels = Vec::<&Label>::new();
        let mut indices = HashMap::<&str, usize>::new();
        let mut edges = Vec::<Vec<usize>>::new();
        for chain in chains {
            let mut previous = None;
            for label in chain {
                self.spend(1)?;
                let point = *indices.entry(&label.text).or_insert_with(|| {
                    labels.push(label);
                    edges.push(Vec::new());
                    labels.len() - 1
                });
                if labels.len() > MAX_SPACE_POINTS {
                    return Err(Error::new(
                        label.span,
                        "space exceeds the limit of 512 points",
                    ));
                }
                if let Some(parent) = previous {
                    let successors: &mut Vec<usize> = &mut edges[parent];
                    successors.push(point);
                }
                previous = Some(point);
            }
        }
        for successors in &mut edges {
            successors.sort_unstable();
            successors.dedup();
        }
        let point_count = labels.len();
        self.allocate(point_count * point_count)?;
        let mut order = Vec::with_capacity(point_count);
        for (source, label) in labels.iter().enumerate() {
            let mut visited = vec![false; point_count];
            let mut stack = vec![source];
            while let Some(point) = stack.pop() {
                self.spend(1 + edges[point].len())?;
                for &next in &edges[point] {
                    if next == source {
                        return Err(Error::new(
                            label.span,
                            format!(
                                "cycle in space '{}' through point '{}'",
                                name.text, label.text
                            ),
                        ));
                    }
                    if !visited[next] {
                        visited[next] = true;
                        stack.push(next);
                    }
                }
            }
            visited[source] = true;
            order.push(visited);
        }
        let names = labels.into_iter().map(|point| point.text.clone()).collect();
        let poset = Poset::new(names, order);
        self.insert(
            Kind::Space(name.text.clone()),
            Some(poset),
            None,
            point_count != 0,
            0,
        );
        Ok(())
    }

    pub(crate) fn product(&mut self, left: TypeId, right: TypeId) -> Result<TypeId> {
        let inhabited = self.types[left.0].inhabited && self.types[right.0].inhabited;
        self.composite(Kind::Product(left, right), left, right, inhabited)
    }

    pub(crate) fn arrow(&mut self, domain: TypeId, range: TypeId) -> Result<TypeId> {
        let inhabited = !self.types[domain.0].inhabited || self.types[range.0].inhabited;
        self.composite(Kind::Arrow(domain, range), domain, range, inhabited)
    }

    pub(crate) fn cohomology_type(&mut self, base: TypeId, degree: usize) -> Result<TypeId> {
        self.inhabited_type(Kind::Cohomology(base, degree), self.types[base.0].depth + 1)
    }

    pub(crate) fn operation_type(&mut self, degree: usize) -> Result<TypeId> {
        self.inhabited_type(Kind::Operations(degree), 0)
    }

    pub(crate) fn milnor_basis(&mut self, degree: usize) -> Result<Rc<MilnorBasis>> {
        if let Some(data) = self.milnor_bases.get(&degree) {
            return Ok(data.clone());
        }
        let data = Rc::new(MilnorBasis::build(self, degree)?);
        self.milnor_bases.insert(degree, data.clone());
        Ok(data)
    }

    pub(crate) fn cobar_type(&mut self, s: usize, t: usize) -> Result<TypeId> {
        self.profile_cobar_type(Profile::Full, s, t)
    }

    pub(crate) fn odd_ext_type(&mut self, p: Prime, s: usize, t: usize) -> Result<TypeId> {
        self.inhabited_type(Kind::PExt(p, s, t), 0)
    }

    pub(crate) fn pcobar_type(&mut self, p: Prime, s: usize, t: usize) -> Result<TypeId> {
        self.inhabited_type(Kind::PCobar(p, s, t), 0)
    }

    pub(crate) fn profile_cobar_type(
        &mut self,
        profile: Profile,
        s: usize,
        t: usize,
    ) -> Result<TypeId> {
        let kind = match profile {
            Profile::Full => Kind::Cobar(s, t),
            Profile::An(n) => Kind::CobarAn(n, s, t),
        };
        self.inhabited_type(kind, 0)
    }

    pub(crate) fn ext_type(&mut self, s: usize, t: usize) -> Result<TypeId> {
        self.profile_ext_type(Profile::Full, s, t)
    }

    pub(crate) fn profile_ext_type(
        &mut self,
        profile: Profile,
        s: usize,
        t: usize,
    ) -> Result<TypeId> {
        let kind = match profile {
            Profile::Full => Kind::Ext(s, t),
            Profile::An(n) => Kind::ExtAn(n, s, t),
        };
        self.inhabited_type(kind, 0)
    }

    pub(crate) fn xext_type(&mut self, base: TypeId, s: usize, t: usize) -> Result<TypeId> {
        self.inhabited_type(Kind::XExt(base, s, t), self.types[base.0].depth + 1)
    }

    pub(crate) fn map_may_type(
        &mut self,
        source: TypeId,
        target: TypeId,
        r: usize,
        s: usize,
        t: i64,
        f: usize,
    ) -> Result<TypeId> {
        let depth = self.types[source.0].depth.max(self.types[target.0].depth) + 1;
        self.inhabited_type(Kind::MapMay(source, target, r, s, t, f), depth)
    }

    pub(crate) fn map_ext_type(
        &mut self,
        source: TypeId,
        target: TypeId,
        s: usize,
        t: i64,
    ) -> Result<TypeId> {
        let depth = self.types[source.0].depth.max(self.types[target.0].depth) + 1;
        self.inhabited_type(Kind::MapExt(source, target, s, t), depth)
    }

    pub(crate) fn map_cobar_type(
        &mut self,
        source: TypeId,
        target: TypeId,
        s: usize,
        t: i64,
    ) -> Result<TypeId> {
        let depth = self.types[source.0].depth.max(self.types[target.0].depth) + 1;
        self.inhabited_type(Kind::MapCobar(source, target, s, t), depth)
    }

    pub(crate) fn may_type(&mut self, r: usize, s: usize, t: usize, f: usize) -> Result<TypeId> {
        self.profile_may_type(Profile::Full, r, s, t, f)
    }

    pub(crate) fn profile_may_type(
        &mut self,
        profile: Profile,
        r: usize,
        s: usize,
        t: usize,
        f: usize,
    ) -> Result<TypeId> {
        let kind = match profile {
            Profile::Full => Kind::May(r, s, t, f),
            Profile::An(n) => Kind::MayAn(n, r, s, t, f),
        };
        self.inhabited_type(kind, 0)
    }

    pub(crate) fn cyclic_type(&mut self, modulus: usize, cohomology: bool) -> Result<TypeId> {
        if modulus < 2 {
            return Err(Error::plain("cyclic modulus must be at least 2"));
        }
        let kind = if cohomology {
            Kind::CyclicH2(modulus)
        } else {
            Kind::Bar(modulus)
        };
        self.inhabited_type(kind, 0)
    }

    pub(crate) fn suspension_type(&mut self, base: TypeId) -> Result<TypeId> {
        self.inhabited_type(Kind::Suspension(base), self.types[base.0].depth + 1)
    }

    pub(crate) fn stable_type(
        &mut self,
        base: TypeId,
        shift: usize,
        degree: usize,
    ) -> Result<TypeId> {
        self.inhabited_type(
            Kind::Stable(base, shift, degree),
            self.types[base.0].depth + 1,
        )
    }

    fn inhabited_type(&mut self, kind: Kind, depth: usize) -> Result<TypeId> {
        if let Some(&id) = self.type_ids.get(&kind) {
            return Ok(id);
        }
        self.check_type_limit(depth)?;
        Ok(self.insert(kind, None, None, true, depth))
    }

    pub(crate) fn cohomology(&mut self, ty: TypeId) -> Result<Rc<Cohomology>> {
        if let Some(data) = &self.types[ty.0].cohomology {
            return Ok(data.clone());
        }
        let Kind::Cohomology(base, degree) = *self.kind(ty) else {
            unreachable!()
        };
        let data = Rc::new(Cohomology::build(self, base, degree)?);
        self.types[ty.0].cohomology = Some(data.clone());
        Ok(data)
    }

    pub(crate) fn odd_ext(&mut self, ty: TypeId) -> Result<Rc<OddExt>> {
        if let Some(data) = &self.types[ty.0].odd_ext {
            return Ok(data.clone());
        }
        let Kind::PExt(prime, s, t) = *self.kind(ty) else {
            unreachable!()
        };
        let data = Rc::new(OddExt::build(self, prime, s, t)?);
        self.types[ty.0].odd_ext = Some(data.clone());
        Ok(data)
    }

    pub(crate) fn adams(&mut self, ty: TypeId) -> Result<Rc<Adams>> {
        if let Some(data) = &self.types[ty.0].adams {
            return Ok(data.clone());
        }
        let (profile, s, t) = self.kind(ty).ext_grade().unwrap();
        let data = Rc::new(Adams::build(self, profile, s, t)?);
        self.types[ty.0].adams = Some(data.clone());
        Ok(data)
    }

    pub(crate) fn xext(&mut self, ty: TypeId) -> Result<Rc<comodule::Adams>> {
        if let Some(data) = &self.types[ty.0].xext {
            return Ok(data.clone());
        }
        let Kind::XExt(base, s, t) = *self.kind(ty) else {
            unreachable!()
        };
        let data = Rc::new(comodule::Adams::build(self, base, s, t)?);
        self.types[ty.0].xext = Some(data.clone());
        Ok(data)
    }

    pub(crate) fn map_may(&mut self, ty: TypeId) -> Result<Rc<MapMayPage>> {
        if let Some(page) = &self.types[ty.0].map_may {
            return Ok(page.clone());
        }
        let Kind::MapMay(source, target, r, s, t, f) = *self.kind(ty) else {
            unreachable!()
        };
        let page = Rc::new(MapMayPage::build(self, source, target, r, s, t, f)?);
        self.types[ty.0].map_may = Some(page.clone());
        Ok(page)
    }

    pub(crate) fn map_ext(&mut self, ty: TypeId) -> Result<Rc<MapAdams>> {
        if let Some(data) = &self.types[ty.0].map_ext {
            return Ok(data.clone());
        }
        let Kind::MapExt(source, target, s, t) = *self.kind(ty) else {
            unreachable!()
        };
        let data = Rc::new(MapAdams::build(self, source, target, s, t)?);
        self.types[ty.0].map_ext = Some(data.clone());
        Ok(data)
    }

    pub(crate) fn may_page(&mut self, ty: TypeId) -> Result<Rc<MayPage>> {
        if let Some(data) = &self.types[ty.0].may_page {
            return Ok(data.clone());
        }
        let (profile, r, s, t, f) = self.kind(ty).may_grade().unwrap();
        let data = Rc::new(MayPage::build(self, profile, r, s, t, f)?);
        self.types[ty.0].may_page = Some(data.clone());
        Ok(data)
    }

    pub(crate) fn page_transition(&mut self, ty: TypeId) -> Result<Rc<Transition>> {
        if let Some(data) = &self.types[ty.0].page_transition {
            return Ok(data.clone());
        }
        let data = Rc::new(Transition::build(self, ty)?);
        self.types[ty.0].page_transition = Some(data.clone());
        Ok(data)
    }

    pub(crate) fn height(&mut self, ty: TypeId, exact: bool) -> Result<usize> {
        if !self.types[ty.0].inhabited {
            return Ok(0);
        }
        let index = usize::from(
            exact
                && matches!(
                    self.kind(ty),
                    Kind::Product(..) | Kind::Arrow(..) | Kind::Suspension(..)
                ),
        );
        if let Some(height) = self.types[ty.0].heights[index] {
            return Ok(height);
        }
        let height = match *self.kind(ty) {
            Kind::Space(_) | Kind::Star(..) => {
                let space = self.poset(ty).unwrap();
                let size = space.names.len();
                self.spend(2 * size * size)?;
                let mut predecessors = vec![0; size];
                for upper in &space.order {
                    for (point, &above) in upper.iter().enumerate() {
                        predecessors[point] += usize::from(above);
                    }
                }
                // Strictly comparable points have strictly increasing downset sizes.
                let mut points: Vec<_> = (0..size).collect();
                points.sort_unstable_by_key(|&point| predecessors[point]);
                let mut levels = vec![1; size];
                for point in points {
                    for (larger, &above) in space.order[point].iter().enumerate() {
                        if above && point != larger {
                            levels[larger] = levels[larger].max(levels[point] + 1);
                        }
                    }
                }
                levels.into_iter().max().unwrap_or(0)
            }
            Kind::Product(left, right) => {
                (self.height(left, exact)? - 1).saturating_add(self.height(right, exact)?)
            }
            Kind::Arrow(domain, range) => {
                if self.types[domain.0].inhabited {
                    let range_height = self.height(range, exact)?;
                    if range_height == 1 {
                        1
                    } else {
                        let points = if exact {
                            self.elements(domain)?.len()
                        } else {
                            self.types[domain.0].cardinality_bound
                        };
                        // Each coordinate can rise h-1 times. Raising points in reverse
                        // topological order attains the bound using a chain in the range.
                        points.saturating_mul(range_height - 1).saturating_add(1)
                    }
                } else {
                    1
                }
            }
            Kind::Suspension(base) => self.height(base, exact)?.saturating_add(1),
            Kind::Cohomology(..)
            | Kind::Operations(..)
            | Kind::Cobar(..)
            | Kind::PCobar(..)
            | Kind::PExt(..)
            | Kind::CobarAn(..)
            | Kind::Ext(..)
            | Kind::ExtAn(..)
            | Kind::XExt(..)
            | Kind::MapCobar(..)
            | Kind::MapExt(..)
            | Kind::MapMay(..)
            | Kind::May(..)
            | Kind::MayAn(..)
            | Kind::Stable(..)
            | Kind::Bar(..)
            | Kind::CyclicH2(..) => 1,
        };
        self.types[ty.0].heights[index] = Some(height);
        Ok(height)
    }

    fn composite(
        &mut self,
        kind: Kind,
        left: TypeId,
        right: TypeId,
        inhabited: bool,
    ) -> Result<TypeId> {
        if let Some(&id) = self.type_ids.get(&kind) {
            return Ok(id);
        }
        let depth = self.types[left.0].depth.max(self.types[right.0].depth) + 1;
        self.check_type_limit(depth)?;
        Ok(self.insert(kind, None, None, inhabited, depth))
    }

    pub(crate) fn poset(&self, id: TypeId) -> Option<Rc<Poset>> {
        self.types[id.0].poset.clone()
    }

    pub(crate) fn point(&self, id: TypeId, label: &Label) -> Result<usize> {
        let name = || self.type_name(id);
        let Some(poset) = &self.types[id.0].poset else {
            return Err(Error::new(
                label.span,
                format!(
                    "point labels require a named or Star space, found {}",
                    name()
                ),
            ));
        };
        poset.lookup.get(&label.text).copied().ok_or_else(|| {
            Error::new(
                label.span,
                format!("unknown point '{}' in {}", label.text, name()),
            )
        })
    }

    pub(crate) fn star(&mut self, base: TypeId, label: &Label) -> Result<TypeId> {
        let center = self.point(base, label)?;
        let kind = Kind::Star(base, center);
        if let Some(&id) = self.type_ids.get(&kind) {
            return Ok(id);
        }
        let depth = self.types[base.0].depth + 1;
        self.check_type_limit(depth)?;
        let base_poset = self.poset(base).unwrap();
        let points: Vec<usize> = base_poset.order[center]
            .iter()
            .enumerate()
            .filter_map(|(index, &in_star)| in_star.then_some(index))
            .collect();
        self.allocate(points.len() * points.len())?;
        self.spend(points.len() * points.len())?;
        let names = points
            .iter()
            .map(|&point| base_poset.names[point].clone())
            .collect();
        let order = points
            .iter()
            .map(|&left| {
                points
                    .iter()
                    .map(|&right| base_poset.order[left][right])
                    .collect()
            })
            .collect();
        Ok(self.insert(
            kind,
            Some(Poset::new(names, order)),
            Some(points.into()),
            true,
            depth,
        ))
    }

    pub(crate) fn parent_points(&self, ty: TypeId) -> Rc<[usize]> {
        self.types[ty.0].parent_points.as_ref().unwrap().clone()
    }

    pub(crate) fn type_name(&self, ty: TypeId) -> String {
        let mut output = Text::new(2048);
        match self.write_type(&mut output, ty, 0) {
            Ok(()) => output.string,
            Err(_) => "<oversized type>".into(),
        }
    }

    fn write_type(&self, out: &mut Text, ty: TypeId, precedence: u8) -> Result<()> {
        match self.kind(ty) {
            Kind::Space(name) => out.push(name),
            &Kind::Star(base, center) => {
                out.push("Star(")?;
                self.write_type(out, base, 0)?;
                out.push(", ")?;
                out.push(&self.types[base.0].poset.as_ref().unwrap().names[center])?;
                out.push(")")
            }
            &Kind::Cohomology(base, degree) => {
                out.push("H(")?;
                self.write_type(out, base, 0)?;
                out.push(&format!(", {degree})"))
            }
            &Kind::Operations(degree) => out.push(&format!("A({degree})")),
            &Kind::Cobar(s, t) => out.push(&format!("Cobar({s}, {t})")),
            &Kind::PCobar(p, s, t) => out.push(&format!("PCobar({}, {s}, {t})", p.get())),
            &Kind::PExt(p, s, t) => out.push(&format!("PExt({}, {s}, {t})", p.get())),
            &Kind::CobarAn(n, s, t) => out.push(&format!("CobarAn({n}, {s}, {t})")),
            &Kind::Ext(s, t) => out.push(&format!("Ext({s}, {t})")),
            &Kind::ExtAn(n, s, t) => out.push(&format!("ExtAn({n}, {s}, {t})")),
            &Kind::XExt(base, s, t) => {
                out.push("XExt(")?;
                self.write_type(out, base, 0)?;
                out.push(&format!(", {s}, {t})"))
            }
            &Kind::MapMay(source, target, r, s, t, f) => {
                out.push("MapMay(")?;
                self.write_type(out, source, 0)?;
                out.push(", ")?;
                self.write_type(out, target, 0)?;
                out.push(&format!(", {r}, {s}, {t}, {f})"))
            }
            &Kind::MapExt(source, target, s, t) | &Kind::MapCobar(source, target, s, t) => {
                out.push(if matches!(self.kind(ty), Kind::MapCobar(..)) {
                    "MapCobar("
                } else {
                    "MapExt("
                })?;
                self.write_type(out, source, 0)?;
                out.push(", ")?;
                self.write_type(out, target, 0)?;
                out.push(&format!(", {s}, {t})"))
            }
            &Kind::May(r, s, t, f) => out.push(&format!("May({r}, {s}, {t}, {f})")),
            &Kind::MayAn(n, r, s, t, f) => out.push(&format!("MayAn({n}, {r}, {s}, {t}, {f})")),
            &Kind::Bar(modulus) => out.push(&format!("Bar({modulus})")),
            &Kind::CyclicH2(modulus) => out.push(&format!("H2({modulus})")),
            &Kind::Suspension(base) => {
                out.push("Susp(")?;
                self.write_type(out, base, 0)?;
                out.push(")")
            }
            &Kind::Stable(base, shift, degree) => {
                out.push("Stable(")?;
                self.write_type(out, base, 0)?;
                out.push(&format!(", {shift}, {degree})"))
            }
            &Kind::Product(left, right) => {
                let parens = precedence > 1;
                if parens {
                    out.push("(")?;
                }
                self.write_type(out, left, 1)?;
                out.push(" * ")?;
                self.write_type(out, right, 2)?;
                if parens {
                    out.push(")")?;
                }
                Ok(())
            }
            &Kind::Arrow(domain, range) => {
                let parens = precedence > 0;
                if parens {
                    out.push("(")?;
                }
                self.write_type(out, domain, 1)?;
                out.push(" -> ")?;
                self.write_type(out, range, 0)?;
                if parens {
                    out.push(")")?;
                }
                Ok(())
            }
        }
    }

    pub(crate) fn leq(&mut self, ty: TypeId, left: &Value, right: &Value) -> Result<bool> {
        self.spend(1)?;
        match (left, right) {
            (&Value::Point(left), &Value::Point(right)) => {
                Ok(self.types[ty.0].poset.as_ref().unwrap().order[left][right])
            }
            (Value::Class(left), Value::Class(right)) | (Value::Bar(left), Value::Bar(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(true);
                }
                if left.len() != right.len() {
                    return Ok(false);
                }
                self.spend(left.len())?;
                Ok(left == right)
            }
            (Value::CyclicClass(left), Value::CyclicClass(right)) => Ok(left == right),
            (Value::Operation(left), Value::Operation(right)) => {
                Ok(Rc::ptr_eq(left, right)
                    || steenrod::compare(self, left, right)? == Ordering::Equal)
            }
            (Value::Cobar(left), Value::Cobar(right)) => Ok(
                Rc::ptr_eq(left, right) || cobar::compare(self, left, right)? == Ordering::Equal
            ),
            (Value::OddClass(left), Value::OddClass(right)) => {
                Ok(Rc::ptr_eq(left, right)
                    || odd_ext::compare(self, left, right)? == Ordering::Equal)
            }
            (Value::PCobar(left), Value::PCobar(right)) => Ok(Rc::ptr_eq(left, right)
                || odd_cobar::compare(self, left, right)? == Ordering::Equal),
            (Value::MapCobar(left), Value::MapCobar(right)) => Ok(Rc::ptr_eq(left, right)
                || map_cochain::compare(self, left, right)? == Ordering::Equal),
            (Value::Suspension(left), Value::Suspension(right)) => {
                use SuspensionPoint::{Equator, North, South};
                match (left, right) {
                    (Equator(left), Equator(right)) => {
                        let Kind::Suspension(base) = *self.kind(ty) else {
                            unreachable!()
                        };
                        Ok(Rc::ptr_eq(left, right) || self.leq(base, left, right)?)
                    }
                    (Equator(_), _) | (North, North) | (South, South) => Ok(true),
                    _ => Ok(false),
                }
            }
            (Value::Pair(left), Value::Pair(right)) => {
                let Kind::Product(a, b) = *self.kind(ty) else {
                    unreachable!()
                };
                Ok(Rc::ptr_eq(left, right)
                    || (self.leq(a, &left.0, &right.0)? && self.leq(b, &left.1, &right.1)?))
            }
            (Value::Map(left), Value::Map(right)) => {
                let Kind::Arrow(_, range) = *self.kind(ty) else {
                    unreachable!()
                };
                if !Rc::ptr_eq(left, right) {
                    for (left, right) in left.iter().zip(right.iter()) {
                        if !self.leq(range, left, right)? {
                            return Ok(false);
                        }
                    }
                }
                Ok(true)
            }
            _ => unreachable!(),
        }
    }

    fn compare(&mut self, left: &Value, right: &Value) -> Result<Ordering> {
        self.spend(1)?;
        match (left, right) {
            (Value::Point(left), Value::Point(right)) => Ok(left.cmp(right)),
            (Value::Class(left), Value::Class(right)) | (Value::Bar(left), Value::Bar(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(Ordering::Equal);
                }
                self.spend(left.len().min(right.len()))?;
                // Compare from the most significant set bit to match mask enumeration.
                Ok(left.iter().rev().cmp(right.iter().rev()))
            }
            (Value::CyclicClass(left), Value::CyclicClass(right)) => Ok(left.cmp(right)),
            (Value::Operation(left), Value::Operation(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    steenrod::compare(self, left, right)
                }
            }
            (Value::Cobar(left), Value::Cobar(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    cobar::compare(self, left, right)
                }
            }
            (Value::OddClass(left), Value::OddClass(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(Ordering::Equal);
                }
                odd_ext::compare(self, left, right)
            }
            (Value::PCobar(left), Value::PCobar(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(Ordering::Equal);
                }
                odd_cobar::compare(self, left, right)
            }
            (Value::MapCobar(left), Value::MapCobar(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    map_cochain::compare(self, left, right)
                }
            }
            (Value::Suspension(left), Value::Suspension(right)) => {
                use SuspensionPoint::{Equator, North, South};
                match (left, right) {
                    (Equator(left), Equator(right)) => {
                        if Rc::ptr_eq(left, right) {
                            Ok(Ordering::Equal)
                        } else {
                            self.compare(left, right)
                        }
                    }
                    (Equator(_), _) | (North, South) => Ok(Ordering::Less),
                    (_, Equator(_)) | (South, North) => Ok(Ordering::Greater),
                    (North, North) | (South, South) => Ok(Ordering::Equal),
                }
            }
            (Value::Pair(left), Value::Pair(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(Ordering::Equal);
                }
                let first = self.compare(&left.0, &right.0)?;
                if first == Ordering::Equal {
                    self.compare(&left.1, &right.1)
                } else {
                    Ok(first)
                }
            }
            (Value::Map(left), Value::Map(right)) => {
                if !Rc::ptr_eq(left, right) {
                    for (left, right) in left.iter().zip(right.iter()) {
                        let comparison = self.compare(left, right)?;
                        if comparison != Ordering::Equal {
                            return Ok(comparison);
                        }
                    }
                }
                Ok(Ordering::Equal)
            }
            _ => unreachable!(),
        }
    }

    pub(crate) fn equal(&mut self, left: &Value, right: &Value) -> Result<bool> {
        Ok(self.compare(left, right)? == Ordering::Equal)
    }

    pub(crate) fn index(&mut self, ty: TypeId, value: &Value) -> Result<usize> {
        if let Value::Point(index) = *value {
            return Ok(index);
        }
        let values = self.elements(ty)?;
        let mut start = 0;
        let mut end = values.len();
        while start < end {
            let mid = start + (end - start) / 2;
            match self.compare(&values[mid], value)? {
                Ordering::Less => start = mid + 1,
                Ordering::Greater => end = mid,
                Ordering::Equal => return Ok(mid),
            }
        }
        Err(Error::plain("value has no index in its finite type"))
    }

    pub(crate) fn elements(&mut self, ty: TypeId) -> Result<Rc<[Value]>> {
        let info = &self.types[ty.0];
        if let Some(elements) = &info.elements {
            return Ok(elements.clone());
        }
        let elements = if !info.inhabited {
            Vec::new()
        } else {
            match *self.kind(ty) {
                Kind::Space(_) | Kind::Star(_, _) => {
                    let size = self.types[ty.0].poset.as_ref().unwrap().names.len();
                    self.allocate(size)?;
                    (0..size).map(Value::Point).collect()
                }
                Kind::Product(left, right) => {
                    let left = self.elements(left)?;
                    let right = self.elements(right)?;
                    let size = left.len() * right.len();
                    self.check_element_limit(size)?;
                    self.allocate(3 * size)?;
                    self.spend(size)?;
                    let mut points = Vec::with_capacity(size);
                    for a in left.iter() {
                        for b in right.iter() {
                            points.push(Value::Pair(Rc::new((a.clone(), b.clone()))));
                        }
                    }
                    points
                }
                Kind::Arrow(domain, range) => self.enumerate_maps(domain, range)?,
                Kind::Cohomology(..) => {
                    let rank = self.cohomology(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::Ext(..) | Kind::ExtAn(..) => {
                    let rank = self.adams(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::XExt(..) => {
                    let rank = self.xext(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::MapMay(..) => {
                    let rank = self.map_may(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::MapExt(..) => {
                    let rank = self.map_ext(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::May(..) | Kind::MayAn(..) => {
                    let rank = self.may_page(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::Operations(degree) => steenrod::elements(self, degree, MAX_ELEMENTS)?,
                Kind::Cobar(..) | Kind::CobarAn(..) => {
                    let (profile, s, t) = self.kind(ty).cobar_grade().unwrap();
                    let basis = profile.basis(self, s, t, MAX_ELEMENTS.ilog2() as usize)?;
                    self.enumerate_supports(0, basis.len(), |support| {
                        Value::Cobar(support.iter().map(|&index| basis[index].clone()).collect())
                    })?
                }
                Kind::PCobar(p, s, t) => odd_cobar::elements(self, p, s, t, MAX_ELEMENTS)?,
                Kind::PExt(p, ..) => odd_ext::elements(self, ty, p, MAX_ELEMENTS)?,
                Kind::MapCobar(source, target, s, t) => {
                    let basis =
                        map_ext::basis(self, source, target, s, t, MAX_ELEMENTS.ilog2() as usize)?;
                    let rank = basis.len();
                    // In addition to the shared word references, every
                    // support entry stores two degree/index pairs.
                    self.allocate(
                        map_ext::Term::GENERATOR_CELLS * (support_space_size(rank) / 2) * rank,
                    )?;
                    self.enumerate_supports(0, rank, |support| {
                        Value::MapCobar(support.iter().map(|&index| basis[index].clone()).collect())
                    })?
                }
                Kind::Bar(modulus) => self.enumerate_supports(1, modulus - 1, Value::Bar)?,
                Kind::CyclicH2(modulus) => {
                    let count = 2 - (modulus & 1);
                    self.allocate(count)?;
                    (0..count)
                        .map(|index| Value::CyclicClass(index == 1))
                        .collect()
                }
                Kind::Suspension(base) => {
                    let equator = self.elements(base)?;
                    let size = equator.len() + 2;
                    self.check_element_limit(size)?;
                    self.allocate(size + equator.len())?;
                    self.spend(size)?;
                    let mut points = Vec::with_capacity(size);
                    points.extend(equator.iter().map(|point| {
                        Value::Suspension(SuspensionPoint::Equator(Rc::new(point.clone())))
                    }));
                    points.push(Value::Suspension(SuspensionPoint::North));
                    points.push(Value::Suspension(SuspensionPoint::South));
                    points
                }
                Kind::Stable(base, shift, degree) => {
                    let rank = if let Some(reduced_degree) = degree.checked_sub(shift) {
                        let underlying = self.cohomology_type(base, reduced_degree)?;
                        self.cohomology(underlying)?.rank()
                    } else {
                        0
                    };
                    self.enumerate_supports(0, rank, Value::Class)?
                }
            }
        };
        let elements: Rc<[Value]> = elements.into();
        self.types[ty.0].elements = Some(elements.clone());
        Ok(elements)
    }

    fn check_element_limit(&self, size: usize) -> Result<()> {
        if size > MAX_ELEMENTS {
            return Err(Error::plain(
                "type exceeds the enumeration limit of 4096 elements",
            ));
        }
        Ok(())
    }

    fn enumerate_supports(
        &mut self,
        first: usize,
        rank: usize,
        value: impl Fn(Rc<[usize]>) -> Value,
    ) -> Result<Vec<Value>> {
        let count = support_space_size(rank);
        self.check_element_limit(count)?;
        self.spend(count * (rank + 1))?;
        self.allocate(count + (count / 2) * rank)?;
        Ok((0..count)
            .map(|mask| {
                value(
                    (0..rank)
                        .filter(|&i| ((mask >> i) & 1) != 0)
                        .map(|i| i + first)
                        .collect(),
                )
            })
            .collect())
    }

    fn enumerate_maps(&mut self, domain: TypeId, range: TypeId) -> Result<Vec<Value>> {
        let inputs = self.elements(domain)?;
        if inputs.is_empty() {
            self.allocate(1)?;
            return Ok(vec![Value::Map(Rc::new([]))]);
        }
        let outputs = self.elements(range)?;
        if outputs.len() == 1 {
            self.allocate(1 + inputs.len())?;
            return Ok(vec![Value::Map(
                vec![outputs[0].clone(); inputs.len()].into(),
            )]);
        }

        // Rows follow the canonical order of the domain's elements. For a new
        // row, only inequalities against earlier rows need to be checked.
        let mut constraints = vec![Vec::<(usize, bool)>::new(); inputs.len()];
        let mut constraint_count = 0;
        for (current, predecessors) in constraints.iter_mut().enumerate() {
            for previous in 0..current {
                if self.leq(domain, &inputs[previous], &inputs[current])? {
                    predecessors.push((previous, true));
                    constraint_count += 1;
                } else if self.leq(domain, &inputs[current], &inputs[previous])? {
                    predecessors.push((previous, false));
                    constraint_count += 1;
                }
                if constraint_count > MAX_STORAGE {
                    return Err(Error::plain(
                        "map enumeration exceeds the relation storage limit",
                    ));
                }
            }
        }

        // A nonrecursive depth-first traversal keeps domain size off the stack.
        let mut maps = Vec::new();
        let mut row = 0;
        let mut chosen = vec![0; inputs.len()];
        let mut next = vec![0; inputs.len()];
        loop {
            self.spend(1)?;
            if next[row] == outputs.len() {
                next[row] = 0;
                if row == 0 {
                    break;
                }
                row -= 1;
                continue;
            }
            let candidate = next[row];
            next[row] += 1;
            let mut allowed = true;
            for &(earlier, earlier_leq_current) in &constraints[row] {
                let before = &outputs[chosen[earlier]];
                let current = &outputs[candidate];
                let (left, right) = if earlier_leq_current {
                    (before, current)
                } else {
                    (current, before)
                };
                if !self.leq(range, left, right)? {
                    allowed = false;
                    break;
                }
            }
            if allowed {
                chosen[row] = candidate;
                if row + 1 < inputs.len() {
                    row += 1;
                } else {
                    self.check_element_limit(maps.len() + 1)?;
                    self.allocate(1 + inputs.len())?;
                    self.spend(inputs.len())?;
                    let rows = chosen.iter().map(|&index| outputs[index].clone()).collect();
                    maps.push(Value::Map(rows));
                }
            }
        }
        Ok(maps)
    }

    pub(crate) fn show(&mut self, ty: TypeId, value: &Value) -> Result<String> {
        let mut out = Text::new(MAX_SHOWN_BYTES);
        self.write_value(&mut out, ty, value)?;
        Ok(out.string)
    }

    pub(crate) fn describe(&mut self, ty: TypeId, value: &Value) -> String {
        let mut out = Text::new(1024);
        match self.write_value(&mut out, ty, value) {
            Ok(()) => out.string,
            Err(_) => format!("<value of type {}>", self.type_name(ty)),
        }
    }

    fn write_value(&mut self, out: &mut Text, ty: TypeId, value: &Value) -> Result<()> {
        match value {
            Value::Bar(support) => {
                let Kind::Bar(modulus) = *self.kind(ty) else {
                    unreachable!()
                };
                out.push(&format!("bar[{modulus}]{{"))?;
                for (index, residue) in support.iter().enumerate() {
                    if index != 0 {
                        out.push(", ")?;
                    }
                    out.push(&residue.to_string())?;
                }
                out.push("}")
            }
            Value::CyclicClass(nonzero) => {
                self.write_type(out, ty, 0)?;
                out.push(if *nonzero { "{carry}" } else { "{0}" })
            }
            Value::OddClass(support) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if support.is_empty() {
                    out.push("0")?;
                }
                for (index, &(generator, coefficient)) in support.iter().enumerate() {
                    if index != 0 {
                        out.push(" + ")?;
                    }
                    if coefficient != 1 {
                        out.push(&format!("{coefficient}*"))?;
                    }
                    out.push(&format!("g{generator}"))?;
                }
                out.push("}")
            }
            Value::Class(support) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if support.is_empty() {
                    out.push("0")?;
                }
                for (index, generator) in support.iter().enumerate() {
                    if index != 0 {
                        out.push(" + ")?;
                    }
                    out.push(&format!("g{generator}"))?;
                }
                out.push("}")
            }
            Value::Operation(words) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if words.is_empty() {
                    out.push("0")?;
                }
                for (index, word) in words.iter().rev().enumerate() {
                    if index > 0 {
                        out.push(" + ")?;
                    }
                    if word.is_empty() {
                        out.push("1")?;
                    }
                    for (index, degree) in word.iter().enumerate() {
                        if index > 0 {
                            out.push(" ")?;
                        }
                        out.push(&format!("Sq^{degree}"))?;
                    }
                }
                out.push("}")
            }
            Value::Cobar(words) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if words.is_empty() {
                    out.push("0")?;
                }
                for (index, word) in words.iter().rev().enumerate() {
                    if index > 0 {
                        out.push(" + ")?;
                    }
                    write_cobar_word(out, word)?;
                }
                out.push("}")
            }
            Value::PCobar(terms) => {
                let Kind::PCobar(p, ..) = *self.kind(ty) else {
                    unreachable!()
                };
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                for (index, (word, coefficient)) in terms.iter().rev().enumerate() {
                    if index > 0 {
                        out.push(" + ")?;
                    }
                    if *coefficient != 1 {
                        out.push(&format!("{coefficient} "))?;
                    }
                    write_odd_word(out, p, word)?;
                }
                out.push("}")
            }
            Value::MapCobar(terms) => {
                let Kind::MapCobar(source, target, ..) = *self.kind(ty) else {
                    unreachable!()
                };
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                for (index, term) in terms.iter().rev().enumerate() {
                    if index > 0 {
                        out.push(" + ")?;
                    }
                    out.push("mword[")?;
                    self.write_type(out, source, 0)?;
                    out.push(&format!(
                        ", {}, {}, ",
                        term.source.degree, term.source.index
                    ))?;
                    self.write_type(out, target, 0)?;
                    out.push(&format!(
                        ", {}, {}](",
                        term.target.degree, term.target.index
                    ))?;
                    write_cobar_word(out, &term.word)?;
                    out.push(")")?;
                }
                out.push("}")
            }
            Value::Suspension(point) => {
                let Kind::Suspension(base) = *self.kind(ty) else {
                    unreachable!()
                };
                out.push(match point {
                    SuspensionPoint::Equator(_) => "equator[",
                    SuspensionPoint::North => "north[",
                    SuspensionPoint::South => "south[",
                })?;
                self.write_type(out, base, 0)?;
                out.push("]")?;
                if let SuspensionPoint::Equator(point) = point {
                    out.push("(")?;
                    self.write_value(out, base, point)?;
                    out.push(")")?;
                }
                Ok(())
            }
            Value::Point(point) => {
                self.write_type(out, ty, 0)?;
                out.push(".")?;
                out.push(&self.types[ty.0].poset.as_ref().unwrap().names[*point])
            }
            Value::Pair(pair) => {
                let Kind::Product(left, right) = *self.kind(ty) else {
                    unreachable!()
                };
                out.push("(")?;
                self.write_value(out, left, &pair.0)?;
                out.push(", ")?;
                self.write_value(out, right, &pair.1)?;
                out.push(")")
            }
            Value::Map(rows) => {
                let Kind::Arrow(domain, range) = *self.kind(ty) else {
                    unreachable!()
                };
                let inputs = self.elements(domain)?;
                out.push("{")?;
                for (index, (input, image)) in inputs.iter().zip(rows.iter()).enumerate() {
                    if index != 0 {
                        out.push(", ")?;
                    }
                    self.write_value(out, domain, input)?;
                    out.push(" => ")?;
                    self.write_value(out, range, image)?;
                }
                out.push("}")
            }
        }
    }

    pub(crate) fn check_map(
        &mut self,
        domain: TypeId,
        range: TypeId,
        rows: &[Value],
    ) -> Result<()> {
        if rows.len() < 2 || self.intrinsically_discrete(domain) || self.is_constant(rows)? {
            return Ok(());
        }
        let inputs = self.elements(domain)?;
        for (left, x) in inputs.iter().enumerate() {
            for (right, y) in inputs.iter().enumerate() {
                if left != right
                    && self.leq(domain, x, y)?
                    && !self.leq(range, &rows[left], &rows[right])?
                {
                    let x = self.describe(domain, x);
                    let y = self.describe(domain, y);
                    let fx = self.describe(range, &rows[left]);
                    let fy = self.describe(range, &rows[right]);
                    return Err(Error::plain(format!(
                        "nonmonotone map: {x} <= {y}, but {fx} is not <= {fy}"
                    )));
                }
            }
        }
        Ok(())
    }

    fn is_constant(&mut self, rows: &[Value]) -> Result<bool> {
        for other in &rows[1..] {
            if !self.equal(&rows[0], other)? {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

fn write_cobar_word(out: &mut Text, word: &cobar::Word) -> Result<()> {
    out.push("cword{")?;
    for (index, factor) in word.iter().enumerate() {
        if index > 0 {
            out.push(", ")?;
        }
        out.push("[")?;
        for (index, exponent) in factor.iter().enumerate() {
            if index > 0 {
                out.push(", ")?;
            }
            out.push(&exponent.to_string())?;
        }
        out.push("]")?;
    }
    out.push("}")
}

fn write_odd_word(out: &mut Text, p: Prime, word: &odd_cobar::Word) -> Result<()> {
    out.push(&format!("pword[{}]{{", p.get()))?;
    for (index, factor) in word.iter().enumerate() {
        if index > 0 {
            out.push(", ")?;
        }
        out.push("(")?;
        write_number_list(out, &factor.exponents)?;
        out.push(", ")?;
        write_number_list(out, &factor.taus)?;
        out.push(")")?;
    }
    out.push("}")
}

fn write_number_list(out: &mut Text, numbers: &[usize]) -> Result<()> {
    out.push("[")?;
    for (index, number) in numbers.iter().enumerate() {
        if index > 0 {
            out.push(", ")?;
        }
        out.push(&number.to_string())?;
    }
    out.push("]")
}

fn support_space_size(rank: usize) -> usize {
    if rank < usize::BITS as usize {
        1usize << rank
    } else {
        usize::MAX
    }
}

struct Text {
    string: String,
    limit: usize,
}

impl Text {
    fn new(limit: usize) -> Self {
        Self {
            string: String::new(),
            limit,
        }
    }

    fn push(&mut self, text: &str) -> Result<()> {
        if self.string.len() + text.len() > self.limit {
            return Err(Error::plain("value exceeds the display size limit"));
        }
        self.string.push_str(text);
        Ok(())
    }
}
