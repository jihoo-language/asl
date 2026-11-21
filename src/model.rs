use std::cmp::Ordering;
use std::collections::HashMap;
use std::rc::Rc;

use crate::adams::Adams;
use crate::binary_field::BinaryField;
use crate::binary_linear;
use crate::bp_classes;
use crate::bp_height;
use crate::bp_page;
use crate::bp_rel;
use crate::bp_triple;
use crate::ce::CePage;
use crate::cobar;
use crate::cohomology::Cohomology;
use crate::comodule;
use crate::error::{Error, Result};
use crate::jet_hom;
use crate::jet_page;
use crate::map_cochain;
use crate::map_ext::{self, MapAdams};
use crate::map_may::MapMayPage;
use crate::may::MayPage;
use crate::milnor::MilnorBasis;
use crate::odd_cobar::{self, Prime};
use crate::odd_ext::{self, OddExt};
use crate::odd_may::{Filtration, OddPage};
use crate::page_transition::Transition;
use crate::prime_linear::Coordinate;
use crate::profile::Profile;
use crate::q_koszul::{self, Grade, QExt};
use crate::q_module::{Module, Presentation};
use crate::rank_jets::Path;
use crate::rv_as_cover::Cover as RvAsCover;
use crate::rv_bar;
use crate::rv_bar_dual;
use crate::rv_bar_page;
use crate::rv_barcobar;
use crate::rv_curve::Curve as RvCurve;
use crate::rv_hom::{self, Hom};
use crate::rv_ideal::Ideal as RvIdeal;
use crate::rv_iwhecke;
use crate::rv_koszul::{self, Grade as RvGrade, Tor};
use crate::rv_laurent;
use crate::rv_poly::{self, Minors, Polynomial};
use crate::rv_witt_cover::Cover as RvWittCover;
use crate::steenrod::{self, Word};
use crate::syntax::Label;
use crate::v_bars;
use crate::v_koszul;
use crate::vbock::{self, BocksteinPage};
use crate::{bp, bp_polynomial};

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
    QModule(String),
    Star(TypeId, usize),
    Product(TypeId, TypeId),
    Arrow(TypeId, TypeId),
    Cohomology(TypeId, usize),
    Operations(usize),
    Cobar(usize, usize),
    PCobar(Prime, usize, usize),
    PCobarAn(Prime, usize, usize, usize),
    PExt(Prime, usize, usize),
    PExtAn(Prime, usize, usize, usize),
    PMay(Prime, usize, usize, usize, usize),
    PMayAn(Prime, usize, usize, usize, usize, usize),
    PCeAn(Prime, usize, usize, usize, usize, usize),
    CobarAn(usize, usize, usize),
    Ext(usize, usize),
    ExtAn(usize, usize, usize),
    XExt(TypeId, usize, usize),
    MapCobar(TypeId, TypeId, usize, i64),
    MapExt(TypeId, TypeId, usize, i64),
    QKoszul(Grade),
    QExt(Grade),
    RVHom(rv_hom::Grade),
    ASDiff(Rc<str>, BinaryField),
    RVLaurent(rv_hom::Grade),
    RVIwHecke(usize),
    JetHom(jet_hom::Grade),
    JetPage(jet_page::Grade),
    VKoszul(v_koszul::Grade),
    VExt(v_koszul::Grade),
    VBock(vbock::Grade),
    RVBarE(rv_bar_page::Grade),
    RVBarCobar(rv_bar::Grade),
    RVBarDual(rv_bar::Grade),
    RVPathExt(rv_bar::Grade),
    RVBar(rv_bar::Grade),
    RVKoszul(RvGrade),
    RVTor(RvGrade),
    RVPoly(q_koszul::Indices, usize),
    BP(bp::Shape),
    BPRel(bp_rel::Shape),
    BPRelExt(bp_rel::Grade),
    BPHeight(bp_height::ClassGrade),
    BPPage(bp_page::Grade),
    BPExt(bp_classes::Grade),
    BPTriple(bp_classes::Grade),
    MapMay(TypeId, TypeId, usize, usize, i64, usize),
    May(usize, usize, usize, usize),
    MayAn(usize, usize, usize, usize, usize),
    CeAn(usize, usize, usize, usize, usize),
    Suspension(TypeId),
    Stable(TypeId, usize, usize),
    Bar(usize),
    CyclicH2(usize),
}

impl Kind {
    pub(crate) fn pce_grade(&self) -> Option<(Prime, usize, usize, usize, usize, usize)> {
        match *self {
            Self::PCeAn(p, n, r, s, t, f) => Some((p, n, r, s, t, f)),
            _ => None,
        }
    }

    pub(crate) fn pmay_grade(&self) -> Option<(Prime, Profile, usize, usize, usize, usize)> {
        match *self {
            Self::PMay(p, r, s, t, f) => Some((p, Profile::Full, r, s, t, f)),
            Self::PMayAn(p, n, r, s, t, f) => Some((p, Profile::An(n), r, s, t, f)),
            _ => None,
        }
    }

    pub(crate) fn pcobar_grade(&self) -> Option<(Prime, Profile, usize, usize)> {
        match *self {
            Self::PCobar(p, s, t) => Some((p, Profile::Full, s, t)),
            Self::PCobarAn(p, n, s, t) => Some((p, Profile::An(n), s, t)),
            _ => None,
        }
    }

    pub(crate) fn pext_grade(&self) -> Option<(Prime, Profile, usize, usize)> {
        match *self {
            Self::PExt(p, s, t) => Some((p, Profile::Full, s, t)),
            Self::PExtAn(p, n, s, t) => Some((p, Profile::An(n), s, t)),
            _ => None,
        }
    }

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
    QKoszul(q_koszul::Support),
    RVHom(rv_hom::Support),
    ASDiff(rv_hom::Support),
    RVLaurent(rv_laurent::Support),
    RVIwHecke(rv_iwhecke::Support),
    JetHom(jet_hom::Support),
    RVBarCobar(rv_barcobar::Support),
    RVBar(rv_bar::Support),
    RVKoszul(rv_koszul::Support),
    JetClass(jet_page::Coordinates),
    RVPoly(Polynomial),
    BP(Rc<bp_polynomial::Polynomial>),
    BPRel(Rc<bp_polynomial::Polynomial>),
    BPExt(bp_classes::Support),
    BPTriple(Rc<bp_triple::Coset>),
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
    bp_ext: Option<Rc<bp_classes::Cohomology>>,
    odd_ext_rank: Option<usize>,
    odd_page: Option<Rc<OddPage>>,
    bp_page: Option<Rc<bp_page::Page>>,
    xext: Option<Rc<comodule::Adams>>,
    map_ext: Option<Rc<MapAdams>>,
    q_ext: Option<Rc<QExt>>,
    rv_hom: Option<Rc<Hom>>,
    q_module: Option<Rc<Module>>,
    vbock: Option<Rc<BocksteinPage>>,
    jet_page: Option<Rc<jet_page::Page>>,
    map_may: Option<Rc<MapMayPage>>,
    may_page: Option<Rc<MayPage>>,
    ce_page: Option<Rc<CePage>>,
    page_transition: Option<Rc<Transition>>,
    cardinality_bound: usize,
    heights: [Option<usize>; 2],
    inhabited: bool,
    depth: usize,
}

type MinorsKey = (TypeId, q_koszul::Indices, usize);

pub(crate) struct Model {
    types: Vec<TypeInfo>,
    type_ids: HashMap<Kind, TypeId>,
    milnor_bases: HashMap<usize, Rc<MilnorBasis>>,
    v_bars: HashMap<v_bars::Family, Rc<v_bars::Decomposition>>,
    rv_ideals: HashMap<Rc<str>, Rc<RvIdeal>>,
    rv_curves: HashMap<String, Rc<RvCurve>>,
    rv_as_covers: HashMap<String, Rc<RvAsCover>>,
    rv_witt_covers: HashMap<String, Rc<RvWittCover>>,
    rv_tor: HashMap<TypeId, Rc<Tor>>,
    minor_bases: HashMap<MinorsKey, Rc<Minors>>,
    rv_bar_paths: Vec<(String, Rc<rv_bar::Path>)>,
    rv_path_ext: HashMap<TypeId, Rc<rv_bar_dual::Ext>>,
    rv_bar_pages: HashMap<TypeId, Rc<rv_bar_page::Page>>,
    rv_bar_ids: HashMap<String, rv_bar::PathId>,
    jet_paths: Vec<(String, Rc<Path>)>,
    jet_path_ids: HashMap<String, jet_hom::PathId>,
    jet_complexes: HashMap<jet_hom::Diagram, Rc<jet_hom::Complex>>,
    work: usize,
    storage: usize,
}

impl Model {
    pub(crate) fn new() -> Self {
        let mut model = Self {
            types: Vec::new(),
            type_ids: HashMap::new(),
            milnor_bases: HashMap::new(),
            v_bars: HashMap::new(),
            rv_ideals: HashMap::new(),
            rv_curves: HashMap::new(),
            rv_as_covers: HashMap::new(),
            rv_witt_covers: HashMap::new(),
            rv_tor: HashMap::new(),
            minor_bases: HashMap::new(),
            rv_bar_paths: Vec::new(),
            rv_path_ext: HashMap::new(),
            rv_bar_pages: HashMap::new(),
            rv_bar_ids: HashMap::new(),
            jet_paths: Vec::new(),
            jet_path_ids: HashMap::new(),
            jet_complexes: HashMap::new(),
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
            Kind::QModule(_) => 1,
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
            Kind::PCobar(_, s, t)
            | Kind::PCobarAn(_, _, s, t)
            | Kind::PExt(_, s, t)
            | Kind::PExtAn(_, _, s, t)
                if s > t || (s == 0 && t != 0) =>
            {
                1
            }
            Kind::PCobar(p, s, t)
            | Kind::PCobarAn(p, _, s, t)
            | Kind::PExt(p, s, t)
            | Kind::PExtAn(p, _, s, t)
                if s == t =>
            {
                p.get() as usize
            }
            Kind::PMay(_, _, s, t, f) | Kind::PMayAn(_, _, _, s, t, f)
                if s > t || f > t - s || (s == 0 && t != 0) =>
            {
                1
            }
            Kind::PMay(prime, _, s, t, _) | Kind::PMayAn(prime, _, _, s, t, _) if s == t => {
                prime.get() as usize
            }
            Kind::PCeAn(_, _, _, s, t, f) if s > t || f > s || (s == 0 && t != 0) => 1,
            Kind::PCeAn(prime, _, _, s, t, _) if s == t => prime.get() as usize,
            Kind::PCobar(..)
            | Kind::PCobarAn(..)
            | Kind::PExt(..)
            | Kind::PExtAn(..)
            | Kind::PMay(..)
            | Kind::PMayAn(..)
            | Kind::PCeAn(..)
            | Kind::BPPage(..) => usize::MAX,
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
            Kind::CeAn(n, _, s, t, f)
                if f > s || s > t || (n == 0 && f != 0) || (s == 0 && t != 0) =>
            {
                1
            }
            Kind::CeAn(_, _, s, t, _) if s == t => 2,
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
            | Kind::MayAn(..)
            | Kind::CeAn(..) => usize::MAX,
            Kind::QKoszul(ref grade) | Kind::QExt(ref grade)
                if !self.types[grade.source.0].inhabited
                    || !self.types[grade.target.0].inhabited
                    || (grade.selected.is_empty() && grade.s != 0) =>
            {
                1
            }
            Kind::QKoszul(..)
            | Kind::QExt(..)
            | Kind::JetHom(..)
            | Kind::JetPage(..)
            | Kind::RVHom(..)
            | Kind::RVLaurent(..)
            | Kind::RVIwHecke(..) => usize::MAX,
            Kind::VKoszul(ref grade) | Kind::VExt(ref grade)
                if !self.types[grade.source.0].inhabited
                    || !self.types[grade.target.0].inhabited =>
            {
                1
            }
            Kind::VKoszul(..) | Kind::VExt(..) | Kind::VBock(..) => usize::MAX,
            Kind::RVKoszul(ref grade) | Kind::RVTor(ref grade)
                if grade.s > grade.t || grade.s > self.rv_ideals[&grade.ideal].selected.len() =>
            {
                1
            }
            Kind::RVBarCobar(..) => usize::MAX,
            Kind::RVBar(..) | Kind::RVBarE(..) | Kind::RVBarDual(..) | Kind::RVPathExt(..) => {
                usize::MAX
            }
            Kind::RVKoszul(..) | Kind::RVTor(..) => usize::MAX,
            Kind::RVPoly(ref selected, d) if selected.is_empty() && d != 0 => 1,
            Kind::RVPoly(ref selected, d) if selected.len() < 2 || d == 0 => 2,
            Kind::ASDiff(ref name, field) => {
                let genus = self.rv_as_covers[name.as_ref()].genus();
                field
                    .order()
                    .saturating_pow(u32::try_from(genus).unwrap_or(u32::MAX))
            }
            Kind::RVPoly(..)
            | Kind::BP(..)
            | Kind::BPRel(..)
            | Kind::BPExt(..)
            | Kind::BPTriple(..)
            | Kind::BPRelExt(..)
            | Kind::BPHeight(..) => usize::MAX,
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
            bp_ext: None,
            odd_ext_rank: None,
            odd_page: None,
            bp_page: None,
            xext: None,
            map_ext: None,
            q_ext: None,
            rv_hom: None,
            q_module: None,
            vbock: None,
            jet_page: None,
            map_may: None,
            may_page: None,
            ce_page: None,
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
                | Kind::PCobarAn(..)
                | Kind::PExt(..)
                | Kind::PExtAn(..)
                | Kind::PMay(..)
                | Kind::PMayAn(..)
                | Kind::PCeAn(..)
                | Kind::BPPage(..)
                | Kind::CobarAn(..)
                | Kind::Ext(..)
                | Kind::ExtAn(..)
                | Kind::XExt(..)
                | Kind::MapCobar(..)
                | Kind::MapExt(..)
                | Kind::QKoszul(..)
                | Kind::JetHom(..)
                | Kind::JetPage(..)
                | Kind::QExt(..)
                | Kind::RVHom(..)
                | Kind::ASDiff(..)
                | Kind::RVLaurent(..)
                | Kind::RVIwHecke(..)
                | Kind::VKoszul(..)
                | Kind::VExt(..)
                | Kind::VBock(..)
                | Kind::RVBarE(..)
                | Kind::RVBarCobar(..)
                | Kind::RVBarDual(..)
                | Kind::RVPathExt(..)
                | Kind::RVBar(..)
                | Kind::RVKoszul(..)
                | Kind::RVTor(..)
                | Kind::RVPoly(..)
                | Kind::BP(..)
                | Kind::BPRel(..)
                | Kind::BPExt(..)
                | Kind::BPTriple(..)
                | Kind::BPRelExt(..)
                | Kind::BPHeight(..)
                | Kind::MapMay(..)
                | Kind::May(..)
                | Kind::MayAn(..)
                | Kind::CeAn(..)
                | Kind::Stable(..)
                | Kind::Bar(..)
                | Kind::CyclicH2(..)
        )
    }

    pub(crate) fn q_endpoint_named(&self, name: &str) -> Option<TypeId> {
        self.type_ids
            .get(&Kind::QModule(name.into()))
            .copied()
            .or_else(|| self.named(name))
    }

    pub(crate) fn q_module(&self, ty: TypeId) -> Option<Rc<Module>> {
        self.types[ty.0].q_module.clone()
    }

    pub(crate) fn define_q_module(
        &mut self,
        name: &Label,
        selected: Vec<usize>,
        presentation: Presentation,
    ) -> Result<()> {
        if self.q_endpoint_named(&name.text).is_some() {
            return Err(Error::new(
                name.span,
                format!("endpoint '{}' is already defined", name.text),
            ));
        }
        self.check_type_limit(0)?;
        let module = Rc::new(Module::build(self, selected, presentation)?);
        let id = self.insert(
            Kind::QModule(name.text.clone()),
            None,
            None,
            module.bound().is_some(),
            0,
        );
        self.types[id.0].q_module = Some(module);
        Ok(())
    }

    pub(crate) fn check_q_endpoints(
        &self,
        source: TypeId,
        target: TypeId,
        selected: &[usize],
    ) -> Result<()> {
        for endpoint in [source, target] {
            if let Some(module) = &self.types[endpoint.0].q_module {
                module.check_selection(selected)?;
            }
        }
        Ok(())
    }

    pub(crate) fn define_space(&mut self, name: &Label, chains: &[Vec<Label>]) -> Result<()> {
        if self.q_endpoint_named(&name.text).is_some() {
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
        self.profile_odd_ext_type(p, Profile::Full, s, t)
    }

    pub(crate) fn profile_odd_ext_type(
        &mut self,
        p: Prime,
        profile: Profile,
        s: usize,
        t: usize,
    ) -> Result<TypeId> {
        let kind = match profile {
            Profile::Full => Kind::PExt(p, s, t),
            Profile::An(n) => Kind::PExtAn(p, n, s, t),
        };
        self.inhabited_type(kind, 0)
    }

    pub(crate) fn pmay_type(
        &mut self,
        prime: Prime,
        r: usize,
        s: usize,
        t: usize,
        f: usize,
    ) -> Result<TypeId> {
        self.profile_pmay_type(prime, Profile::Full, r, s, t, f)
    }

    pub(crate) fn profile_pmay_type(
        &mut self,
        prime: Prime,
        profile: Profile,
        r: usize,
        s: usize,
        t: usize,
        f: usize,
    ) -> Result<TypeId> {
        let kind = match profile {
            Profile::Full => Kind::PMay(prime, r, s, t, f),
            Profile::An(n) => Kind::PMayAn(prime, n, r, s, t, f),
        };
        self.inhabited_type(kind, 0)
    }

    pub(crate) fn pce_type(
        &mut self,
        prime: Prime,
        n: usize,
        r: usize,
        s: usize,
        t: usize,
        f: usize,
    ) -> Result<TypeId> {
        self.inhabited_type(Kind::PCeAn(prime, n, r, s, t, f), 0)
    }

    pub(crate) fn pcobar_type(&mut self, p: Prime, s: usize, t: usize) -> Result<TypeId> {
        self.profile_pcobar_type(p, Profile::Full, s, t)
    }

    pub(crate) fn profile_pcobar_type(
        &mut self,
        p: Prime,
        profile: Profile,
        s: usize,
        t: usize,
    ) -> Result<TypeId> {
        let kind = match profile {
            Profile::Full => Kind::PCobar(p, s, t),
            Profile::An(n) => Kind::PCobarAn(p, n, s, t),
        };
        self.inhabited_type(kind, 0)
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

    pub(crate) fn define_rv_bar_path(
        &mut self,
        name: Label,
        ideal: Label,
        steps: Vec<Vec<usize>>,
    ) -> Result<()> {
        if self.rv_bar_ids.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("bar path '{}' is already defined", name.text),
            ));
        }
        let path = rv_bar::Path::build(self, ideal, steps)?;
        self.allocate(2 + 2 * name.text.len())?;
        let id = rv_bar::PathId(self.rv_bar_paths.len());
        self.rv_bar_ids.insert(name.text.clone(), id);
        self.rv_bar_paths.push((name.text, Rc::new(path)));
        Ok(())
    }

    pub(crate) fn rv_bar_path_id(&self, name: &Label) -> Result<rv_bar::PathId> {
        self.rv_bar_ids
            .get(&name.text)
            .copied()
            .ok_or_else(|| Error::new(name.span, format!("unknown bar path '{}'", name.text)))
    }

    pub(crate) fn rv_bar_path(&self, id: rv_bar::PathId) -> Rc<rv_bar::Path> {
        self.rv_bar_paths[id.0].1.clone()
    }

    pub(crate) fn rv_bar_path_name(&self, id: rv_bar::PathId) -> &str {
        &self.rv_bar_paths[id.0].0
    }

    pub(crate) fn rv_bar_type(&mut self, grade: rv_bar::Grade) -> Result<TypeId> {
        self.rv_bar_path(grade.path)
            .check_interval(grade.start, grade.end)?;
        self.inhabited_type(Kind::RVBar(grade), 0)
    }

    pub(crate) fn rv_barcobar_type(&mut self, grade: rv_bar::Grade) -> Result<TypeId> {
        self.rv_bar_path(grade.path)
            .check_interval(grade.start, grade.end)?;
        self.inhabited_type(Kind::RVBarCobar(grade), 0)
    }

    pub(crate) fn rv_path_type(&mut self, grade: rv_bar::Grade, homology: bool) -> Result<TypeId> {
        self.rv_bar_path(grade.path)
            .check_interval(grade.start, grade.end)?;
        let kind = if homology {
            Kind::RVPathExt(grade)
        } else {
            Kind::RVBarDual(grade)
        };
        self.inhabited_type(kind, 0)
    }

    pub(crate) fn rv_path_ext(&mut self, ty: TypeId) -> Result<Rc<rv_bar_dual::Ext>> {
        if let Some(ext) = self.rv_path_ext.get(&ty) {
            return Ok(ext.clone());
        }
        let Kind::RVPathExt(grade) = *self.kind(ty) else {
            unreachable!()
        };
        let ext = Rc::new(rv_bar_dual::Ext::build(self, grade)?);
        self.rv_path_ext.insert(ty, ext.clone());
        Ok(ext)
    }

    pub(crate) fn rv_bar_page_type(&mut self, grade: rv_bar_page::Grade) -> Result<TypeId> {
        self.rv_bar_path(grade.chain.path)
            .check_interval(grade.chain.start, grade.chain.end)?;
        self.inhabited_type(Kind::RVBarE(grade), 0)
    }

    pub(crate) fn rv_bar_page(&mut self, ty: TypeId) -> Result<Rc<rv_bar_page::Page>> {
        if let Some(page) = self.rv_bar_pages.get(&ty) {
            return Ok(page.clone());
        }
        let Kind::RVBarE(grade) = *self.kind(ty) else {
            unreachable!();
        };
        let page = Rc::new(rv_bar_page::Page::build(self, grade)?);
        self.rv_bar_pages.insert(ty, page.clone());
        Ok(page)
    }

    pub(crate) fn define_jet_path(&mut self, name: Label, path: Path) -> Result<()> {
        if self.jet_path_ids.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Margolis arc '{}' is already defined", name.text),
            ));
        }
        self.allocate(2 + path.cells() + 2 * name.text.len())?;
        let id = jet_hom::PathId(self.jet_paths.len());
        self.jet_path_ids.insert(name.text.clone(), id);
        self.jet_paths.push((name.text, Rc::new(path)));
        Ok(())
    }

    pub(crate) fn jet_path_id(&self, label: &Label) -> Result<jet_hom::PathId> {
        self.jet_path_ids
            .get(&label.text)
            .copied()
            .ok_or_else(|| Error::new(label.span, format!("unknown Margolis arc '{}'", label.text)))
    }

    pub(crate) fn jet_path(&self, id: jet_hom::PathId) -> Rc<Path> {
        self.jet_paths[id.0].1.clone()
    }

    pub(crate) fn jet_complex(
        &mut self,
        diagram: jet_hom::Diagram,
    ) -> Result<Rc<jet_hom::Complex>> {
        if let Some(complex) = self.jet_complexes.get(&diagram) {
            return Ok(complex.clone());
        }
        let complex = Rc::new(jet_hom::Complex::build(self, diagram)?);
        self.jet_complexes.insert(diagram, complex.clone());
        Ok(complex)
    }

    pub(crate) fn jet_page_type(&mut self, grade: jet_page::Grade) -> Result<TypeId> {
        let chain = self.jet_hom_type(grade.chain)?;
        self.inhabited_type(Kind::JetPage(grade), self.types[chain.0].depth)
    }

    pub(crate) fn jet_page(&mut self, ty: TypeId) -> Result<Rc<jet_page::Page>> {
        if let Some(page) = &self.types[ty.0].jet_page {
            return Ok(page.clone());
        }
        let Kind::JetPage(grade) = self.types[ty.0].kind else {
            unreachable!()
        };
        let page = Rc::new(jet_page::Page::build(self, grade)?);
        self.types[ty.0].jet_page = Some(page.clone());
        Ok(page)
    }

    pub(crate) fn jet_hom_type(&mut self, grade: jet_hom::Grade) -> Result<TypeId> {
        if grade.parity > 1 {
            return Err(Error::plain("JetHom parity must be zero or one"));
        }
        let diagram = grade.diagram;
        let path = self.jet_path(diagram.path);
        self.check_q_endpoints(diagram.source, diagram.target, &path.selected)?;
        let depth = 1 + self.types[diagram.source.0]
            .depth
            .max(self.types[diagram.target.0].depth);
        self.inhabited_type(Kind::JetHom(grade), depth)
    }

    pub(crate) fn q_koszul_type(&mut self, grade: Grade) -> Result<TypeId> {
        self.check_q_endpoints(grade.source, grade.target, &grade.selected)?;
        let depth = self.koszul_depth(&grade);
        self.inhabited_type(Kind::QKoszul(grade), depth)
    }

    pub(crate) fn rv_hom_type(&mut self, grade: rv_hom::Grade) -> Result<TypeId> {
        let depth = self.types[grade.source.0]
            .depth
            .max(self.types[grade.target.0].depth)
            + 1;
        self.inhabited_type(Kind::RVHom(grade), depth)
    }

    pub(crate) fn iw_hecke_type(&mut self, rank: usize) -> Result<TypeId> {
        self.inhabited_type(Kind::RVIwHecke(rank), 0)
    }

    pub(crate) fn rv_laurent_type(&mut self, grade: rv_hom::Grade) -> Result<TypeId> {
        let depth = self.types[grade.source.0]
            .depth
            .max(self.types[grade.target.0].depth)
            + 1;
        self.inhabited_type(Kind::RVLaurent(grade), depth)
    }

    pub(crate) fn q_ext_type(&mut self, grade: Grade) -> Result<TypeId> {
        self.check_q_endpoints(grade.source, grade.target, &grade.selected)?;
        let depth = self.koszul_depth(&grade);
        self.inhabited_type(Kind::QExt(grade), depth)
    }

    pub(crate) fn v_koszul_type(&mut self, grade: v_koszul::Grade) -> Result<TypeId> {
        self.check_q_endpoints(grade.source, grade.target, &grade.selected)?;
        v_koszul::top(&grade.selected)?;
        let depth = self.koszul_depth(&grade);
        self.inhabited_type(Kind::VKoszul(grade), depth)
    }

    pub(crate) fn v_ext_type(&mut self, grade: v_koszul::Grade) -> Result<TypeId> {
        self.check_q_endpoints(grade.source, grade.target, &grade.selected)?;
        v_koszul::top(&grade.selected)?;
        let depth = self.koszul_depth(&grade);
        self.inhabited_type(Kind::VExt(grade), depth)
    }

    pub(crate) fn bp_page_type(&mut self, grade: bp_page::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::BPPage(grade), 0)
    }

    pub(crate) fn bp_triple_type(&mut self, grade: bp_classes::Grade) -> Result<TypeId> {
        if grade.shape.arrows < 2 {
            return Err(Error::plain(
                "BPTriple requires target cobar degree at least two",
            ));
        }
        self.inhabited_type(Kind::BPTriple(grade), 0)
    }

    pub(crate) fn bp_ext_type(&mut self, grade: bp_classes::Grade) -> Result<TypeId> {
        grade.shape.require_integral("BPExt")?;
        self.inhabited_type(Kind::BPExt(grade), 0)
    }

    pub(crate) fn bp_height_type(&mut self, grade: bp_height::ClassGrade) -> Result<TypeId> {
        self.inhabited_type(Kind::BPHeight(grade.validate()?), 0)
    }

    pub(crate) fn bp_rel_ext_type(&mut self, grade: bp_rel::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::BPRelExt(grade), 0)
    }

    pub(crate) fn bp_ext(&mut self, ty: TypeId) -> Result<Rc<bp_classes::Cohomology>> {
        if let Some(group) = &self.types[ty.0].bp_ext {
            return Ok(group.clone());
        }
        let group = Rc::new(match self.types[ty.0].kind {
            Kind::BPExt(grade) => bp_classes::Cohomology::build(self, grade)?,
            Kind::BPRelExt(grade) => bp_classes::Cohomology::build_relative(self, grade)?,
            Kind::BPHeight(grade) => bp_classes::Cohomology::build_height(self, grade)?,
            _ => unreachable!(),
        });
        self.types[ty.0].bp_ext = Some(group.clone());
        Ok(group)
    }

    pub(crate) fn bp_rel_type(&mut self, shape: bp_rel::Shape) -> Result<TypeId> {
        self.inhabited_type(Kind::BPRel(shape), 0)
    }

    pub(crate) fn bp_type(&mut self, shape: bp::Shape) -> Result<TypeId> {
        self.inhabited_type(Kind::BP(shape), 0)
    }

    pub(crate) fn rv_ideal_named(&self, name: &str) -> Option<Rc<RvIdeal>> {
        self.rv_ideals.get(name).cloned()
    }

    pub(crate) fn rv_ideal(&self, name: &Label) -> Result<Rc<RvIdeal>> {
        self.rv_ideal_named(&name.text)
            .ok_or_else(|| Error::new(name.span, format!("unknown rank ideal '{}'", name.text)))
    }

    pub(crate) fn rv_witt_cover(&self, name: &Label) -> Result<Rc<RvWittCover>> {
        self.rv_witt_covers
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown Witt cover '{}'", name.text)))
    }

    pub(crate) fn check_new_rv_witt_cover(&self, name: &Label) -> Result<()> {
        if self.rv_witt_covers.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Witt cover '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn define_rv_witt_cover(&mut self, name: String, cover: RvWittCover) {
        self.rv_witt_covers.insert(name, Rc::new(cover));
    }

    pub(crate) fn rv_as_diff_type(&mut self, name: &Label, modulus: usize) -> Result<TypeId> {
        self.rv_as_cover(name)?;
        let field = BinaryField::new(self, modulus)?;
        self.inhabited_type(Kind::ASDiff(name.text.as_str().into(), field), 0)
    }

    pub(crate) fn rv_as_cover_named(&self, name: &str) -> Option<Rc<RvAsCover>> {
        self.rv_as_covers.get(name).cloned()
    }

    pub(crate) fn rv_as_cover(&self, name: &Label) -> Result<Rc<RvAsCover>> {
        self.rv_as_cover_named(&name.text).ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown Artin–Schreier cover '{}'", name.text),
            )
        })
    }

    pub(crate) fn check_new_rv_as_cover(&self, name: &Label) -> Result<()> {
        if self.rv_as_covers.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Artin–Schreier cover '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn define_rv_as_cover(&mut self, name: String, cover: RvAsCover) {
        self.rv_as_covers.insert(name, Rc::new(cover));
    }

    pub(crate) fn rv_curve(&self, name: &Label) -> Result<Rc<RvCurve>> {
        self.rv_curves
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown rational curve '{}'", name.text)))
    }

    pub(crate) fn check_new_rv_curve(&self, name: &Label) -> Result<()> {
        if self.rv_curves.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("rational curve '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn define_rv_curve(&mut self, name: String, curve: RvCurve) {
        self.rv_curves.insert(name, Rc::new(curve));
    }

    pub(crate) fn define_rv_ideal(&mut self, name: String, ideal: RvIdeal) {
        self.rv_ideals.insert(name.into(), Rc::new(ideal));
    }

    pub(crate) fn rv_grade(&self, name: Label, s: usize, t: usize) -> Result<RvGrade> {
        self.rv_ideal(&name)?;
        Ok(RvGrade::total(name.text.into(), s, t))
    }

    pub(crate) fn rv_complex_type(&mut self, grade: RvGrade, homology: bool) -> Result<TypeId> {
        let kind = if homology {
            Kind::RVTor(grade)
        } else {
            Kind::RVKoszul(grade)
        };
        self.inhabited_type(kind, 0)
    }

    pub(crate) fn rv_tor(&mut self, ty: TypeId) -> Result<Rc<Tor>> {
        if let Some(tor) = self.rv_tor.get(&ty) {
            return Ok(tor.clone());
        }
        let Kind::RVTor(grade) = self.kind(ty) else {
            unreachable!();
        };
        let grade = grade.clone();
        let ideal = self.rv_ideals[&grade.ideal].clone();
        let tor = Rc::new(Tor::build(self, &ideal, &grade)?);
        self.rv_tor.insert(ty, tor.clone());
        Ok(tor)
    }

    pub(crate) fn rv_poly_type(
        &mut self,
        selected: q_koszul::Indices,
        degree: usize,
    ) -> Result<TypeId> {
        self.inhabited_type(Kind::RVPoly(selected, degree), 0)
    }

    pub(crate) fn rv_minors(
        &mut self,
        base: TypeId,
        selected: q_koszul::Indices,
        size: usize,
    ) -> Result<Rc<Minors>> {
        let key = (base, selected, size);
        if let Some(basis) = self.minor_bases.get(&key) {
            return Ok(basis.clone());
        }
        let basis = Rc::new(Minors::build(self, key.0, &key.1, key.2)?);
        self.allocate(3 + basis.storage())?;
        self.minor_bases.insert(key, basis.clone());
        Ok(basis)
    }

    pub(crate) fn show_rv_minors(
        &mut self,
        base: TypeId,
        selected: q_koszul::Indices,
        size: usize,
    ) -> Result<String> {
        let basis = self.rv_minors(base, selected.clone(), size)?;
        let prefix = format!("rvminors[{}, {selected:?}, {size}]", self.type_name(base));
        self.show_rv_polynomials(&prefix, &selected, basis.polynomials())
    }

    pub(crate) fn show_rv_polynomials<'a>(
        &mut self,
        prefix: &str,
        selected: &[usize],
        polynomials: impl Iterator<Item = &'a Polynomial>,
    ) -> Result<String> {
        let mut out = Text::new(MAX_SHOWN_BYTES);
        out.push(prefix)?;
        out.push(" = [")?;
        for (i, polynomial) in polynomials.enumerate() {
            if i != 0 {
                out.push(", ")?;
            }
            self.spend(rv_poly::cells(polynomial))?;
            write_rv_polynomial(&mut out, selected, polynomial)?;
        }
        out.push("]")?;
        Ok(out.string)
    }

    pub(crate) fn vbock_type(&mut self, grade: vbock::Grade) -> Result<TypeId> {
        self.check_q_endpoints(
            grade.localized.source,
            grade.localized.target,
            &grade.localized.selected,
        )?;
        v_koszul::top(&grade.localized.selected)?;
        let depth = self.koszul_depth(&grade.localized);
        self.inhabited_type(Kind::VBock(grade), depth)
    }

    fn koszul_depth<S>(&self, grade: &Grade<S>) -> usize {
        self.types[grade.source.0]
            .depth
            .max(self.types[grade.target.0].depth)
            + 1
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

    pub(crate) fn ce_type(
        &mut self,
        n: usize,
        r: usize,
        s: usize,
        t: usize,
        f: usize,
    ) -> Result<TypeId> {
        self.inhabited_type(Kind::CeAn(n, r, s, t, f), 0)
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

    pub(crate) fn cached_odd_ext(&self, ty: TypeId) -> Option<Rc<OddExt>> {
        self.types[ty.0].odd_ext.clone()
    }

    pub(crate) fn v_bars(&mut self, family: v_bars::Family) -> Result<Rc<v_bars::Decomposition>> {
        if let Some(decomposition) = self.v_bars.get(&family) {
            return Ok(decomposition.clone());
        }
        let decomposition = Rc::new(v_bars::Decomposition::build(self, family)?);
        self.v_bars.insert(family, decomposition.clone());
        Ok(decomposition)
    }

    pub(crate) fn odd_ext_rank(&mut self, ty: TypeId) -> Result<usize> {
        if let Some(group) = &self.types[ty.0].odd_ext {
            return Ok(group.rank());
        }
        if let Some(rank) = self.types[ty.0].odd_ext_rank {
            return Ok(rank);
        }
        let (prime, profile, s, t) = self.kind(ty).pext_grade().unwrap();
        let rank = odd_ext::rank(self, prime, profile, s, t)?;
        self.allocate(1)?;
        self.types[ty.0].odd_ext_rank = Some(rank);
        Ok(rank)
    }

    pub(crate) fn odd_ext(&mut self, ty: TypeId) -> Result<Rc<OddExt>> {
        if let Some(data) = &self.types[ty.0].odd_ext {
            return Ok(data.clone());
        }
        let (prime, profile, s, t) = self.kind(ty).pext_grade().unwrap();
        let data = Rc::new(OddExt::build(self, prime, profile, s, t)?);
        self.types[ty.0].odd_ext = Some(data.clone());
        Ok(data)
    }

    pub(crate) fn bp_page(&mut self, ty: TypeId) -> Result<Rc<bp_page::Page>> {
        if let Some(data) = &self.types[ty.0].bp_page {
            return Ok(data.clone());
        }
        let Kind::BPPage(grade) = *self.kind(ty) else {
            unreachable!()
        };
        let data = Rc::new(bp_page::Page::build(self, grade)?);
        self.types[ty.0].bp_page = Some(data.clone());
        Ok(data)
    }

    pub(crate) fn odd_page(&mut self, ty: TypeId) -> Result<Rc<OddPage>> {
        if let Some(data) = &self.types[ty.0].odd_page {
            return Ok(data.clone());
        }
        let (prime, filter, r, s, t, f) =
            if let Some((prime, n, r, s, t, f)) = self.kind(ty).pce_grade() {
                (prime, Filtration::Ce(n), r, s, t, f)
            } else {
                let (prime, profile, r, s, t, f) = self.kind(ty).pmay_grade().unwrap();
                (prime, Filtration::May(profile), r, s, t, f)
            };
        let data = Rc::new(OddPage::build(self, prime, filter, r, s, t, f)?);
        self.types[ty.0].odd_page = Some(data.clone());
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

    pub(crate) fn rv_hom(&mut self, ty: TypeId) -> Result<Rc<Hom>> {
        if let Some(group) = &self.types[ty.0].rv_hom {
            return Ok(group.clone());
        }
        let Kind::RVHom(grade) = self.kind(ty) else {
            unreachable!()
        };
        let grade = grade.clone();
        let group = Rc::new(Hom::build(self, &grade)?);
        self.types[ty.0].rv_hom = Some(group.clone());
        Ok(group)
    }

    pub(crate) fn q_ext(&mut self, ty: TypeId) -> Result<Rc<QExt>> {
        if let Some(group) = &self.types[ty.0].q_ext {
            return Ok(group.clone());
        }
        let group = match self.kind(ty).clone() {
            Kind::QExt(grade) => QExt::build(
                self,
                grade.source,
                grade.target,
                &grade.selected,
                grade.s,
                grade.t,
            )?,
            Kind::VExt(grade) => QExt::build_localized(
                self,
                grade.source,
                grade.target,
                &grade.selected,
                grade.s,
                grade.t,
            )?,
            _ => unreachable!(),
        };
        let group = Rc::new(group);
        self.types[ty.0].q_ext = Some(group.clone());
        Ok(group)
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

    pub(crate) fn vbock_page(&mut self, ty: TypeId) -> Result<Rc<BocksteinPage>> {
        if let Some(page) = &self.types[ty.0].vbock {
            return Ok(page.clone());
        }
        let Kind::VBock(grade) = self.kind(ty) else {
            unreachable!()
        };
        let grade = grade.clone();
        let page = Rc::new(BocksteinPage::build(self, &grade)?);
        self.types[ty.0].vbock = Some(page.clone());
        Ok(page)
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

    pub(crate) fn ce_page(&mut self, ty: TypeId) -> Result<Rc<CePage>> {
        if let Some(page) = &self.types[ty.0].ce_page {
            return Ok(page.clone());
        }
        let Kind::CeAn(n, r, s, t, f) = *self.kind(ty) else {
            unreachable!()
        };
        let page = Rc::new(CePage::build(self, n, r, s, t, f)?);
        self.types[ty.0].ce_page = Some(page.clone());
        Ok(page)
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
            Kind::QModule(_) => {
                return Err(Error::plain(
                    "coefficient modules are allowed only as Koszul endpoints",
                ))
            }
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
            | Kind::PCobarAn(..)
            | Kind::PExt(..)
            | Kind::PExtAn(..)
            | Kind::PMay(..)
            | Kind::PMayAn(..)
            | Kind::PCeAn(..)
            | Kind::BPPage(..)
            | Kind::CobarAn(..)
            | Kind::Ext(..)
            | Kind::ExtAn(..)
            | Kind::XExt(..)
            | Kind::MapCobar(..)
            | Kind::MapExt(..)
            | Kind::QKoszul(..)
            | Kind::JetHom(..)
            | Kind::JetPage(..)
            | Kind::QExt(..)
            | Kind::RVHom(..)
            | Kind::ASDiff(..)
            | Kind::RVLaurent(..)
            | Kind::RVIwHecke(..)
            | Kind::VKoszul(..)
            | Kind::VExt(..)
            | Kind::VBock(..)
            | Kind::RVBarE(..)
            | Kind::RVBarCobar(..)
            | Kind::RVBarDual(..)
            | Kind::RVPathExt(..)
            | Kind::RVBar(..)
            | Kind::RVKoszul(..)
            | Kind::RVTor(..)
            | Kind::RVPoly(..)
            | Kind::BP(..)
            | Kind::BPRel(..)
            | Kind::BPExt(..)
            | Kind::BPTriple(..)
            | Kind::BPRelExt(..)
            | Kind::BPHeight(..)
            | Kind::MapMay(..)
            | Kind::May(..)
            | Kind::MayAn(..)
            | Kind::CeAn(..)
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
            Kind::Space(name) | Kind::QModule(name) => out.push(name),
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
            &Kind::PCobarAn(p, n, s, t) => {
                out.push(&format!("PCobarAn({}, {n}, {s}, {t})", p.get()))
            }
            &Kind::PExt(p, s, t) => out.push(&format!("PExt({}, {s}, {t})", p.get())),
            &Kind::PExtAn(p, n, s, t) => out.push(&format!("PExtAn({}, {n}, {s}, {t})", p.get())),
            &Kind::PMay(p, r, s, t, f) => {
                out.push(&format!("PMay({}, {r}, {s}, {t}, {f})", p.get()))
            }
            &Kind::PMayAn(p, n, r, s, t, f) => {
                out.push(&format!("PMayAn({}, {n}, {r}, {s}, {t}, {f})", p.get()))
            }
            &Kind::PCeAn(p, n, r, s, t, f) => {
                out.push(&format!("PCEAn({}, {n}, {r}, {s}, {t}, {f})", p.get()))
            }
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
            Kind::RVBarE(grade) => out.push(&format!(
                "RVBarE({}, {}, {}, {}, {}, {})",
                self.rv_bar_path_name(grade.chain.path),
                grade.chain.start,
                grade.chain.end,
                grade.r,
                grade.chain.b,
                grade.length,
            )),
            Kind::RVBar(grade)
            | Kind::RVBarDual(grade)
            | Kind::RVPathExt(grade)
            | Kind::RVBarCobar(grade) => {
                let family = match self.kind(ty) {
                    Kind::RVBar(_) => "RVBar",
                    Kind::RVBarCobar(_) => "RVBarCobar",
                    Kind::RVBarDual(_) => "RVBarDual",
                    _ => "RVPathExt",
                };
                out.push(&format!(
                    "{family}({}, {}, {}, {})",
                    self.rv_bar_path_name(grade.path),
                    grade.start,
                    grade.end,
                    grade.b
                ))
            }
            Kind::RVKoszul(grade) | Kind::RVTor(grade) => {
                let family = if matches!(self.kind(ty), Kind::RVKoszul(_)) {
                    "RVKoszul"
                } else {
                    "RVTor"
                };
                if let Some(multidegree) = &grade.multidegree {
                    let family = if matches!(self.kind(ty), Kind::RVKoszul(_)) {
                        "RVMKoszul"
                    } else {
                        "RVMTor"
                    };
                    out.push(&format!(
                        "{family}({}, {}, {multidegree:?})",
                        grade.ideal, grade.s
                    ))
                } else {
                    out.push(&format!(
                        "{family}({}, {}, {})",
                        grade.ideal, grade.s, grade.t
                    ))
                }
            }
            Kind::ASDiff(name, field) => out.push(&format!("ASDiff({name}, {})", field.modulus())),
            Kind::RVPoly(selected, degree) => out.push(&format!("RVPoly({selected:?}, {degree})")),
            Kind::BP(shape) => out.push(&shape.to_string()),
            Kind::BPHeight(grade) => out.push(&format!(
                "BPHeight({}, {}, {}, {}, {}, {})",
                grade.top.prime, grade.top.height, grade.r, grade.h, grade.top.arrows, grade.t,
            )),
            Kind::BPRelExt(grade) => out.push(&format!(
                "BPRelExt({}, {}, {}, {}, {})",
                grade.shape.upper.prime,
                grade.shape.lower,
                grade.shape.upper.height,
                grade.shape.upper.arrows,
                grade.t,
            )),
            Kind::BPRel(shape) => out.push(&format!(
                "BPRel({}, {}, {}, {})",
                shape.upper.prime, shape.lower, shape.upper.height, shape.upper.arrows,
            )),
            Kind::BPPage(grade) => out.push(&grade.to_string()),
            Kind::BPTriple(grade) => out.push(&format!(
                "BPTriple({}, {}, {}, {})",
                grade.shape.prime, grade.shape.height, grade.shape.arrows, grade.t
            )),
            Kind::BPExt(grade) => out.push(&format!(
                "BPExt({}, {}, {}, {})",
                grade.shape.prime, grade.shape.height, grade.shape.arrows, grade.t,
            )),
            Kind::JetPage(grade) => {
                out.push("JetPage(")?;
                self.write_jet_diagram(out, grade.chain.diagram)?;
                out.push(&format!(
                    ", {}, {}, {})",
                    grade.r, grade.chain.parity, grade.f
                ))
            }
            Kind::JetHom(grade) => {
                out.push("JetHom(")?;
                self.write_jet_diagram(out, grade.diagram)?;
                out.push(&format!(", {})", grade.parity))
            }
            Kind::QKoszul(grade) => self.write_koszul_type(out, "QKoszul", grade),
            Kind::QExt(grade) => self.write_koszul_type(out, "QExt", grade),
            Kind::RVIwHecke(rank) => out.push(&format!("RVIwHecke({rank})")),
            Kind::RVHom(grade) | Kind::RVLaurent(grade) => {
                let name = if matches!(self.kind(ty), Kind::RVHom(..)) {
                    "RVHom("
                } else {
                    "RVLaurent("
                };
                out.push(name)?;
                self.write_type(out, grade.source, 0)?;
                out.push(", ")?;
                self.write_type(out, grade.target, 0)?;
                out.push(&format!(
                    ", {:?}, {}, {:?})",
                    grade.selected,
                    grade.field.modulus(),
                    grade.point
                ))
            }
            Kind::VKoszul(grade) => self.write_koszul_type(out, "VKoszul", grade),
            Kind::VExt(grade) => self.write_koszul_type(out, "VExt", grade),
            Kind::VBock(grade) => {
                let localized = &grade.localized;
                out.push("VBock(")?;
                self.write_type(out, localized.source, 0)?;
                out.push(", ")?;
                self.write_type(out, localized.target, 0)?;
                out.push(&format!(
                    ", {:?}, {}, {}, {}, {})",
                    localized.selected, grade.r, localized.s, localized.t, grade.f,
                ))
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
            &Kind::CeAn(n, r, s, t, f) => out.push(&format!("CEAn({n}, {r}, {s}, {t}, {f})")),
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
            (Value::BPTriple(left), Value::BPTriple(right)) => {
                Ok(Rc::ptr_eq(left, right) || left.compare(self, right)?.is_eq())
            }
            (Value::BPExt(left), Value::BPExt(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(true);
                }
                self.spend(left.len().min(right.len()))?;
                Ok(left == right)
            }
            (Value::BP(left), Value::BP(right)) | (Value::BPRel(left), Value::BPRel(right)) => {
                Ok(Rc::ptr_eq(left, right) || left.compare(self, right)? == Ordering::Equal)
            }
            (Value::RVBarCobar(left), Value::RVBarCobar(right)) => Ok(Rc::ptr_eq(left, right)
                || rv_barcobar::compare(self, left, right)? == Ordering::Equal),
            (Value::RVBar(left), Value::RVBar(right)) => {
                Ok(Rc::ptr_eq(left, right)
                    || rv_bar::compare(self, left, right)? == Ordering::Equal)
            }
            (Value::RVKoszul(left), Value::RVKoszul(right)) => Ok(Rc::ptr_eq(left, right)
                || rv_koszul::compare(self, left, right)? == Ordering::Equal),
            (Value::RVPoly(left), Value::RVPoly(right)) => {
                Ok(Rc::ptr_eq(left, right)
                    || rv_poly::compare(self, left, right)? == Ordering::Equal)
            }
            (Value::RVIwHecke(left), Value::RVIwHecke(right)) => Ok(Rc::ptr_eq(left, right)
                || rv_iwhecke::compare(self, left, right)? == Ordering::Equal),
            (Value::RVLaurent(left), Value::RVLaurent(right)) => Ok(Rc::ptr_eq(left, right)
                || rv_laurent::compare(self, left, right)? == Ordering::Equal),
            (Value::RVHom(left), Value::RVHom(right))
            | (Value::ASDiff(left), Value::ASDiff(right))
            | (Value::JetClass(left), Value::JetClass(right)) => {
                Ok(Rc::ptr_eq(left, right)
                    || rv_hom::compare(self, left, right)? == Ordering::Equal)
            }
            (Value::JetHom(left), Value::JetHom(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(true);
                }
                self.spend(6 * left.len().min(right.len()))?;
                Ok(left == right)
            }
            (Value::QKoszul(left), Value::QKoszul(right)) => {
                Ok(Rc::ptr_eq(left, right)
                    || q_koszul::compare(self, left, right)? == Ordering::Equal)
            }
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
            (Value::BPTriple(left), Value::BPTriple(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    left.compare(self, right)
                }
            }
            (Value::BPExt(left), Value::BPExt(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(Ordering::Equal);
                }
                self.spend(left.len().min(right.len()))?;
                Ok(left.iter().rev().cmp(right.iter().rev()))
            }
            (Value::BP(left), Value::BP(right)) | (Value::BPRel(left), Value::BPRel(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    left.compare(self, right)
                }
            }
            (Value::RVBarCobar(left), Value::RVBarCobar(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    rv_barcobar::compare(self, left, right)
                }
            }
            (Value::RVBar(left), Value::RVBar(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    rv_bar::compare(self, left, right)
                }
            }
            (Value::RVKoszul(left), Value::RVKoszul(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    rv_koszul::compare(self, left, right)
                }
            }
            (Value::RVPoly(left), Value::RVPoly(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    rv_poly::compare(self, left, right)
                }
            }
            (Value::RVIwHecke(left), Value::RVIwHecke(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(Ordering::Equal);
                }
                rv_iwhecke::compare(self, left, right)
            }
            (Value::RVLaurent(left), Value::RVLaurent(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(Ordering::Equal);
                }
                rv_laurent::compare(self, left, right)
            }
            (Value::RVHom(left), Value::RVHom(right))
            | (Value::ASDiff(left), Value::ASDiff(right))
            | (Value::JetClass(left), Value::JetClass(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(Ordering::Equal);
                }
                rv_hom::compare(self, left, right)
            }
            (Value::JetHom(left), Value::JetHom(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(Ordering::Equal);
                }
                self.spend(6 * left.len().min(right.len()))?;
                Ok(left.iter().rev().cmp(right.iter().rev()))
            }
            (Value::QKoszul(left), Value::QKoszul(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(Ordering::Equal);
                }
                q_koszul::compare(self, left, right)
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
                Kind::QModule(_) => {
                    return Err(Error::plain(
                        "coefficient modules are allowed only as Koszul endpoints",
                    ))
                }
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
                Kind::RVIwHecke(..) => {
                    return Err(Error::plain("cannot enumerate an affine Hecke algebra"))
                }
                Kind::RVLaurent(..) => {
                    let (source, target) = rv_laurent::dimensions(self, ty)?;
                    if source != 0 && target != 0 {
                        return Err(Error::plain(
                            "cannot enumerate a positive-dimensional Laurent Margolis space",
                        ));
                    }
                    vec![rv_laurent::zero()]
                }
                Kind::ASDiff(ref name, field) => {
                    let genus = self.rv_as_covers[name.as_ref()].genus();
                    binary_linear::elements(self, field, genus, MAX_ELEMENTS, Value::ASDiff)?
                }
                Kind::RVHom(ref grade) => {
                    let field = grade.field;
                    rv_hom::elements(self, ty, field, MAX_ELEMENTS)?
                }
                Kind::QExt(..) | Kind::VExt(..) => {
                    let rank = self.q_ext(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::VBock(..) => {
                    let rank = self.vbock_page(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::BP(..) => return Err(Error::plain("cannot enumerate a BP polynomial ring")),
                Kind::BPRel(..) => {
                    return Err(Error::plain(
                        "cannot enumerate an integral BP relative cochain space",
                    ))
                }
                Kind::BPTriple(..) => {
                    return Err(Error::plain(
                        "BPTriple types do not enumerate ambient cosets",
                    ))
                }
                Kind::BPExt(..) | Kind::BPRelExt(..) | Kind::BPHeight(..) => {
                    self.bp_ext(ty)?.elements(self, MAX_ELEMENTS)?
                }
                Kind::RVBarCobar(grade) => {
                    let basis = rv_barcobar::basis(self, grade, MAX_ELEMENTS.ilog2() as usize)?;
                    self.allocate(basis.iter().map(|word| rv_barcobar::word_cells(word)).sum())?;
                    self.enumerate_supports(0, basis.len(), |class| {
                        Value::RVBarCobar(class.iter().map(|&item| basis[item].clone()).collect())
                    })?
                }
                Kind::RVPathExt(..) => {
                    let rank = self.rv_path_ext(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::RVBarE(..) => {
                    let rank = self.rv_bar_page(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::RVBar(grade) | Kind::RVBarDual(grade) => {
                    let basis = rv_bar::basis(
                        self,
                        grade,
                        grade.max_length(),
                        MAX_ELEMENTS.ilog2() as usize,
                    )?;
                    self.allocate(basis.iter().map(|word| rv_bar::word_cells(word)).sum())?;
                    self.enumerate_supports(0, basis.len(), |class| {
                        Value::RVBar(class.iter().map(|&item| basis[item].clone()).collect())
                    })?
                }
                Kind::RVKoszul(ref grade) => {
                    let grade = grade.clone();
                    let ideal = self.rv_ideals[&grade.ideal].clone();
                    let basis = rv_koszul::Slice::build_grade(
                        self,
                        &ideal,
                        &grade,
                        MAX_ELEMENTS.ilog2() as usize,
                    )?;
                    let rank = basis.len();
                    let terms: Vec<_> = (0..rank).map(|index| basis.term(index)).collect();
                    let retained = if grade.multidegree.is_some() {
                        rv_koszul::cells(&terms)
                    } else {
                        basis.cells()
                    };
                    self.allocate(retained + 2 * (support_space_size(rank) / 2) * rank)?;
                    self.enumerate_supports(0, rank, |support| {
                        Value::RVKoszul(support.iter().map(|&index| terms[index].clone()).collect())
                    })?
                }
                Kind::RVTor(..) => {
                    let rank = self.rv_tor(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::RVPoly(ref selected, degree) => {
                    let basis = rv_poly::basis(
                        self,
                        selected.len(),
                        degree,
                        MAX_ELEMENTS.ilog2() as usize,
                    )?;
                    self.enumerate_supports(0, basis.len(), |support| {
                        Value::RVPoly(support.iter().map(|&index| basis[index].clone()).collect())
                    })?
                }
                Kind::JetPage(grade) => {
                    let rank = self.jet_page(ty)?.rank();
                    let field = self.jet_path(grade.chain.diagram.path).field;
                    jet_page::elements(self, field, rank, MAX_ELEMENTS)?
                }
                Kind::JetHom(grade) => jet_hom::elements(self, grade, MAX_ELEMENTS)?,
                Kind::QKoszul(ref grade) => {
                    let grade = grade.clone();
                    let basis = q_koszul::basis(
                        self,
                        grade.source,
                        grade.target,
                        &grade.selected,
                        grade.s,
                        grade.t,
                        MAX_ELEMENTS.ilog2() as usize,
                    )?;
                    self.koszul_elements(&basis)?
                }
                Kind::VKoszul(ref grade) => {
                    let grade = grade.clone();
                    let basis = q_koszul::basis_localized(
                        self,
                        grade.source,
                        grade.target,
                        &grade.selected,
                        i128::from(grade.s),
                        grade.t,
                        MAX_ELEMENTS.ilog2() as usize,
                    )?;
                    self.koszul_elements(&basis)?
                }
                Kind::MapExt(..) => {
                    let rank = self.map_ext(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::May(..) | Kind::MayAn(..) => {
                    let rank = self.may_page(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::CeAn(..) => {
                    let rank = self.ce_page(ty)?.rank();
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
                Kind::PCobar(..) | Kind::PCobarAn(..) => {
                    let (p, profile, s, t) = self.kind(ty).pcobar_grade().unwrap();
                    odd_cobar::elements(self, p, profile, s, t, MAX_ELEMENTS)?
                }
                Kind::PExt(p, ..) | Kind::PExtAn(p, ..) => {
                    let rank = self.odd_ext_rank(ty)?;
                    odd_ext::elements(self, p, rank, MAX_ELEMENTS)?
                }
                Kind::BPPage(grade) => {
                    let rank = self.bp_page(ty)?.rank();
                    odd_ext::elements(
                        self,
                        Prime::validated(grade.shape.prime),
                        rank,
                        MAX_ELEMENTS,
                    )?
                }
                Kind::PMay(p, ..) | Kind::PMayAn(p, ..) | Kind::PCeAn(p, ..) => {
                    let rank = self.odd_page(ty)?.rank();
                    odd_ext::elements(self, p, rank, MAX_ELEMENTS)?
                }
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

    pub(crate) fn show_integer_coefficients(
        &self,
        name: &str,
        coefficients: &[i128],
    ) -> Result<String> {
        let mut out = Text::new(MAX_SHOWN_BYTES);
        out.push(name)?;
        out.push(" = [")?;
        for (index, coefficient) in coefficients.iter().enumerate() {
            if index != 0 {
                out.push(", ")?;
            }
            out.push(&coefficient.to_string())?;
        }
        out.push("]")?;
        Ok(out.string)
    }

    pub(crate) fn show(&mut self, ty: TypeId, value: &Value) -> Result<String> {
        let mut out = Text::new(MAX_SHOWN_BYTES);
        self.write_value(&mut out, ty, value)?;
        Ok(out.string)
    }

    fn koszul_elements(&mut self, basis: &[q_koszul::Term]) -> Result<Vec<Value>> {
        let rank = basis.len();
        self.allocate(
            basis.iter().map(q_koszul::Term::cells).sum::<usize>()
                + 4 * (support_space_size(rank) / 2) * rank,
        )?;
        self.enumerate_supports(0, rank, |support| {
            Value::QKoszul(support.iter().map(|&index| basis[index].clone()).collect())
        })
    }

    fn write_jet_diagram(&self, out: &mut Text, diagram: jet_hom::Diagram) -> Result<()> {
        self.write_type(out, diagram.source, 0)?;
        out.push(", ")?;
        self.write_type(out, diagram.target, 0)?;
        out.push(&format!(", {}", self.jet_paths[diagram.path.0].0))
    }

    fn write_koszul_type<S: std::fmt::Display>(
        &self,
        out: &mut Text,
        name: &str,
        grade: &Grade<S>,
    ) -> Result<()> {
        out.push(&format!("{name}("))?;
        self.write_type(out, grade.source, 0)?;
        out.push(", ")?;
        self.write_type(out, grade.target, 0)?;
        out.push(&format!(
            ", {:?}, {}, {})",
            grade.selected, grade.s, grade.t
        ))
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
                let (p, _, _, _) = self.kind(ty).pcobar_grade().unwrap();
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
            Value::BPTriple(coset) => {
                let Kind::BPTriple(grade) = *self.kind(ty) else {
                    unreachable!()
                };
                self.write_type(out, ty, 0)?;
                out.push("{base = ")?;
                let ext_type = self.bp_ext_type(grade)?;
                self.write_value(out, ext_type, &coset.base())?;
                out.push(&format!(", size = {}}}", coset.size()))
            }
            Value::BPExt(support) => {
                self.spend(support.len())?;
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
            Value::BP(polynomial) | Value::BPRel(polynomial) => {
                self.spend(polynomial.cells())?;
                let shape = match *self.kind(ty) {
                    Kind::BP(shape) => shape,
                    Kind::BPRel(shape) => shape.upper,
                    _ => unreachable!(),
                };
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                polynomial.write(|piece| out.push(piece), shape.height)?;
                out.push("}")
            }
            Value::RVBarCobar(words) => {
                let Kind::RVBarCobar(grade) = self.kind(ty) else {
                    unreachable!()
                };
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if words.is_empty() {
                    out.push("0")?;
                }
                for (index, word) in words.iter().enumerate() {
                    if index != 0 {
                        out.push(" + ")?;
                    }
                    let mut start = grade.start;
                    for (outer, block) in word.iter().enumerate() {
                        if outer != 0 {
                            out.push(" * ")?;
                        }
                        start = write_bar_word(out, block, start)?;
                    }
                }
                out.push("}")
            }
            Value::RVBar(words) => {
                let (Kind::RVBar(grade) | Kind::RVBarDual(grade)) = self.kind(ty) else {
                    unreachable!();
                };
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if words.is_empty() {
                    out.push("0")?;
                }
                for (index, word) in words.iter().enumerate() {
                    if index != 0 {
                        out.push(" + ")?;
                    }
                    write_bar_word(out, word, grade.start)?;
                }
                out.push("}")
            }
            Value::RVKoszul(terms) => {
                let Kind::RVKoszul(grade) = self.kind(ty) else {
                    unreachable!();
                };
                let selected = &self.rv_ideals[&grade.ideal].selected;
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                for (index, term) in terms.iter().rev().enumerate() {
                    if index != 0 {
                        out.push(" + ")?;
                    }
                    write_rv_polynomial(out, selected, std::slice::from_ref(&term.powers))?;
                    for (bit, primitive) in selected.iter().enumerate() {
                        if term.wedge & (1 << bit) != 0 {
                            out.push(&format!(" e{primitive}"))?;
                        }
                    }
                }
                out.push("}")
            }
            Value::RVPoly(polynomial) => {
                self.spend(rv_poly::cells(polynomial))?;
                let Kind::RVPoly(selected, degree) = self.kind(ty) else {
                    unreachable!()
                };
                out.push(&format!("RVPoly({selected:?}, {degree}){{"))?;
                write_rv_polynomial(out, selected, polynomial)?;
                out.push("}")
            }
            Value::RVIwHecke(terms) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                let mut first = true;
                for term in terms.iter() {
                    for &(degree, coefficient) in term.coefficient.iter() {
                        self.spend(1 + term.window.len())?;
                        if !first {
                            out.push(if coefficient < 0 { " - " } else { " + " })?;
                        } else if coefficient < 0 {
                            out.push("-")?;
                        }
                        first = false;
                        let magnitude = coefficient.unsigned_abs();
                        if magnitude != 1 {
                            out.push(&format!("{magnitude}*"))?;
                        }
                        match degree {
                            0 => (),
                            1 => out.push("q*")?,
                            degree => out.push(&format!("q^{degree}*"))?,
                        }
                        out.push(&format!("T{:?}", term.window))?;
                    }
                }
                out.push("}")
            }
            Value::RVLaurent(terms) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                self.spend(terms.len())?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                for (i, term) in terms.iter().enumerate() {
                    if i != 0 {
                        out.push(" + ")?;
                    }
                    if term.coefficient != 1 {
                        out.push(&format!("{}*", term.coefficient))?;
                    }
                    match term.power {
                        0 => (),
                        1 => out.push("u*")?,
                        power => out.push(&format!("u^{power}*"))?,
                    }
                    out.push(&format!("g{}", term.coordinate))?;
                }
                out.push("}")
            }
            Value::RVHom(support) | Value::JetClass(support) | Value::ASDiff(support) => {
                let symbol = if matches!(value, Value::ASDiff(_)) {
                    "omega"
                } else {
                    "g"
                };
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if support.is_empty() {
                    out.push("0")?;
                }
                for (i, &(coordinate, coefficient)) in support.iter().enumerate() {
                    if i != 0 {
                        out.push(" + ")?;
                    }
                    if coefficient != 1 {
                        out.push(&format!("{coefficient}*"))?;
                    }
                    out.push(&format!("{symbol}{coordinate}"))?;
                }
                out.push("}")
            }
            Value::JetHom(terms) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                self.spend(6 * terms.len())?;
                for (index, (term, coefficient)) in terms.iter().enumerate() {
                    if index != 0 {
                        out.push(" + ")?;
                    }
                    if *coefficient != 1 {
                        out.push(&format!("{coefficient}*"))?;
                    }
                    out.push(&format!(
                        "u^{}*e({},{}, {},{})",
                        term.power,
                        term.source.degree,
                        term.source.index,
                        term.target.degree,
                        term.target.index,
                    ))?;
                }
                out.push("}")
            }
            Value::QKoszul(terms) => {
                let selected = match self.kind(ty) {
                    Kind::QKoszul(grade) => &grade.selected,
                    Kind::VKoszul(grade) => &grade.selected,
                    _ => unreachable!(),
                };
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                for (index, term) in terms.iter().enumerate() {
                    if index != 0 {
                        out.push(" + ")?;
                    }
                    out.push(&format!(
                        "[{}, {} -> {}, {}]",
                        term.source.degree,
                        term.source.index,
                        term.target.degree,
                        term.target.index
                    ))?;
                    for (&primitive, &power) in selected.iter().zip(term.powers.iter()) {
                        match power {
                            0 => {}
                            1 => out.push(&format!(" v{primitive}"))?,
                            _ => out.push(&format!(" v{primitive}^{power}"))?,
                        }
                    }
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

fn write_bar_word(out: &mut Text, word: &rv_bar::Word, mut start: usize) -> Result<usize> {
    out.push("[")?;
    for (position, factor) in word.iter().enumerate() {
        if position != 0 {
            out.push(" | ")?;
        }
        out.push(&format!(
            "({}, {}, {}):g{}",
            start, factor.end, factor.s, factor.generator
        ))?;
        start = factor.end;
    }
    out.push("]")?;
    Ok(start)
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

fn write_rv_polynomial(
    out: &mut Text,
    selected: &[usize],
    polynomial: &[rv_poly::Monomial],
) -> Result<()> {
    if polynomial.is_empty() {
        return out.push("0");
    }
    for (i, term) in polynomial.iter().enumerate() {
        if i != 0 {
            out.push(" + ")?;
        }
        let mut constant = true;
        for (&primitive, &power) in selected.iter().zip(term.iter()) {
            if power == 0 {
                continue;
            }
            if !constant {
                out.push("*")?;
            }
            out.push(&format!("x{primitive}"))?;
            if power != 1 {
                out.push(&format!("^{power}"))?;
            }
            constant = false;
        }
        if constant {
            out.push("1")?;
        }
    }
    Ok(())
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
