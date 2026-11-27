use std::collections::{HashMap, HashSet};
use std::rc::Rc;

use crate::adams;
use crate::binary_field::{BinaryField, Embedding};
use crate::binary_linear;
use crate::bp_height;
use crate::bp_higher;
use crate::bp_mod_ext;
use crate::bp_page;
use crate::bp_triple;
use crate::bp_vpage;
use crate::bracket;
use crate::ce;
use crate::cobar;
use crate::cobar_squares;
use crate::cohomology;
use crate::comodule::{self, Generator};
use crate::cyclic::{self, Layout};
use crate::dependent_sets::{self, Operation as SetOperation};
use crate::eo;
use crate::error::{Error, Result, Span};
use crate::jet_hom;
use crate::jet_page;
use crate::map_bracket;
use crate::map_cochain;
use crate::map_ext;
use crate::map_may;
use crate::map_system;
use crate::map_transfer;
use crate::may;
use crate::may_lift;
use crate::milnor;
use crate::mixed_cone;
use crate::mixed_hom;
use crate::mixed_homology;
use crate::mixed_minimal;
use crate::mixed_tensor;
use crate::model::{DiffFamily, Kind, Model, Poset, SuspensionPoint, TypeId, Value};
use crate::motivic;
use crate::nakamura;
use crate::observations::{betti, homotopy};
use crate::odd_cobar::{self, Prime};
use crate::odd_ext;
use crate::odd_massey;
use crate::odd_may;
use crate::odd_profile;
use crate::odd_restriction_transfer;
use crate::odd_transfer;
use crate::page_transition::{self, Direction as PageDirection, PageFamily};
use crate::profile::Profile;
use crate::q_koszul::{self, Grade};
use crate::q_module;
use crate::rank_jets;
use crate::rank_variety;
use crate::real_ext;
use crate::real_extensions::{self, Operation as RealExtOp};
use crate::real_massey;
use crate::real_motivic;
use crate::real_rho::{self, Operation as RealRhoOp};
use crate::real_transfer;
use crate::restriction_transfer;
use crate::rho_bock;
use crate::rho_localized;
use crate::rho_minimal::{self, Operation as RhoMinOp};
use crate::rv_arithmetic;
use crate::rv_as_cartier;
use crate::rv_as_cover::Cover as RvAsCover;
use crate::rv_as_derham;
use crate::rv_as_odd::Cover as RvAsOddCover;
use crate::rv_as_odd_cartier;
use crate::rv_as_odd_derham;
use crate::rv_as_odd_tate;
use crate::rv_as_two::Cover as RvAsTwoCover;
use crate::rv_as_two_cartier;
use crate::rv_bar;
use crate::rv_bar_dual;
use crate::rv_bar_page;
use crate::rv_barcobar;
use crate::rv_conformal;
use crate::rv_curve::{self, Curve as RvCurve};
use crate::rv_cyclic_homology;
use crate::rv_cyclotomic;
use crate::rv_descent;
use crate::rv_field;
use crate::rv_fusion;
use crate::rv_gauss;
use crate::rv_hh;
use crate::rv_hh_homology;
use crate::rv_hh_mixed;
use crate::rv_hive;
use crate::rv_hochster;
use crate::rv_hodge;
use crate::rv_hom;
use crate::rv_ideal::Ideal as RvIdeal;
use crate::rv_iwahori;
use crate::rv_iwhecke;
use crate::rv_iwkl;
use crate::rv_koszul;
use crate::rv_laurent;
use crate::rv_lens;
use crate::rv_negative;
use crate::rv_newton;
use crate::rv_ore;
use crate::rv_pcurv;
use crate::rv_period::{self, Grade as PeriodGrade, Level as PeriodLevel};
use crate::rv_poly;
use crate::rv_qconnection;
use crate::rv_qgr;
use crate::rv_qp_quad::{self, Element as QuadElement, Field as QuadField, Unary as QuadUnary};
use crate::rv_qregular;
use crate::rv_rees;
use crate::rv_res_space::{self, Grade as ResGrade, Space as ResSpace};
use crate::rv_satake;
use crate::rv_scheme;
use crate::rv_shintani;
use crate::rv_sphhecke;
use crate::rv_t_bracket;
use crate::rv_t_system;
use crate::rv_t_transfer;
use crate::rv_tl;
use crate::rv_tree;
use crate::rv_verlinde;
use crate::rv_witt_cover::Cover as RvWittCover;
use crate::rv_witt_tower::Tower as RvWittTower;
use crate::rv_wzw;
use crate::set_family;
use crate::set_j;
use crate::set_logic::{self, Operation as SetLogic};
use crate::set_topology;
use crate::sheaf::{self, Cochains, Form as SheafForm, Operation as SheafOp, Sheaf};
use crate::sheaf_cellular::{Cellular, Resolution, Total as SheafTotal};
use crate::sheaf_deform;
use crate::sheaf_ext::{self, Ends as SheafEnds};
use crate::sheaf_higher::{self, Operation as SheafHigherOp};
use crate::sheaf_hyper::{self, NamedSheaf};
use crate::sheaf_image::{Operation as SheafImageOp, Site as SheafSite};
use crate::sheaf_kuranishi;
use crate::sheaf_leray::{Input as LerayInput, InputForm as SheafInputForm, Leray};
use crate::sheaf_page::{
    self, Grade as SheafPageGrade, Lift as SheafPageLift, Operation as SheafPageOp,
};
use crate::sheaf_support::{Grade as SupportGrade, Operation as SupportOp, Supports};
use crate::stable_marks::{Reduction as StableCurve, Tree as StableTree};
use crate::steenrod;
use crate::syntax::BarFamilySpec;
use crate::syntax::MixedSpec;
use crate::syntax::RvHHSpec;
use crate::syntax::RvNCSpec;
use crate::syntax::{
    BinaryOp, BpExtOp, BpHeightOp, BpHeightSpec, BpRelOp, ClassForm, CyclicForm, Entries, ExprForm,
    JetHomSpec, JetOp, JetPageSpec, JetPathArg, Label, MapMaySpec, ParsedExpr, ParsedType, QSpec,
    RealBarOp, RhoLocOp, RhoOp, RvIdealOp, RvSchemeOp, Statement, StatementForm, SuspensionForm,
    TypeForm, VBockSpec, Witness,
};
use crate::syntax::{OreOp, QuantumOp};
use crate::syntax::{RvBarPageSpec, RvBarSpec, RvDegree, RvGradeSpec};
use crate::syntax::{TlOp, WzwOp};
use crate::transfer;
use crate::tropical_matroid::{self, Grade as ValGrade};
use crate::v_bars;
use crate::v_koszul;
use crate::vbock;
use crate::{bp, bp_classes, bp_ext, bp_polynomial::Polynomial as BpPolynomial, bp_rel};
use crate::{mixed, mixed_page};
use crate::{sheaf_serre, sheaf_stalk, sheaf_truncate};

pub(crate) struct Engine {
    model: Model,
    bindings: HashMap<String, Binding>,
    values: Vec<Value>,
    locals: Vec<Value>,
}

#[derive(Clone, Copy)]
struct Binding {
    index: usize,
    ty: TypeId,
}

struct TypedExpr {
    span: Span,
    ty: TypeId,
    kind: TypedForm,
}

enum CobarOp {
    Differential,
    Class,
    Representative,
    Primitive,
    MayClass,
    MayRepresentative,
    MayDifferential,
    CeClass,
    CeRepresentative,
    CeDifferential,
    PageTransition(PageFamily, PageDirection),
    MapDifferential,
    MapClass,
    MapRepresentative,
    RvCyclicD,
    RvCyclicClass,
    RvCyclicRep,
    RvCyclicI,
    RvCyclicB,
    RvNCd,
    RvNCClass,
    RvNCRep,
    MixDown,
    MixUp,
    MixD,
    RvMixedToHH,
    MixHClass,
    MixHRep,
    MixHPrimitive,
    MixMinI,
    MixClass,
    MixRep,
    MixNext,
    MixPrev,
    RvHHClass,
    RvHHRep,
    RvHHd,
    RvConnesB,
    RvBCDesuspend,
    RvBCD,
    RvBCProject,
    RvBCHomotopy,
    RvBarDualize,
    RvBarDualD,
    RvPathClass,
    RvPathRep,
    RvBarClass,
    RvBarRep,
    RvBarNext,
    RvBarPrev,
    RvBarD,
    RvMForget,
    RvKoszulDifferential,
    RvKoszulPrimitive,
    RvTorClass,
    RvTorRepresentative,
    QDifferential,
    QClass,
    QRepresentative,
    VLocal,
    VBockClass,
    VBockDifferential,
    VBockRepresentative,
}

impl CobarOp {
    fn from_name(name: &str) -> Option<Self> {
        match name {
            "cobarD" => Some(Self::Differential),
            "extClass" => Some(Self::Class),
            "extRep" => Some(Self::Representative),
            "extPrimitive" => Some(Self::Primitive),
            "mayRep" => Some(Self::MayRepresentative),
            "mayD" => Some(Self::MayDifferential),
            "mayNext" => Some(Self::PageTransition(PageFamily::May, PageDirection::Next)),
            "mayPrev" => Some(Self::PageTransition(
                PageFamily::May,
                PageDirection::Previous,
            )),
            "ceNext" => Some(Self::PageTransition(PageFamily::Ce, PageDirection::Next)),
            "cePrev" => Some(Self::PageTransition(
                PageFamily::Ce,
                PageDirection::Previous,
            )),
            "ceRep" => Some(Self::CeRepresentative),
            "ceD" => Some(Self::CeDifferential),
            "eD" => Some(Self::MapDifferential),
            "eClass" => Some(Self::MapClass),
            "eRep" => Some(Self::MapRepresentative),
            "rvCyclicD" => Some(Self::RvCyclicD),
            "rvCyclicClass" => Some(Self::RvCyclicClass),
            "rvCyclicRep" => Some(Self::RvCyclicRep),
            "rvCyclicI" => Some(Self::RvCyclicI),
            "rvCyclicB" => Some(Self::RvCyclicB),
            "rvNCd" => Some(Self::RvNCd),
            "rvNCClass" => Some(Self::RvNCClass),
            "rvNCRep" => Some(Self::RvNCRep),
            "rvMixedToHH" => Some(Self::RvMixedToHH),
            "mixb" => Some(Self::MixDown),
            "mixB" => Some(Self::MixUp),
            "mixHClass" => Some(Self::MixHClass),
            "mixHRep" => Some(Self::MixHRep),
            "mixHPrimitive" => Some(Self::MixHPrimitive),
            "mixMinI" => Some(Self::MixMinI),
            "mixD" => Some(Self::MixD),
            "mixPageRep" => Some(Self::MixRep),
            "mixNext" => Some(Self::MixNext),
            "mixPrev" => Some(Self::MixPrev),
            "rvHHClass" => Some(Self::RvHHClass),
            "rvHHRep" => Some(Self::RvHHRep),
            "rvHHd" => Some(Self::RvHHd),
            "rvConnesB" => Some(Self::RvConnesB),
            "rvBCDesuspend" => Some(Self::RvBCDesuspend),
            "rvBCD" => Some(Self::RvBCD),
            "rvBCProject" => Some(Self::RvBCProject),
            "rvBCHomotopy" => Some(Self::RvBCHomotopy),
            "rvBarDualize" => Some(Self::RvBarDualize),
            "rvBarDualD" => Some(Self::RvBarDualD),
            "rvPathClass" => Some(Self::RvPathClass),
            "rvPathRep" => Some(Self::RvPathRep),
            "rvBarRep" => Some(Self::RvBarRep),
            "rvBarNext" => Some(Self::RvBarNext),
            "rvBarPrev" => Some(Self::RvBarPrev),
            "rvBarD" => Some(Self::RvBarD),
            "rvMForget" => Some(Self::RvMForget),
            "rvkD" => Some(Self::RvKoszulDifferential),
            "rvkPrimitive" => Some(Self::RvKoszulPrimitive),
            "rvTorClass" => Some(Self::RvTorClass),
            "rvTorRep" => Some(Self::RvTorRepresentative),
            "qD" => Some(Self::QDifferential),
            "qClass" => Some(Self::QClass),
            "qRep" => Some(Self::QRepresentative),
            "vLocal" => Some(Self::VLocal),
            "vbD" => Some(Self::VBockDifferential),
            "vbRep" => Some(Self::VBockRepresentative),
            "vbNext" => Some(Self::PageTransition(PageFamily::VBock, PageDirection::Next)),
            "vbPrev" => Some(Self::PageTransition(
                PageFamily::VBock,
                PageDirection::Previous,
            )),
            _ => None,
        }
    }
}

struct DefiningEntry {
    from: usize,
    to: usize,
    span: Span,
    expr: TypedExpr,
}

#[derive(Clone, Copy)]
enum ScalarField {
    ModTwo(Profile),
    Odd(Prime, Profile),
}

#[derive(Clone, Copy)]
enum ClassBracketFamily {
    Map,
    Koszul,
}

enum SheafComputation {
    Ordinary(Rc<Sheaf>),
    Total(SheafTotal),
    Hom(Rc<sheaf_ext::Hom>),
}

impl SheafComputation {
    fn data(&self) -> (&dyn Cochains, &sheaf::Cohomology) {
        match self {
            Self::Ordinary(data) => (data.as_ref(), &data.cohomology),
            Self::Total(data) => (data.filtered(), data.cohomology()),
            Self::Hom(data) => (data.as_ref(), &data.cohomology),
        }
    }
}

enum DependentAction {
    PullNatural(Rc<dependent_sets::Pullback>, Rc<set_family::Diagram>),
    Evaluation(Rc<dependent_sets::Product>),
    Transpose(
        Rc<dependent_sets::Product>,
        Rc<set_family::Diagram>,
        Rc<set_family::Diagram>,
    ),
    Untranspose(
        Rc<dependent_sets::Product>,
        Rc<set_family::Diagram>,
        Rc<set_family::Diagram>,
    ),
    SumLift(Rc<dependent_sets::Sum>, Rc<set_family::Diagram>),
    SumRestrict(
        Rc<dependent_sets::Sum>,
        Rc<set_family::Diagram>,
        Rc<set_family::Diagram>,
    ),
}

enum LogicAction {
    Constant(Rc<set_logic::Omega>, Rc<set_family::Diagram>, bool),
    Connect(Rc<set_logic::Omega>, Rc<set_family::Diagram>, SetLogic),
    Classify(
        Rc<set_logic::Omega>,
        Rc<set_family::Diagram>,
        Rc<set_family::Diagram>,
    ),
    Quantify(
        Rc<set_logic::Omega>,
        Rc<set_family::Diagram>,
        Rc<set_family::Diagram>,
        bool,
    ),
}

enum TypedForm {
    SetLogic(LogicAction, Vec<TypedExpr>),
    SetReflectExtend {
        reflection: Rc<set_topology::Reflection>,
        source: Rc<set_family::Diagram>,
        target: Rc<set_family::Diagram>,
        certificate: Rc<set_topology::Certificate>,
        map: Box<TypedExpr>,
    },
    SetJExtend {
        reflection: Rc<set_j::Sheafification>,
        source: Rc<set_family::Diagram>,
        target: Rc<set_family::Diagram>,
        certificate: Rc<set_j::Certificate>,
        map: Box<TypedExpr>,
    },
    SetDependent(DependentAction, Option<Box<TypedExpr>>),
    SetKan(crate::set_kan::Action, Option<Box<TypedExpr>>),
    SetIdentity(Rc<set_family::Diagram>),
    SetSection(Rc<set_family::Elements>, Rc<set_family::Diagram>, usize),
    SetCompose(
        Rc<set_family::Diagram>,
        Rc<set_family::Diagram>,
        Box<[TypedExpr; 2]>,
    ),
    SetTotal(
        Rc<set_family::Elements>,
        Rc<set_family::Diagram>,
        Box<TypedExpr>,
    ),
    Global(usize),
    Local(usize),
    Point(usize),
    Constant(Value),
    Binary(BinaryOp, Box<TypedExpr>, Box<TypedExpr>),
    Cobar(CobarOp, Box<TypedExpr>),
    RestrictAlgebra(Profile, Box<TypedExpr>),
    QRestrict(Box<TypedExpr>),
    RVScale(BinaryField, usize, Box<TypedExpr>),
    RVLaurentMonomial(i64, Box<TypedExpr>),
    HeckeNegate(Box<TypedExpr>),
    QpQuadUnary(QuadField, QuadUnary, Box<TypedExpr>),
    ResSpecializeQuad(ResGrade, QuadField, Vec<i64>, Vec<Vec<TypedExpr>>),
    ResOrthogonal(Box<TypedExpr>),
    ValMatQuad(ValGrade, QuadField, Vec<Vec<TypedExpr>>),
    ValDual(Box<TypedExpr>),
    StableLineQuad(QuadField, Vec<[TypedExpr; 2]>),
    StableProjection(Box<TypedExpr>),
    HeckeProduct(Box<[TypedExpr; 2]>),
    SphericalProduct(Box<[TypedExpr; 2]>),
    SphericalEmbedding(Box<TypedExpr>),
    SatakeTransform(Box<TypedExpr>),
    QuantumProduct(rv_qgr::Grade, Box<[TypedExpr; 2]>),
    PeriodInitial(PeriodGrade),
    PeriodQuery(PeriodGrade, u128),
    PeriodSection(PeriodGrade, u32, Box<TypedExpr>),
    PeriodRead(PeriodGrade, Box<TypedExpr>),
    QuantumNabla(rv_qgr::Grade, Box<TypedExpr>),
    QuantumPNabla(rv_pcurv::ModGrade, Box<TypedExpr>),
    QuantumPUnary(rv_pcurv::ModGrade, QuantumOp, Box<TypedExpr>),
    QuantumReduce(rv_pcurv::ModGrade, Box<TypedExpr>),
    QuantumCurvature(rv_pcurv::ModGrade, bool, Box<TypedExpr>),
    OreProduct(Box<[TypedExpr; 2]>),
    OreUnary(OreOp, Box<TypedExpr>),
    OreApply(rv_qgr::Grade, Box<[TypedExpr; 2]>),
    QuantumUnary(QuantumOp, Box<TypedExpr>),
    SatakeInverse(Box<TypedExpr>),
    FusionProduct(rv_fusion::Grade, Box<[TypedExpr; 2]>),
    FusionScale(rv_fusion::Grade, i128, Box<TypedExpr>),
    FusionDual(rv_fusion::Grade, Box<TypedExpr>),
    CycloProduct(usize, Box<[TypedExpr; 2]>),
    CycloScale(i128, Box<TypedExpr>),
    CycloConjugate(usize, Box<TypedExpr>),
    CycloEmbed(usize, usize, Box<TypedExpr>),
    CycloGalois(usize, i128, Box<TypedExpr>),
    WzwUnary(rv_wzw::Grade, WzwOp, Box<TypedExpr>),
    WzwScale(rv_wzw::Grade, Box<[TypedExpr; 2]>),
    WzwCoordinate(rv_wzw::Grade, Vec<i128>, Box<TypedExpr>),
    CycloRatio(u128, Box<TypedExpr>),
    CycloRatioNegate(Box<TypedExpr>),
    CycloRatioProduct(usize, Box<[TypedExpr; 2]>),
    CycloRatioGalois(usize, i128, Box<TypedExpr>),
    TlOperation(rv_tl::Grade, TlOp, Box<TypedExpr>),
    TlCoordinate(rv_tl::Grade, Vec<usize>, Box<TypedExpr>),
    TlScale(rv_tl::Grade, Box<[TypedExpr; 2]>),
    Verlinde(usize, Vec<i128>, Box<TypedExpr>),
    SatakeProduct(Box<[TypedExpr; 2]>),
    RVLaurentFrobenius(BinaryField, usize, Box<TypedExpr>),
    RVNewton(rv_newton::Operation, Embedding, Box<TypedExpr>),
    RVAdlv(rv_hodge::Boundary, usize, Vec<i64>, Box<[TypedExpr; 2]>),
    RVIwAdlv(rv_iwahori::Boundary, usize, Vec<i128>, Box<[TypedExpr; 2]>),
    RVSpecialize(Box<TypedExpr>),
    RVField(rv_field::Map, Box<TypedExpr>),
    RVInverse(Box<TypedExpr>),
    RVHilbert90(Embedding, Box<TypedExpr>),
    RVShintani(rv_shintani::Operation, Embedding, Box<TypedExpr>),
    RVTwisted(rv_descent::Twist, Embedding, Box<[TypedExpr; 3]>),
    VShift(i64, Box<TypedExpr>),
    VRatio(usize, Box<TypedExpr>),
    BpOperation(bp::Operation, bp::Shape, Box<TypedExpr>),
    BpRelOperation(BpRelOp, bp_rel::Shape, Box<TypedExpr>),
    BpHeightOperation(BpHeightOp, Box<TypedExpr>),
    BpExtOperation(BpExtOp, Box<TypedExpr>),
    BpMassey(bp::Massey, Vec<TypedExpr>),
    BpHigher(bp_higher::Layout, Vec<TypedExpr>, Option<Vec<TypedExpr>>),
    BpTriple(bp_triple::Layout, Box<[TypedExpr; 3]>),
    BpTripleBase(Box<TypedExpr>),
    BpRaise(Box<TypedExpr>),
    BpPageClass(Box<TypedExpr>),
    BpModExtClass(Box<TypedExpr>),
    BpModExtRep(Box<TypedExpr>),
    BpModExtOperation(bp::Operation, Box<TypedExpr>),
    BpPageOperation(bp_page::Operation, Box<TypedExpr>),
    BpVPageClass(Box<TypedExpr>),
    BpVPageOperation(bp_vpage::Operation, Box<TypedExpr>),
    RvNCTerm(i128, Box<TypedExpr>),
    RvNCPart(i128, Box<TypedExpr>),
    RvCyclicShift(rv_cyclic_homology::Shift, Box<TypedExpr>),
    RvNCTransfer(rv_negative::Transfer, Box<TypedExpr>),
    RvHHWord(Vec<TypedExpr>),
    RvHHLength(usize, Box<TypedExpr>),
    RvBarWord(Vec<usize>, Vec<usize>, Vec<TypedExpr>),
    RvBarLength(usize, Box<TypedExpr>),
    RvMComponent(Box<TypedExpr>),
    RvKoszulTerm(usize, Rc<RvIdeal>, Box<TypedExpr>),
    RvHasse(Vec<usize>, Box<TypedExpr>),
    RvVanish(BinaryField, Vec<usize>, Box<TypedExpr>),
    RvInMinors(TypeId, q_koszul::Indices, usize, Box<TypedExpr>),
    RvIdeal(RvIdealOp, Rc<RvIdeal>, Box<TypedExpr>),
    DiffScale(BinaryField, usize, Box<TypedExpr>),
    AsDeRhamOperation(
        rv_as_derham::Operation,
        Rc<RvAsCover>,
        BinaryField,
        Box<TypedExpr>,
    ),
    AsOddTateOperation(
        rv_as_odd_tate::Operation,
        Rc<RvAsOddCover>,
        i64,
        Box<TypedExpr>,
    ),
    AsOddDeRhamOperation(
        rv_as_odd_derham::Operation,
        Rc<RvAsOddCover>,
        Box<TypedExpr>,
    ),
    AsTwoDiffOperation(
        rv_as_two_cartier::Operation,
        Rc<RvAsTwoCover>,
        Box<TypedExpr>,
    ),
    AsOddDiffOperation(
        rv_as_odd_cartier::Operation,
        Rc<RvAsOddCover>,
        Box<TypedExpr>,
    ),
    AsCartier(Rc<RvAsCover>, BinaryField, Box<TypedExpr>),
    WittCartier(Rc<RvWittCover>, BinaryField, Box<TypedExpr>),
    CobarCup(Profile, usize, Box<TypedExpr>, Box<TypedExpr>),
    ExtSquare(usize, Box<TypedExpr>),
    Nakamura(usize, Box<TypedExpr>),
    MayLiftTo(usize, usize, Box<TypedExpr>),
    Coact(usize, usize, Box<TypedExpr>),
    MotivicCast(Box<TypedExpr>),
    RhoShift(usize, Box<TypedExpr>),
    RealComplex(Box<TypedExpr>),
    RhoNull(usize, Box<TypedExpr>),
    RealBar(RealBarOp, usize, i128, usize),
    RealExtension(RealExtOp, motivic::Grade, Box<TypedExpr>),
    RhoProductNull(usize, Box<[TypedExpr; 2]>),
    RealRhoOperation(RealRhoOp, motivic::Grade, Box<TypedExpr>),
    RealMassey(bool, real_massey::Layout, Box<[TypedExpr; 3]>),
    RealMasseyNull(motivic::Grade, Box<[TypedExpr; 2]>),
    RealTripleClass(Option<usize>, Box<TypedExpr>),
    RhoOperation(RhoOp, Box<TypedExpr>),
    RhoLocOperation(RhoLocOp, Box<TypedExpr>),
    RhoMinOperation(RhoMinOp, motivic::Grade, Box<TypedExpr>),
    MapExtInduced(Box<TypedExpr>),
    MapWord(Generator, Generator, Box<TypedExpr>),
    Equator,
    Suspend(Box<TypedExpr>),
    Shift(Box<TypedExpr>),
    SheafPageOperation(SheafPageOp, SheafTotal, SheafPageGrade, Box<TypedExpr>),
    SheafPageLift(SheafPageLift, SheafTotal, SheafPageGrade, Box<TypedExpr>),
    RvHHToMixed(Rc<rv_hh_mixed::Origin>, Box<TypedExpr>),
    MixedConeMap(Rc<mixed_cone::Cone>, mixed_cone::Operation, Box<TypedExpr>),
    MixedConeShear(
        Rc<mixed_cone::Cone>,
        Rc<mixed_cone::Cone>,
        mixed::Id,
        Box<TypedExpr>,
        Box<TypedExpr>,
    ),
    MixedRawDelta(Rc<mixed::Complex>, usize, Box<TypedExpr>),
    MixedTensorProduct(Rc<mixed_tensor::Tensor>, Box<TypedExpr>, Box<TypedExpr>),
    MixedTensorSwap(
        Rc<mixed_tensor::Tensor>,
        Rc<mixed_tensor::Tensor>,
        Box<TypedExpr>,
    ),
    MixedComposition([Rc<mixed_hom::Hom>; 3], Box<TypedExpr>, Box<TypedExpr>),
    MixedMinMap(
        Rc<mixed_minimal::Minimal>,
        mixed_minimal::Operation,
        bool,
        Box<TypedExpr>,
    ),
    MixedDelta(Rc<mixed_minimal::Minimal>, usize, Box<TypedExpr>),
    MixedHShift(Box<TypedExpr>),
    MixedEmbed(Box<TypedExpr>),
    SheafOperation(SheafOp, SheafComputation, usize, Box<TypedExpr>),
    SheafLerayIncl(Rc<Leray>, usize, Box<TypedExpr>),
    SheafLocalClass(Rc<Supports>, SupportGrade, Box<TypedExpr>),
    SheafLocalOperation(Rc<Supports>, SupportGrade, SupportOp, Box<TypedExpr>),
    SheafApply(Rc<sheaf_hyper::Map>, usize, Box<TypedExpr>),
    SheafDefOperation(Rc<sheaf_ext::Hom>, sheaf_deform::Operation, Box<TypedExpr>),
    SheafDefAppend(Rc<sheaf_ext::Hom>, Box<[TypedExpr; 2]>),
    SheafDefGauge(Rc<sheaf_ext::Hom>, usize, bool, Box<[TypedExpr; 2]>),
    SheafCompose(Box<[TypedExpr; 2]>),
    SheafHigher(
        SheafHigherOp,
        sheaf_higher::Signature,
        Vec<TypedExpr>,
        Vec<DefiningEntry>,
        Option<Box<TypedExpr>>,
    ),
    SheafEmbed(Rc<sheaf_hyper::Complex>, usize, Box<TypedExpr>),
    SheafPart(Rc<sheaf_hyper::Complex>, usize, Box<TypedExpr>),
    PScale(Prime, usize, Box<TypedExpr>),
    JetOperation(JetOp, Box<TypedExpr>),
    Bracket(Vec<TypedExpr>),
    Transfer(ScalarField, transfer::Operation, Vec<TypedExpr>),
    RealTransfer(transfer::Operation, Vec<TypedExpr>),
    HigherRestriction(ScalarField, Profile, Vec<TypedExpr>),
    MapTransfer(map_transfer::Operation, Vec<TypedExpr>),
    RvTTransfer(rv_t_transfer::Operation, Vec<TypedExpr>),
    RvTSystem(rv_t_system::Operation, usize, Vec<DefiningEntry>),
    MapSystem(map_system::Operation, usize, Vec<DefiningEntry>),
    OddSystem(odd_massey::Operation, usize, Vec<DefiningEntry>),
    ClassBracket {
        family: ClassBracketFamily,
        target: TypeId,
        inputs: Vec<TypedExpr>,
        query: Option<Box<TypedExpr>>,
    },
    DefiningSystem {
        arity: usize,
        layout: Layout,
        cochains: Vec<TypedExpr>,
    },
    Lambda(Box<TypedExpr>),
    Call(Box<TypedExpr>, Box<TypedExpr>),
    Pair(Box<TypedExpr>, Box<TypedExpr>),
    Fst(Box<TypedExpr>),
    Snd(Box<TypedExpr>),
    Table(Vec<TypedExpr>),
    Restrict(Box<TypedExpr>),
    Sew(Vec<(usize, TypedExpr)>),
}

impl Engine {
    pub(crate) fn new() -> Self {
        Self {
            model: Model::new(),
            bindings: HashMap::new(),
            values: Vec::new(),
            locals: Vec::new(),
        }
    }

    fn checker(&mut self) -> Checker<'_> {
        Checker {
            model: &mut self.model,
            bindings: &self.bindings,
            locals: Vec::new(),
        }
    }

    pub(crate) fn execute(&mut self, statement: Statement) -> Result<Option<String>> {
        self.model.reset_work();
        self.run_statement(statement.kind)
            .map_err(|error| error.at(statement.span))
    }

    fn set_topology_command(&mut self, command: set_topology::Command) -> Result<Option<String>> {
        use set_topology::{Command, Topology};
        match command {
            Command::Subset(name, omega, points) => {
                let classifier = self.model.set_omega(&omega.text)?;
                let topology =
                    Topology::subset(&mut self.model, omega.text.into(), classifier, &points)?;
                self.model.define_set_topology(&name, topology)?;
            }
            Command::FromMap(name, expression) => {
                let expression = self.checker().expr(*expression)?;
                let (source, target) = self.checker().set_nat_ends(expression.ty)?;
                if source != target {
                    return Err(Error::plain(
                        "settopologyMap requires an endomorphism of a set truth classifier",
                    ));
                }
                let omega = self.model.set_omega(&source)?;
                let diagram = self
                    .model
                    .set_family_named(&source)
                    .expect("registered classifier");
                let Value::SetNat(rows) = self.eval(&expression)? else {
                    unreachable!()
                };
                let topology = Topology::from_map(&mut self.model, source, omega, &diagram, rows)?;
                self.model.define_set_topology(&name, topology)?;
            }
            Command::Lattice(name, left, right, meet) => {
                let left = self.model.set_topology(&left.text)?;
                let right = self.model.set_topology(&right.text)?;
                let topology = left.combine(&mut self.model, &right, meet)?;
                self.model.define_set_topology(&name, topology)?;
            }
            Command::Support(name) => {
                let topology = self.model.set_topology(&name.text)?;
                let points = topology.members(&mut self.model)?.join(", ");
                return Ok(Some(format!(
                    "settopologysupport[{}] = {{{points}}}",
                    name.text
                )));
            }
            Command::JoinSteps(left, right) => {
                let a = self.model.set_topology(&left.text)?;
                let b = self.model.set_topology(&right.text)?;
                let omega = self
                    .model
                    .set_family_named(&a.omega_name)
                    .expect("registered classifier");
                let steps = a.join_steps(&mut self.model, &omega, &b)?;
                return Ok(Some(format!(
                    "settopologyjoinsteps[{}, {}] = {steps}",
                    left.text, right.text
                )));
            }
            Command::Reflect(name, topology, family) => {
                self.model.check_new_set_family(&name)?;
                let data = self.model.set_topology(&topology.text)?;
                let source = self.model.set_family(&family)?;
                let (diagram, reflection) = set_topology::Reflection::build(
                    &mut self.model,
                    topology.text.into(),
                    &data,
                    family.text.into(),
                    &source,
                )?;
                self.model
                    .define_set_reflection(&name, diagram, reflection)?;
            }
            Command::Check(topology, family) => {
                self.model.set_family(&family)?;
                let sheaf = self
                    .model
                    .set_reflection_certificate(
                        topology.text.as_str().into(),
                        family.text.as_str().into(),
                    )?
                    .is_some();
                let answer = if sheaf { "yes" } else { "no" };
                return Ok(Some(format!(
                    "setreflectcheck[{}, {}] = {answer}",
                    topology.text, family.text
                )));
            }
        }
        Ok(None)
    }

    fn define_sheaf_kuranishi(
        &mut self,
        name: Label,
        sheaf: Label,
        order: usize,
        directions: Option<Vec<ParsedExpr>>,
    ) -> Result<()> {
        self.model.check_new_sheaf_kuranishi(&name)?;
        let ty = self
            .model
            .sheaf_hom_type(SheafForm::Class, &sheaf, &sheaf, 1)?;
        let Kind::SheafHom(_, ends, _) = self.model.kind(ty) else {
            unreachable!()
        };
        let hom = self.model.sheaf_hom(ends);
        let directions = if let Some(directions) = directions {
            let mut classes = Vec::with_capacity(directions.len());
            for direction in directions {
                let direction = self.checker().expr(direction)?;
                if direction.ty != ty {
                    return Err(Error::new(direction.span, "Kuranishi directions must be degree-one self-Ext classes on the selected sheaf"));
                }
                let Value::OddClass(coordinates) = self.eval(&direction)? else {
                    unreachable!()
                };
                classes.push(coordinates);
            }
            Some(classes)
        } else {
            None
        };
        let chart = sheaf_kuranishi::Chart::build(
            &mut self.model,
            &hom,
            sheaf.text.into(),
            order,
            directions,
        )?;
        self.model.define_sheaf_kuranishi(&name, chart)
    }

    fn run_statement(&mut self, statement: StatementForm) -> Result<Option<String>> {
        let line = match statement {
            StatementForm::Legendre(query) => query.run(&mut self.model)?,
            StatementForm::QModule(name, selected, presentation) => {
                self.model.define_q_module(&name, selected, presentation)?;
                return Ok(None);
            }
            StatementForm::MixedHigher(name, prime, generators, rules) => {
                self.model.check_new_mixed(&name)?;
                let complex =
                    mixed::Complex::build_higher(&mut self.model, prime, generators, rules)?;
                self.model.define_mixed(name.text, complex);
                return Ok(None);
            }
            StatementForm::MixedCone(name, map) => {
                self.model.check_new_mixed(&name)?;
                let map = self.checker().expr(map)?;
                let Kind::MixedChain(
                    mixed::Grade {
                        complex: map_id,
                        b: 0,
                    },
                    true,
                ) = *self.model.kind(map.ty)
                else {
                    return Err(Error::plain(
                        "mixcone requires a degree-zero MixedU map in a mixhom complex",
                    ));
                };
                let Value::OddClass(terms) = self.eval(&map)? else {
                    unreachable!()
                };
                let (cone, complex) = mixed_cone::Cone::build(&mut self.model, map_id, terms)?;
                self.model.define_mixed_cone(name.text, complex, cone);
                return Ok(None);
            }
            StatementForm::MixedTensor(name, left, right) => {
                self.model.check_new_mixed(&name)?;
                let left = self.model.mixed_id(&left)?;
                let right = self.model.mixed_id(&right)?;
                let (tensor, complex) = mixed_tensor::Tensor::build(&mut self.model, left, right)?;
                self.model.define_mixed_tensor(name.text, complex, tensor);
                return Ok(None);
            }
            StatementForm::MixedHom(name, source, target) => {
                self.model.check_new_mixed(&name)?;
                let source = self.model.mixed_id(&source)?;
                let target = self.model.mixed_id(&target)?;
                let (hom, complex) = mixed_hom::Hom::build(&mut self.model, source, target)?;
                self.model.define_mixed_hom(name.text, complex, hom);
                return Ok(None);
            }
            StatementForm::MixedComplex(name, prime, generators, b, big_b) => {
                self.model.check_new_mixed(&name)?;
                let complex = mixed::Complex::build(&mut self.model, prime, generators, b, big_b)?;
                self.model.define_mixed(name.text, complex);
                return Ok(None);
            }
            StatementForm::RvHHMixed(name, ideal, weight) => {
                self.model.check_new_mixed(&name)?;
                let grade = self.checker().rv_hh_grade((ideal, 0, weight))?;
                let (origin, complex) = rv_hh_mixed::Origin::build(&mut self.model, grade)?;
                self.model.define_hh_mixed(name.text, complex, origin);
                return Ok(None);
            }
            StatementForm::MixedHDimension(name, b) => {
                let complex = self.model.mixed_id(&name)?;
                let rank = self
                    .model
                    .mixed_bars(complex)?
                    .rank(&mut self.model, b, b)?;
                format!("mixHdim[{}, {b}] = {rank}", name.text)
            }
            StatementForm::MixedBars(name) => {
                let id = self.model.mixed_id(&name)?;
                self.model.show_mixed_bars(id)?
            }
            StatementForm::MixedBarRank(name, source, target) => {
                let id = self.model.mixed_id(&name)?;
                let bars = self.model.mixed_bars(id)?;
                let rank = bars.rank(&mut self.model, source, target.unwrap_or(source))?;
                if let Some(target) = target {
                    format!("mixbarrank[{}, {source}, {target}] = {rank}", name.text)
                } else {
                    format!("mixbardim[{}, {source}] = {rank}", name.text)
                }
            }
            StatementForm::MixedPeriodicRank(name, degree) => {
                let id = self.model.mixed_id(&name)?;
                let rank = self
                    .model
                    .mixed_bars(id)?
                    .periodic_rank(&mut self.model, degree)?;
                format!("mixbarperiodicdim[{}, {degree}] = {rank}", name.text)
            }
            StatementForm::MixedMinimal(name, source) => {
                self.model.check_new_mixed_minimal(&name)?;
                let source = self.model.mixed_id(&source)?;
                let data = mixed_minimal::Minimal::build(&mut self.model, source)?;
                self.model.define_mixed_minimal(name.text, data);
                return Ok(None);
            }
            StatementForm::MixedMinDimension(name, b) => {
                let id = self.model.mixed_minimal_id(&name)?;
                let dimension = self
                    .model
                    .mixed_minimal(id)
                    .basis(&mut self.model, b, false)?
                    .len();
                format!("mixmindim[{}, {b}] = {dimension}", name.text)
            }
            StatementForm::MixedDimension(spec) => {
                let grade = self.checker().mixed_grade(spec)?;
                let ty = self.model.mixed_page_type(grade)?;
                let rank = self.model.mixed_page(ty)?.rank();
                format!(
                    "mixdim[{}, {}, {}, {}] = {rank}",
                    self.model.mixed_name(grade.chain.complex),
                    grade.r,
                    grade.chain.b,
                    grade.f
                )
            }
            StatementForm::SheafKuranishi(name, sheaf, order, directions) => {
                self.define_sheaf_kuranishi(name, sheaf, order, directions)?;
                return Ok(None);
            }
            StatementForm::SheafKuranishiEquation(name, index) => {
                self.model.show_kuranishi_equation(&name, index)?
            }
            StatementForm::SheafKuranishiDimensions(name) => {
                let chart = self.model.sheaf_kuranishi(&name)?;
                format!(
                    "shkuranishiDims[{}] = tangent={}, obstruction={}, order={}",
                    name.text, chart.directions, chart.obstructions, chart.order
                )
            }
            StatementForm::SheafExtDimension(source, target, degree) => {
                let ty = self
                    .model
                    .sheaf_hom_type(SheafForm::Class, &source, &target, degree)?;
                let Kind::SheafHom(_, ends, _) = self.model.kind(ty) else {
                    unreachable!()
                };
                let data = self.model.sheaf_hom(ends);
                let rank = data
                    .cohomology
                    .rank(&mut self.model, data.as_ref(), degree)?;
                format!(
                    "shExtDim[{}, {}, {degree}] = {rank}",
                    source.text, target.text
                )
            }
            StatementForm::Sheaf(name, base, prime, presentation) => {
                self.model.check_new_sheaf(&name)?;
                let ty = self.model.named(&base.text).ok_or_else(|| {
                    Error::new(
                        base.span,
                        format!("unknown named sheaf base '{}'", base.text),
                    )
                })?;
                let data = Sheaf::build(&mut self.model, ty, prime, presentation)?;
                self.model.define_sheaf(name.text, data);
                return Ok(None);
            }
            StatementForm::SheafSerre(form, name, source) => {
                self.model.check_new_sheaf_complex(&name)?;
                let (shift, components, complex) = match form {
                    SheafInputForm::Sheaf => {
                        let source = self.model.sheaf(&source)?;
                        sheaf_serre::build(&mut self.model, &name, &source)?
                    }
                    SheafInputForm::Complex => {
                        let source = self.model.sheaf_complex(&source)?;
                        sheaf_serre::build_complex(&mut self.model, &name, &source)?
                    }
                };
                for component in components {
                    self.model
                        .define_sheaf_shared(component.name.to_string(), component.sheaf);
                }
                self.model.define_sheaf_complex(name.text.clone(), complex);
                self.model.define_serre_shift(name.text, shift);
                return Ok(None);
            }
            StatementForm::SheafSerreShift(name) => {
                let shift = self.model.serre_shift(&name)?;
                format!("shserreshift[{}] = {shift}", name.text)
            }
            StatementForm::SheafStalkCohomology(name, complex, degree) => {
                self.model.check_new_sheaf(&name)?;
                let complex = self.model.sheaf_complex(&complex)?;
                let sheaf = sheaf_stalk::cohomology(&mut self.model, &complex, degree)?;
                self.model.define_sheaf(name.text, sheaf);
                return Ok(None);
            }
            StatementForm::SheafTruncate(name, original, degree) => {
                self.model.check_new_sheaf_complex(&name)?;
                let original = self.model.sheaf_complex(&original)?;
                let (result, terminal) =
                    sheaf_truncate::truncate(&mut self.model, &name, original, degree)?;
                if let Some(entry) = terminal {
                    self.model
                        .define_sheaf_shared(entry.name.to_string(), entry.sheaf);
                }
                self.model.define_sheaf_complex_shared(name.text, result);
                return Ok(None);
            }
            StatementForm::SheafDeligne(form, name, map, source, degree) => {
                self.model.check_new_sheaf_complex(&name)?;
                let site = self.sheaf_site(&map, form.deligne_name())?;
                let source = LerayInput::resolve(&self.model, &source, form)?;
                let result = sheaf_truncate::deligne(&mut self.model, &name, source, site, degree)?;
                for entry in result.components() {
                    let label = Label {
                        text: entry.name.to_string(),
                        span: name.span,
                    };
                    self.model.check_new_sheaf(&label)?;
                }
                for entry in result.components() {
                    self.model
                        .define_sheaf_shared(entry.name.to_string(), entry.sheaf.clone());
                }
                self.model.define_sheaf_complex_shared(name.text, result);
                return Ok(None);
            }
            StatementForm::SheafExtra(form, name, map, source) => {
                self.model.check_new_sheaf_complex(&name)?;
                let site = self.sheaf_site(&map, form.extra_name())?;
                let source = LerayInput::resolve(&self.model, &source, form)?;
                let (supports, components, complex) =
                    Supports::build(&mut self.model, &name, source, site)?;
                for entry in components {
                    self.model
                        .define_sheaf_shared(entry.name.to_string(), entry.sheaf);
                }
                self.model.define_sheaf_complex(name.text.clone(), complex);
                self.model.define_sheaf_supports(name.text, supports);
                return Ok(None);
            }
            StatementForm::SheafLocalDimension(name, point, mode, degree) => {
                let data = self.model.sheaf_supports(&name)?;
                let index = data.point(&self.model, &point)?;
                let grade = SupportGrade {
                    point: index,
                    mode,
                    degree,
                };
                let rank = data.rank(&mut self.model, grade)?;
                format!(
                    "shlocaldim[{}, {}, {}, {degree}] = {rank}",
                    name.text,
                    point.text,
                    mode.name()
                )
            }
            StatementForm::SheafLeray(form, name, map, source) => {
                self.model.check_new_sheaf_complex(&name)?;
                let site = self.sheaf_site(&map, form.leray_name())?;
                let source = LerayInput::resolve(&self.model, &source, form)?;
                let (leray, components, complex) =
                    Leray::build(&mut self.model, &name, source, site)?;
                for entry in components {
                    self.model
                        .define_sheaf_shared(entry.name.to_string(), entry.sheaf);
                }
                self.model.define_sheaf_complex(name.text.clone(), complex);
                self.model.define_leray_complex(name.text, leray);
                return Ok(None);
            }
            StatementForm::SheafImage(operation, name, function, sheaf) => {
                self.model.check_new_sheaf(&name)?;
                let source = self.model.sheaf(&sheaf)?;
                let site = self.sheaf_site(&function, operation.name())?;
                let image = match operation {
                    SheafImageOp::Pull => site.pull(&mut self.model, &source)?,
                    SheafImageOp::Push(degree) => site.push(&mut self.model, &source, degree)?,
                };
                self.model.define_sheaf(name.text, image);
                return Ok(None);
            }
            StatementForm::SheafMap(name, source, target, matrices) => {
                self.model.check_new_sheaf_map(&name)?;
                let source = NamedSheaf::resolve(&self.model, &source)?;
                let target = NamedSheaf::resolve(&self.model, &target)?;
                let map = sheaf_hyper::Map::build(&mut self.model, source, target, matrices)?;
                self.model.define_sheaf_map(name.text, map);
                return Ok(None);
            }
            StatementForm::SheafResolution(name, base, characteristic, presentation) => {
                self.model.check_new_sheaf_resolution(&name)?;
                let ty = self.model.named(&base.text).ok_or_else(|| {
                    Error::new(
                        base.span,
                        format!("unknown named resolution base '{}'", base.text),
                    )
                })?;
                let prime = sheaf::field(&mut self.model, characteristic)?;
                let resolution = Resolution::build(&mut self.model, ty, prime, presentation)?;
                self.model.define_sheaf_resolution(&name, resolution);
                return Ok(None);
            }
            StatementForm::SheafCellular(name, resolution, complex) => {
                self.model.check_new_sheaf_complex(&name)?;
                let resolution = self.model.sheaf_resolution(&resolution)?;
                let complex = self.model.sheaf_complex(&complex)?;
                let total = Cellular::new(&mut self.model, resolution, complex)?;
                self.model.define_sheaf_cellular(name.text, total);
                return Ok(None);
            }
            StatementForm::SheafComplex(name, entries) => {
                self.model.check_new_sheaf_complex(&name)?;
                if entries.len() % 2 == 0 {
                    return Err(Error::new(
                        name.span,
                        "shcomplex requires a sheaf followed by alternating maps and sheaves",
                    ));
                }
                let components = entries
                    .iter()
                    .step_by(2)
                    .map(|label| NamedSheaf::resolve(&self.model, label))
                    .collect::<Result<_>>()?;
                let maps = entries
                    .iter()
                    .skip(1)
                    .step_by(2)
                    .map(|label| Ok((label.clone(), self.model.sheaf_map(label)?)))
                    .collect::<Result<_>>()?;
                let complex = sheaf_hyper::Complex::build(&mut self.model, components, maps)?;
                self.model.define_sheaf_complex(name.text, complex);
                return Ok(None);
            }
            StatementForm::SheafPageDimension(name, grade) => {
                let complex = self.model.sheaf_total(&name)?;
                let rank = complex.filtered().page(&mut self.model, grade)?.rank();
                format!(
                    "shedim[{}, {}, {}, {}] = {rank}",
                    name.text, grade.r, grade.p, grade.q
                )
            }
            StatementForm::SheafHyperDimension(name, degree) => {
                let data = self.model.sheaf_total(&name)?;
                let rank = data.rank(&mut self.model, degree)?;
                format!("shhyperdim[{}, {degree}] = {rank}", name.text)
            }
            StatementForm::SheafDimension(name, degree) => {
                let data = self.model.sheaf(&name)?;
                let rank = data.rank(&mut self.model, degree)?;
                format!("shdim[{}, {degree}] = {rank}", name.text)
            }
            StatementForm::SetTopology(command) => return self.set_topology_command(command),
            StatementForm::SetOmega(name, site) => {
                self.model.check_new_set_family(&name)?;
                let site = self.model.named(&site.text).ok_or_else(|| {
                    Error::new(site.span, format!("unknown named site '{}'", site.text))
                })?;
                let (diagram, omega) = set_logic::Omega::build(&mut self.model, site)?;
                self.model.define_set_logic(
                    &name,
                    diagram,
                    set_logic::Origin::Omega(Rc::new(omega)),
                )?;
                return Ok(None);
            }
            StatementForm::SetSubobject(name, predicate) => {
                self.model.check_new_set_family(&name)?;
                let predicate = self.checker().expr(predicate)?;
                let (family, _, omega) = self.checker().set_predicate(predicate.ty)?;
                let domain = self
                    .model
                    .set_family_named(&family)
                    .expect("predicate domain");
                let Value::SetNat(rows) = self.eval(&predicate)? else {
                    unreachable!()
                };
                let (diagram, inclusion) =
                    set_logic::Subobject::build(&mut self.model, family, &domain, &omega, &rows)?;
                self.model.define_set_logic(
                    &name,
                    diagram,
                    set_logic::Origin::Subobject(Rc::new(inclusion)),
                )?;
                return Ok(None);
            }
            StatementForm::SetJ(name, family) => {
                self.model.check_new_set_family(&name)?;
                let source = self.model.set_family(&family)?;
                let (diagram, reflection) =
                    set_j::Sheafification::build(&mut self.model, family.text.into(), &source)?;
                self.model
                    .define_set_sheafification(&name, diagram, reflection)?;
                return Ok(None);
            }
            StatementForm::SetJCheck(family) => {
                self.model.set_family(&family)?;
                let is_sheaf = self
                    .model
                    .set_sheaf_certificate(Rc::from(family.text.as_str()))?
                    .is_some();
                let answer = if is_sheaf { "yes" } else { "no" };
                format!("setjcheck[{}] = {answer}", family.text)
            }
            StatementForm::SetSieve(name, point, index, dense) => {
                let omega = self.model.set_omega(&name.text)?;
                let x = self.model.point(omega.base, &point)?;
                if dense {
                    let answer = if omega.dense(&mut self.model, x, index)? {
                        "yes"
                    } else {
                        "no"
                    };
                    format!(
                        "setdense[{}, {}, {index}] = {answer}",
                        name.text, point.text
                    )
                } else {
                    let elements = omega.members(&mut self.model, x, index)?.join(", ");
                    format!(
                        "setsieve[{}, {}, {index}] = {{{elements}}}",
                        name.text, point.text
                    )
                }
            }
            StatementForm::SetPullback(name, family, map) => {
                self.model.check_new_set_family(&name)?;
                let source = self.model.set_family(&family)?;
                let map = self.checker().expr(map)?;
                let Kind::Arrow(domain, range) = *self.model.kind(map.ty) else {
                    return Err(Error::new(
                        map.span,
                        "setpullback requires a monotone map of named sites",
                    ));
                };
                if !matches!(self.model.kind(domain), Kind::Space(_)) || range != source.base {
                    return Err(Error::new(
                        map.span,
                        "setpullback requires a named domain and the family's site as codomain",
                    ));
                }
                let Value::Map(rows) = self.eval(&map)? else {
                    unreachable!()
                };
                let (diagram, origin) = dependent_sets::Pullback::build(
                    &mut self.model,
                    domain,
                    family.text.into(),
                    &source,
                    &rows,
                )?;
                self.model.define_dependent_set(
                    &name,
                    diagram,
                    dependent_sets::Origin::Pullback(Rc::new(origin)),
                )?;
                return Ok(None);
            }
            StatementForm::SetKan(name, family, map, right) => {
                self.model.check_new_set_family(&name)?;
                let body = self.model.set_family(&family)?;
                let map = self.checker().expr(map)?;
                let operation = if right { "setran" } else { "setlan" };
                let Kind::Arrow(domain, range) = *self.model.kind(map.ty) else {
                    return Err(Error::new(
                        map.span,
                        format!("{operation} requires a monotone map of named sites"),
                    ));
                };
                if domain != body.base || !matches!(self.model.kind(range), Kind::Space(_)) {
                    return Err(Error::new(map.span, format!(
                        "{operation} requires the input family's site as domain and a named codomain"
                    )));
                }
                let Value::Map(rows) = self.eval(&map)? else {
                    unreachable!()
                };
                let along = crate::set_kan::Along::new(
                    &mut self.model,
                    family.text.into(),
                    domain,
                    range,
                    &rows,
                )?;
                let (diagram, origin) = if right {
                    let (diagram, origin) =
                        crate::set_kan::Ran::build(&mut self.model, along, &body)?;
                    (diagram, dependent_sets::Origin::Ran(Rc::new(origin)))
                } else {
                    let (diagram, origin) =
                        crate::set_kan::Lan::build(&mut self.model, along, &body)?;
                    (diagram, dependent_sets::Origin::Lan(Rc::new(origin)))
                };
                self.model.define_dependent_set(&name, diagram, origin)?;
                return Ok(None);
            }
            StatementForm::SetDependent(name, total, family, product) => {
                self.model.check_new_set_family(&name)?;
                let (total, elements) = self.model.set_elements(&total)?;
                let body = self.model.set_family(&family)?;
                if body.base != total {
                    return Err(Error::new(
                        family.span,
                        "dependent family must live on the declared elements space",
                    ));
                }
                let base = self
                    .model
                    .set_family_named(&elements.family)
                    .expect("elements family");
                let projection = dependent_sets::Projection {
                    total,
                    elements,
                    body: family.text.into(),
                };
                let (diagram, origin) = if product {
                    let (diagram, origin) =
                        dependent_sets::Product::build(&mut self.model, projection, &base, &body)?;
                    (diagram, dependent_sets::Origin::Product(Rc::new(origin)))
                } else {
                    let (diagram, origin) =
                        dependent_sets::Sum::build(&mut self.model, projection, &base, &body)?;
                    (diagram, dependent_sets::Origin::Sum(Rc::new(origin)))
                };
                self.model.define_dependent_set(&name, diagram, origin)?;
                return Ok(None);
            }
            StatementForm::SetFamily(family, site, fibers, transitions) => {
                self.model.check_new_set_family(&family)?;
                let site = self.model.named(&site.text).ok_or_else(|| {
                    Error::new(site.span, format!("unknown named site '{}'", site.text))
                })?;
                let diagram =
                    set_family::Diagram::declare(&mut self.model, site, fibers, transitions)?;
                self.model.define_set_family(&family, diagram)?;
                return Ok(None);
            }
            StatementForm::SetElements(space, family) => {
                self.model.define_set_elements(&space, &family)?;
                return Ok(None);
            }
            StatementForm::SetNatural(name, source, target, components) => {
                self.check_new_binding(&name)?;
                let a = self.model.set_family(&source)?;
                let b = self.model.set_family(&target)?;
                if a.base != b.base {
                    return Err(Error::plain(
                        "SetNat requires two families on the same named site",
                    ));
                }
                let value = a.natural(&mut self.model, &b, components)?;
                let ty = self.model.register_set_nat_type(
                    source.text.into(),
                    target.text.into(),
                    true,
                )?;
                self.bind(name.text, ty, value)?;
                return Ok(None);
            }
            StatementForm::SetSections(family) => {
                let diagram = self.model.set_family(&family)?;
                let count = diagram.section_count(&mut self.model)?;
                format!("setsections[{}] = {count}", family.text)
            }
            StatementForm::SetMonodromy(name, source, target, point) => {
                self.model.check_new_set_monodromy(&name)?;
                let action = crate::set_monodromy::Monodromy::build(
                    &mut self.model,
                    &source,
                    &target,
                    &point,
                )?;
                self.model.define_set_monodromy(&name, action)?;
                return Ok(None);
            }
            StatementForm::SetMonoQuery(name, query) => {
                use crate::set_monodromy::Query;
                let action = self.model.set_monodromy(&name)?;
                match query {
                    Query::Size => format!("setmonosize[{}] = {}", name.text, action.size()),
                    Query::Characters => {
                        let histogram = action.characters(&mut self.model)?;
                        format!("setmonochars[{}] = {histogram:?}", name.text)
                    }
                    Query::Element(index) => {
                        let (source, target) = action.element(index)?;
                        self.model.spend(source.len() + target.len())?;
                        format!(
                            "setmonoat[{}, {index}] = {source:?} / {target:?}",
                            name.text
                        )
                    }
                    Query::Stabilizer(index) => {
                        let (size, source, target) = action.stabilizer(&mut self.model, index)?;
                        format!(
                            "setmonostabilizer[{}, {index}] = [{size}, {source}, {target}]",
                            name.text
                        )
                    }
                    Query::Maps => {
                        let count = action.natural_count(&mut self.model)?;
                        format!("setmonomaps[{}] = {count}", name.text)
                    }
                    Query::Orbitals => {
                        let sizes = action.orbital_sizes();
                        self.model.spend(sizes.len())?;
                        format!("setmonoorbitals[{}] = {sizes:?}", name.text)
                    }
                }
            }
            StatementForm::SetLinearize(name, family, prime) => {
                self.model.check_new_sheaf(&name)?;
                let diagram = self.model.set_family(&family)?;
                let prime = crate::sheaf::field(&mut self.model, prime)?;
                let sheaf = crate::set_monodromy::linearize(&mut self.model, &diagram, prime)?;
                self.model.define_set_linearization(name, &family, sheaf)?;
                return Ok(None);
            }
            StatementForm::SetOrbitalMap(name, action, source, target, coefficients) => {
                self.model.check_new_sheaf_map(&name)?;
                let action = self.model.set_monodromy(&action)?;
                self.model
                    .require_set_linearization(&source, &action.source)?;
                self.model
                    .require_set_linearization(&target, &action.target)?;
                let source = NamedSheaf::resolve(&self.model, &source)?;
                let target = NamedSheaf::resolve(&self.model, &target)?;
                let map = action.orbital_map(&mut self.model, source, target, coefficients)?;
                self.model.define_sheaf_map(name.text, map);
                return Ok(None);
            }
            StatementForm::SheafMapRank(name, point) => {
                let map = self.model.sheaf_map(&name)?;
                let index = self.model.point(map.source.sheaf.base(), &point)?;
                let rank = crate::set_monodromy::map_rank(&mut self.model, &map, index)?;
                format!("shmaprank[{}, {}] = {rank}", name.text, point.text)
            }
            StatementForm::Space(name, chains) => {
                self.model.define_space(&name, &chains)?;
                return Ok(None);
            }
            StatementForm::Let(name, expr) => {
                self.check_new_binding(&name)?;
                let expr = self.checker().expr(expr)?;
                let value = self.eval(&expr)?;
                self.bind(name.text, expr.ty, value)?;
                return Ok(None);
            }
            StatementForm::RVTree(report, modulus, root, radius, input) => {
                let input = self.checker().expr(input)?;
                let plan =
                    rv_tree::plan(&mut self.model, input.ty, modulus, root, radius, &report)?;
                let Value::RVLaurent(value) = self.eval(&input)? else {
                    unreachable!()
                };
                rv_tree::evaluate(
                    &mut self.model,
                    plan,
                    &value,
                    &report,
                    modulus,
                    root,
                    radius,
                )?
            }
            StatementForm::CycloCoefficient(power, input) => {
                let input = self.checker().expr(input)?;
                if !matches!(self.model.kind(input.ty), Kind::RVCyclo(..)) {
                    return Err(Error::new(
                        input.span,
                        "rvCycloCoeff requires a cyclotomic integer",
                    ));
                }
                let Value::RVCyclo(terms) = self.eval(&input)? else {
                    unreachable!()
                };
                format!("rvCycloCoeff = {}", terms.get(power).copied().unwrap_or(0))
            }
            StatementForm::ConformalSeries(report, rank, level, punctures) => {
                let grade = rv_fusion::Grade::checked(rank, level)?;
                rv_conformal::report(&mut self.model, grade, punctures, report)?
            }
            StatementForm::ConformalModulo(rank, level, genus, modulus, punctures) => {
                let grade = rv_fusion::Grade::checked(rank, level)?;
                let value =
                    rv_conformal::modular(&mut self.model, grade, genus, modulus, punctures)?;
                format!("rvConformalMod = {value}")
            }
            StatementForm::ConformalExact(rank, level, genus, punctures) => {
                let grade = rv_fusion::Grade::checked(rank, level)?;
                let value = rv_conformal::exact(&mut self.model, grade, genus, punctures)?;
                format!("rvConformalExact = {value}")
            }
            StatementForm::ConformalDimension(rank, level, genus, punctures) => {
                let grade = rv_fusion::Grade::checked(rank, level)?;
                let dimension =
                    rv_fusion::conformal_dimension(&mut self.model, grade, genus, punctures)?;
                format!("rvConformalDim = {dimension}")
            }
            StatementForm::FusionCoefficient(weight, input) => {
                let input = self.checker().expr(input)?;
                let Kind::RVFusion(grade) = *self.model.kind(input.ty) else {
                    return Err(Error::new(
                        input.span,
                        "rvFusionCoeff requires a fusion element",
                    ));
                };
                let Value::RVFusion(terms) = self.eval(&input)? else {
                    unreachable!()
                };
                let coefficient = rv_fusion::coefficient(&mut self.model, grade, &terms, weight)?;
                format!("rvFusionCoeff = {coefficient}")
            }
            StatementForm::WzCheck(mirror, certificate, input) => {
                let operation = if mirror.is_some() {
                    "rvCTWZCheck"
                } else {
                    "rvWZCheck"
                };
                let input = self.checker().expr(input)?;
                if !matches!(self.model.kind(input.ty), Kind::RVOre) {
                    return Err(Error::new(
                        input.span,
                        format!("{operation} requires an Ore operator"),
                    ));
                }
                let Value::RVOre(terms) = self.eval(&input)? else {
                    unreachable!()
                };
                let label = if let Some(mirror) = mirror {
                    let period = self.model.rv_period(&mirror.period)?;
                    crate::rv_period::wz::check(
                        &mut self.model,
                        &period,
                        &certificate,
                        mirror.columns.as_deref(),
                    )?;
                    "certified Laurent identity and recurrence"
                } else {
                    "certified"
                };
                crate::rv_wz::verify(&mut self.model, &terms, &certificate)?;
                let start = certificate.start;
                let description = if start == 1 {
                    "coefficient"
                } else {
                    "coefficients"
                };
                format!("{operation} = {label} ({start} initial {description} checked)")
            }
            StatementForm::OreCoefficient(q, theta, input) => {
                let input = self.checker().expr(input)?;
                if !matches!(self.model.kind(input.ty), Kind::RVOre) {
                    return Err(Error::new(
                        input.span,
                        "rvOreCoeff requires an Ore operator",
                    ));
                }
                let Value::RVOre(terms) = self.eval(&input)? else {
                    unreachable!()
                };
                let coefficient = rv_ore::coefficient(&mut self.model, &terms, q, theta)?;
                format!("rvOreCoeff = {coefficient}")
            }
            StatementForm::QuantumRegularized(rank, ambient, degree) => {
                let grade = rv_qgr::Grade::checked(rank, ambient)?;
                let period = rv_qconnection::regularized_period(&mut self.model, grade, degree)?;
                format!("rvQGrRegularized = {period}")
            }
            StatementForm::QuantumRegularizedRhs(rank, ambient) => {
                let grade = rv_qgr::Grade::checked(rank, ambient)?;
                let rhs = rv_qregular::rhs(&mut self.model, grade)?;
                format!("rvQGrRegularizedRhs = {rhs}")
            }
            StatementForm::QuantumPeriod(rank, ambient, degree) => {
                let grade = rv_qgr::Grade::checked(rank, ambient)?;
                let period = rv_qconnection::period(&mut self.model, grade, degree)?;
                format!("rvQGrPeriod = {period}")
            }
            StatementForm::QuantumCurvature(report, prime, rank, ambient) => {
                let grade = rv_pcurv::ModGrade::checked(&mut self.model, prime, rank, ambient)?;
                rv_pcurv::report(&mut self.model, grade, report)?
            }
            StatementForm::QuantumCoefficient(degree, shape, input) => {
                let input = self.checker().expr(input)?;
                let kind = self.model.kind(input.ty).clone();
                let coefficient = match kind {
                    Kind::RVQGr(grade) => {
                        let Value::RVQGr(terms) = self.eval(&input)? else {
                            unreachable!()
                        };
                        rv_qgr::coefficient(&mut self.model, grade, &terms, degree, shape)?
                    }
                    Kind::RVQGrP(grade) => {
                        let Value::RVQGrP(terms) = self.eval(&input)? else {
                            unreachable!()
                        };
                        rv_pcurv::coefficient(&mut self.model, grade, &terms, degree, shape)?.into()
                    }
                    _ => {
                        return Err(Error::new(
                            input.span,
                            "rvQGrCoeff requires a quantum Grassmannian element",
                        ))
                    }
                };
                format!("rvQGrCoeff = {coefficient}")
            }
            StatementForm::GromovWitten(rank, ambient, degree, shapes) => {
                let grade = rv_qgr::Grade::checked(rank, ambient)?;
                let coefficient = rv_qgr::gw(&mut self.model, grade, degree, shapes)?;
                format!("rvQGrGW = {coefficient}")
            }
            StatementForm::HiveReport(weights, horn) => {
                let boundary = rv_hive::Boundary::checked(&mut self.model, weights)?;
                if horn {
                    match rv_hive::horn(&mut self.model, &boundary)? {
                        rv_hive::Horn::Admissible => "rvHorn = admissible".into(),
                        rv_hive::Horn::Total { expected, actual } => {
                            format!("rvHorn = total mismatch ({expected} != {actual})")
                        }
                        rv_hive::Horn::Violation {
                            indices: [i, j, k],
                            left,
                            right,
                        } => {
                            format!("rvHorn = {i:?}, {j:?}, {k:?} gives {left} < {right}")
                        }
                    }
                } else {
                    let count = rv_hive::count(&mut self.model, &boundary)?;
                    format!("rvHiveCount = {count}")
                }
            }
            StatementForm::HiveVertex(i, j, input) => {
                let input = self.checker().expr(input)?;
                let Kind::RVHive(boundary) = self.model.kind(input.ty) else {
                    return Err(Error::new(input.span, "rvHiveAt requires a hive"));
                };
                let boundary = boundary.clone();
                let Value::RVHive(interior) = self.eval(&input)? else {
                    unreachable!()
                };
                let value = rv_hive::vertex(&mut self.model, &boundary, &interior, i, j)?;
                format!("rvHiveAt[{i}, {j}] = {value}")
            }
            StatementForm::SatakeQuery(v, coordinates, input) => {
                let input = self.checker().expr(input)?;
                let name = if v.is_some() {
                    "rvSatakeAt"
                } else {
                    "rvSatakeCoeff"
                };
                let Kind::RVSatake(rank) = self.model.kind(input.ty) else {
                    return Err(Error::plain(format!(
                        "{name} requires a Satake Laurent polynomial"
                    )));
                };
                if coordinates.len() != *rank {
                    return Err(Error::plain(
                        "Satake coordinates must match the polynomial rank",
                    ));
                }
                let Value::RVSatake(terms) = self.eval(&input)? else {
                    unreachable!()
                };
                if let Some(v) = v {
                    let numeric = rv_satake::evaluate(&mut self.model, &terms, v, &coordinates)?;
                    format!("{name} = {numeric}")
                } else {
                    let coefficient =
                        rv_satake::coefficient(&mut self.model, &terms, &coordinates)?;
                    self.model.show_satake_coefficient(coefficient)?
                }
            }
            StatementForm::HeckeCoefficient(spherical, parameter, index, input) => {
                let input = self.checker().expr(input)?;
                let name = match (spherical, parameter.is_some()) {
                    (true, true) => "rvSphAt",
                    (true, false) => "rvSphCoeff",
                    (false, true) => "rvHeckeAt",
                    (false, false) => "rvHeckeCoeff",
                };
                let rank = match (spherical, self.model.kind(input.ty)) {
                    (true, Kind::RVSphHecke(rank)) | (false, Kind::RVIwHecke(rank)) => *rank,
                    (true, _) => {
                        return Err(Error::plain(format!(
                            "{name} requires a spherical Hecke element"
                        )))
                    }
                    (false, _) => {
                        return Err(Error::plain(format!(
                            "{name} requires an affine Hecke element"
                        )))
                    }
                };
                if index.len() != rank {
                    let label = if spherical { "coweight" } else { "window" };
                    return Err(Error::plain(format!(
                        "Hecke coefficient {label} rank must match the algebra"
                    )));
                }
                if spherical {
                    rv_sphhecke::validate(&mut self.model, &index)?;
                } else {
                    rv_iwahori::validate_window(&mut self.model, &index)?;
                }
                let value = self.eval(&input)?;
                let coefficient = match &value {
                    Value::RVIwHecke(support) => {
                        rv_iwhecke::coefficient(&mut self.model, support, &index)?
                    }
                    Value::RVSphHecke(support) => {
                        rv_sphhecke::coefficient(&mut self.model, support, &index)?
                    }
                    _ => unreachable!(),
                };
                if let Some(parameter) = parameter {
                    let result = rv_iwhecke::evaluate(&mut self.model, coefficient, parameter)?;
                    format!("{name} = {result}")
                } else {
                    let result = rv_iwhecke::coefficients(&mut self.model, coefficient)?;
                    self.model.show_integer_coefficients(name, &result)?
                }
            }
            StatementForm::RVIwKL(lower, upper) => {
                let polynomial = rv_iwkl::polynomial(&mut self.model, &lower, &upper)?;
                format!("rvIwKL = {polynomial:?}")
            }
            StatementForm::RVIwSchubertKL(upper, input) => {
                let input = self.checker().expr(input)?;
                let name = "rvIwSchubertKL";
                rv_hodge::square(&self.model, input.ty, name)?;
                rv_iwahori::check_window(&mut self.model, input.ty, &upper)?;
                let Value::RVLaurent(value) = self.eval(&input)? else {
                    unreachable!()
                };
                let lower = rv_iwahori::standard(&mut self.model, input.ty, &value, name)?;
                let polynomial = rv_iwkl::polynomial(&mut self.model, &lower, &upper)?;
                format!("{name} = {polynomial:?}")
            }
            StatementForm::RVIw(length, input) => {
                let name = if length { "rvIwLength" } else { "rvIw" };
                let input = self.checker().expr(input)?;
                rv_hodge::square(&self.model, input.ty, name)?;
                let Value::RVLaurent(value) = self.eval(&input)? else {
                    unreachable!()
                };
                let window = rv_iwahori::standard(&mut self.model, input.ty, &value, name)?;
                if length {
                    let length = rv_iwahori::length(&mut self.model, &window)?;
                    format!("{name} = {length}")
                } else {
                    format!("{name} = {window:?}")
                }
            }
            StatementForm::RVIwPosition(arguments) => {
                let [x, y] = *arguments;
                let x = self.checker().expr(x)?;
                let y = self.checker().expr(y)?;
                rv_hodge::compatible(&self.model, x.ty, y.ty, "rvIwPosition")?;
                let values = [self.eval(&x)?, self.eval(&y)?];
                let [Value::RVLaurent(frame), Value::RVLaurent(other)] = values else {
                    unreachable!()
                };
                let window = rv_iwahori::position(&mut self.model, x.ty, &frame, &other)?;
                format!("rvIwPosition = {window:?}")
            }
            StatementForm::RVIwAt(modulus, root, arguments) => {
                let [a, x] = *arguments;
                let a = self.checker().expr(a)?;
                let x = self.checker().expr(x)?;
                let step =
                    rv_hodge::relative_step(&mut self.model, a.ty, x.ty, modulus, root, "rvIwAt")?;
                let values = [self.eval(&a)?, self.eval(&x)?];
                let [Value::RVLaurent(action), Value::RVLaurent(frame)] = values else {
                    unreachable!()
                };
                let window =
                    rv_iwahori::at(&mut self.model, a.ty, &action, &frame, step, "rvIwAt")?;
                format!("rvIwAt[{modulus}, {root}] = {window:?}")
            }
            StatementForm::TlLinking(strands, word, framings) => {
                let matrix = rv_tl::linking(&mut self.model, strands, &word, &framings)?;
                format!("rvTLLinking = {matrix:?}")
            }
            StatementForm::RVHodge(input) => {
                let input = self.checker().expr(input)?;
                rv_hodge::square(&self.model, input.ty, "rvHodge")?;
                let Value::RVLaurent(value) = self.eval(&input)? else {
                    unreachable!()
                };
                let exponents = rv_hodge::hodge(&mut self.model, input.ty, &value)?;
                format!("rvHodge = {exponents:?}")
            }
            StatementForm::RVPosition(arguments) => {
                let [x, y] = *arguments;
                let x = self.checker().expr(x)?;
                let y = self.checker().expr(y)?;
                rv_hodge::compatible(&self.model, x.ty, y.ty, "rvPosition")?;
                let values = [self.eval(&x)?, self.eval(&y)?];
                let [Value::RVLaurent(frame), Value::RVLaurent(other)] = values else {
                    unreachable!()
                };
                let exponents = rv_hodge::relative(&mut self.model, x.ty, &frame, &other)?;
                format!("rvPosition = {exponents:?}")
            }
            StatementForm::RVMuAt(modulus, root, arguments) => {
                let [a, x] = *arguments;
                let a = self.checker().expr(a)?;
                let x = self.checker().expr(x)?;
                let step =
                    rv_hodge::relative_step(&mut self.model, a.ty, x.ty, modulus, root, "rvMuAt")?;
                let values = [self.eval(&a)?, self.eval(&x)?];
                let [Value::RVLaurent(action), Value::RVLaurent(frame)] = values else {
                    unreachable!()
                };
                let exponents =
                    rv_hodge::mu_at(&mut self.model, a.ty, &action, &frame, step, "rvMuAt")?;
                format!("rvMuAt[{modulus}, {root}] = {exponents:?}")
            }
            StatementForm::RVNewton(modulus, root, input) => {
                let input = self.checker().expr(input)?;
                let (_, embedding) = rv_newton::plan(
                    &mut self.model,
                    input.ty,
                    modulus,
                    root,
                    rv_newton::Operation::Slopes,
                )?;
                let Value::RVLaurent(support) = self.eval(&input)? else {
                    unreachable!()
                };
                let slopes = rv_newton::slopes(&mut self.model, input.ty, &embedding, &support)?;
                format!("rvNewton[{modulus}, {root}] = {{{slopes}}}")
            }
            StatementForm::RVInvariantFactors(modulus, root, input) => {
                let input = self.checker().expr(input)?;
                let (_, embedding) = rv_shintani::plan(
                    &mut self.model,
                    input.ty,
                    modulus,
                    root,
                    rv_shintani::Operation::Factors,
                )?;
                let Value::RVHom(support) = self.eval(&input)? else {
                    unreachable!()
                };
                let factors = rv_shintani::invariant_factors(
                    &mut self.model,
                    input.ty,
                    &embedding,
                    &support,
                )?;
                format!("rvInvariantFactors[{modulus}, {root}] = {factors:?}")
            }
            StatementForm::Show(expr) => {
                let expr = self.checker().expr(expr)?;
                let value = self.eval(&expr)?;
                self.model
                    .show(expr.ty, &value)
                    .map_err(|error| error.at(expr.span))?
            }
            StatementForm::Card(expr) => {
                let expr = self.checker().expr(expr)?;
                let sierpinski = self.model.named("S").unwrap();
                if !matches!(self.model.kind(expr.ty), Kind::Arrow(_, range) if *range == sierpinski)
                {
                    return Err(Error::new(
                        expr.span,
                        format!(
                            "card requires a predicate of type X -> S, found {}",
                            self.model.type_name(expr.ty),
                        ),
                    ));
                }
                let Value::Map(rows) = self.eval(&expr)? else {
                    unreachable!()
                };
                self.model.spend(rows.len())?;
                let count = rows
                    .iter()
                    .filter(|row| matches!(row, Value::Point(1)))
                    .count();
                format!("card = {count} / {}", rows.len())
            }
            StatementForm::Assert(left, right) => {
                let (left_expr, right_expr) = self.comparable_exprs(left, right, "assertion")?;
                let left = self.eval(&left_expr)?;
                let right = self.eval(&right_expr)?;
                if !self.model.equal(&left, &right)? {
                    let left = self.model.describe(left_expr.ty, &left);
                    let right = self.model.describe(right_expr.ty, &right);
                    return Err(Error::plain(format!("assertion failed: {left} != {right}")));
                }
                return Ok(None);
            }
            StatementForm::Betti(ty) => {
                let ty = self.checker().ty(ty)?;
                let ranks = betti(&mut self.model, ty)?;
                format!("betti({}) = {ranks:?}", self.model.type_name(ty))
            }
            StatementForm::BpRelExt(spec, t) => {
                let shape = self.checker().bp_relative_shape(spec)?;
                let group = shape.domain().group(&mut self.model, t)?;
                let (p, n, m, a) = (
                    shape.upper.prime,
                    shape.lower,
                    shape.upper.height,
                    shape.upper.arrows,
                );
                format!(
                    "bprelext[{p}, {n}, {m}, {a}, {t}] = free {}, torsion {:?}",
                    group.free, group.torsion
                )
            }
            StatementForm::BpExt(crate::syntax::BpSpec(p, n, a), t) => {
                let shape = bp::Shape::new(&mut self.model, p, n, a)?;
                let group = bp_ext::group(&mut self.model, shape, t)?;
                format!(
                    "bpext[{p}, {n}, {a}, {t}] = free {}, torsion {:?}",
                    group.free, group.torsion,
                )
            }
            StatementForm::BpTripleDefined(inputs) => {
                let (layout, inputs) = self
                    .checker()
                    .bp_triple_inputs(*inputs, "bpTripleDefined")?;
                let values = self.eval_bp_triple_inputs(&inputs)?;
                let answer =
                    layout.defined(&mut self.model, values.each_ref().map(|(ty, x)| (*ty, x)))?;
                format!("bpTripleDefined = {}", if answer { "yes" } else { "no" })
            }
            StatementForm::BpTripleContains(pair) => {
                let [bracket, candidate] = *pair;
                let mut checker = self.checker();
                let bracket = checker.expr(bracket)?;
                let grade = checker.bp_triple_grade(&bracket, "bpTripleContains")?;
                let candidate = checker.expr(candidate)?;
                let ty = checker.model.bp_ext_type(grade)?;
                checker.require_type(ty, &candidate, "bpTripleContains")?;
                let Value::BPTriple(coset) = self.eval(&bracket)? else {
                    unreachable!()
                };
                let Value::BPExt(candidate) = self.eval(&candidate)? else {
                    unreachable!()
                };
                let answer = coset.contains(&mut self.model, &candidate)?;
                format!("bpTripleContains = {}", if answer { "yes" } else { "no" })
            }
            StatementForm::BpTripleSize(bracket) => {
                let mut checker = self.checker();
                let bracket = checker.expr(bracket)?;
                checker.bp_triple_grade(&bracket, "bpTripleSize")?;
                let Value::BPTriple(coset) = self.eval(&bracket)? else {
                    unreachable!()
                };
                format!("bpTripleSize = {}", coset.size())
            }
            StatementForm::BpOrder(t, input) => {
                let input = self.checker().expr(input)?;
                let relative = match *self.model.kind(input.ty) {
                    Kind::BPRel(shape) => Some(shape),
                    _ => None,
                };
                let order = if let Some(shape) = relative {
                    let Value::BPRel(cycle) = self.eval(&input)? else {
                        unreachable!()
                    };
                    shape.domain().order(&mut self.model, t, &cycle)?
                } else {
                    let shape = self.checker().bp_shape(&input, "bporder")?;
                    let Value::BP(cycle) = self.eval(&input)? else {
                        unreachable!()
                    };
                    bp_ext::order(&mut self.model, shape, t, &cycle)?
                };
                match order {
                    Some(order) => format!("bporder[{t}] = {order}"),
                    None => format!("bporder[{t}] = infinite"),
                }
            }
            StatementForm::PExtDimension(p, s, t) => {
                let prime = Prime::new(&mut self.model, p)?;
                let ty = self.model.odd_ext_type(prime, s, t)?;
                let rank = self.model.odd_ext_rank(ty)?;
                format!("pextdim[{p}, {s}, {t}] = {rank}")
            }
            StatementForm::BpHeightPage(image, (p, n, r, h, a, t)) => {
                let top = bp::Shape::new(&mut self.model, p, n, a)?;
                let grade = bp_height::Grade { top, r, h, t }.validate()?;
                if image {
                    let label = format!("bpHeightImage[{p}, {n}, {r}, {h}, {a}, {t}]");
                    match bp_height::image(&mut self.model, grade)? {
                        bp_height::Image::Finite(order) => format!("{label} = order {order}"),
                        bp_height::Image::Free(rank) => {
                            format!("{label} = free {rank}, order infinite")
                        }
                    }
                } else {
                    let group = bp_height::page(&mut self.model, grade)?;
                    format!(
                        "bpHeightPage[{p}, {n}, {r}, {h}, {a}, {t}] = free {}, torsion {:?}",
                        group.free, group.torsion
                    )
                }
            }
            StatementForm::BpModExtDimension(spec) => {
                let grade = self.checker().bp_mod_ext_grade(spec)?;
                let ty = self.model.bp_mod_ext_type(grade)?;
                let rank = self.model.bp_mod_ext(ty)?.rank();
                let shape = grade.shape;
                let (p, n, a, t) = (shape.prime, shape.height, shape.arrows, grade.t);
                if let Some(j) = shape.thickness {
                    let h = shape.depth + 1;
                    format!("bpThickExtdim[{p}, {n}, {a}, {h}, {j}, {t}] = {rank}")
                } else {
                    let h = shape.depth;
                    format!("bpIExtdim[{p}, {n}, {a}, {h}, {t}] = {rank}")
                }
            }
            StatementForm::BpVPageDimension(spec @ (p, n, r, a, h, t, f)) => {
                let grade = self.checker().bp_vpage_grade(spec)?;
                let ty = self.model.bp_vpage_type(grade)?;
                let rank = self.model.bp_vpage(ty)?.rank();
                format!("bpVdim[{p}, {n}, {r}, {a}, {h}, {t}, {f}] = {rank}")
            }
            StatementForm::BpPageDimension(mode, (p, n, r, a, t, f)) => {
                let shape = bp::Shape::new(&mut self.model, p, n, a)?;
                let grade = bp_page::Grade {
                    mode,
                    shape,
                    page: r,
                    t,
                    filtration: f,
                };
                let rank = bp_page::rank(&mut self.model, grade)?;
                let name = mode.dimension_name();
                format!("{name}[{p}, {n}, {r}, {a}, {t}, {f}] = {rank}")
            }
            StatementForm::PExtAnDimension(p, n, s, t) => {
                let prime = Prime::new(&mut self.model, p)?;
                let ty = self
                    .model
                    .profile_odd_ext_type(prime, Profile::An(n), s, t)?;
                let rank = self.model.odd_ext_rank(ty)?;
                format!("pextandim[{p}, {n}, {s}, {t}] = {rank}")
            }
            StatementForm::PMayDimension(p, r, s, t, f) => {
                let prime = Prime::new(&mut self.model, p)?;
                let ty = self.model.pmay_type(prime, r, s, t, f)?;
                let rank = self.model.odd_page(ty)?.rank();
                format!("pmaydim[{p}, {r}, {s}, {t}, {f}] = {rank}")
            }
            StatementForm::PMayAnDimension(p, n, r, s, t, f) => {
                let prime = Prime::new(&mut self.model, p)?;
                let ty = self
                    .model
                    .profile_pmay_type(prime, Profile::An(n), r, s, t, f)?;
                let rank = self.model.odd_page(ty)?.rank();
                format!("pmayandim[{p}, {n}, {r}, {s}, {t}, {f}] = {rank}")
            }
            StatementForm::PCeAnDimension(p, n, r, s, t, f) => {
                let prime = Prime::new(&mut self.model, p)?;
                let ty = self.model.pce_type(prime, n, r, s, t, f)?;
                let rank = self.model.odd_page(ty)?.rank();
                format!("pceandim[{p}, {n}, {r}, {s}, {t}, {f}] = {rank}")
            }
            StatementForm::RhoDimension(grade) => {
                let ty = self.model.rho_page_type(grade)?;
                let rank = self.model.rho_page(ty)?.rank();
                let rho_bock::Grade {
                    chain: motivic::Grade { s, t, w },
                    r,
                    f,
                } = grade;
                format!("rhobockdim[{r}, {s}, {t}, {w}, {f}] = {rank}")
            }
            StatementForm::RealBars(s, k) => {
                let group = self.model.real_bars(s, k)?;
                let bars = group.describe(&mut self.model)?;
                format!("rmbars[{s}, {k}] = {bars}")
            }
            StatementForm::RealBarRank(s, k, t, power) => {
                let group = self.model.real_bars(s, k)?;
                let rank = group.shift_rank(&mut self.model, t, power)?;
                format!("rmbarRank[{s}, {k}, {t}, {power}] = {rank}")
            }
            StatementForm::RhoOrder(input) => {
                let input = self.checker().expr(*input)?;
                let grade = self.checker().real_grade(&input, "rhoorder")?;
                let (terms, known_cycle) = self.real_terms(&input)?;
                match real_ext::order(&mut self.model, grade, &terms, known_cycle)? {
                    Some(exponent) => format!("rhoorder = {exponent}"),
                    None => "rhoorder = infinite".to_string(),
                }
            }
            StatementForm::RhoDepth(input) => {
                let input = self.checker().expr(*input)?;
                let kind = self.model.kind(input.ty);
                let (grade, is_bar) = match *kind {
                    Kind::RMExt(grade) => (grade, false),
                    Kind::RMBar(grade) => (grade, true),
                    _ => {
                        return Err(self.checker().invalid_operand(
                            &input,
                            "rhodepth",
                            "an RMExt class or RMBar coordinates",
                        ));
                    }
                };
                let Value::Class(coords) = self.eval(&input)? else {
                    unreachable!()
                };
                let depth = if is_bar {
                    real_extensions::depth_of_support(&mut self.model, grade, &coords)?
                } else {
                    real_extensions::depth(&mut self.model, (input.ty, &coords))?
                };
                match depth {
                    Some(depth) => format!("rhodepth = {depth}"),
                    None => "rhodepth = infinite".to_owned(),
                }
            }
            StatementForm::RhoHidden(pair) => {
                let [a, b] = *pair;
                let mut checker = self.checker();
                let a = checker.expr(a)?;
                let b = checker.expr(b)?;
                checker.rmext_grade(&a, "rhohidden")?;
                checker.rmext_grade(&b, "rhohidden")?;
                let (Value::Class(x), Value::Class(y)) = (self.eval(&a)?, self.eval(&b)?) else {
                    unreachable!()
                };
                real_extensions::hidden(&mut self.model, (a.ty, &x), (b.ty, &y))?
            }
            StatementForm::RealTripleDim(bracket) => {
                let bracket = self.checker().expr(*bracket)?;
                self.checker().rm_triple_grade(&bracket, "rmTripleDim")?;
                let Value::RMTriple(coset) = self.eval(&bracket)? else {
                    unreachable!();
                };
                format!("rmTripleDim = {}", coset.dimension())
            }
            StatementForm::RealTripleContains(pair) => {
                let [bracket, candidate] = *pair;
                let mut checker = self.checker();
                let bracket = checker.expr(bracket)?;
                let grade = checker.rm_triple_grade(&bracket, "rmTripleContains")?;
                let candidate = checker.expr(candidate)?;
                let candidate_grade = checker.real_grade(&candidate, "rmTripleContains")?;
                if candidate_grade != grade {
                    let motivic::Grade { s, t, w } = grade;
                    return Err(checker.invalid_operand(
                        &candidate,
                        "rmTripleContains",
                        &format!("an RMCobar cochain or RMExt class in degree ({s}, {t}, {w})"),
                    ));
                }
                let Value::RMTriple(coset) = self.eval(&bracket)? else {
                    unreachable!();
                };
                let coords = match self.eval(&candidate)? {
                    Value::Class(coords) => coords,
                    Value::RMCobar(terms) => {
                        let ty = self.model.rmext_type(grade)?;
                        let Value::Class(coords) = real_ext::class(&mut self.model, ty, &terms)?
                        else {
                            unreachable!();
                        };
                        coords
                    }
                    _ => unreachable!(),
                };
                let answer = coset.contains(&mut self.model, &coords)?;
                format!("rmTripleContains = {}", if answer { "yes" } else { "no" })
            }
            StatementForm::RMExtDimension(grade) => {
                let ty = self.model.rmext_type(grade)?;
                let rank = self.model.real_ext(ty)?.rank();
                let motivic::Grade { s, t, w } = grade;
                format!("rmextdim[{s}, {t}, {w}] = {rank}")
            }
            StatementForm::RealRhoDimension(grade) => {
                let ty = self.model.real_rho_ext_type(grade)?;
                let rank = self.model.real_ext(ty)?.rank();
                let real_rho::Grade {
                    exponent,
                    chain: motivic::Grade { s, t, w },
                } = grade;
                format!("rmrhoextdim[{exponent}, {s}, {t}, {w}] = {rank}")
            }
            StatementForm::MExtDimension(grade) => {
                let ty = self.model.mext_type(grade)?;
                let rank = self.model.adams(ty)?.rank();
                let motivic::Grade { s, t, w } = grade;
                format!("mextdim[{s}, {t}, {w}] = {rank}")
            }
            StatementForm::ExtDimension(s, t) => {
                let ty = self.model.ext_type(s, t)?;
                let rank = self.model.adams(ty)?.rank();
                format!("extdim[{s}, {t}] = {rank}")
            }
            StatementForm::ExtAnDimension(n, s, t) => {
                let ty = self.model.profile_ext_type(Profile::An(n), s, t)?;
                let rank = self.model.adams(ty)?.rank();
                format!("extandim[{n}, {s}, {t}] = {rank}")
            }
            StatementForm::XExtDimension(base, s, t) => {
                let base = self.checker().ty(base)?;
                let ty = self.model.xext_type(base, s, t)?;
                let rank = self.model.xext(ty)?.rank();
                format!("xextdim[{}, {s}, {t}] = {rank}", self.model.type_name(base))
            }
            StatementForm::MapExtDimension(source, target, s, t) => {
                let source = self.checker().ty(source)?;
                let target = self.checker().ty(target)?;
                let ty = self.model.map_ext_type(source, target, s, t)?;
                let rank = self.model.map_ext(ty)?.rank();
                format!(
                    "mapextdim[{}, {}, {s}, {t}] = {rank}",
                    self.model.type_name(source),
                    self.model.type_name(target),
                )
            }
            StatementForm::RvIdeal(name, selected, expressions) => {
                check_new_ideal(&self.model, &name)?;
                let selected = q_koszul::indices(&mut self.model, selected)?;
                let generators = self.rv_generator_values(expressions, &selected, "rvideal")?;
                let ideal = RvIdeal::build(&mut self.model, selected, generators)?;
                self.model.define_rv_ideal(name.text, ideal);
                return Ok(None);
            }
            StatementForm::RvWittCover(name, variables, expressions) => {
                self.model.check_new_rv_witt_cover(&name)?;
                let forms =
                    self.rv_witt_coordinates(variables, expressions, &name, "rvwittcover")?;
                let cover = RvWittCover::build(&mut self.model, forms)?;
                self.model.define_rv_witt_cover(name.text, cover);
                return Ok(None);
            }
            StatementForm::RvWittQuery(query, name) => {
                let cover = self.model.rv_witt_cover(&name)?;
                cover.report(&mut self.model, &name.text, query)?
            }
            StatementForm::RvWittTower(name, variables, expressions) => {
                self.model.check_new_rv_witt_tower(&name)?;
                let coordinates =
                    self.rv_witt_coordinates(variables, expressions, &name, "rvwitttower")?;
                let tower = RvWittTower::build(&mut self.model, coordinates)?;
                self.model.define_rv_witt_tower(name.text, tower);
                return Ok(None);
            }
            StatementForm::RvWittTowerQuery(query, name) => {
                let tower = self.model.rv_witt_tower(&name)?;
                tower.report(&mut self.model, &name.text, query)?
            }
            StatementForm::RvAsCover(name, variables, expressions) => {
                self.model.check_new_rv_as_cover(&name)?;
                let selected = q_koszul::indices(&mut self.model, variables)?;
                let forms = self.rv_generator_values(expressions, &selected, "rvascover")?;
                let cover = RvAsCover::build(&mut self.model, selected.len(), forms)?;
                self.model.define_rv_as_cover(name.text, cover);
                return Ok(None);
            }
            StatementForm::RvCrystal(name, prime, genus, entries) => {
                self.model
                    .register_rv_crystal(name, prime, genus, entries)?;
                return Ok(None);
            }
            StatementForm::RvCrystalQuery(name, query) => {
                self.model
                    .rv_crystal(&name)?
                    .report(&mut self.model, &name.text, query)?
            }
            StatementForm::EOReport(report, expression) => {
                let expression = self.checker().expr(expression)?;
                self.checker().require_eo(&expression, report.name())?;
                let Value::EOType(sequence) = self.eval(&expression)? else {
                    unreachable!()
                };
                report.evaluate(&mut self.model, &sequence)?
            }
            StatementForm::EOWitness(lower, upper) => {
                let lower = self.checker().expr(lower)?;
                let upper = self.checker().expr(upper)?;
                self.checker().require_eo(&lower, "eoClosureWitness")?;
                self.checker()
                    .require_type(lower.ty, &upper, "eoClosureWitness")?;
                let (Value::EOType(lower), Value::EOType(upper)) =
                    (self.eval(&lower)?, self.eval(&upper)?)
                else {
                    unreachable!()
                };
                if let Some(mut permutation) = eo::closure_witness(&mut self.model, &lower, &upper)?
                {
                    for image in &mut permutation {
                        *image += 1;
                    }
                    format!("eoClosureWitness = {permutation:?}")
                } else {
                    "eoClosureWitness = none".to_owned()
                }
            }
            StatementForm::RvAsWordQuery(query, name, words) => {
                let cover = self.model.rv_as_cover(&name)?;
                rv_as_derham::word_report(&mut self.model, &cover, &name.text, query, &words)?
            }
            StatementForm::RvAsTwo(name, prime, terms) => {
                self.model.register_rv_as_two_cover(name, prime, terms)?;
                return Ok(None);
            }
            StatementForm::RvAsTwoQuery(query, name) => {
                self.model
                    .rv_as_two_cover(&name)?
                    .report(&mut self.model, &name.text, query)?
            }
            StatementForm::RvPeriod(name, variables, terms) => {
                self.model.define_rv_period(name, variables, terms)?;
                return Ok(None);
            }
            StatementForm::RVCTExact(name, index) => {
                let period = self.model.rv_period(&name)?;
                let value = period.exact(&mut self.model, index)?;
                format!("rvCTExact[{}, {index}] = {value}", name.text)
            }
            StatementForm::RvAsOdd(name, prime, terms) => {
                self.model.check_new_rv_as_odd_cover(&name)?;
                let cover = RvAsOddCover::build(&mut self.model, prime, terms)?;
                self.model.define_rv_as_odd_cover(name.text, cover);
                return Ok(None);
            }
            StatementForm::RvAsOddQuery(query, name) => {
                self.model
                    .rv_as_odd_cover(&name)?
                    .report(&mut self.model, &name.text, query)?
            }
            StatementForm::RvQpQuad(name, p, d, r) => {
                self.model.register_qp_quad(name, p, d, r)?;
                return Ok(None);
            }
            StatementForm::ResReport(query, input) => {
                let input = self.checker().expr(input)?;
                self.checker().require_res_space(&input, query.name())?;
                let Value::RVResSpace(space) = self.eval(&input)? else {
                    unreachable!()
                };
                query.report(&mut self.model, &space)?
            }
            StatementForm::ResLeq(lower, upper) => {
                let lower = self.checker().expr(lower)?;
                let upper = self.checker().expr(upper)?;
                self.checker().require_res_space(&lower, "rvResLeq")?;
                self.checker().require_type(lower.ty, &upper, "rvResLeq")?;
                let (Value::RVResSpace(lower), Value::RVResSpace(upper)) =
                    (self.eval(&lower)?, self.eval(&upper)?)
                else {
                    unreachable!()
                };
                format!("rvResLeq = {}", lower.is_subspace(&mut self.model, &upper)?)
            }
            StatementForm::StableReport(query, input) => {
                let input = self.checker().expr(input)?;
                let (accepts, expected) = if query.requires_tree() {
                    (
                        matches!(self.model.kind(input.ty), Kind::RVStableTree(..)),
                        "an RVStableTree",
                    )
                } else {
                    (
                        matches!(self.model.kind(input.ty), Kind::RVStableMarks(..)),
                        "an RVStableMarks",
                    )
                };
                if !accepts {
                    return Err(Error::new(
                        input.span,
                        format!("stable reduction query requires {expected}"),
                    ));
                }
                match self.eval(&input)? {
                    Value::RVStableTree(tree) => query.tree(&mut self.model, &tree)?,
                    Value::RVStableMarks(curve) => query.curve(&mut self.model, &curve)?,
                    _ => unreachable!(),
                }
            }
            StatementForm::ValMatReport(query, input) => {
                let input = self.checker().expr(input)?;
                if !matches!(self.model.kind(input.ty), Kind::RVValMat(..)) {
                    return Err(Error::new(
                        input.span,
                        "tropical Plücker query requires an RVValMat",
                    ));
                }
                let Value::RVValMat(value) = self.eval(&input)? else {
                    unreachable!()
                };
                query.report(&mut self.model, &value)?
            }
            StatementForm::RvQpQuadReport(precision, input) => {
                let input = self.checker().expr(input)?;
                let operation = if precision.is_some() {
                    "rvQpQuadResidue"
                } else {
                    "rvQpQuadVal"
                };
                let field = self.checker().qp_quad_field(&input, operation)?;
                let Value::RVQpQuad(value) = self.eval(&input)? else {
                    unreachable!()
                };
                match precision {
                    Some(s) => format!(
                        "rvQpQuadResidue = {}",
                        field.residue(&mut self.model, value, s)?
                    ),
                    None => match field.valuation(&mut self.model, value)? {
                        Some(order) => format!("rvQpQuadVal = {order}"),
                        None => "rvQpQuadVal = infinity".to_owned(),
                    },
                }
            }
            StatementForm::RvPHodgeQuad(name, field, spec) => {
                let ty = self.model.qp_quad_type(&field)?;
                let field = self.model.qp_quad(&field)?;
                let mut read = |input| -> Result<QuadElement> {
                    let input = self.checker().expr(input)?;
                    self.checker().require_type(ty, &input, "rvphodgequad")?;
                    let Value::RVQpQuad(value) = self.eval(&input)? else {
                        unreachable!()
                    };
                    Ok(value)
                };
                let eigenvalues = spec
                    .eigenvalues
                    .into_iter()
                    .map(&mut read)
                    .collect::<Result<Vec<_>>>()?;
                let monodromy = spec
                    .monodromy
                    .into_iter()
                    .map(|(i, j, value)| Ok((i, j, read(value)?)))
                    .collect::<Result<Vec<_>>>()?;
                let filtration = spec
                    .filtration
                    .into_iter()
                    .map(|(weight, column)| {
                        Ok((
                            weight,
                            column
                                .into_iter()
                                .map(&mut read)
                                .collect::<Result<Vec<_>>>()?,
                        ))
                    })
                    .collect::<Result<Vec<_>>>()?;
                let spec = crate::rv_phodge::QuadraticSpec {
                    eigenvalues,
                    monodromy,
                    filtration,
                };
                self.model.register_rv_phodge_quad(name, field, spec)?;
                return Ok(None);
            }
            StatementForm::RvPHodge(name, spec) => {
                self.model.register_rv_phodge(name, spec)?;
                return Ok(None);
            }
            StatementForm::RvPHodgeQuery(name, query) => {
                self.model
                    .rv_phodge(&name)?
                    .report(&mut self.model, &name.text, query)?
            }
            StatementForm::RvAsQuery(query, name) => {
                let cover = self.model.rv_as_cover(&name)?;
                cover.report(&mut self.model, &name.text, query)?
            }
            StatementForm::RvCurve(name, variables, expressions) => {
                self.model.check_new_rv_curve(&name)?;
                let selected = q_koszul::indices(&mut self.model, variables)?;
                let forms = self.rv_generator_values(expressions, &selected, "rvcurve")?;
                let curve = RvCurve::build(&mut self.model, selected, forms)?;
                self.model.define_rv_curve(name.text, curve);
                return Ok(None);
            }
            StatementForm::RvFrobeniusOrders(name, exponent) => {
                let curve = self.model.rv_curve(&name)?;
                let orders = curve.frobenius_orders(&mut self.model, exponent)?;
                format!("rvFrobOrders[{}, {exponent}] = {orders:?}", name.text)
            }
            StatementForm::RvOrderAt(modulus, point, polynomial) => {
                let input = self.checker().expr(polynomial)?;
                let (selected, _) = self.checker().rv_poly_grade(&input, "rvOrdAt")?;
                if selected.len() != 2 {
                    return Err(Error::new(input.span, "rvOrdAt requires a binary form"));
                }
                let Value::RVPoly(polynomial) = self.eval(&input)? else {
                    unreachable!()
                };
                let order = rv_curve::order_at(&mut self.model, &polynomial, modulus, &point)?;
                let order = order.map_or_else(|| "infinite".into(), |order| order.to_string());
                format!("rvOrdAt[{modulus}, {point:?}] = {order}")
            }
            StatementForm::RvOrders(name) => {
                let curve = self.model.rv_curve(&name)?;
                format!("rvorders[{}] = {:?}", name.text, curve.orders)
            }
            StatementForm::RvOsculating(name, curve, k, variables) => {
                check_new_ideal(&self.model, &name)?;
                let curve = self.model.rv_curve(&curve)?;
                let selected = q_koszul::indices(&mut self.model, variables)?;
                let ideal = curve.image(&mut self.model, k, selected)?;
                self.model.define_rv_ideal(name.text, ideal);
                return Ok(None);
            }
            StatementForm::RvGauss(operation, name, dual, input) => {
                check_new_ideal(&self.model, &name)?;
                let input = self.checker().expr(input)?;
                let (selected, degree) = self.checker().rv_poly_grade(&input, operation.name())?;
                let Value::RVPoly(polynomial) = self.eval(&input)? else {
                    unreachable!()
                };
                let ideal = rv_gauss::build(
                    &mut self.model,
                    operation,
                    selected,
                    degree,
                    &polynomial,
                    dual,
                )?;
                self.model.define_rv_ideal(name.text, ideal);
                return Ok(None);
            }
            StatementForm::RvRees(operation, name, original, variables, expressions) => {
                check_new_ideal(&self.model, &name)?;
                let original = self.model.rv_ideal(&original)?;
                let generators =
                    self.rv_generator_values(expressions, &original.selected, operation.name())?;
                let ideal =
                    rv_rees::build(&mut self.model, operation, &original, variables, generators)?;
                self.model.define_rv_ideal(name.text, ideal);
                return Ok(None);
            }
            StatementForm::RvFitting(name, spec, size) => {
                check_new_ideal(&self.model, &name)?;
                let base = self.checker().q_endpoint(*spec.base)?;
                let selected = q_koszul::indices(&mut self.model, spec.selected)?;
                let minors = self.model.rv_minors(base, selected.clone(), size)?;
                let generators = minors.polynomials().map(|p| (size, p.clone())).collect();
                let ideal = RvIdeal::build(&mut self.model, selected, generators)?;
                self.model.define_rv_ideal(name.text, ideal);
                return Ok(None);
            }
            StatementForm::RvIntersect(name, left, right) => {
                check_new_ideal(&self.model, &name)?;
                let left = self.model.rv_ideal(&left)?;
                let right = self.model.rv_ideal(&right)?;
                let ideal = rv_scheme::intersect(&mut self.model, &left, &right)?;
                self.model.define_rv_ideal(name.text, ideal);
                return Ok(None);
            }
            StatementForm::RvScheme(operation, name, original, input) => {
                check_new_ideal(&self.model, &name)?;
                let original = self.model.rv_ideal(&original)?;
                let name_of_operation = match operation {
                    RvSchemeOp::Colon => "rvcolon",
                    RvSchemeOp::Saturate => "rvsaturate",
                };
                let input = self.checker().expr(input)?;
                let (selected, _) = self.checker().rv_poly_grade(&input, name_of_operation)?;
                if selected != original.selected {
                    return Err(Error::new(
                        input.span,
                        format!(
                            "{name_of_operation} polynomial uses a different primitive selection"
                        ),
                    ));
                }
                let Value::RVPoly(polynomial) = self.eval(&input)? else {
                    unreachable!()
                };
                let ideal = match operation {
                    RvSchemeOp::Colon => rv_scheme::colon(&mut self.model, &original, polynomial)?,
                    RvSchemeOp::Saturate => {
                        rv_scheme::saturate(&mut self.model, &original, &polynomial)?
                    }
                };
                self.model.define_rv_ideal(name.text, ideal);
                return Ok(None);
            }
            StatementForm::RvProjective(name, original) => {
                check_new_ideal(&self.model, &name)?;
                let original = self.model.rv_ideal(&original)?;
                let ideal = rv_scheme::projective(&mut self.model, &original)?;
                self.model.define_rv_ideal(name.text, ideal);
                return Ok(None);
            }
            StatementForm::RvGroebner(name) => {
                let ideal = self.model.rv_ideal(&name)?;
                let prefix = format!("rvgroebner[{}]", name.text);
                self.model
                    .show_rv_polynomials(&prefix, &ideal.selected, ideal.polynomials())?
            }
            StatementForm::RvHilbert(name, degree) => {
                let ideal = self.model.rv_ideal(&name)?;
                let count = ideal.hilbert(&mut self.model, degree)?;
                format!("rvhilbert[{}, {degree}] = {count}", name.text)
            }
            StatementForm::RvBarPath(name, ideal, steps) => {
                self.model.define_rv_bar_path(name, ideal, steps)?;
                return Ok(None);
            }
            StatementForm::RvCyclicDimension(spec, mode) => {
                let ordinary = self.checker().rv_hh_grade(spec)?;
                let ty = self.model.rv_cyclic_type(
                    rv_cyclic_homology::Grade {
                        ordinary: ordinary.clone(),
                        mode,
                    },
                    true,
                )?;
                let rank = self.model.rv_cyclic_homology(ty)?.rank();
                format!(
                    "{}[{}, {}, {:?}] = {rank}",
                    mode.dimension_name(),
                    ordinary.ideal,
                    ordinary.b,
                    ordinary.weight
                )
            }
            StatementForm::RvNCDimension(spec) => {
                let grade = self.checker().rv_negative_grade(spec)?;
                let ty = self.model.rv_negative_type(grade.clone(), true)?;
                let rank = self.model.rv_negative_homology(ty)?.rank();
                format!(
                    "rvNCdim[{}, {}, {}, {:?}] = {rank}",
                    grade.ordinary.ideal, grade.length, grade.ordinary.b, grade.ordinary.weight
                )
            }
            StatementForm::RvHHDimension(spec) => {
                let grade = self.checker().rv_hh_grade(spec)?;
                let ty = self.model.rv_hh_class_type(grade.clone())?;
                let rank = self.model.rv_hh_homology(ty)?.rank();
                format!(
                    "rvHHdim[{}, {}, {:?}] = {rank}",
                    grade.ideal, grade.b, grade.weight
                )
            }
            StatementForm::RvPathDimension(spec) => {
                let grade = self.checker().rv_bar_grade(spec)?;
                let ty = self.model.rv_path_type(grade, true)?;
                let rank = self.model.rv_path_ext(ty)?.rank();
                format!(
                    "rvpathextdim[{}, {}, {}, {}] = {rank}",
                    self.model.rv_bar_path_name(grade.path),
                    grade.start,
                    grade.end,
                    grade.b
                )
            }
            StatementForm::RvBarPageDimension(spec) => {
                let grade = self.checker().rv_bar_page_grade(spec)?;
                let ty = self.model.rv_bar_page_type(grade)?;
                let rank = self.model.rv_bar_page(ty)?.rank();
                format!(
                    "rvbaredim[{}, {}, {}, {}, {}, {}] = {rank}",
                    self.model.rv_bar_path_name(grade.chain.path),
                    grade.chain.start,
                    grade.chain.end,
                    grade.r,
                    grade.chain.b,
                    grade.length
                )
            }
            StatementForm::RvBarDimension(spec, length) => {
                let grade = self.checker().rv_bar_grade(spec)?;
                let rank = rv_bar::dimension(&mut self.model, grade, length)?;
                let name = self.model.rv_bar_path_name(grade.path);
                let rv_bar::Grade { start, end, b, .. } = grade;
                if let Some(length) = length {
                    format!("rvbardim[{name}, {start}, {end}, {b}, {length}] = {rank}")
                } else {
                    format!("rvbardim[{name}, {start}, {end}, {b}] = {rank}")
                }
            }
            StatementForm::RvHochster(name, vertices, degree) => {
                let rank = rv_hochster::rank(&mut self.model, &name, &vertices, degree)?;
                format!("rvhochster[{}, {vertices:?}, {degree}] = {rank}", name.text)
            }
            StatementForm::RvBetti(spec) => {
                let grade = self.checker().rv_koszul_grade(spec)?;
                let ty = self.model.rv_complex_type(grade.clone(), true)?;
                let rank = self.model.rv_tor(ty)?.rank();
                if let Some(multidegree) = grade.multidegree {
                    format!(
                        "rvmbetti[{}, {}, {multidegree:?}] = {rank}",
                        grade.ideal, grade.s
                    )
                } else {
                    format!(
                        "rvbetti[{}, {}, {}] = {rank}",
                        grade.ideal, grade.s, grade.t
                    )
                }
            }
            StatementForm::RvMinors(spec, size) => {
                let base = self.checker().q_endpoint(*spec.base)?;
                let selected = q_koszul::indices(&mut self.model, spec.selected)?;
                self.model.show_rv_minors(base, selected, size)?
            }
            StatementForm::JetDimension(spec) => {
                let name = spec.chain.path.text.clone();
                let grade = self.checker().jet_page_grade(spec)?;
                let ty = self.model.jet_page_type(grade)?;
                let rank = self.model.jet_page(ty)?.rank();
                let source = self.model.type_name(grade.chain.diagram.source);
                let target = self.model.type_name(grade.chain.diagram.target);
                format!(
                    "jetdim[{source}, {target}, {name}, {}, {}, {}] = {rank}",
                    grade.r, grade.chain.parity, grade.f
                )
            }
            StatementForm::JetPath(name, spec) => {
                let selected = q_koszul::indices(&mut self.model, spec.selected)?;
                let path = rank_jets::Path::new(
                    &mut self.model,
                    &selected,
                    spec.modulus,
                    spec.order,
                    &spec.coefficients,
                )?;
                self.model.define_jet_path(name, path)?;
                return Ok(None);
            }
            StatementForm::RankJet(spec) => {
                let source = self.checker().q_endpoint(*spec.source)?;
                let target = self.checker().q_endpoint(*spec.target)?;
                let from = self.model.type_name(source);
                let to = self.model.type_name(target);
                let (label, report) = match spec.path {
                    JetPathArg::Inline(spec) => {
                        let selected = q_koszul::indices(&mut self.model, spec.selected)?;
                        let modulus = spec.modulus;
                        let order = spec.order;
                        let coefficients = spec.coefficients;
                        let report = rank_jets::compute(
                            &mut self.model,
                            source,
                            target,
                            &selected,
                            modulus,
                            order,
                            &coefficients,
                        )?;
                        (
                            format!("{selected:?}, {modulus}, {order}, {coefficients:?}"),
                            report,
                        )
                    }
                    JetPathArg::Named(label) => {
                        let id = self.model.jet_path_id(&label)?;
                        let report =
                            self.model
                                .jet_path(id)
                                .report(&mut self.model, source, target)?;
                        (label.text, report)
                    }
                };
                let rank_jets::Report { dimension, blocks } = report;
                format!("rvjet[{from}, {to}, {label}] = dimension {dimension}, blocks {blocks:?}")
            }
            StatementForm::RankOrbit(modulus, point) => {
                let orbit = rv_arithmetic::orbit(&mut self.model, modulus, &point)?;
                format!(
                    "rvorbit[{modulus}, {point:?}] = degree {}, orbit {orbit:?}",
                    orbit.len(),
                )
            }
            StatementForm::RankArithmetic(series, spec) => {
                let source = self.checker().q_endpoint(*spec.source)?;
                let target = self.checker().q_endpoint(*spec.target)?;
                let selected = q_koszul::indices(&mut self.model, spec.selected)?;
                let max_degree = spec.max_degree;
                let counts = rv_arithmetic::evaluate(
                    &mut self.model,
                    source,
                    target,
                    &selected,
                    max_degree,
                    series,
                )?;
                format!(
                    "{}[{}, {}, {selected:?}, {max_degree}] = {counts:?}",
                    series.name(),
                    self.model.type_name(source),
                    self.model.type_name(target),
                )
            }
            StatementForm::RankVariety(spec) => {
                let source = self.checker().q_endpoint(*spec.source)?;
                let target = self.checker().q_endpoint(*spec.target)?;
                let selected = q_koszul::indices(&mut self.model, spec.selected)?;
                let from = self.model.type_name(source);
                let to = self.model.type_name(target);
                let modulus = spec.modulus;
                if let Some(point) = spec.point {
                    let dimension = rank_variety::dimension(
                        &mut self.model,
                        source,
                        target,
                        &selected,
                        modulus,
                        &point,
                    )?;
                    format!("rvdim[{from}, {to}, {selected:?}, {modulus}, {point:?}] = {dimension}")
                } else {
                    let counts = rank_variety::spectrum(
                        &mut self.model,
                        source,
                        target,
                        &selected,
                        modulus,
                    )?;
                    format!("rvspectrum[{from}, {to}, {selected:?}, {modulus}] = {counts:?}")
                }
            }
            StatementForm::QExtDimension(spec) => {
                let grade = self.checker().q_grade(spec)?;
                let ty = self.model.q_ext_type(grade.clone())?;
                let rank = self.model.q_ext(ty)?.rank();
                format!(
                    "qextdim[{}, {}, {:?}, {}, {}] = {rank}",
                    self.model.type_name(grade.source),
                    self.model.type_name(grade.target),
                    grade.selected,
                    grade.s,
                    grade.t
                )
            }
            StatementForm::VExtDimension(spec) => {
                let grade = self.checker().q_grade(spec)?;
                let ty = self.model.v_ext_type(grade.clone())?;
                let rank = self.model.q_ext(ty)?.rank();
                format!(
                    "vextdim[{}, {}, {:?}, {}, {}] = {rank}",
                    self.model.type_name(grade.source),
                    self.model.type_name(grade.target),
                    grade.selected,
                    grade.s,
                    grade.t
                )
            }
            StatementForm::VBockDimension(spec) => {
                let ty = self.checker().vbock_type(spec)?;
                let Kind::VBock(grade) = self.model.kind(ty).clone() else {
                    unreachable!()
                };
                let rank = self.model.vbock_page(ty)?.rank();
                let localized = grade.localized;
                format!(
                    "vbockdim[{}, {}, {:?}, {}, {}, {}, {}] = {rank}",
                    self.model.type_name(localized.source),
                    self.model.type_name(localized.target),
                    localized.selected,
                    grade.r,
                    localized.s,
                    localized.t,
                    grade.f,
                )
            }
            StatementForm::VBars(spec) => {
                let family = self.checker().vbar_family(spec)?;
                let bars = self.model.v_bars(family)?;
                let description = bars.report(&mut self.model)?;
                format!(
                    "vbars[{}, {}, {:?}] = {description}",
                    self.model.type_name(family.source),
                    self.model.type_name(family.target),
                    family.indices,
                )
            }
            StatementForm::VBarDimension(spec) => {
                let grade = self.checker().q_grade(spec)?;
                let family = v_bars::Family::new(grade.source, grade.target, &grade.selected)?;
                let bars = self.model.v_bars(family)?;
                let rank = bars.rank(&mut self.model, family, grade.s, grade.t)?;
                format!(
                    "vbardim[{}, {}, {:?}, {}, {}] = {rank}",
                    self.model.type_name(family.source),
                    self.model.type_name(family.target),
                    family.indices,
                    grade.s,
                    grade.t,
                )
            }
            StatementForm::MapMayDimension(spec) => {
                let ty = self.checker().map_may_type(spec)?;
                let Kind::MapMay(source, target, r, s, t, f) = *self.model.kind(ty) else {
                    unreachable!()
                };
                let rank = self.model.map_may(ty)?.rank();
                format!(
                    "mapmaydim[{}, {}, {r}, {s}, {t}, {f}] = {rank}",
                    self.model.type_name(source),
                    self.model.type_name(target),
                )
            }
            StatementForm::MayAnDimension(n, r, s, t, f) => {
                let ty = self.model.profile_may_type(Profile::An(n), r, s, t, f)?;
                let rank = self.model.may_page(ty)?.rank();
                format!("mayandim[{n}, {r}, {s}, {t}, {f}] = {rank}")
            }
            StatementForm::CeDimension(n, r, s, t, f) => {
                let ty = self.model.ce_type(n, r, s, t, f)?;
                let rank = self.model.ce_page(ty)?.rank();
                format!("ceandim[{n}, {r}, {s}, {t}, {f}] = {rank}")
            }
            StatementForm::MayDimension(r, s, t, f) => {
                let ty = self.model.may_type(r, s, t, f)?;
                let rank = self.model.may_page(ty)?.rank();
                format!("maydim[{r}, {s}, {t}, {f}] = {rank}")
            }
            StatementForm::Homotopy(left, right) => {
                let (left_expr, right_expr) = self.comparable_exprs(left, right, "homotopy")?;
                if !matches!(self.model.kind(left_expr.ty), Kind::Arrow(..)) {
                    return Err(Error::new(
                        left_expr.span,
                        format!(
                            "homotopy requires functions, found {}",
                            self.model.type_name(left_expr.ty)
                        ),
                    ));
                }
                let left = self.eval(&left_expr)?;
                let right = self.eval(&right_expr)?;
                match homotopy(&mut self.model, left_expr.ty, &left, &right)? {
                    Some(1) => "homotopy = yes (1 step)".into(),
                    Some(length) => format!("homotopy = yes ({length} steps)"),
                    None => "homotopy = no".into(),
                }
            }
        };
        Ok(Some(line))
    }

    fn rv_koszul_regrade(&mut self, source: TypeId, target: TypeId, value: Value) -> Result<Value> {
        let (Kind::RVKoszul(grade) | Kind::RVTor(grade)) = self.model.kind(target) else {
            unreachable!();
        };
        let grade = grade.clone();
        let ideal = self.model.rv_ideal_named(&grade.ideal).unwrap();
        let filter = |model: &mut Model, terms: &[rv_koszul::Term], mu: &[usize]| {
            model.spend(rv_koszul::cells(terms))?;
            Ok::<_, Error>(
                terms
                    .iter()
                    .filter(|term| rv_koszul::matches_multidegree(term, mu))
                    .cloned()
                    .collect::<Vec<_>>(),
            )
        };
        match value {
            Value::RVKoszul(terms) => match &grade.multidegree {
                None => Ok(Value::RVKoszul(terms)),
                Some(mu) => {
                    let selected = filter(&mut self.model, &terms, mu)?;
                    Ok(Value::RVKoszul(rv_koszul::persist(
                        &mut self.model,
                        selected,
                    )?))
                }
            },
            Value::Class(class) if class.is_empty() => Ok(cohomology::zero()),
            Value::Class(class) => {
                let mut representative = self
                    .model
                    .rv_tor(source)?
                    .representative_terms(&mut self.model, &class)?;
                if let Some(mu) = &grade.multidegree {
                    representative = filter(&mut self.model, &representative, mu)?;
                }
                self.rv_tor_class(target, &ideal, &representative)
            }
            _ => unreachable!(),
        }
    }

    fn rv_tor_class(
        &mut self,
        ty: TypeId,
        ideal: &RvIdeal,
        terms: &[rv_koszul::Term],
    ) -> Result<Value> {
        if terms.is_empty() {
            return Ok(cohomology::zero());
        }
        let coordinates = self
            .model
            .rv_tor(ty)?
            .project(&mut self.model, ideal, terms)?;
        self.model.allocate(coordinates.len())?;
        Ok(Value::Class(coordinates.into()))
    }

    fn quadratic_rows(&mut self, rows: &[Vec<TypedExpr>]) -> Result<Vec<Vec<QuadElement>>> {
        rows.iter()
            .map(|row| {
                row.iter()
                    .map(|entry| {
                        let Value::RVQpQuad(value) = self.eval(entry)? else {
                            unreachable!()
                        };
                        Ok(value)
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()
    }

    fn rv_witt_coordinates(
        &mut self,
        variables: Vec<usize>,
        expressions: Vec<ParsedExpr>,
        name: &Label,
        operation: &str,
    ) -> Result<Vec<rv_poly::Polynomial>> {
        let selected = q_koszul::indices(&mut self.model, variables)?;
        let &[variable] = selected.as_ref() else {
            return Err(Error::new(
                name.span,
                format!("{operation} requires exactly one parameter variable"),
            ));
        };
        expressions
            .into_iter()
            .map(|expression| self.rv_witt_polynomial(expression, variable, operation))
            .collect()
    }

    fn rv_witt_polynomial(
        &mut self,
        expression: ParsedExpr,
        variable: usize,
        operation: &str,
    ) -> Result<rv_poly::Polynomial> {
        let input = self.checker().expr(expression)?;
        let (variables, degree) = self.checker().rv_poly_grade(&input, operation)?;
        let Value::RVPoly(polynomial) = self.eval(&input)? else {
            unreachable!()
        };
        self.model.spend(rv_poly::cells(&polynomial))?;
        match variables.as_ref() {
            [index] if *index == variable => Ok(polynomial),
            [index, _] if *index == variable => Ok(polynomial.iter().map(|term| [term[0]].into()).collect()),
            [_, index] if *index == variable => Ok(rv_curve::affine_form(&polynomial)),
            [] if degree == 0 => Ok(polynomial.iter().map(|_| [0].into()).collect()),
            _ => Err(Error::new(input.span, format!("{operation} requires a univariate polynomial or a binary homogenization in the chosen variable"))),
        }
    }

    fn rv_generator_values(
        &mut self,
        expressions: Vec<ParsedExpr>,
        selected: &[usize],
        operation: &str,
    ) -> Result<Vec<(usize, rv_poly::Polynomial)>> {
        let mut generators = Vec::with_capacity(expressions.len());
        for expression in expressions {
            let input = self.checker().expr(expression)?;
            let (variables, degree) = self.checker().rv_poly_grade(&input, operation)?;
            if variables.as_ref() != selected {
                return Err(Error::new(
                    input.span,
                    format!("{operation} generator uses a different primitive selection"),
                ));
            }
            let Value::RVPoly(polynomial) = self.eval(&input)? else {
                unreachable!()
            };
            generators.push((degree, polynomial));
        }
        Ok(generators)
    }

    fn comparable_exprs(
        &mut self,
        left: ParsedExpr,
        right: ParsedExpr,
        context: &str,
    ) -> Result<(TypedExpr, TypedExpr)> {
        let mut checker = self.checker();
        let left = checker.expr(left)?;
        let right = checker.expr(right)?;
        checker.require_type(left.ty, &right, context)?;
        Ok((left, right))
    }

    fn eval_bp_triple_inputs(
        &mut self,
        [a, b, c]: &[TypedExpr; 3],
    ) -> Result<[(TypeId, bp_classes::Support); 3]> {
        let mut evaluate = |input: &TypedExpr| {
            let Value::BPExt(class) = self.eval(input)? else {
                unreachable!()
            };
            Ok((input.ty, class))
        };
        Ok([evaluate(a)?, evaluate(b)?, evaluate(c)?])
    }

    fn persist_bp(&mut self, ty: TypeId, mut value: BpPolynomial) -> Result<Value> {
        let Kind::BP(shape) = *self.model.kind(ty) else {
            unreachable!()
        };
        shape.normalize(&mut self.model, &mut value)?;
        bp::persist(&mut self.model, value)
    }

    fn sheaf_site(&mut self, name: &Label, operation: &str) -> Result<SheafSite> {
        let binding = self
            .bindings
            .get(&name.text)
            .ok_or_else(|| Error::new(name.span, format!("unknown named map '{}'", name.text)))?;
        let invalid = || {
            Error::new(
                name.span,
                format!("{operation} requires a named monotone map between named finite spaces"),
            )
        };
        let Kind::Arrow(source, target) = *self.model.kind(binding.ty) else {
            return Err(invalid());
        };
        if !matches!(self.model.kind(source), Kind::Space(_))
            || !matches!(self.model.kind(target), Kind::Space(_))
        {
            return Err(invalid());
        }
        let Value::Map(table) = &self.values[binding.index] else {
            unreachable!()
        };
        SheafSite::new(&mut self.model, source, target, table)
    }

    fn sheaf_prime(&self, ty: TypeId) -> Prime {
        match self.model.kind(ty) {
            Kind::Sheaf(_, name, _) => {
                self.model
                    .sheaf_named(name)
                    .expect("registered sheaf")
                    .prime
            }
            Kind::SheafHom(_, ends, _) => self.model.sheaf_hom(ends).prime(),
            Kind::SheafComplex(_, name, _)
            | Kind::SheafPage(name, _)
            | Kind::SheafLocal(name, _) => self
                .model
                .sheaf_total_named(name)
                .expect("registered sheaf complex")
                .filtered()
                .prime(),
            _ => unreachable!(),
        }
    }

    fn real_terms(&mut self, input: &TypedExpr) -> Result<(real_motivic::Cochain, bool)> {
        match self.eval(input)? {
            Value::RMCobar(terms) => Ok((terms, false)),
            Value::Class(coords) if coords.is_empty() => Ok(([].into(), true)),
            Value::Class(coords) => {
                let terms = self
                    .model
                    .real_ext(input.ty)?
                    .representative(&mut self.model, &coords)?;
                Ok((terms.into(), true))
            }
            _ => unreachable!(),
        }
    }

    fn check_new_binding(&self, name: &Label) -> Result<()> {
        if self.bindings.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("value '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    fn bind(&mut self, name: String, ty: TypeId, value: Value) -> Result<()> {
        self.model.allocate(1)?;
        self.bindings.insert(
            name,
            Binding {
                index: self.values.len(),
                ty,
            },
        );
        self.values.push(value);
        Ok(())
    }

    fn eval_mixed_induce(
        &mut self,
        types: [TypeId; 3],
        map: &[crate::prime_linear::Coordinate],
        class: &[crate::prime_linear::Coordinate],
    ) -> Result<Value> {
        let Kind::MixedChain(map_grade, true) = *self.model.kind(types[0]) else {
            unreachable!()
        };
        let Kind::MixedHomology(source) = *self.model.kind(types[1]) else {
            unreachable!()
        };
        let Kind::MixedHomology(target) = *self.model.kind(types[2]) else {
            unreachable!()
        };
        let image = mixed_homology::induce(
            &mut self.model,
            map_grade.complex,
            map,
            source,
            target,
            class,
        )?;
        odd_ext::persist(&mut self.model, image)
    }

    fn eval_mixed_h_unary(
        &mut self,
        op: &CobarOp,
        ty: TypeId,
        terms: Rc<[crate::prime_linear::Coordinate]>,
    ) -> Result<Value> {
        match op {
            CobarOp::MixHClass | CobarOp::MixHPrimitive => {
                if terms.is_empty() {
                    return Ok(Value::OddClass(terms));
                }
                let Kind::MixedChain(grade, true) = *self.model.kind(ty) else {
                    unreachable!()
                };
                let support = if matches!(op, CobarOp::MixHClass) {
                    self.model
                        .mixed_homology(grade)?
                        .project(&mut self.model, &terms)?
                } else {
                    mixed_homology::primitive(&mut self.model, grade, &terms)?
                };
                odd_ext::persist(&mut self.model, support)
            }
            CobarOp::MixHRep => {
                if terms.is_empty() {
                    return Ok(Value::OddClass(terms));
                }
                let Kind::MixedHomology(grade) = *self.model.kind(ty) else {
                    unreachable!()
                };
                let support = self
                    .model
                    .mixed_homology(grade)?
                    .representative(&mut self.model, &terms)?;
                odd_ext::persist(&mut self.model, support)
            }
            _ => unreachable!(),
        }
    }

    fn eval_mixed_cone_map(
        &mut self,
        cone: &mixed_cone::Cone,
        operation: &mixed_cone::Operation,
        input: &TypedExpr,
        output: TypeId,
    ) -> Result<Value> {
        let Value::OddClass(support) = self.eval(input)? else {
            unreachable!()
        };
        let image = match (self.model.kind(input.ty), self.model.kind(output)) {
            (&Kind::MixedHomology(source), &Kind::MixedHomology(target)) => {
                mixed_homology::map(&mut self.model, source, target, &support, |model, chain| {
                    cone.apply(model, *operation, chain)
                })?
            }
            _ => cone.apply(&mut self.model, *operation, &support)?,
        };
        odd_ext::persist(&mut self.model, image)
    }

    fn eval_mixed_h_shift(&mut self, input: &TypedExpr, output: TypeId) -> Result<Value> {
        let Value::OddClass(class) = self.eval(input)? else {
            unreachable!()
        };
        if class.is_empty() {
            return Ok(Value::OddClass(class));
        }
        let Kind::MixedHomology(source) = *self.model.kind(input.ty) else {
            unreachable!()
        };
        let Kind::MixedHomology(target) = *self.model.kind(output) else {
            unreachable!()
        };
        let source = self.model.mixed_homology(source)?;
        let target = self.model.mixed_homology(target)?;
        let support = source.shift(&mut self.model, &target, &class)?;
        odd_ext::persist(&mut self.model, support)
    }

    fn eval_sheaf_deformation(&mut self, form: &TypedForm) -> Result<Value> {
        match form {
            TypedForm::SheafDefOperation(hom, operation, input) => {
                let Value::SheafDeformation(jet) = self.eval(input)? else {
                    unreachable!()
                };
                match operation {
                    sheaf_deform::Operation::Term(0) => {
                        sheaf_deform::constant(&mut self.model, hom)
                    }
                    sheaf_deform::Operation::Term(k) => {
                        Ok(Value::SheafC(jet.coefficients()[k - 1].clone()))
                    }
                    sheaf_deform::Operation::Curvature => {
                        Ok(Value::SheafC(jet.curvature(&mut self.model, hom)?))
                    }
                    sheaf_deform::Operation::Obstruction => {
                        let curvature = jet.curvature(&mut self.model, hom)?;
                        hom.cohomology
                            .class(&mut self.model, hom.as_ref(), 2, &curvature)
                    }
                    sheaf_deform::Operation::Lift => {
                        Ok(Value::SheafDeformation(jet.lift(&mut self.model, hom)?))
                    }
                }
            }
            TypedForm::SheafDefAppend(hom, args) => {
                let Value::SheafDeformation(jet) = self.eval(&args[0])? else {
                    unreachable!()
                };
                let Value::SheafC(terms) = self.eval(&args[1])? else {
                    unreachable!()
                };
                Ok(Value::SheafDeformation(jet.append(
                    &mut self.model,
                    hom,
                    terms,
                )?))
            }
            TypedForm::SheafDefGauge(hom, exponent, witness, args) => {
                let Value::SheafDeformation(jet) = self.eval(&args[0])? else {
                    unreachable!()
                };
                let Value::SheafC(terms) = self.eval(&args[1])? else {
                    unreachable!()
                };
                jet.gauge(&mut self.model, hom, *exponent, &terms, *witness)
            }
            _ => unreachable!(),
        }
    }

    fn eval(&mut self, expr: &TypedExpr) -> Result<Value> {
        self.eval_inner(expr).map_err(|error| error.at(expr.span))
    }

    fn eval_inner(&mut self, expr: &TypedExpr) -> Result<Value> {
        self.model.spend(1)?;
        match &expr.kind {
            &TypedForm::Global(index) => Ok(self.values[index].clone()),
            &TypedForm::Local(index) => Ok(self.locals[index].clone()),
            &TypedForm::Point(index) => Ok(Value::Point(index)),
            TypedForm::Constant(value) => Ok(value.clone()),
            TypedForm::SetReflectExtend {
                reflection,
                source,
                target,
                certificate,
                map,
            } => {
                let Value::SetNat(map) = self.eval(map)? else {
                    unreachable!()
                };
                reflection.extend(&mut self.model, source, target, certificate, &map)
            }
            TypedForm::SetJExtend {
                reflection,
                source,
                target,
                certificate,
                map,
            } => {
                let Value::SetNat(rows) = self.eval(map)? else {
                    unreachable!()
                };
                reflection.extend(&mut self.model, source, target, certificate, &rows)
            }
            TypedForm::SetLogic(action, inputs) => {
                let rows = inputs
                    .iter()
                    .map(|expr| {
                        let Value::SetNat(value) = self.eval(expr)? else {
                            unreachable!()
                        };
                        Ok(value)
                    })
                    .collect::<Result<Vec<_>>>()?;
                match action {
                    LogicAction::Constant(omega, source, top) => {
                        omega.constant(&mut self.model, source, *top)
                    }
                    LogicAction::Connect(omega, source, operation) => omega.connect(
                        &mut self.model,
                        source,
                        *operation,
                        &rows[0],
                        rows.get(1).map(AsRef::as_ref),
                    ),
                    LogicAction::Classify(omega, source, target) => {
                        omega.classify(&mut self.model, source, target, &rows[0])
                    }
                    LogicAction::Quantify(omega, source, target, universal) => omega.quantify(
                        &mut self.model,
                        source,
                        target,
                        &rows[0],
                        &rows[1],
                        *universal,
                    ),
                }
            }
            TypedForm::SetKan(action, input) => {
                let rows = input.as_ref().map(|input| self.eval(input)).transpose()?;
                let rows = rows.as_ref().map(|value| {
                    let Value::SetNat(rows) = value else {
                        unreachable!()
                    };
                    rows.as_ref()
                });
                action.evaluate(&mut self.model, rows)
            }
            TypedForm::SetDependent(action, input) => {
                let rows = if let Some(input) = input {
                    let Value::SetNat(rows) = self.eval(input)? else {
                        unreachable!()
                    };
                    rows
                } else {
                    Rc::new([])
                };
                match action {
                    DependentAction::PullNatural(pullback, source) => {
                        pullback.natural(&mut self.model, source, &rows)
                    }
                    DependentAction::Evaluation(product) => product.evaluate(&mut self.model),
                    DependentAction::Transpose(product, source, lifted) => {
                        product.transpose(&mut self.model, source, lifted, &rows)
                    }
                    DependentAction::Untranspose(product, source, lifted) => {
                        product.untranspose(&mut self.model, source, lifted, &rows)
                    }
                    DependentAction::SumLift(sum, body) => sum.lift(&mut self.model, body, &rows),
                    DependentAction::SumRestrict(sum, source, body) => {
                        sum.restrict(&mut self.model, source, body, &rows)
                    }
                }
            }
            TypedForm::SetIdentity(diagram) => diagram.identity(&mut self.model),
            TypedForm::SetSection(elements, diagram, index) => {
                elements.section(&mut self.model, diagram, *index)
            }
            TypedForm::SetCompose(source, middle, inputs) => {
                let Value::SetNat(outer) = self.eval(&inputs[0])? else {
                    unreachable!()
                };
                let Value::SetNat(inner) = self.eval(&inputs[1])? else {
                    unreachable!()
                };
                source.compose(&mut self.model, middle, &outer, &inner)
            }
            TypedForm::SetTotal(elements, target, input) => {
                let Value::SetNat(rows) = self.eval(input)? else {
                    unreachable!()
                };
                elements.map(&mut self.model, &rows, target)
            }
            TypedForm::Binary(op, left, right) => {
                let first = self.eval(left)?;
                let second = self.eval(right)?;
                match (op, first, second) {
                    (BinaryOp::Add, Value::Class(a), Value::Class(b)) => {
                        cohomology::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::Operation(a), Value::Operation(b)) => {
                        steenrod::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::BPExt(a), Value::BPExt(b)) => {
                        bp_classes::add(&mut self.model, expr.ty, &a, &b)
                    }
                    (BinaryOp::Add, Value::BPRel(a), Value::BPRel(b)) => {
                        self.model.spend(a.cells())?;
                        let mut result = (*a).clone();
                        result.add_assign(&mut self.model, &b, false)?;
                        bp_rel::persist(&mut self.model, result)
                    }
                    (BinaryOp::Add, Value::BP(a), Value::BP(b)) => {
                        self.model.spend(a.cells())?;
                        let mut result = (*a).clone();
                        result.add_assign(&mut self.model, &b, false)?;
                        self.persist_bp(expr.ty, result)
                    }
                    (BinaryOp::BpCup, Value::BPExt(a), Value::BPExt(b)) => {
                        bp_classes::cup(&mut self.model, (left.ty, &a), (right.ty, &b), expr.ty)
                    }
                    (BinaryOp::BpCup, Value::BP(a), Value::BP(b)) => {
                        let (&Kind::BP(l), &Kind::BP(r)) =
                            (self.model.kind(left.ty), self.model.kind(right.ty))
                        else {
                            unreachable!()
                        };
                        let result = bp::cup(&mut self.model, l, r, &a, &b)?;
                        bp::persist(&mut self.model, result)
                    }
                    (BinaryOp::BpCup, Value::OddClass(a), Value::OddClass(b)) => {
                        if matches!(self.model.kind(expr.ty), Kind::BPModExt(..)) {
                            return bp_mod_ext::cup(
                                &mut self.model,
                                (left.ty, &a),
                                (right.ty, &b),
                                expr.ty,
                            );
                        }
                        if a.is_empty() || b.is_empty() {
                            return Ok(odd_ext::zero());
                        }
                        if matches!(self.model.kind(expr.ty), Kind::BPVPage(..)) {
                            return bp_vpage::cup(
                                &mut self.model,
                                (left.ty, &a),
                                (right.ty, &b),
                                expr.ty,
                            );
                        }
                        let source = self.model.bp_page(left.ty)?;
                        let other = self.model.bp_page(right.ty)?;
                        let target = self.model.bp_page(expr.ty)?;
                        let product = source.cup(&mut self.model, &other, &target, &a, &b)?;
                        odd_ext::persist(&mut self.model, product)
                    }
                    (BinaryOp::BpMultiply, Value::BP(a), Value::BP(b)) => {
                        let Kind::BP(shape) = *self.model.kind(expr.ty) else {
                            unreachable!()
                        };
                        let result = a.multiply_in(&mut self.model, &b, shape.reduction())?;
                        bp::persist(&mut self.model, result)
                    }
                    (BinaryOp::Add, Value::RVPoly(a), Value::RVPoly(b)) => {
                        let polynomial = rv_poly::add(&mut self.model, &a, &b)?;
                        Ok(Value::RVPoly(rv_poly::persist(
                            &mut self.model,
                            polynomial,
                        )?))
                    }
                    (BinaryOp::RvBCMul, Value::RVBarCobar(a), Value::RVBarCobar(b)) => Ok(
                        Value::RVBarCobar(rv_barcobar::product(&mut self.model, &a, &b)?),
                    ),
                    (BinaryOp::Add, Value::RVNCChain(a), Value::RVNCChain(b)) => {
                        Ok(Value::RVNCChain(rv_negative::add(&mut self.model, &a, &b)?))
                    }
                    (BinaryOp::Add, Value::RVHHChain(a), Value::RVHHChain(b)) => {
                        Ok(Value::RVHHChain(rv_hh::add(&mut self.model, &a, &b)?))
                    }
                    (BinaryOp::Add, Value::RVBarCobar(a), Value::RVBarCobar(b)) => Ok(
                        Value::RVBarCobar(rv_barcobar::add(&mut self.model, &a, &b)?),
                    ),
                    (BinaryOp::MixInduce, Value::OddClass(map), Value::OddClass(class)) => {
                        self.eval_mixed_induce([left.ty, right.ty, expr.ty], &map, &class)
                    }
                    (BinaryOp::MixEval, Value::OddClass(map), Value::OddClass(argument)) => {
                        let Kind::MixedChain(grade, _) = *self.model.kind(left.ty) else {
                            unreachable!()
                        };
                        let hom = self.model.mixed_hom(grade.complex)?;
                        let value = hom.evaluate(&mut self.model, &map, &argument)?;
                        odd_ext::persist(&mut self.model, value)
                    }
                    (BinaryOp::RvBarPair, Value::RVBar(dual), Value::RVBar(raw)) => {
                        Ok(Value::Point(usize::from(rv_bar_dual::pair(
                            &mut self.model,
                            &dual,
                            &raw,
                        )?)))
                    }
                    (BinaryOp::RvPathCup, Value::RVBar(a), Value::RVBar(b)) => {
                        let terms = rv_bar_dual::cup_terms(&mut self.model, &a, &b)?;
                        Ok(Value::RVBar(rv_bar::persist(&mut self.model, terms)?))
                    }
                    (BinaryOp::RvPathCup, Value::Class(a), Value::Class(b)) => {
                        rv_bar_dual::cup(&mut self.model, (left.ty, &a), (right.ty, &b), expr.ty)
                    }
                    (BinaryOp::Add, Value::RVBar(a), Value::RVBar(b)) => {
                        Ok(Value::RVBar(rv_bar::add(&mut self.model, &a, &b)?))
                    }
                    (BinaryOp::Add, Value::RVKoszul(a), Value::RVKoszul(b)) => {
                        Ok(Value::RVKoszul(rv_koszul::add(&mut self.model, &a, &b)?))
                    }
                    (BinaryOp::RvKoszulProduct, Value::RVKoszul(a), Value::RVKoszul(b)) => {
                        let Kind::RVKoszul(grade) = self.model.kind(expr.ty) else {
                            unreachable!();
                        };
                        let ideal = self.model.rv_ideal_named(&grade.ideal).unwrap();
                        Ok(Value::RVKoszul(rv_koszul::product(
                            &mut self.model,
                            &ideal,
                            &a,
                            &b,
                        )?))
                    }
                    (BinaryOp::RvTorProduct, Value::Class(a), Value::Class(b)) => {
                        if a.is_empty() || b.is_empty() {
                            return Ok(cohomology::zero());
                        }
                        let Kind::RVTor(grade) = self.model.kind(expr.ty) else {
                            unreachable!();
                        };
                        let ideal = self.model.rv_ideal_named(&grade.ideal).unwrap();
                        let a = self
                            .model
                            .rv_tor(left.ty)?
                            .representative(&mut self.model, &a)?;
                        let b = self
                            .model
                            .rv_tor(right.ty)?
                            .representative(&mut self.model, &b)?;
                        let terms = rv_koszul::product(&mut self.model, &ideal, &a, &b)?;
                        self.rv_tor_class(expr.ty, &ideal, &terms)
                    }
                    (BinaryOp::RvMultiply, Value::RVPoly(a), Value::RVPoly(b)) => {
                        let polynomial = rv_poly::multiply(&mut self.model, &a, &b)?;
                        Ok(Value::RVPoly(rv_poly::persist(
                            &mut self.model,
                            polynomial,
                        )?))
                    }
                    (BinaryOp::Add, Value::RhoLoc(a), Value::RhoLoc(b)) => {
                        rho_localized::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::RhoMin(a), Value::RhoMin(b)) => {
                        rho_minimal::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::RMCobar(a), Value::RMCobar(b)) => {
                        let terms = real_motivic::add(&mut self.model, &a, &b)?;
                        real_motivic::persist(&mut self.model, terms)
                    }
                    (BinaryOp::Add, Value::Cobar(a), Value::Cobar(b)) => {
                        cobar::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::SheafC(a), Value::SheafC(b)) => {
                        let prime = self.sheaf_prime(expr.ty);
                        sheaf::add(&mut self.model, prime, &a, &b)
                    }
                    (BinaryOp::AsOddTatePair, Value::OddClass(a), Value::OddClass(b)) => {
                        if a.is_empty() || b.is_empty() {
                            return Ok(Value::PrimeScalar(0));
                        }
                        let Kind::ASOddTate(name, _, first) = self.model.kind(left.ty) else {
                            unreachable!()
                        };
                        let Kind::ASOddTate(_, _, second) = self.model.kind(right.ty) else {
                            unreachable!()
                        };
                        let (first, second) = (*first, *second);
                        let cover = self
                            .model
                            .rv_as_odd_cover_named(name)
                            .expect("registered odd Artin–Schreier cover");
                        cover.tate(&mut self.model)?.pair(
                            &mut self.model,
                            &cover,
                            (first, &a),
                            (second, &b),
                        )
                    }
                    (BinaryOp::AsOddDeRhamPair, Value::OddClass(a), Value::OddClass(b)) => {
                        let Kind::ASOddDeRham(name, prime) = self.model.kind(left.ty) else {
                            unreachable!()
                        };
                        if a.is_empty() || b.is_empty() {
                            return Ok(Value::PrimeScalar(0));
                        }
                        let cover = self
                            .model
                            .rv_as_odd_cover_named(name)
                            .expect("registered odd Artin–Schreier cover");
                        debug_assert!(*prime == cover.prime());
                        cover
                            .derham(&mut self.model)?
                            .pair(&mut self.model, &cover, &a, &b)
                    }
                    (BinaryOp::Add, Value::PrimeScalar(a), Value::PrimeScalar(b)) => {
                        let Kind::PrimeScalar(prime) = *self.model.kind(expr.ty) else {
                            unreachable!()
                        };
                        Ok(Value::PrimeScalar(prime.add(a, b)))
                    }
                    (BinaryOp::Add, Value::OddClass(a), Value::OddClass(b)) => {
                        let p = match *self.model.kind(expr.ty) {
                            Kind::PExt(p, ..)
                            | Kind::PExtAn(p, ..)
                            | Kind::PMay(p, ..)
                            | Kind::PMayAn(p, ..)
                            | Kind::PCeAn(p, ..)
                            | Kind::ASOddDiff(_, p)
                            | Kind::ASTwoDiff(_, p)
                            | Kind::ASOddDeRham(_, p)
                            | Kind::ASOddTate(_, p, _) => p,
                            Kind::Sheaf(..)
                            | Kind::SheafComplex(..)
                            | Kind::SheafPage(..)
                            | Kind::SheafLocal(..)
                            | Kind::SheafHom(..) => self.sheaf_prime(expr.ty),
                            Kind::MixedChain(grade, _) | Kind::MixedHomology(grade) => {
                                self.model.mixed_complex(grade.complex).prime
                            }
                            Kind::MixedMinimum(grade, _) => {
                                self.model.mixed_minimal(grade.minimal).prime()
                            }
                            Kind::MixedPage(grade) => {
                                self.model.mixed_complex(grade.chain.complex).prime
                            }
                            Kind::BPPage(grade) => Prime::validated(grade.shape.prime),
                            Kind::BPVPage(grade) => grade.prime(),
                            Kind::BPModExt(grade) => grade.prime(),
                            _ => unreachable!(),
                        };
                        odd_ext::add(&mut self.model, p, &a, &b)
                    }
                    (BinaryOp::Add, Value::PCobar(a), Value::PCobar(b)) => {
                        let (p, _, _, _) = self.model.kind(expr.ty).pcobar_grade().unwrap();
                        odd_cobar::add(&mut self.model, p, &a, &b)
                    }
                    (BinaryOp::Add, Value::RVWzw(a), Value::RVWzw(b)) => {
                        rv_wzw::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::RvTl(a), Value::RvTl(b)) => {
                        rv_tl::add(&mut self.model, &a, &b)
                    }
                    (
                        operation @ (BinaryOp::ResSum | BinaryOp::ResIntersection),
                        Value::RVResSpace(a),
                        Value::RVResSpace(b),
                    ) => {
                        let value = match operation {
                            BinaryOp::ResSum => a.sum(&mut self.model, &b)?,
                            _ => a.intersection(&mut self.model, &b)?,
                        };
                        value.persist(&mut self.model)
                    }
                    (
                        operation @ (BinaryOp::Add | BinaryOp::QpQuadMultiply),
                        Value::RVQpQuad(a),
                        Value::RVQpQuad(b),
                    ) => {
                        let Kind::RVQpQuad(name) = self.model.kind(expr.ty) else {
                            unreachable!()
                        };
                        let field = self
                            .model
                            .qp_quad_named(name)
                            .expect("registered quadratic field");
                        let value = match operation {
                            BinaryOp::Add => field.add(&mut self.model, a, b)?,
                            _ => field.multiply(&mut self.model, a, b)?,
                        };
                        rv_qp_quad::persist(&mut self.model, value)
                    }
                    (BinaryOp::Add, Value::RVCycloQ(a), Value::RVCycloQ(b)) => {
                        rv_lens::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::RVCyclo(a), Value::RVCyclo(b)) => {
                        rv_cyclotomic::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::RVFusion(a), Value::RVFusion(b)) => {
                        let Kind::RVFusion(grade) = *self.model.kind(expr.ty) else {
                            unreachable!()
                        };
                        rv_fusion::add(&mut self.model, grade, &a, &b)
                    }
                    (BinaryOp::Add, Value::PeriodState(a), Value::PeriodState(b)) => {
                        let Kind::RVCTState(grade) = self.model.kind(expr.ty) else {
                            unreachable!()
                        };
                        let level = grade.level;
                        Ok(Value::PeriodState(rv_period::add(
                            &mut self.model,
                            level,
                            &a,
                            &b,
                        )?))
                    }
                    (BinaryOp::Add, Value::PeriodResidue(a), Value::PeriodResidue(b)) => {
                        let Kind::RVCTResidue(level) = *self.model.kind(expr.ty) else {
                            unreachable!()
                        };
                        Ok(Value::PeriodResidue(level.add(a, b)))
                    }
                    (BinaryOp::Add, Value::RVOre(a), Value::RVOre(b)) => {
                        rv_ore::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::RVQGrP(a), Value::RVQGrP(b)) => {
                        let Kind::RVQGrP(grade) = *self.model.kind(expr.ty) else {
                            unreachable!()
                        };
                        rv_pcurv::add(&mut self.model, grade.prime, &a, &b)
                    }
                    (BinaryOp::Add, Value::RVQGr(a), Value::RVQGr(b)) => {
                        rv_qgr::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::RVSatake(a), Value::RVSatake(b)) => {
                        rv_satake::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::RVSphHecke(a), Value::RVSphHecke(b)) => {
                        rv_sphhecke::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::RVIwHecke(a), Value::RVIwHecke(b)) => {
                        rv_iwhecke::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::RVLaurent(a), Value::RVLaurent(b)) => {
                        rv_laurent::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::EOBruhat, Value::EOType(a), Value::EOType(b)) => Ok(Value::Point(
                        usize::from(eo::bruhat(&mut self.model, &a, &b)?),
                    )),
                    (BinaryOp::EOClosure, Value::EOType(a), Value::EOType(b)) => {
                        let contained = eo::closure_witness(&mut self.model, &a, &b)?.is_some();
                        Ok(Value::Point(usize::from(contained)))
                    }
                    (BinaryOp::AsDeRhamPair, Value::ASDeRham(a), Value::ASDeRham(b)) => {
                        let Kind::ASDeRham(name, field) = self.model.kind(left.ty) else {
                            unreachable!()
                        };
                        let field = *field;
                        let cover = self
                            .model
                            .rv_as_cover_named(name)
                            .expect("registered Artin–Schreier cover");
                        rv_as_derham::pair(&mut self.model, &cover, field, &a, &b)
                    }
                    (BinaryOp::Add, Value::ASDeRham(a), Value::ASDeRham(b)) => {
                        let Kind::ASDeRham(_, field) = *self.model.kind(expr.ty) else {
                            unreachable!()
                        };
                        let sum = binary_linear::add_multiple(&mut self.model, field, &a, &b, 1)?;
                        Ok(Value::ASDeRham(binary_linear::persist(
                            &mut self.model,
                            sum,
                        )?))
                    }
                    (BinaryOp::Add, Value::BinaryScalar(a), Value::BinaryScalar(b)) => {
                        Ok(Value::BinaryScalar(a ^ b))
                    }
                    (BinaryOp::Add, Value::CurveDiff(a), Value::CurveDiff(b)) => {
                        let Kind::CurveDiff(_, _, field) = *self.model.kind(expr.ty) else {
                            unreachable!()
                        };
                        let coordinates =
                            binary_linear::add_multiple(&mut self.model, field, &a, &b, 1)?;
                        Ok(Value::CurveDiff(binary_linear::persist(
                            &mut self.model,
                            coordinates,
                        )?))
                    }
                    (BinaryOp::Add, Value::RVHom(a), Value::RVHom(b)) => {
                        let Kind::RVHom(grade) = self.model.kind(expr.ty) else {
                            unreachable!()
                        };
                        let field = grade.field;
                        rv_hom::add(&mut self.model, field, &a, &b)
                    }
                    (BinaryOp::Add, Value::JetClass(a), Value::JetClass(b)) => {
                        let Kind::JetPage(grade) = *self.model.kind(expr.ty) else {
                            unreachable!()
                        };
                        let field = self.model.jet_path(grade.chain.diagram.path).field;
                        let sum =
                            crate::binary_linear::add_multiple(&mut self.model, field, &a, &b, 1)?;
                        jet_page::persist(&mut self.model, sum)
                    }
                    (BinaryOp::Add, Value::JetHom(a), Value::JetHom(b)) => {
                        jet_hom::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::QKoszul(a), Value::QKoszul(b)) => {
                        Ok(Value::QKoszul(q_koszul::add(&mut self.model, &a, &b)?))
                    }
                    (BinaryOp::Add, Value::MapCobar(a), Value::MapCobar(b)) => {
                        Ok(Value::MapCobar(map_cochain::add(&mut self.model, &a, &b)?))
                    }
                    (BinaryOp::Concat, Value::RhoLoc(a), Value::RhoLoc(b)) => {
                        rho_localized::concat(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Concat, Value::RMCobar(a), Value::RMCobar(b)) => {
                        if let Kind::RMRhoCobar(grade) = *self.model.kind(expr.ty) {
                            return real_rho::concat(&mut self.model, grade, &a, &b);
                        }
                        let terms = real_motivic::concat_support(&mut self.model, &a, &b)?;
                        real_motivic::persist(&mut self.model, terms)
                    }
                    (BinaryOp::Concat, Value::Cobar(a), Value::Cobar(b)) => {
                        cobar::concat(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Concat, Value::PCobar(a), Value::PCobar(b)) => {
                        let (p, _, _, _) = self.model.kind(expr.ty).pcobar_grade().unwrap();
                        odd_cobar::concat(&mut self.model, p, &a, &b)
                    }
                    (BinaryOp::Yoneda, Value::OddClass(a), Value::OddClass(b)) => {
                        odd_ext::yoneda(&mut self.model, left.ty, &a, right.ty, &b, expr.ty)
                    }
                    (BinaryOp::Yoneda, Value::Class(a), Value::Class(b)) => {
                        if let Kind::RMRhoExt(grade) = *self.model.kind(expr.ty) {
                            real_rho::yoneda(
                                &mut self.model,
                                (left.ty, &a),
                                (right.ty, &b),
                                expr.ty,
                                grade,
                            )
                        } else if matches!(self.model.kind(expr.ty), Kind::RMExt(..)) {
                            real_ext::yoneda(
                                &mut self.model,
                                (left.ty, &a),
                                (right.ty, &b),
                                expr.ty,
                            )
                        } else {
                            adams::yoneda(&mut self.model, left.ty, &a, right.ty, &b, expr.ty)
                        }
                    }
                    (BinaryOp::XAction, Value::Class(a), Value::Class(b)) => {
                        comodule::action(&mut self.model, left.ty, &a, right.ty, &b, expr.ty)
                    }
                    (BinaryOp::XMap, Value::Map(f), Value::Class(b)) => {
                        comodule::map(&mut self.model, left.ty, &f, right.ty, &b, expr.ty)
                    }
                    (BinaryOp::MapExtCompose, Value::RVLaurent(g), Value::RVLaurent(f)) => {
                        rv_laurent::compose(&mut self.model, left.ty, &g, right.ty, &f)
                    }
                    (BinaryOp::MapExtCompose, Value::RVHom(g), Value::RVHom(f)) => {
                        rv_hom::compose(&mut self.model, left.ty, &g, right.ty, &f)
                    }
                    (BinaryOp::MapExtCompose, Value::JetClass(g), Value::JetClass(f)) => {
                        jet_page::compose(&mut self.model, (left.ty, &g), (right.ty, &f), expr.ty)
                    }
                    (BinaryOp::MapExtCompose, Value::JetHom(g), Value::JetHom(f)) => {
                        let Kind::JetHom(grade) = *self.model.kind(expr.ty) else {
                            unreachable!()
                        };
                        let path = self.model.jet_path(grade.diagram.path);
                        let terms = jet_hom::compose(&mut self.model, &path, &g, &f)?;
                        jet_hom::persist(&mut self.model, terms)
                    }
                    (BinaryOp::MapExtCompose, Value::QKoszul(g), Value::QKoszul(f)) => {
                        Ok(Value::QKoszul(q_koszul::compose(&mut self.model, &g, &f)?))
                    }
                    (BinaryOp::MapExtCompose, Value::Class(g), Value::Class(f)) => {
                        if matches!(self.model.kind(expr.ty), Kind::VBock(..)) {
                            return vbock::compose(
                                &mut self.model,
                                (left.ty, &g),
                                (right.ty, &f),
                                expr.ty,
                            );
                        }
                        if matches!(self.model.kind(expr.ty), Kind::QExt(..) | Kind::VExt(..)) {
                            return q_koszul::compose_classes(
                                &mut self.model,
                                (left.ty, &g),
                                (right.ty, &f),
                                expr.ty,
                            );
                        }
                        let compose = if matches!(self.model.kind(expr.ty), Kind::MapMay(..)) {
                            map_may::compose
                        } else {
                            map_ext::compose
                        };
                        compose(&mut self.model, (left.ty, &g), (right.ty, &f), expr.ty)
                    }
                    (BinaryOp::MapExtCompose, Value::MapCobar(g), Value::MapCobar(f)) => Ok(
                        Value::MapCobar(map_cochain::compose(&mut self.model, &g, &f)?),
                    ),
                    (BinaryOp::RhoMul, Value::Class(a), Value::Class(b)) => {
                        rho_bock::product(&mut self.model, (left.ty, &a), (right.ty, &b), expr.ty)
                    }
                    (BinaryOp::MayMul, Value::Class(a), Value::Class(b)) => {
                        may::multiply(&mut self.model, left.ty, &a, right.ty, &b, expr.ty)
                    }
                    (
                        BinaryOp::MayMul | BinaryOp::CeMul,
                        Value::OddClass(a),
                        Value::OddClass(b),
                    ) => odd_may::multiply(&mut self.model, left.ty, &a, right.ty, &b, expr.ty),
                    (BinaryOp::CeMul, Value::Class(a), Value::Class(b)) => {
                        ce::multiply(&mut self.model, left.ty, &a, right.ty, &b, expr.ty)
                    }
                    (BinaryOp::Add, Value::Bar(a), Value::Bar(b)) => {
                        cyclic::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::CyclicClass(a), Value::CyclicClass(b)) => {
                        Ok(Value::CyclicClass(a ^ b))
                    }
                    (BinaryOp::Mul, Value::Operation(a), Value::Operation(b)) => {
                        steenrod::multiply(&mut self.model, &a, &b)
                    }
                    (BinaryOp::MilnorPair, Value::Operation(a), Value::Cobar(b)) => {
                        let paired = milnor::pair(&mut self.model, &a, &b)?;
                        Ok(Value::Point(usize::from(paired)))
                    }
                    (BinaryOp::Act, Value::Operation(operation), Value::Class(class)) => {
                        let (base, degree) = match *self.model.kind(right.ty) {
                            Kind::Cohomology(base, degree) => (base, degree),
                            Kind::Stable(base, shift, degree) => {
                                let Some(reduced_degree) = degree.checked_sub(shift) else {
                                    return Ok(cohomology::zero());
                                };
                                (base, reduced_degree)
                            }
                            _ => unreachable!(),
                        };
                        cohomology::act(&mut self.model, base, degree, &operation, class)
                    }
                    (BinaryOp::Cup, Value::Class(a), Value::Class(b)) => {
                        cohomology::cup(&mut self.model, left.ty, &a, right.ty, &b, expr.ty)
                    }
                    (BinaryOp::Pull, Value::Map(f), Value::Class(b)) => {
                        cohomology::pull(&mut self.model, left.ty, &f, right.ty, &b, expr.ty)
                    }
                    _ => unreachable!(),
                }
            }
            TypedForm::Coact(n, index, input) => {
                let Kind::Cohomology(base, q) = *self.model.kind(input.ty) else {
                    unreachable!()
                };
                let Value::Class(class) = self.eval(input)? else {
                    unreachable!()
                };
                let support = milnor::coact(&mut self.model, base, *n, *index, q, class)?;
                Ok(Value::Cobar(support.into()))
            }
            TypedForm::RhoLocOperation(operation, input) => match (operation, self.eval(input)?) {
                (RhoLocOp::Localize, Value::RMCobar(terms)) => {
                    rho_localized::localize(&mut self.model, &terms)
                }
                (RhoLocOp::Section(_), Value::Cobar(words)) => {
                    rho_localized::section(&mut self.model, &words)
                }
                (RhoLocOp::Project, Value::RhoLoc(terms)) => {
                    rho_localized::project(&mut self.model, &terms)
                }
                (RhoLocOp::Project, Value::Class(coords)) => {
                    real_ext::project_classical(&mut self.model, input.ty, &coords, expr.ty)
                }
                (RhoLocOp::Homotopy, Value::RhoLoc(terms)) => {
                    rho_localized::homotopy(&mut self.model, &terms)
                }
                _ => unreachable!(),
            },
            TypedForm::RhoOperation(operation, input) => match (operation, self.eval(input)?) {
                (RhoOp::Class(..), Value::RMCobar(terms)) => {
                    rho_bock::class(&mut self.model, expr.ty, &terms)
                }
                (RhoOp::Representative, Value::Class(coords)) => {
                    rho_bock::representative(&mut self.model, input.ty, &coords)
                }
                (RhoOp::Differential, Value::Class(coords)) => {
                    rho_bock::differential(&mut self.model, input.ty, expr.ty, &coords)
                }
                (RhoOp::Next, Value::Class(coords)) => {
                    page_transition::next(&mut self.model, input.ty, &coords)
                }
                (RhoOp::Previous, Value::Class(coords)) => {
                    rho_bock::previous(&mut self.model, input.ty, expr.ty, &coords)
                }
                _ => unreachable!(),
            },
            TypedForm::RhoShift(power, input) => match self.eval(input)? {
                Value::RMCobar(terms) => {
                    if let Kind::RMRhoCobar(grade) = *self.model.kind(input.ty) {
                        if *power == 0 {
                            return Ok(Value::RMCobar(terms));
                        }
                        let terms =
                            real_rho::shift(&mut self.model, &terms, grade.exponent, *power)?;
                        return real_motivic::persist(&mut self.model, terms);
                    }
                    real_motivic::rho_shift(&mut self.model, &terms, *power)
                }
                Value::RhoMin(terms) => rho_minimal::shift(&mut self.model, &terms, *power),
                Value::Class(coords) if *power == 0 => Ok(Value::Class(coords)),
                Value::Class(coords) => {
                    if let Kind::RMRhoExt(grade) = *self.model.kind(input.ty) {
                        return real_rho::shift_class(
                            &mut self.model,
                            input.ty,
                            &coords,
                            expr.ty,
                            grade.exponent,
                            *power,
                        );
                    }
                    real_ext::shift(&mut self.model, input.ty, &coords, expr.ty, *power)
                }
                _ => unreachable!(),
            },
            TypedForm::RhoMinOperation(op, grade, input) => match (op, self.eval(input)?) {
                (RhoMinOp::Inject(f), Value::Class(coords)) => {
                    rho_minimal::inject(&mut self.model, *f, &coords)
                }
                (RhoMinOp::Coefficient(f), Value::RhoMin(terms)) => {
                    rho_minimal::coefficient(&mut self.model, &terms, *f)
                }
                (RhoMinOp::Differential, Value::RhoMin(terms)) => {
                    rho_minimal::differential(&mut self.model, *grade, &terms)
                }
                (RhoMinOp::Lift, Value::RhoMin(terms)) => {
                    rho_minimal::lift(&mut self.model, *grade, &terms)
                }
                (RhoMinOp::Project, Value::RMCobar(terms)) => {
                    rho_minimal::project(&mut self.model, *grade, &terms)
                }
                (RhoMinOp::Homotopy, Value::RMCobar(terms)) => {
                    rho_minimal::homotopy(&mut self.model, *grade, &terms)
                }
                _ => unreachable!(),
            },
            TypedForm::RealComplex(input) => match self.eval(input)? {
                Value::RMCobar(terms) => real_motivic::specialize(&mut self.model, &terms),
                Value::Class(coords) => {
                    real_ext::specialize(&mut self.model, input.ty, &coords, expr.ty)
                }
                _ => unreachable!(),
            },
            TypedForm::RealRhoOperation(op, grade, input) => {
                let value = self.eval(input)?;
                if matches!(op, RealRhoOp::Reduce(_)) && expr.ty == input.ty {
                    return Ok(value);
                }
                match (op, value) {
                    (RealRhoOp::Reduce(e), Value::RMCobar(terms)) => {
                        real_rho::reduce(&mut self.model, &terms, *e)
                    }
                    (RealRhoOp::Reduce(_), Value::Cobar(words)) => {
                        let terms = real_rho::from_complex(&mut self.model, *grade, &words)?;
                        real_motivic::persist(&mut self.model, terms)
                    }
                    (RealRhoOp::Reduce(_), Value::Class(coords)) if coords.is_empty() => {
                        Ok(Value::Class(coords))
                    }
                    (RealRhoOp::Reduce(e), Value::Class(coords)) => {
                        if matches!(self.model.kind(input.ty), Kind::MExt(..)) {
                            let words = self
                                .model
                                .adams(input.ty)?
                                .representative(&mut self.model, &coords)?;
                            let terms = real_rho::from_complex(&mut self.model, *grade, &words)?;
                            real_ext::class(&mut self.model, expr.ty, &terms)
                        } else {
                            real_rho::reduce_class(&mut self.model, input.ty, &coords, expr.ty, *e)
                        }
                    }
                    (RealRhoOp::Lift, Value::RMCobar(terms)) => Ok(Value::RMCobar(terms)),
                    (RealRhoOp::Lift, Value::Class(coords)) => {
                        real_ext::representative(&mut self.model, input.ty, &coords)
                    }
                    (
                        RealRhoOp::Connect
                        | RealRhoOp::ConnectRep
                        | RealRhoOp::Beta
                        | RealRhoOp::BetaRep
                        | RealRhoOp::BetaNull,
                        Value::Class(coords),
                    ) => real_rho::connecting(&mut self.model, *op, input.ty, &coords, expr.ty),
                    _ => unreachable!(),
                }
            }
            TypedForm::RealBar(op, s, k, index) => {
                let group = self.model.real_bars(*s, *k)?;
                let bar = group.get(*index)?;
                match op {
                    RealBarOp::Generator => bar.generator(&mut self.model, expr.ty),
                    RealBarOp::Representative => Ok(bar.representative()),
                    RealBarOp::Witness => bar.witness(&mut self.model),
                }
            }
            TypedForm::RealExtension(op, grade, input) => {
                let Value::Class(coords) = self.eval(input)? else {
                    unreachable!()
                };
                match op {
                    RealExtOp::Coordinates => {
                        real_extensions::decompose(&mut self.model, (input.ty, &coords))
                    }
                    RealExtOp::Class => {
                        real_extensions::reconstruct(&mut self.model, *grade, &coords, expr.ty)
                    }
                    RealExtOp::Lift => real_extensions::lift(&mut self.model, *grade, &coords),
                    RealExtOp::Divide(power) => real_extensions::divide(
                        &mut self.model,
                        (input.ty, &coords),
                        *power,
                        expr.ty,
                    ),
                    RealExtOp::Edge(_, f) => {
                        real_extensions::edge(&mut self.model, (input.ty, &coords), expr.ty, *f)
                    }
                }
            }
            TypedForm::RhoProductNull(power, inputs) => {
                let [a, b] = inputs.as_ref();
                let (Value::Class(x), Value::Class(y)) = (self.eval(a)?, self.eval(b)?) else {
                    unreachable!()
                };
                real_extensions::product_null(&mut self.model, (a.ty, &x), (b.ty, &y), *power)
            }
            TypedForm::RhoNull(power, input) => {
                let grade = self.model.kind(input.ty).real_grade().unwrap();
                let (terms, known_cycle) = self.real_terms(input)?;
                real_ext::null(&mut self.model, grade, &terms, *power, known_cycle)
            }
            TypedForm::RealMassey(coset, layout, inputs) => {
                let [a, b, c] = inputs.as_ref();
                let values = [
                    self.real_terms(a)?,
                    self.real_terms(b)?,
                    self.real_terms(c)?,
                ];
                let inputs = values
                    .each_ref()
                    .map(|(terms, known)| (terms.as_ref(), *known));
                if *coset {
                    let coset = layout.coset(&mut self.model, inputs)?;
                    Ok(Value::RMTriple(Rc::new(coset)))
                } else {
                    layout.representative(&mut self.model, inputs)
                }
            }
            TypedForm::RealMasseyNull(grade, inputs) => {
                let [left, right] = inputs.as_ref();
                let (a, a_closed) = self.real_terms(left)?;
                let (b, b_closed) = self.real_terms(right)?;
                real_massey::adjacent_null(&mut self.model, *grade, (&a, a_closed), (&b, b_closed))
            }
            TypedForm::RealTripleClass(index, input) => {
                let Value::RMTriple(coset) = self.eval(input)? else {
                    unreachable!();
                };
                match index {
                    Some(index) => coset.direction(&mut self.model, *index),
                    None => Ok(coset.base()),
                }
            }
            TypedForm::MotivicCast(input) => {
                let value = self.eval(input)?;
                let Value::Class(support) = value else {
                    return Ok(value);
                };
                if support.is_empty() || input.ty == expr.ty {
                    return Ok(Value::Class(support));
                }
                let group = self.model.adams(input.ty)?;
                let representative = group.representative(&mut self.model, &support)?;
                adams::class(&mut self.model, expr.ty, &representative)
            }
            TypedForm::Cobar(op, input) => match (op, self.eval(input)?) {
                (CobarOp::Differential, Value::RhoLoc(terms)) => {
                    rho_localized::differential(&mut self.model, &terms)
                }
                (CobarOp::Differential, Value::RMCobar(terms)) => {
                    if let Kind::RMRhoCobar(grade) = *self.model.kind(input.ty) {
                        return real_rho::differential(&mut self.model, grade, &terms);
                    }
                    let terms = real_motivic::differential_support(&mut self.model, &terms)?;
                    real_motivic::persist(&mut self.model, terms)
                }
                (CobarOp::Differential, Value::Cobar(words)) => {
                    let kind = self.model.kind(input.ty);
                    let profile = if matches!(kind, Kind::MCobar(..)) {
                        Profile::Full
                    } else {
                        kind.cobar_grade().unwrap().0
                    };
                    cobar::differential(&mut self.model, profile, &words)
                }
                (CobarOp::Differential, Value::PCobar(terms)) => {
                    let (p, profile, _, _) = self.model.kind(input.ty).pcobar_grade().unwrap();
                    odd_cobar::differential(&mut self.model, p, profile, &terms)
                }
                (CobarOp::Class, Value::RMCobar(terms)) => {
                    real_ext::class(&mut self.model, expr.ty, &terms)
                }
                (CobarOp::Class, Value::PCobar(terms)) => {
                    odd_ext::class(&mut self.model, expr.ty, &terms)
                }
                (CobarOp::Class, Value::Cobar(words)) => {
                    adams::class(&mut self.model, expr.ty, &words)
                }
                (CobarOp::QDifferential, Value::QKoszul(terms)) => {
                    let (source, target, selected) = match self.model.kind(input.ty) {
                        Kind::QKoszul(grade) => {
                            (grade.source, grade.target, grade.selected.clone())
                        }
                        Kind::VKoszul(grade) => {
                            (grade.source, grade.target, grade.selected.clone())
                        }
                        _ => unreachable!(),
                    };
                    let terms =
                        q_koszul::differential(&mut self.model, source, target, &selected, &terms)?;
                    Ok(Value::QKoszul(terms))
                }
                (CobarOp::RvCyclicD, Value::RVNCChain(terms)) => {
                    let Kind::RvCyclicChain(grade) = self.model.kind(input.ty).clone() else {
                        unreachable!()
                    };
                    Ok(Value::RVNCChain(rv_cyclic_homology::differential(
                        &mut self.model,
                        &grade,
                        &terms,
                    )?))
                }
                (CobarOp::RvCyclicClass, Value::RVNCChain(terms)) => {
                    rv_cyclic_homology::class(&mut self.model, expr.ty, &terms)
                }
                (CobarOp::RvCyclicRep, Value::Class(class)) => Ok(Value::RVNCChain(
                    rv_cyclic_homology::representative(&mut self.model, input.ty, &class)?,
                )),
                (CobarOp::RvCyclicI, Value::RVHHChain(words)) => Ok(Value::RVNCChain(
                    rv_negative::term(&mut self.model, 0, &words)?,
                )),
                (CobarOp::RvCyclicI, Value::Class(class)) => {
                    rv_cyclic_homology::include_class(&mut self.model, input.ty, expr.ty, &class)
                }
                (CobarOp::RvCyclicB, Value::RVNCChain(terms)) => {
                    let Kind::RvCyclicChain(grade) = self.model.kind(input.ty).clone() else {
                        unreachable!()
                    };
                    Ok(Value::RVHHChain(rv_cyclic_homology::connes(
                        &mut self.model,
                        &grade,
                        &terms,
                    )?))
                }
                (CobarOp::RvCyclicB, Value::Class(class)) => {
                    rv_cyclic_homology::connes_class(&mut self.model, input.ty, expr.ty, &class)
                }
                (CobarOp::RvNCd, Value::RVNCChain(terms)) => {
                    let Kind::RVNCChain(grade) = self.model.kind(input.ty).clone() else {
                        unreachable!()
                    };
                    Ok(Value::RVNCChain(rv_negative::differential(
                        &mut self.model,
                        &grade,
                        &terms,
                    )?))
                }
                (CobarOp::RvNCClass, Value::RVNCChain(terms)) => {
                    rv_negative::class(&mut self.model, expr.ty, &terms)
                }
                (CobarOp::RvNCRep, Value::Class(class)) => Ok(Value::RVNCChain(
                    rv_negative::representative(&mut self.model, input.ty, &class)?,
                )),
                (CobarOp::MixDown | CobarOp::MixUp | CobarOp::MixD, Value::OddClass(terms)) => {
                    match *self.model.kind(input.ty) {
                        Kind::MixedChain(grade, _) => {
                            let complex = self.model.mixed_complex(grade.complex);
                            let differential = match op {
                                CobarOp::MixDown => mixed::Differential::Down,
                                CobarOp::MixUp => mixed::Differential::Up,
                                _ => mixed::Differential::Total,
                            };
                            let image = complex.apply(&mut self.model, differential, &terms)?;
                            odd_ext::persist(&mut self.model, image)
                        }
                        Kind::MixedMinimum(grade, true) => {
                            let data = self.model.mixed_minimal(grade.minimal);
                            let image = data.differential(&mut self.model, &terms)?;
                            odd_ext::persist(&mut self.model, image)
                        }
                        Kind::MixedPage(_) => {
                            mixed_page::differential(&mut self.model, input.ty, expr.ty, &terms)
                        }
                        _ => unreachable!(),
                    }
                }
                (
                    CobarOp::MixHClass | CobarOp::MixHPrimitive | CobarOp::MixHRep,
                    Value::OddClass(terms),
                ) => self.eval_mixed_h_unary(op, input.ty, terms),
                (CobarOp::MixMinI, Value::OddClass(terms)) => {
                    let Kind::MixedMinimum(grade, total) = *self.model.kind(input.ty) else {
                        unreachable!()
                    };
                    let data = self.model.mixed_minimal(grade.minimal);
                    let image = data.inclusion(&mut self.model, &terms, total)?;
                    odd_ext::persist(&mut self.model, image)
                }
                (CobarOp::RvMixedToHH, Value::OddClass(terms)) => {
                    let (id, is_chain) = match *self.model.kind(input.ty) {
                        Kind::MixedChain(grade, false) => (grade.complex, true),
                        Kind::MixedPage(grade) => (grade.chain.complex, false),
                        _ => unreachable!(),
                    };
                    let origin = self.model.mixed_hh_origin(id)?;
                    if is_chain {
                        origin.chain_to_hh(&mut self.model, &terms)
                    } else {
                        origin.class_to_hh(&mut self.model, input.ty, expr.ty, &terms)
                    }
                }
                (CobarOp::MixClass, Value::OddClass(terms)) => {
                    mixed_page::class(&mut self.model, expr.ty, &terms)
                }
                (CobarOp::MixRep, Value::OddClass(terms)) => {
                    mixed_page::representative(&mut self.model, input.ty, &terms)
                }
                (CobarOp::MixNext, Value::OddClass(terms)) => {
                    mixed_page::next(&mut self.model, input.ty, expr.ty, &terms)
                }
                (CobarOp::MixPrev, Value::OddClass(terms)) => {
                    mixed_page::previous(&mut self.model, input.ty, expr.ty, &terms)
                }
                (CobarOp::RvHHClass, Value::RVHHChain(terms)) => {
                    rv_hh_homology::class(&mut self.model, expr.ty, &terms)
                }
                (CobarOp::RvHHRep, Value::Class(class)) => Ok(Value::RVHHChain(
                    rv_hh_homology::representative(&mut self.model, input.ty, &class)?,
                )),
                (CobarOp::RvConnesB, Value::Class(class)) => {
                    rv_hh_homology::connes(&mut self.model, input.ty, &class, expr.ty)
                }
                (CobarOp::RvHHd | CobarOp::RvConnesB, Value::RVHHChain(terms)) => {
                    let Kind::RVHHChain(grade) = self.model.kind(input.ty) else {
                        unreachable!()
                    };
                    let grade = grade.clone();
                    let words = if matches!(op, CobarOp::RvHHd) {
                        rv_hh::differential(&mut self.model, &grade, &terms)?
                    } else {
                        rv_hh::connes(&mut self.model, &grade, &terms)?
                    };
                    Ok(Value::RVHHChain(rv_hh::persist(&mut self.model, words)?))
                }
                (CobarOp::RvBCDesuspend, Value::RVBar(terms)) => Ok(Value::RVBarCobar(
                    rv_barcobar::desuspend(&mut self.model, &terms)?,
                )),
                (CobarOp::RvBCD, Value::RVBarCobar(terms)) => {
                    let Kind::RVBarCobar(grade) = *self.model.kind(input.ty) else {
                        unreachable!()
                    };
                    Ok(Value::RVBarCobar(rv_barcobar::differential(
                        &mut self.model,
                        grade,
                        &terms,
                    )?))
                }
                (CobarOp::RvBCProject | CobarOp::RvBCHomotopy, Value::RVBarCobar(terms)) => {
                    let Kind::RVBarCobar(grade) = *self.model.kind(input.ty) else {
                        unreachable!()
                    };
                    let operation = if matches!(op, CobarOp::RvBCProject) {
                        rv_barcobar::Contraction::Projection
                    } else {
                        rv_barcobar::Contraction::Homotopy
                    };
                    rv_barcobar::contract(&mut self.model, grade, &terms, operation)
                }
                (CobarOp::RvBarDualize, Value::RVBar(terms)) => Ok(Value::RVBar(terms)),
                (CobarOp::RvBarDualD, Value::RVBar(terms)) => {
                    let Kind::RVBarDual(grade) = *self.model.kind(input.ty) else {
                        unreachable!()
                    };
                    Ok(Value::RVBar(rv_bar_dual::differential(
                        &mut self.model,
                        grade,
                        &terms,
                    )?))
                }
                (CobarOp::RvPathClass, Value::RVBar(terms)) => {
                    rv_bar_dual::class(&mut self.model, expr.ty, &terms)
                }
                (CobarOp::RvPathRep, Value::Class(class)) => Ok(Value::RVBar(
                    rv_bar_dual::representative(&mut self.model, input.ty, &class)?,
                )),
                (CobarOp::RvBarD, Value::Class(class)) => {
                    rv_bar_page::differential(&mut self.model, input.ty, expr.ty, &class)
                }
                (CobarOp::RvBarClass, Value::RVBar(terms)) => {
                    rv_bar_page::class(&mut self.model, expr.ty, &terms)
                }
                (CobarOp::RvBarRep, Value::Class(class)) => Ok(Value::RVBar(
                    rv_bar_page::representative(&mut self.model, input.ty, &class)?,
                )),
                (CobarOp::RvBarNext, Value::Class(class)) => {
                    page_transition::next(&mut self.model, input.ty, &class)
                }
                (CobarOp::RvBarPrev, Value::Class(class)) => {
                    rv_bar_page::previous(&mut self.model, input.ty, expr.ty, &class)
                }
                (CobarOp::RvBarD, Value::RVBar(terms)) => {
                    let Kind::RVBar(grade) = *self.model.kind(input.ty) else {
                        unreachable!();
                    };
                    let terms = rv_bar::differential(&mut self.model, grade, &terms)?;
                    Ok(Value::RVBar(rv_bar::persist(&mut self.model, terms)?))
                }
                (CobarOp::RvMForget, value) => self.rv_koszul_regrade(input.ty, expr.ty, value),
                (CobarOp::RvKoszulDifferential, Value::RVKoszul(terms)) => {
                    let Kind::RVKoszul(grade) = self.model.kind(input.ty) else {
                        unreachable!();
                    };
                    let ideal = self.model.rv_ideal_named(&grade.ideal).unwrap();
                    Ok(Value::RVKoszul(rv_koszul::differential(
                        &mut self.model,
                        &ideal,
                        &terms,
                    )?))
                }
                (CobarOp::RvKoszulPrimitive, Value::RVKoszul(terms)) => {
                    if terms.is_empty() {
                        return Ok(Value::RVKoszul([].into()));
                    }
                    let Kind::RVKoszul(grade) = self.model.kind(input.ty) else {
                        unreachable!();
                    };
                    let grade = grade.clone();
                    let ideal = self.model.rv_ideal_named(&grade.ideal).unwrap();
                    let group = self.model.rv_complex_type(grade.clone(), true)?;
                    let primitive = self
                        .model
                        .rv_tor(group)?
                        .primitive(&mut self.model, &ideal, &grade, &terms)?
                        .ok_or_else(|| Error::plain("rvkPrimitive requires a Koszul boundary"))?;
                    Ok(Value::RVKoszul(rv_koszul::persist(
                        &mut self.model,
                        primitive,
                    )?))
                }
                (CobarOp::RvTorClass, Value::RVKoszul(terms)) => {
                    let Kind::RVTor(grade) = self.model.kind(expr.ty) else {
                        unreachable!();
                    };
                    let ideal = self.model.rv_ideal_named(&grade.ideal).unwrap();
                    self.rv_tor_class(expr.ty, &ideal, &terms)
                }
                (CobarOp::RvTorRepresentative, Value::Class(class)) => {
                    if class.is_empty() {
                        return Ok(Value::RVKoszul([].into()));
                    }
                    let tor = self.model.rv_tor(input.ty)?;
                    Ok(Value::RVKoszul(
                        tor.representative(&mut self.model, &class)?,
                    ))
                }
                (CobarOp::QClass, Value::QKoszul(terms)) => {
                    q_koszul::class(&mut self.model, expr.ty, &terms)
                }
                (CobarOp::QRepresentative, Value::Class(class)) => {
                    let terms = q_koszul::representative(&mut self.model, input.ty, &class)?;
                    Ok(Value::QKoszul(terms))
                }
                (CobarOp::VLocal, value) => {
                    v_koszul::localize(&mut self.model, input.ty, expr.ty, value)
                }
                (CobarOp::VBockClass, Value::QKoszul(terms)) => {
                    vbock::class(&mut self.model, expr.ty, &terms)
                }
                (CobarOp::VBockDifferential, Value::Class(class)) => {
                    vbock::differential(&mut self.model, input.ty, expr.ty, &class)
                }
                (CobarOp::VBockRepresentative, Value::Class(class)) => {
                    let terms = vbock::representative(&mut self.model, input.ty, &class)?;
                    Ok(Value::QKoszul(terms))
                }
                (CobarOp::MapDifferential, Value::MapCobar(terms)) => {
                    let Kind::MapCobar(source, target, ..) = *self.model.kind(input.ty) else {
                        unreachable!()
                    };
                    let terms = map_cochain::differential(&mut self.model, source, target, &terms)?;
                    Ok(Value::MapCobar(terms))
                }
                (CobarOp::MapClass, Value::MapCobar(terms)) => {
                    map_cochain::class(&mut self.model, expr.ty, &terms)
                }
                (CobarOp::MapRepresentative, Value::Class(class)) => {
                    let terms = map_cochain::representative(&mut self.model, input.ty, &class)?;
                    Ok(Value::MapCobar(terms))
                }
                (CobarOp::MayClass, Value::Cobar(words)) => {
                    may::class(&mut self.model, expr.ty, &words)
                }
                (CobarOp::MayClass, Value::MapCobar(terms)) => {
                    map_may::class(&mut self.model, expr.ty, &terms)
                }
                (CobarOp::MayClass | CobarOp::CeClass, Value::PCobar(terms)) => {
                    odd_may::class(&mut self.model, expr.ty, &terms)
                }
                (CobarOp::MayDifferential, Value::Class(class)) => {
                    if matches!(self.model.kind(input.ty), Kind::MapMay(..)) {
                        map_may::differential(&mut self.model, input.ty, &class, expr.ty)
                    } else {
                        may::differential(&mut self.model, input.ty, &class, expr.ty)
                    }
                }
                (CobarOp::CeClass, Value::Cobar(words)) => {
                    ce::class(&mut self.model, expr.ty, &words)
                }
                (CobarOp::CeDifferential, Value::Class(class)) => {
                    ce::differential(&mut self.model, input.ty, &class, expr.ty)
                }
                (CobarOp::CeRepresentative, Value::Class(class)) => {
                    if class.is_empty() {
                        return Ok(cobar::zero());
                    }
                    let page = self.model.ce_page(input.ty)?;
                    let words = page.representative(&mut self.model, &class)?;
                    self.model.allocate(words.len())?;
                    Ok(Value::Cobar(words.into()))
                }
                (CobarOp::PageTransition(_, PageDirection::Next), Value::Class(class)) => {
                    page_transition::next(&mut self.model, input.ty, &class)
                }
                (CobarOp::PageTransition(_, PageDirection::Previous), Value::Class(class)) => {
                    page_transition::previous(&mut self.model, input.ty, expr.ty, &class)
                }
                (CobarOp::MayDifferential | CobarOp::CeDifferential, Value::OddClass(class)) => {
                    odd_may::differential(&mut self.model, input.ty, &class, expr.ty)
                }
                (CobarOp::PageTransition(_, PageDirection::Next), Value::OddClass(class)) => {
                    odd_may::next(&mut self.model, input.ty, &class, expr.ty)
                }
                (CobarOp::PageTransition(_, PageDirection::Previous), Value::OddClass(class)) => {
                    odd_may::previous(&mut self.model, input.ty, &class, expr.ty)
                }
                (
                    CobarOp::MayRepresentative | CobarOp::CeRepresentative,
                    Value::OddClass(class),
                ) => {
                    if class.is_empty() {
                        return Ok(odd_cobar::zero());
                    }
                    let page = self.model.odd_page(input.ty)?;
                    let terms = page.representative(&mut self.model, &class)?;
                    odd_cobar::persist(&mut self.model, terms)
                }
                (CobarOp::MayRepresentative, Value::Class(class)) => {
                    if matches!(self.model.kind(input.ty), Kind::MapMay(..)) {
                        let terms = map_may::representative(&mut self.model, input.ty, &class)?;
                        return Ok(Value::MapCobar(terms));
                    }
                    if class.is_empty() {
                        return Ok(cobar::zero());
                    }
                    let page = self.model.may_page(input.ty)?;
                    let words = page.representative(&mut self.model, &class)?;
                    self.model.allocate(words.len())?;
                    Ok(Value::Cobar(words.into()))
                }
                (CobarOp::Representative, Value::OddClass(class)) => {
                    if class.is_empty() {
                        return Ok(odd_cobar::zero());
                    }
                    let group = self.model.odd_ext(input.ty)?;
                    let terms = group.representative(&mut self.model, &class)?;
                    odd_cobar::persist(&mut self.model, terms)
                }
                (CobarOp::Representative, Value::Class(class)) => {
                    if matches!(
                        self.model.kind(input.ty),
                        Kind::RMExt(..) | Kind::RMRhoExt(..)
                    ) {
                        return real_ext::representative(&mut self.model, input.ty, &class);
                    }
                    if class.is_empty() {
                        return Ok(cobar::zero());
                    }
                    let adams = self.model.adams(input.ty)?;
                    let words = adams.representative(&mut self.model, &class)?;
                    self.model.allocate(words.len())?;
                    Ok(Value::Cobar(words.into()))
                }
                (CobarOp::Primitive, Value::RMCobar(terms)) => {
                    if let Kind::RMRhoCobar(grade) = *self.model.kind(input.ty) {
                        real_rho::primitive(&mut self.model, grade, &terms, "extPrimitive")
                    } else {
                        let grade = self.model.kind(input.ty).real_grade().unwrap();
                        real_ext::primitive(&mut self.model, grade, &terms)
                    }
                }
                (CobarOp::Primitive, Value::PCobar(terms)) => {
                    if terms.is_empty() {
                        return Ok(odd_cobar::zero());
                    }
                    let (p, profile, s, t) = self.model.kind(input.ty).pcobar_grade().unwrap();
                    let ty = self.model.profile_odd_ext_type(p, profile, s, t)?;
                    let group = self.model.odd_ext(ty)?;
                    let primitive = group.primitive(&mut self.model, &terms)?.ok_or_else(|| {
                        Error::plain("extPrimitive requires an odd-prime cobar boundary")
                    })?;
                    odd_cobar::persist(&mut self.model, primitive)
                }
                (CobarOp::Primitive, Value::Cobar(words)) => {
                    if words.is_empty() {
                        return Ok(cobar::zero());
                    }
                    let ty = if let Kind::MCobar(grade) = *self.model.kind(input.ty) {
                        self.model.mext_type(grade)?
                    } else {
                        let (profile, s, t) = self.model.kind(input.ty).cobar_grade().unwrap();
                        self.model.profile_ext_type(profile, s, t)?
                    };
                    let adams = self.model.adams(ty)?;
                    let primitive = adams
                        .primitive(&mut self.model, &words)?
                        .ok_or_else(|| Error::plain("extPrimitive requires a cobar boundary"))?;
                    self.model.allocate(primitive.len())?;
                    Ok(Value::Cobar(primitive.into()))
                }
                _ => unreachable!(),
            },
            TypedForm::ResSpecializeQuad(grade, field, weights, rows) => {
                let rows = self.quadratic_rows(rows)?;
                rv_res_space::quadratic(&mut self.model, *grade, weights, *field, rows)?
                    .persist(&mut self.model)
            }
            TypedForm::ResOrthogonal(input) => {
                let Value::RVResSpace(space) = self.eval(input)? else {
                    unreachable!()
                };
                space.orthogonal(&mut self.model)?.persist(&mut self.model)
            }
            TypedForm::ValMatQuad(grade, field, rows) => {
                let rows = self.quadratic_rows(rows)?;
                Ok(Value::RVValMat(tropical_matroid::quadratic(
                    &mut self.model,
                    *grade,
                    *field,
                    &rows,
                )?))
            }
            TypedForm::StableLineQuad(field, marks) => {
                let mut points = Vec::with_capacity(marks.len());
                for [first, second] in marks {
                    let Value::RVQpQuad(a) = self.eval(first)? else {
                        unreachable!()
                    };
                    let Value::RVQpQuad(b) = self.eval(second)? else {
                        unreachable!()
                    };
                    points.push([a, b]);
                }
                Ok(Value::RVStableMarks(StableCurve::quadratic(
                    &mut self.model,
                    *field,
                    &points,
                )?))
            }
            TypedForm::StableProjection(input) => {
                let tree = match self.eval(input)? {
                    Value::RVValMat(matroid) => {
                        StableTree::from_matroid(&mut self.model, &matroid)?
                    }
                    Value::RVStableMarks(curve) => curve.tree.clone(),
                    _ => unreachable!(),
                };
                Ok(Value::RVStableTree(tree))
            }
            TypedForm::ValDual(input) => {
                let Value::RVValMat(value) = self.eval(input)? else {
                    unreachable!()
                };
                Ok(Value::RVValMat(value.dual(&mut self.model)?))
            }
            TypedForm::QpQuadUnary(field, operation, input) => {
                let Value::RVQpQuad(value) = self.eval(input)? else {
                    unreachable!()
                };
                let value = operation.apply(&mut self.model, *field, value)?;
                rv_qp_quad::persist(&mut self.model, value)
            }
            TypedForm::HeckeNegate(input) => match self.eval(input)? {
                Value::RVIwHecke(value) => rv_iwhecke::neg(&mut self.model, &value),
                Value::RVSphHecke(value) => rv_sphhecke::neg(&mut self.model, &value),
                Value::RVOre(value) => rv_ore::negate(&mut self.model, &value),
                Value::RVSatake(value) => rv_satake::neg(&mut self.model, &value),
                _ => unreachable!(),
            },
            TypedForm::SphericalProduct(arguments) => {
                let Kind::RVSphHecke(rank) = *self.model.kind(expr.ty) else {
                    unreachable!()
                };
                let [left, right] = arguments.as_ref();
                let values = [self.eval(left)?, self.eval(right)?];
                let [Value::RVSphHecke(left), Value::RVSphHecke(right)] = values else {
                    unreachable!()
                };
                rv_sphhecke::multiply(&mut self.model, rank, &left, &right)
            }
            TypedForm::CycloProduct(order, arguments) => {
                let values = [self.eval(&arguments[0])?, self.eval(&arguments[1])?];
                let [Value::RVCyclo(a), Value::RVCyclo(b)] = values else {
                    unreachable!()
                };
                rv_cyclotomic::multiply(&mut self.model, *order, &a, &b)
            }
            TypedForm::CycloScale(factor, input) => {
                let Value::RVCyclo(terms) = self.eval(input)? else {
                    unreachable!()
                };
                rv_cyclotomic::scale(&mut self.model, &terms, *factor)
            }
            TypedForm::CycloGalois(order, power, input) => {
                let Value::RVCyclo(coefficients) = self.eval(input)? else {
                    unreachable!()
                };
                rv_cyclotomic::galois(&mut self.model, *order, *power, &coefficients)
            }
            TypedForm::CycloEmbed(source, target, input) => {
                let Value::RVCyclo(terms) = self.eval(input)? else {
                    unreachable!()
                };
                rv_cyclotomic::embed(&mut self.model, *source, *target, &terms)
            }
            TypedForm::WzwUnary(grade, operation, input) => {
                let Value::RVWzw(terms) = self.eval(input)? else {
                    unreachable!()
                };
                match *operation {
                    WzwOp::S(adjoint) => rv_wzw::s_action(&mut self.model, *grade, &terms, adjoint),
                    WzwOp::T(power) => rv_wzw::t_action(&mut self.model, *grade, &terms, power),
                    WzwOp::Charge => rv_wzw::charge(&mut self.model, &terms),
                    WzwOp::Galois(power) => rv_wzw::galois(&mut self.model, *grade, &terms, power),
                    WzwOp::CoeffGalois(power) => {
                        rv_wzw::coeff_galois(&mut self.model, *grade, &terms, power)
                    }
                    WzwOp::Negate => rv_wzw::negate(&mut self.model, &terms),
                }
            }
            TypedForm::WzwScale(grade, arguments) => {
                let values = [self.eval(&arguments[0])?, self.eval(&arguments[1])?];
                let [Value::RVCyclo(scalar), Value::RVWzw(state)] = values else {
                    unreachable!()
                };
                rv_wzw::scale(&mut self.model, *grade, &state, &scalar)
            }
            TypedForm::TlOperation(grade, operation, input) => {
                let Value::RvTl(state) = self.eval(input)? else {
                    unreachable!()
                };
                match operation {
                    TlOp::Generator(index) => rv_tl::tl(&mut self.model, *grade, &state, *index),
                    TlOp::Braid(word) => rv_tl::braid(&mut self.model, *grade, &state, word),
                    TlOp::Project(first, length) => {
                        rv_tl::project(&mut self.model, *grade, &state, *first, *length)
                    }
                    TlOp::Negate => rv_tl::negate(&mut self.model, &state),
                }
            }
            TypedForm::TlCoordinate(grade, path, input) => {
                let Value::RvTl(state) = self.eval(input)? else {
                    unreachable!()
                };
                rv_tl::coordinate(&mut self.model, *grade, &state, path)
            }
            TypedForm::TlScale(grade, arguments) => {
                let values = [self.eval(&arguments[0])?, self.eval(&arguments[1])?];
                let [Value::RVCycloQ(scalar), Value::RvTl(state)] = values else {
                    unreachable!()
                };
                rv_tl::scale(&mut self.model, *grade, &state, &scalar)
            }
            TypedForm::CycloRatioProduct(order, arguments) => {
                let values = [self.eval(&arguments[0])?, self.eval(&arguments[1])?];
                let [Value::RVCycloQ(a), Value::RVCycloQ(b)] = values else {
                    unreachable!()
                };
                let value = rv_lens::multiply(&mut self.model, *order, &a, &b)?;
                rv_lens::rational_value(&mut self.model, value)
            }
            TypedForm::CycloRatioGalois(order, power, input) => {
                let Value::RVCycloQ(value) = self.eval(input)? else {
                    unreachable!()
                };
                rv_lens::galois(&mut self.model, *order, &value, *power)
            }
            TypedForm::CycloRatio(denominator, input) => {
                let Value::RVCyclo(numerator) = self.eval(input)? else {
                    unreachable!()
                };
                rv_lens::ratio(&mut self.model, &numerator, *denominator)
            }
            TypedForm::CycloRatioNegate(input) => {
                let Value::RVCycloQ(ratio) = self.eval(input)? else {
                    unreachable!()
                };
                rv_lens::negate(&mut self.model, &ratio)
            }
            TypedForm::WzwCoordinate(grade, weight, input) => {
                let Value::RVWzw(terms) = self.eval(input)? else {
                    unreachable!()
                };
                rv_wzw::coordinate(&mut self.model, *grade, &terms, weight.clone())
            }
            TypedForm::CycloConjugate(order, input) => {
                let Value::RVCyclo(terms) = self.eval(input)? else {
                    unreachable!()
                };
                rv_cyclotomic::conjugate(&mut self.model, *order, &terms)
            }
            TypedForm::Verlinde(order, point, input) => {
                let Value::RVFusion(terms) = self.eval(input)? else {
                    unreachable!()
                };
                rv_verlinde::character(&mut self.model, *order, point, &terms)
            }
            TypedForm::FusionProduct(grade, arguments) => {
                let values = [self.eval(&arguments[0])?, self.eval(&arguments[1])?];
                let [Value::RVFusion(a), Value::RVFusion(b)] = values else {
                    unreachable!()
                };
                rv_fusion::multiply(&mut self.model, *grade, &a, &b)
            }
            TypedForm::FusionScale(grade, factor, input) => {
                let Value::RVFusion(terms) = self.eval(input)? else {
                    unreachable!()
                };
                rv_fusion::scale(&mut self.model, *grade, &terms, *factor)
            }
            TypedForm::FusionDual(grade, input) => {
                let Value::RVFusion(terms) = self.eval(input)? else {
                    unreachable!()
                };
                rv_fusion::dual(&mut self.model, *grade, &terms)
            }
            TypedForm::SatakeInverse(input) => {
                let Value::RVSatake(terms) = self.eval(input)? else {
                    unreachable!()
                };
                rv_satake::inverse(&mut self.model, &terms)
            }
            TypedForm::OreUnary(operation, input) => {
                let Value::RVOre(operator) = self.eval(input)? else {
                    unreachable!()
                };
                match operation {
                    OreOp::Factorial(height) => {
                        rv_qregular::factorial(&mut self.model, *height, &operator)
                    }
                    OreOp::EulerFactor => rv_qregular::euler_factor(&mut self.model, &operator),
                    OreOp::EulerReduce => rv_qregular::euler_reduce(&mut self.model, &operator),
                }
            }
            TypedForm::OreProduct(arguments) => {
                let values = [self.eval(&arguments[0])?, self.eval(&arguments[1])?];
                let [Value::RVOre(a), Value::RVOre(b)] = values else {
                    unreachable!()
                };
                rv_ore::multiply(&mut self.model, &a, &b)
            }
            TypedForm::OreApply(grade, arguments) => {
                let values = [self.eval(&arguments[0])?, self.eval(&arguments[1])?];
                let [Value::RVOre(a), Value::RVQGr(b)] = values else {
                    unreachable!()
                };
                rv_qconnection::apply(&mut self.model, *grade, &a, &b)
            }
            TypedForm::QuantumReduce(grade, input) => {
                let Value::RVQGr(terms) = self.eval(input)? else {
                    unreachable!()
                };
                rv_pcurv::reduction(&mut self.model, *grade, &terms)
            }
            TypedForm::QuantumPNabla(grade, input) => {
                let Value::RVQGrP(terms) = self.eval(input)? else {
                    unreachable!()
                };
                rv_pcurv::nabla(&mut self.model, *grade, &terms)
            }
            TypedForm::QuantumPUnary(grade, operation, input) => {
                let Value::RVQGrP(terms) = self.eval(input)? else {
                    unreachable!()
                };
                match *operation {
                    QuantumOp::Scale(factor) => {
                        rv_pcurv::scale(&mut self.model, grade.prime, &terms, factor)
                    }
                    QuantumOp::Shift(power) => rv_pcurv::shift(&mut self.model, &terms, power),
                    QuantumOp::Transpose => {
                        rv_pcurv::transpose(&mut self.model, grade.prime, &terms)
                    }
                    QuantumOp::Classical => rv_pcurv::classical(&mut self.model, &terms),
                }
            }
            TypedForm::QuantumCurvature(grade, satake, input) => {
                let Value::RVQGrP(terms) = self.eval(input)? else {
                    unreachable!()
                };
                rv_pcurv::apply(&mut self.model, *grade, &terms, *satake)
            }
            TypedForm::QuantumNabla(grade, input) => {
                let Value::RVQGr(v) = self.eval(input)? else {
                    unreachable!()
                };
                rv_qconnection::nabla(&mut self.model, *grade, &v)
            }
            TypedForm::PeriodInitial(grade) => {
                let period = self.model.rv_period_for(grade);
                Ok(Value::PeriodState(
                    period.initial(&mut self.model, grade.level)?,
                ))
            }
            TypedForm::PeriodQuery(grade, index) => {
                let period = self.model.rv_period_for(grade);
                Ok(Value::PeriodResidue(period.query(
                    &mut self.model,
                    grade.level,
                    *index,
                )?))
            }
            TypedForm::PeriodSection(grade, digit, input) => {
                let Value::PeriodState(state) = self.eval(input)? else {
                    unreachable!()
                };
                let period = self.model.rv_period_for(grade);
                let image = period.section(&mut self.model, grade.level, *digit, &state)?;
                Ok(Value::PeriodState(rv_period::persist(
                    &mut self.model,
                    image,
                )?))
            }
            TypedForm::PeriodRead(grade, input) => {
                let Value::PeriodState(state) = self.eval(input)? else {
                    unreachable!()
                };
                let period = self.model.rv_period_for(grade);
                Ok(Value::PeriodResidue(period.read(&mut self.model, &state)?))
            }
            TypedForm::QuantumProduct(grade, arguments) => {
                let values = [self.eval(&arguments[0])?, self.eval(&arguments[1])?];
                let [Value::RVQGr(a), Value::RVQGr(b)] = values else {
                    unreachable!()
                };
                rv_qgr::multiply(&mut self.model, *grade, &a, &b)
            }
            TypedForm::QuantumUnary(operation, input) => {
                let Value::RVQGr(terms) = self.eval(input)? else {
                    unreachable!()
                };
                match *operation {
                    QuantumOp::Scale(factor) => rv_qgr::scale(&mut self.model, &terms, factor),
                    QuantumOp::Shift(power) => rv_qgr::shift(&mut self.model, &terms, power),
                    QuantumOp::Transpose => rv_qgr::transpose(&mut self.model, &terms),
                    QuantumOp::Classical => rv_qgr::classical(&mut self.model, &terms),
                }
            }
            TypedForm::SatakeTransform(input) => {
                let Value::RVSphHecke(terms) = self.eval(input)? else {
                    unreachable!()
                };
                rv_satake::transform(&mut self.model, &terms)
            }
            TypedForm::SatakeProduct(arguments) => {
                let [left, right] = arguments.as_ref();
                let values = [self.eval(left)?, self.eval(right)?];
                let [Value::RVSatake(left), Value::RVSatake(right)] = values else {
                    unreachable!()
                };
                rv_satake::multiply(&mut self.model, &left, &right)
            }
            TypedForm::SphericalEmbedding(input) => {
                let Kind::RVIwHecke(rank) = *self.model.kind(expr.ty) else {
                    unreachable!()
                };
                let Value::RVSphHecke(support) = self.eval(input)? else {
                    unreachable!()
                };
                rv_sphhecke::embed(&mut self.model, rank, &support)
            }
            TypedForm::HeckeProduct(arguments) => {
                let [left, right] = arguments.as_ref();
                let values = [self.eval(left)?, self.eval(right)?];
                let [Value::RVIwHecke(left), Value::RVIwHecke(right)] = values else {
                    unreachable!()
                };
                rv_iwhecke::multiply(&mut self.model, &left, &right)
            }
            TypedForm::RVLaurentMonomial(power, input) => {
                let Value::RVHom(value) = self.eval(input)? else {
                    unreachable!()
                };
                rv_laurent::monomial(&mut self.model, &value, *power)
            }
            TypedForm::RVLaurentFrobenius(field, powers, input) => {
                let Value::RVLaurent(value) = self.eval(input)? else {
                    unreachable!()
                };
                let image = rv_laurent::frobenius(&mut self.model, *field, &value, *powers)?;
                rv_laurent::persist(&mut self.model, image)
            }
            TypedForm::RVIwAdlv(boundary, step, target, arguments) => {
                let [a, x] = arguments.as_ref();
                let values = [self.eval(a)?, self.eval(x)?];
                let [Value::RVLaurent(action), Value::RVLaurent(frame)] = values else {
                    unreachable!()
                };
                let actual = rv_iwahori::at(
                    &mut self.model,
                    a.ty,
                    &action,
                    &frame,
                    *step,
                    boundary.name(),
                )?;
                let contains = boundary.contains(&mut self.model, &actual, target)?;
                Ok(Value::Point(usize::from(contains)))
            }
            TypedForm::RVAdlv(boundary, step, mu, arguments) => {
                let [a, x] = arguments.as_ref();
                let values = [self.eval(a)?, self.eval(x)?];
                let [Value::RVLaurent(action), Value::RVLaurent(frame)] = values else {
                    unreachable!()
                };
                let position = rv_hodge::mu_at(
                    &mut self.model,
                    a.ty,
                    &action,
                    &frame,
                    *step,
                    boundary.name(),
                )?;
                let contains = boundary.contains(&mut self.model, &position, mu)?;
                Ok(Value::Point(usize::from(contains)))
            }
            TypedForm::RVNewton(operation, embedding, input) => {
                let Value::RVLaurent(value) = self.eval(input)? else {
                    unreachable!()
                };
                rv_newton::evaluate(&mut self.model, input.ty, embedding, &value, *operation)
            }
            TypedForm::RVScale(field, scalar, input) => match self.eval(input)? {
                Value::RVLaurent(value) => {
                    rv_laurent::scale(&mut self.model, *field, &value, *scalar)
                }
                Value::RVHom(value) => rv_hom::scale(&mut self.model, *field, &value, *scalar),
                Value::JetHom(value) => jet_hom::scale(&mut self.model, *field, &value, *scalar),
                Value::JetClass(value) => {
                    let scaled = crate::binary_linear::add_multiple(
                        &mut self.model,
                        *field,
                        &[],
                        &value,
                        *scalar,
                    )?;
                    jet_page::persist(&mut self.model, scaled)
                }
                _ => unreachable!(),
            },
            TypedForm::RVField(map, input) => {
                let Value::RVHom(value) = self.eval(input)? else {
                    unreachable!()
                };
                map.apply(&mut self.model, &value)
            }
            TypedForm::RVInverse(input) => {
                let Value::RVHom(value) = self.eval(input)? else {
                    unreachable!()
                };
                rv_descent::inverse(&mut self.model, input.ty, &value)
            }
            TypedForm::RVHilbert90(embedding, input) => {
                let Value::RVHom(value) = self.eval(input)? else {
                    unreachable!()
                };
                rv_descent::hilbert(&mut self.model, embedding, input.ty, &value)
            }
            TypedForm::RVTwisted(operation, embedding, arguments) => {
                let [left, input, right] = arguments.as_ref();
                let types = [left.ty, input.ty, right.ty];
                let values = [self.eval(left)?, self.eval(input)?, self.eval(right)?];
                let [Value::RVHom(left), Value::RVHom(input), Value::RVHom(right)] = values else {
                    unreachable!()
                };
                rv_descent::twisted(
                    &mut self.model,
                    *operation,
                    embedding,
                    types,
                    [&left, &input, &right],
                )
            }
            TypedForm::RVShintani(operation, embedding, input) => {
                let Value::RVHom(value) = self.eval(input)? else {
                    unreachable!()
                };
                rv_shintani::evaluate(&mut self.model, input.ty, embedding, &value, *operation)
            }
            TypedForm::RVSpecialize(input) => {
                let value = self.eval(input)?;
                rv_hom::specialize(&mut self.model, input.ty, value, expr.ty)
            }
            TypedForm::QRestrict(input) => {
                let value = self.eval(input)?;
                q_koszul::restrict(&mut self.model, input.ty, value, expr.ty)
            }
            TypedForm::BpRaise(input) => {
                let value = self.eval(input)?;
                if input.ty == expr.ty {
                    return Ok(value);
                }
                match value {
                    Value::BP(polynomial) => {
                        let (&Kind::BP(from), &Kind::BP(to)) =
                            (self.model.kind(input.ty), self.model.kind(expr.ty))
                        else {
                            unreachable!()
                        };
                        let output = bp::reindex_height(
                            &mut self.model,
                            from.height,
                            to.height,
                            &polynomial,
                        )?;
                        bp::persist(&mut self.model, output)
                    }
                    Value::OddClass(class) => {
                        if class.is_empty() {
                            return Ok(odd_ext::zero());
                        }
                        let source = self.model.bp_page(input.ty)?;
                        let target = self.model.bp_page(expr.ty)?;
                        let class = source.raise(&mut self.model, &target, &class)?;
                        odd_ext::persist(&mut self.model, class)
                    }
                    Value::BPExt(class) => {
                        if class.is_empty() {
                            return Ok(bp_classes::zero());
                        }
                        let (&Kind::BPExt(from), &Kind::BPExt(to)) =
                            (self.model.kind(input.ty), self.model.kind(expr.ty))
                        else {
                            unreachable!()
                        };
                        let cycle = self
                            .model
                            .bp_ext(input.ty)?
                            .representative(&mut self.model, &class)?;
                        let output = bp::reindex_height(
                            &mut self.model,
                            from.shape.height,
                            to.shape.height,
                            &cycle,
                        )?;
                        bp_classes::class(&mut self.model, expr.ty, &output)
                    }
                    _ => unreachable!(),
                }
            }
            TypedForm::BpTriple(layout, inputs) => {
                let values = self.eval_bp_triple_inputs(inputs)?;
                let coset = layout.build(
                    &mut self.model,
                    values.each_ref().map(|(ty, class)| (*ty, class)),
                )?;
                Ok(Value::BPTriple(Rc::new(coset)))
            }
            TypedForm::BpTripleBase(input) => {
                let Value::BPTriple(coset) = self.eval(input)? else {
                    unreachable!()
                };
                Ok(coset.base())
            }
            TypedForm::BpHigher(layout, inputs, witnesses) => {
                let inputs = inputs
                    .iter()
                    .map(|input| {
                        let Value::BPExt(class) = self.eval(input)? else {
                            unreachable!()
                        };
                        Ok((input.ty, class))
                    })
                    .collect::<Result<Vec<_>>>()?;
                let witnesses = witnesses
                    .as_ref()
                    .map(|witnesses| {
                        witnesses
                            .iter()
                            .map(|witness| {
                                let Value::BP(polynomial) = self.eval(witness)? else {
                                    unreachable!()
                                };
                                Ok((witness.span, polynomial))
                            })
                            .collect::<Result<Vec<_>>>()
                    })
                    .transpose()?;
                layout.apply(&mut self.model, inputs, witnesses, expr.ty)
            }
            TypedForm::BpMassey(layout, inputs) => {
                let values = inputs
                    .iter()
                    .map(|input| self.eval(input))
                    .map(|value| {
                        let Value::BP(polynomial) = value? else {
                            unreachable!()
                        };
                        Ok(polynomial)
                    })
                    .collect::<Result<Vec<_>>>()?;
                let values: [_; 5] = values.try_into().unwrap_or_else(|_| unreachable!());
                let result = layout.apply(&mut self.model, values.each_ref().map(Rc::as_ref))?;
                bp::persist(&mut self.model, result)
            }
            TypedForm::BpModExtClass(input) => {
                let Value::BP(polynomial) = self.eval(input)? else {
                    unreachable!()
                };
                bp_mod_ext::class(&mut self.model, expr.ty, &polynomial)
            }
            TypedForm::BpModExtRep(input) => {
                let Value::OddClass(class) = self.eval(input)? else {
                    unreachable!()
                };
                if class.is_empty() {
                    return Ok(bp::zero());
                }
                let representative = self
                    .model
                    .bp_mod_ext(input.ty)?
                    .representative(&mut self.model, &class)?;
                bp::persist(&mut self.model, representative)
            }
            TypedForm::BpModExtOperation(operation, input) => {
                let Value::OddClass(class) = self.eval(input)? else {
                    unreachable!()
                };
                bp_mod_ext::apply(&mut self.model, *operation, input.ty, expr.ty, &class)
            }
            TypedForm::BpVPageClass(input) => {
                let Value::BP(polynomial) = self.eval(input)? else {
                    unreachable!()
                };
                bp_vpage::class(&mut self.model, expr.ty, &polynomial)
            }
            TypedForm::BpVPageOperation(operation, input) => {
                let Value::OddClass(coordinates) = self.eval(input)? else {
                    unreachable!()
                };
                bp_vpage::apply(&mut self.model, *operation, input.ty, expr.ty, &coordinates)
            }
            TypedForm::BpPageClass(input) => {
                let Value::BP(polynomial) = self.eval(input)? else {
                    unreachable!()
                };
                if polynomial.is_zero() {
                    return Ok(odd_ext::zero());
                }
                let coordinates = self
                    .model
                    .bp_page(expr.ty)?
                    .class(&mut self.model, &polynomial)?;
                odd_ext::persist(&mut self.model, coordinates)
            }
            TypedForm::BpPageOperation(operation, input) => {
                use bp_page::Operation::{Differential, Next, Previous, Representative};
                let Value::OddClass(coordinates) = self.eval(input)? else {
                    unreachable!()
                };
                if coordinates.is_empty() {
                    return Ok(match operation {
                        Representative => bp::zero(),
                        _ => odd_ext::zero(),
                    });
                }
                let source = self.model.bp_page(input.ty)?;
                if matches!(operation, Representative) {
                    let polynomial = source.representative(&mut self.model, &coordinates, true)?;
                    return bp::persist(&mut self.model, polynomial);
                }
                let target = self.model.bp_page(expr.ty)?;
                let output = match operation {
                    Differential => source.differential(&mut self.model, &target, &coordinates)?,
                    Next => source.next(&mut self.model, &target, &coordinates)?,
                    Previous => source.previous(&mut self.model, &target, &coordinates)?,
                    Representative => unreachable!(),
                };
                odd_ext::persist(&mut self.model, output)
            }
            TypedForm::BpHeightOperation(operation, input) => {
                use BpHeightOp::{Class, Differential, Next, Previous, Representative};
                let value = self.eval(input)?;
                if matches!(operation, Class(..)) {
                    let Value::BP(cycle) = value else {
                        unreachable!()
                    };
                    return bp_classes::class(&mut self.model, expr.ty, &cycle);
                }
                let Value::BPExt(support) = value else {
                    unreachable!()
                };
                if matches!(operation, Representative) {
                    return bp_classes::representative(&mut self.model, input.ty, &support);
                }
                let Kind::BPHeight(grade) = *self.model.kind(input.ty) else {
                    unreachable!();
                };
                if support.is_empty()
                    || matches!(operation, Differential) && grade.h < grade.r as i128
                {
                    return Ok(bp_classes::zero());
                }
                let group = self.model.bp_ext(input.ty)?;
                let cycle = group.representative(&mut self.model, &support)?;
                let result = match operation {
                    Differential => bp::differential(&mut self.model, grade.top, &cycle)?,
                    Next => grade.correct_next(&mut self.model, cycle)?,
                    Previous => cycle,
                    _ => unreachable!(),
                };
                bp_classes::class(&mut self.model, expr.ty, &result)
            }
            TypedForm::BpRelOperation(operation, shape, input) => {
                let value = self.eval(input)?;
                match (*operation, value) {
                    (BpRelOp::Quotient(_), Value::BP(cochain)) => {
                        let quotient = bp_rel::quotient(&mut self.model, *shape, &cochain)?;
                        bp_rel::persist(&mut self.model, quotient)
                    }
                    (BpRelOp::Quotient(_), Value::BPExt(class)) => {
                        if class.is_empty() {
                            return Ok(bp_classes::zero());
                        }
                        let representative = self
                            .model
                            .bp_ext(input.ty)?
                            .representative(&mut self.model, &class)?;
                        let quotient = shape.domain().project(&mut self.model, &representative)?;
                        bp_classes::class(&mut self.model, expr.ty, &quotient)
                    }
                    (BpRelOp::Class(_), Value::BPRel(cycle)) => {
                        bp_classes::class(&mut self.model, expr.ty, &cycle)
                    }
                    (BpRelOp::Representative, Value::BPExt(class)) => {
                        bp_classes::representative(&mut self.model, input.ty, &class)
                    }
                    (BpRelOp::HeightDelta, Value::BPRel(cycle)) => {
                        let image = bp_rel::connecting(&mut self.model, *shape, &cycle)?;
                        bp::persist(&mut self.model, image)
                    }
                    (BpRelOp::HeightDelta, Value::BPExt(class)) => {
                        if class.is_empty() {
                            return Ok(bp_classes::zero());
                        }
                        let representative = self
                            .model
                            .bp_ext(input.ty)?
                            .representative(&mut self.model, &class)?;
                        let image = bp_rel::connecting(&mut self.model, *shape, &representative)?;
                        bp_classes::class(&mut self.model, expr.ty, &image)
                    }
                    (BpRelOp::Lift, Value::BPRel(cochain)) => Ok(Value::BP(cochain)),
                    (BpRelOp::Differential, Value::BPRel(cochain)) => {
                        let result = shape.domain().differential(&mut self.model, &cochain)?;
                        bp_rel::persist(&mut self.model, result)
                    }
                    (BpRelOp::Null(t), Value::BPRel(cycle)) => {
                        let witness = shape.domain().null(&mut self.model, t, &cycle)?;
                        bp_rel::persist(&mut self.model, witness)
                    }
                    _ => unreachable!(),
                }
            }
            TypedForm::BpExtOperation(operation, input) => match (*operation, self.eval(input)?) {
                (BpExtOp::Class(_), Value::BP(cycle)) => {
                    bp_classes::class(&mut self.model, expr.ty, &cycle)
                }
                (BpExtOp::Representative, Value::BPExt(support)) => {
                    bp_classes::representative(&mut self.model, input.ty, &support)
                }
                (BpExtOp::Scale(scalar), Value::BPRel(cochain)) => {
                    self.model.spend(cochain.cells())?;
                    let mut result = (*cochain).clone();
                    result.scale(&mut self.model, scalar)?;
                    bp_rel::persist(&mut self.model, result)
                }
                (BpExtOp::Scale(scalar), Value::BPExt(support)) => {
                    bp_classes::scale(&mut self.model, input.ty, &support, scalar)
                }
                _ => unreachable!(),
            },
            TypedForm::BpOperation(operation, shape, input) => {
                let Value::BP(polynomial) = self.eval(input)? else {
                    unreachable!()
                };
                let result = operation.apply(&mut self.model, *shape, &polynomial)?;
                bp::persist(&mut self.model, result)
            }
            TypedForm::RvHasse(orders, input) => {
                let Value::RVPoly(polynomial) = self.eval(input)? else {
                    unreachable!()
                };
                let result = rv_poly::hasse(&mut self.model, &polynomial, orders)?;
                Ok(Value::RVPoly(rv_poly::persist(&mut self.model, result)?))
            }
            TypedForm::RvVanish(field, point, input) => {
                let Value::RVPoly(polynomial) = self.eval(input)? else {
                    unreachable!()
                };
                let vanishes = rv_poly::vanishes(&mut self.model, &polynomial, *field, point)?;
                Ok(Value::Point(usize::from(vanishes)))
            }
            TypedForm::RvNCTerm(power, input) => {
                let Value::RVHHChain(words) = self.eval(input)? else {
                    unreachable!()
                };
                Ok(Value::RVNCChain(rv_negative::term(
                    &mut self.model,
                    *power,
                    &words,
                )?))
            }
            TypedForm::RvNCPart(power, input) => {
                let Value::RVNCChain(terms) = self.eval(input)? else {
                    unreachable!()
                };
                Ok(Value::RVHHChain(rv_negative::part(
                    &mut self.model,
                    *power,
                    &terms,
                )?))
            }
            TypedForm::RvCyclicShift(shift, input) => match self.eval(input)? {
                Value::RVNCChain(terms) => {
                    let Kind::RvCyclicChain(grade) = self.model.kind(input.ty) else {
                        unreachable!()
                    };
                    let mode = grade.mode;
                    let terms = shift.terms(&mut self.model, mode, &terms)?;
                    Ok(Value::RVNCChain(rv_negative::persist(
                        &mut self.model,
                        terms,
                    )?))
                }
                Value::Class(class) => rv_cyclic_homology::shift_class(
                    &mut self.model,
                    *shift,
                    input.ty,
                    expr.ty,
                    &class,
                ),
                _ => unreachable!(),
            },
            TypedForm::RvNCTransfer(operation, input) => match self.eval(input)? {
                Value::RVNCChain(terms) => {
                    let Kind::RVNCChain(grade) = self.model.kind(expr.ty).clone() else {
                        unreachable!()
                    };
                    let support = operation.terms(&mut self.model, &grade, &terms)?;
                    Ok(Value::RVNCChain(rv_negative::persist(
                        &mut self.model,
                        support,
                    )?))
                }
                Value::Class(class) => rv_negative::transfer_class(
                    &mut self.model,
                    *operation,
                    input.ty,
                    expr.ty,
                    &class,
                ),
                _ => unreachable!(),
            },
            TypedForm::RvHHWord(inputs) => {
                let classes = inputs
                    .iter()
                    .map(|item| {
                        let Value::Class(class) = self.eval(item)? else {
                            unreachable!()
                        };
                        Ok((item.ty, class))
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(Value::RVHHChain(rv_hh::word(&mut self.model, &classes)?))
            }
            TypedForm::RvHHLength(length, input) => {
                let Value::RVHHChain(terms) = self.eval(input)? else {
                    unreachable!()
                };
                Ok(Value::RVHHChain(rv_hh::length_part(
                    &mut self.model,
                    *length,
                    &terms,
                )?))
            }
            TypedForm::RvBarWord(ends, degrees, inputs) => {
                let classes = inputs
                    .iter()
                    .map(|input| {
                        let Value::Class(class) = self.eval(input)? else {
                            unreachable!();
                        };
                        Ok(class)
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(Value::RVBar(rv_bar::word(
                    &mut self.model,
                    ends,
                    degrees,
                    &classes,
                )?))
            }
            TypedForm::RvBarLength(length, input) => {
                let Value::RVBar(terms) = self.eval(input)? else {
                    unreachable!();
                };
                Ok(Value::RVBar(rv_bar::length_part(
                    &mut self.model,
                    *length,
                    &terms,
                )?))
            }
            TypedForm::RvMComponent(input) => {
                let value = self.eval(input)?;
                self.rv_koszul_regrade(input.ty, expr.ty, value)
            }
            TypedForm::RvKoszulTerm(wedge, ideal, input) => {
                let Value::RVPoly(coefficient) = self.eval(input)? else {
                    unreachable!();
                };
                Ok(Value::RVKoszul(rv_koszul::term(
                    &mut self.model,
                    ideal,
                    *wedge,
                    coefficient,
                )?))
            }
            TypedForm::DiffScale(field, scalar, input) => {
                let Value::CurveDiff(coordinates) = self.eval(input)? else {
                    unreachable!()
                };
                let scaled = binary_linear::scaled(&mut self.model, *field, &coordinates, *scalar)?;
                Ok(Value::CurveDiff(binary_linear::persist(
                    &mut self.model,
                    scaled,
                )?))
            }
            TypedForm::AsDeRhamOperation(operation, cover, field, input) => {
                let value = self.eval(input)?;
                rv_as_derham::apply(&mut self.model, cover, *field, *operation, value)
            }
            TypedForm::AsOddTateOperation(operation, cover, degree, input) => {
                let Value::OddClass(coordinates) = self.eval(input)? else {
                    unreachable!()
                };
                if coordinates.is_empty() {
                    return Ok(odd_ext::zero());
                }
                cover.tate(&mut self.model)?.apply(
                    &mut self.model,
                    cover,
                    *degree,
                    *operation,
                    coordinates,
                )
            }
            TypedForm::AsOddDeRhamOperation(operation, cover, input) => {
                let Value::OddClass(coordinates) = self.eval(input)? else {
                    unreachable!()
                };
                if coordinates.is_empty() {
                    return Ok(odd_ext::zero());
                }
                cover.derham(&mut self.model)?.apply(
                    &mut self.model,
                    cover,
                    *operation,
                    coordinates,
                )
            }
            TypedForm::AsTwoDiffOperation(operation, cover, input) => {
                let Value::OddClass(coordinates) = self.eval(input)? else {
                    unreachable!()
                };
                if coordinates.is_empty() {
                    return Ok(odd_ext::zero());
                }
                cover
                    .cartier(&mut self.model)?
                    .apply(&mut self.model, *operation, coordinates)
            }
            TypedForm::AsOddDiffOperation(operation, cover, input) => {
                let Value::OddClass(coordinates) = self.eval(input)? else {
                    unreachable!()
                };
                if coordinates.is_empty() {
                    return Ok(odd_ext::zero());
                }
                cover.cartier(&mut self.model)?.apply(
                    &mut self.model,
                    cover,
                    *operation,
                    coordinates,
                )
            }
            TypedForm::AsCartier(cover, field, input) => {
                let Value::CurveDiff(coordinates) = self.eval(input)? else {
                    unreachable!()
                };
                rv_as_cartier::apply(&mut self.model, cover, *field, &coordinates)
            }
            TypedForm::WittCartier(cover, field, input) => {
                let Value::CurveDiff(coordinates) = self.eval(input)? else {
                    unreachable!()
                };
                cover.apply_cartier(&mut self.model, *field, &coordinates)
            }
            TypedForm::RvIdeal(op, ideal, input) => {
                let Value::RVPoly(original) = self.eval(input)? else {
                    unreachable!()
                };
                let polynomial = match op {
                    RvIdealOp::Normal => ideal.normal(&mut self.model, original)?,
                    RvIdealOp::Contains => {
                        let contains = ideal.normal(&mut self.model, original)?.is_empty();
                        return Ok(Value::Point(usize::from(contains)));
                    }
                    RvIdealOp::Witness(index) => {
                        ideal.witness(&mut self.model, *index, original)?
                    }
                };
                Ok(Value::RVPoly(rv_poly::persist(
                    &mut self.model,
                    polynomial,
                )?))
            }
            TypedForm::RvInMinors(base, selected, size, input) => {
                let Value::RVPoly(polynomial) = self.eval(input)? else {
                    unreachable!()
                };
                if polynomial.is_empty() {
                    return Ok(Value::Point(1));
                }
                let basis = self.model.rv_minors(*base, selected.clone(), *size)?;
                let contained = basis.contains(&mut self.model, polynomial)?;
                Ok(Value::Point(usize::from(contained)))
            }
            TypedForm::VShift(k, input) => {
                let value = self.eval(input)?;
                if matches!(self.model.kind(input.ty), Kind::VBock(..)) && *k != 0 {
                    let Value::Class(class) = value else {
                        unreachable!()
                    };
                    return vbock::shift(&mut self.model, input.ty, expr.ty, *k, &class);
                }
                v_koszul::shift(&mut self.model, input.ty, expr.ty, *k, value)
            }
            TypedForm::VRatio(k, input) => {
                let value = self.eval(input)?;
                v_bars::ratio(&mut self.model, input.ty, expr.ty, *k, value)
            }
            TypedForm::RestrictAlgebra(profile, input) => {
                let value = self.eval(input)?;
                if input.ty == expr.ty {
                    return Ok(value);
                }
                match value {
                    Value::Cobar(words) => {
                        let words = profile.restrict(&mut self.model, &words)?;
                        self.model.allocate(words.len())?;
                        Ok(Value::Cobar(words.into()))
                    }
                    Value::Class(class) if self.model.kind(input.ty).may_grade().is_none() => {
                        adams::restrict(&mut self.model, input.ty, &class, expr.ty)
                    }
                    Value::PCobar(terms) => {
                        let Kind::PCobarAn(p, n, ..) = *self.model.kind(expr.ty) else {
                            unreachable!()
                        };
                        let terms = odd_profile::restrict(&mut self.model, p, n, &terms)?;
                        odd_cobar::persist(&mut self.model, terms)
                    }
                    Value::OddClass(class) => {
                        if self.model.kind(input.ty).pmay_grade().is_some() {
                            odd_may::restrict(&mut self.model, input.ty, &class, expr.ty)
                        } else {
                            odd_ext::restrict(&mut self.model, input.ty, &class, expr.ty)
                        }
                    }
                    Value::Class(class) => {
                        if class.is_empty() {
                            return Ok(cohomology::zero());
                        }
                        let representative = self
                            .model
                            .may_page(input.ty)?
                            .representative(&mut self.model, &class)?;
                        let words = profile.restrict(&mut self.model, &representative)?;
                        may::class(&mut self.model, expr.ty, &words)
                    }
                    _ => unreachable!(),
                }
            }
            TypedForm::CobarCup(profile, i, left, right) => {
                let Value::Cobar(a) = self.eval(left)? else {
                    unreachable!()
                };
                let Value::Cobar(b) = self.eval(right)? else {
                    unreachable!()
                };
                cobar_squares::cup(&mut self.model, *profile, *i, &a, &b)
            }
            TypedForm::ExtSquare(k, input) => {
                let Value::Class(class) = self.eval(input)? else {
                    unreachable!()
                };
                cobar_squares::square(&mut self.model, input.ty, expr.ty, *k, &class)
            }
            TypedForm::Nakamura(k, input) => {
                let (profile, s, _) = self.model.kind(input.ty).cobar_grade().unwrap();
                let Value::Cobar(words) = self.eval(input)? else {
                    unreachable!()
                };
                nakamura::apply(&mut self.model, profile, *k, s, &words)
            }
            TypedForm::MayLiftTo(r, f, input) => {
                let cochain = self.eval(input)?;
                may_lift::evaluate(&mut self.model, input.ty, *r, *f, cochain)
            }
            TypedForm::MapExtInduced(function) => {
                let Value::Map(table) = self.eval(function)? else {
                    unreachable!()
                };
                map_ext::induced(&mut self.model, function.ty, &table, expr.ty)
            }
            TypedForm::MapWord(source, target, input) => {
                let Value::Cobar(words) = self.eval(input)? else {
                    unreachable!()
                };
                let terms = map_cochain::elementary(&mut self.model, *source, *target, &words)?;
                Ok(Value::MapCobar(terms))
            }
            TypedForm::Equator => {
                let Kind::Arrow(base, _) = *self.model.kind(expr.ty) else {
                    unreachable!()
                };
                let points = self.model.elements(base)?;
                self.model.allocate(1 + 2 * points.len())?;
                self.model.spend(points.len())?;
                let rows = points
                    .iter()
                    .map(|point| {
                        Value::Suspension(SuspensionPoint::Equator(Rc::new(point.clone())))
                    })
                    .collect();
                Ok(Value::Map(rows))
            }
            TypedForm::Suspend(value) => {
                let Value::Class(class) = self.eval(value)? else {
                    unreachable!()
                };
                cohomology::suspend(&mut self.model, value.ty, &class, expr.ty)
            }
            TypedForm::Shift(value) => self.eval(value),
            TypedForm::JetOperation(op, input) => match (*op, self.eval(input)?) {
                (JetOp::Class(..), Value::JetHom(terms)) => {
                    jet_page::class(&mut self.model, expr.ty, &terms)
                }
                (JetOp::Representative, Value::JetClass(class)) => {
                    let terms = jet_page::representative(&mut self.model, input.ty, &class)?;
                    jet_hom::persist(&mut self.model, terms)
                }
                (JetOp::Differential, Value::JetHom(terms)) => {
                    if terms.is_empty() {
                        return Ok(jet_hom::zero());
                    }
                    let Kind::JetHom(grade) = *self.model.kind(input.ty) else {
                        unreachable!()
                    };
                    let path = self.model.jet_path(grade.diagram.path);
                    let terms = self.model.jet_complex(grade.diagram)?.differential(
                        &mut self.model,
                        &path,
                        &terms,
                    )?;
                    jet_hom::persist(&mut self.model, terms)
                }
                (JetOp::Differential, Value::JetClass(class)) => {
                    jet_page::differential(&mut self.model, input.ty, expr.ty, &class)
                }
                (JetOp::Next, Value::JetClass(class)) => {
                    jet_page::next(&mut self.model, input.ty, expr.ty, &class)
                }
                (JetOp::Prev, Value::JetClass(class)) => {
                    jet_page::prev(&mut self.model, input.ty, expr.ty, &class)
                }
                _ => unreachable!(),
            },
            TypedForm::SheafPageOperation(operation, complex, grade, input) => {
                let Value::OddClass(coordinates) = self.eval(input)? else {
                    unreachable!()
                };
                sheaf_page::apply(
                    &mut self.model,
                    complex.filtered(),
                    *grade,
                    *operation,
                    &coordinates,
                )
            }
            TypedForm::SheafPageLift(operation, complex, grade, input) => {
                let Value::SheafC(terms) = self.eval(input)? else {
                    unreachable!()
                };
                match operation {
                    SheafPageLift::Class => {
                        sheaf_page::class(&mut self.model, complex.filtered(), *grade, &terms)
                    }
                    SheafPageLift::Correct => {
                        sheaf_page::correct(&mut self.model, complex.filtered(), *grade, &terms)
                    }
                }
            }
            TypedForm::RvHHToMixed(origin, input) => match self.eval(input)? {
                Value::RVHHChain(terms) => origin.chain_to_mixed(&mut self.model, &terms),
                Value::Class(class) => {
                    origin.class_to_mixed(&mut self.model, input.ty, expr.ty, &class)
                }
                _ => unreachable!(),
            },
            TypedForm::MixedConeMap(cone, operation, input) => {
                self.eval_mixed_cone_map(cone, operation, input, expr.ty)
            }
            TypedForm::MixedConeShear(target, source, map, homotopy, input) => {
                let Value::OddClass(homotopy) = self.eval(homotopy)? else {
                    unreachable!()
                };
                let Value::OddClass(input) = self.eval(input)? else {
                    unreachable!()
                };
                let image = target.shear(&mut self.model, source, *map, &homotopy, &input)?;
                odd_ext::persist(&mut self.model, image)
            }
            TypedForm::MixedRawDelta(data, order, input) => {
                let Value::OddClass(support) = self.eval(input)? else {
                    unreachable!()
                };
                let image = data.delta(&mut self.model, *order, &support)?;
                odd_ext::persist(&mut self.model, image)
            }
            TypedForm::MixedTensorProduct(tensor, left, right) => {
                let Value::OddClass(left) = self.eval(left)? else {
                    unreachable!()
                };
                let Value::OddClass(right) = self.eval(right)? else {
                    unreachable!()
                };
                let support = tensor.product(&mut self.model, &left, &right)?;
                odd_ext::persist(&mut self.model, support)
            }
            TypedForm::MixedTensorSwap(source, target, input) => {
                let Value::OddClass(terms) = self.eval(input)? else {
                    unreachable!()
                };
                let support = source.transpose(&mut self.model, target, &terms)?;
                odd_ext::persist(&mut self.model, support)
            }
            TypedForm::MixedComposition(homs, outer, inner) => {
                let Value::OddClass(outer) = self.eval(outer)? else {
                    unreachable!()
                };
                let Value::OddClass(inner) = self.eval(inner)? else {
                    unreachable!()
                };
                let [output, outer_hom, inner_hom] = homs;
                let value =
                    output.compose(&mut self.model, outer_hom, inner_hom, &outer, &inner)?;
                odd_ext::persist(&mut self.model, value)
            }
            TypedForm::MixedMinMap(data, operation, total, input) => {
                let Value::OddClass(support) = self.eval(input)? else {
                    unreachable!()
                };
                let image = data.apply(&mut self.model, &support, *operation, *total)?;
                odd_ext::persist(&mut self.model, image)
            }
            TypedForm::MixedDelta(data, order, input) => {
                let Value::OddClass(support) = self.eval(input)? else {
                    unreachable!()
                };
                let image = data.delta(&mut self.model, *order, &support)?;
                odd_ext::persist(&mut self.model, image)
            }
            TypedForm::MixedHShift(input) => self.eval_mixed_h_shift(input, expr.ty),
            TypedForm::MixedEmbed(input) => self.eval(input),
            TypedForm::SheafOperation(operation, data, degree, input) => {
                let value = self.eval(input)?;
                let (cochains, cohomology) = data.data();
                match (*operation, value) {
                    (SheafOp::Differential, Value::SheafC(raw)) => {
                        cochains.differential(&mut self.model, &raw)
                    }
                    (SheafOp::Class, Value::SheafC(raw)) => {
                        cohomology.class(&mut self.model, cochains, *degree, &raw)
                    }
                    (SheafOp::Representative, Value::OddClass(class)) => {
                        cohomology.representative(&mut self.model, cochains, *degree, &class)
                    }
                    (SheafOp::Scale(scalar), Value::SheafC(raw)) => {
                        let prime = cochains.prime();
                        sheaf::scale(&mut self.model, prime, sheaf::residue(prime, scalar), &raw)
                    }
                    (SheafOp::Scale(scalar), Value::OddClass(class)) => {
                        let prime = cochains.prime();
                        let scalar = sheaf::residue(prime, scalar) as usize;
                        odd_ext::scale(&mut self.model, prime, scalar, &class)
                    }
                    _ => unreachable!(),
                }
            }
            TypedForm::SheafLocalClass(data, grade, input) => {
                let Value::SheafC(raw) = self.eval(input)? else {
                    unreachable!()
                };
                data.class(&mut self.model, *grade, &raw)
            }
            TypedForm::SheafLocalOperation(data, grade, operation, input) => {
                let Value::OddClass(coordinates) = self.eval(input)? else {
                    unreachable!()
                };
                data.apply(&mut self.model, *grade, *operation, &coordinates)
            }
            TypedForm::SheafLerayIncl(leray, degree, input) => {
                let Value::SheafC(raw) = self.eval(input)? else {
                    unreachable!()
                };
                leray.include(&mut self.model, *degree, &raw)
            }
            TypedForm::SheafHigher(operation, signature, args, witnesses, target) => {
                let values = args
                    .iter()
                    .map(|arg| self.eval(arg))
                    .collect::<Result<Vec<_>>>()?;
                let expected = if let Some(target) = target {
                    let Value::OddClass(value) = self.eval(target)? else {
                        unreachable!()
                    };
                    Some(value)
                } else {
                    None
                };
                let witnesses = witnesses
                    .iter()
                    .map(|entry| {
                        Ok(sheaf_higher::Witness {
                            from: entry.from,
                            to: entry.to,
                            span: entry.span,
                            value: self.eval(&entry.expr)?,
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                sheaf_higher::evaluate(
                    &mut self.model,
                    signature,
                    *operation,
                    values,
                    witnesses,
                    expected.as_deref(),
                )
            }
            form @ (TypedForm::SheafDefOperation(..)
            | TypedForm::SheafDefAppend(..)
            | TypedForm::SheafDefGauge(..)) => self.eval_sheaf_deformation(form),
            TypedForm::SheafCompose(arguments) => {
                let [after, before] = arguments.as_ref();
                let after_value = self.eval(after)?;
                let before_value = self.eval(before)?;
                let Kind::SheafHom(_, a, m) = self.model.kind(after.ty) else {
                    unreachable!()
                };
                let Kind::SheafHom(_, b, n) = self.model.kind(before.ty) else {
                    unreachable!()
                };
                let (a, m, b, n) = (self.model.sheaf_hom(a), *m, self.model.sheaf_hom(b), *n);
                match (after_value, before_value) {
                    (Value::SheafC(x), Value::SheafC(y)) => {
                        let result = sheaf_ext::compose(&mut self.model, (&a, m, &x), (&b, n, &y))?;
                        sheaf::persist(&mut self.model, result)
                    }
                    (Value::OddClass(x), Value::OddClass(y)) => {
                        let x = a.cohomology.representative_support(
                            &mut self.model,
                            a.as_ref(),
                            m,
                            &x,
                        )?;
                        let y = b.cohomology.representative_support(
                            &mut self.model,
                            b.as_ref(),
                            n,
                            &y,
                        )?;
                        let result = sheaf_ext::compose(&mut self.model, (&a, m, &x), (&b, n, &y))?;
                        let Kind::SheafHom(_, ends, degree) = self.model.kind(expr.ty) else {
                            unreachable!()
                        };
                        let (data, degree) = (self.model.sheaf_hom(ends), *degree);
                        data.cohomology
                            .class(&mut self.model, data.as_ref(), degree, &result)
                    }
                    _ => unreachable!(),
                }
            }
            TypedForm::SheafApply(map, degree, input) => match self.eval(input)? {
                Value::SheafC(raw) => map.apply_cochain(&mut self.model, &raw),
                Value::OddClass(class) => map.apply_class(&mut self.model, *degree, &class),
                _ => unreachable!(),
            },
            TypedForm::SheafEmbed(data, degree, input) => {
                let Value::SheafC(raw) = self.eval(input)? else {
                    unreachable!()
                };
                data.embed(&mut self.model, *degree, &raw)
            }
            TypedForm::SheafPart(data, degree, input) => {
                let Value::SheafC(total) = self.eval(input)? else {
                    unreachable!()
                };
                data.part(&mut self.model, *degree, &total)
            }
            TypedForm::PScale(prime, scalar, value) => match self.eval(value)? {
                Value::PCobar(terms) => odd_cobar::scale(&mut self.model, *prime, *scalar, &terms),
                Value::OddClass(class) => odd_ext::scale(&mut self.model, *prime, *scalar, &class),
                _ => unreachable!(),
            },
            TypedForm::Bracket(inputs) => {
                let Kind::Arrow(target, _) = *self.model.kind(expr.ty) else {
                    unreachable!()
                };
                let inputs = inputs
                    .iter()
                    .map(|term| {
                        let Value::Class(coordinates) = self.eval(term)? else {
                            unreachable!()
                        };
                        Ok((term.ty, coordinates))
                    })
                    .collect::<Result<Vec<_>>>()?;
                bracket::evaluate(&mut self.model, target, &inputs)
            }
            TypedForm::Transfer(ScalarField::ModTwo(profile), operation, inputs) => {
                let inputs = self.scalar_arguments(inputs)?;
                transfer::evaluate(&mut self.model, *profile, &inputs, *operation)
            }
            TypedForm::Transfer(ScalarField::Odd(prime, profile), operation, inputs) => {
                let inputs = self.odd_arguments(inputs)?;
                odd_transfer::evaluate(&mut self.model, *prime, *profile, &inputs, *operation)
            }
            TypedForm::RealTransfer(operation, inputs) => {
                let inputs = inputs
                    .iter()
                    .map(|term| {
                        let Value::Class(coordinates) = self.eval(term)? else {
                            unreachable!()
                        };
                        Ok((term.ty, coordinates))
                    })
                    .collect::<Result<Vec<_>>>()?;
                real_transfer::evaluate(&mut self.model, *operation, expr.ty, &inputs)
            }
            TypedForm::HigherRestriction(ScalarField::ModTwo(source), target, inputs) => {
                let inputs = self.scalar_arguments(inputs)?;
                restriction_transfer::evaluate(&mut self.model, *source, *target, &inputs)
            }
            TypedForm::HigherRestriction(ScalarField::Odd(prime, source), target, inputs) => {
                let inputs = self.odd_arguments(inputs)?;
                odd_restriction_transfer::evaluate(
                    &mut self.model,
                    *prime,
                    *source,
                    *target,
                    &inputs,
                )
            }
            TypedForm::MapTransfer(operation, inputs) => {
                let inputs = inputs
                    .iter()
                    .map(|term| {
                        let Value::Class(coordinates) = self.eval(term)? else {
                            unreachable!()
                        };
                        Ok((term.ty, coordinates))
                    })
                    .collect::<Result<Vec<_>>>()?;
                map_transfer::evaluate(&mut self.model, &inputs, *operation)
            }
            TypedForm::RvTTransfer(operation, terms) => {
                let inputs = terms
                    .iter()
                    .map(|term| {
                        let Value::Class(class) = self.eval(term)? else {
                            unreachable!();
                        };
                        Ok((term.ty, class))
                    })
                    .collect::<Result<Vec<_>>>()?;
                rv_t_transfer::evaluate(&mut self.model, *operation, expr.ty, &inputs)
            }
            TypedForm::RvTSystem(operation, arity, entries) => {
                let entries = entries
                    .iter()
                    .map(|entry| {
                        Ok(rv_t_system::Entry {
                            from: entry.from,
                            to: entry.to,
                            ty: entry.expr.ty,
                            span: entry.span,
                            value: self.eval(&entry.expr)?,
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                rv_t_system::evaluate(&mut self.model, *arity, *operation, expr.ty, &entries)
            }
            TypedForm::MapSystem(operation, arity, entries) => {
                let entries = entries
                    .iter()
                    .map(|entry| {
                        Ok(map_system::Entry {
                            from: entry.from,
                            to: entry.to,
                            ty: entry.expr.ty,
                            span: entry.span,
                            value: self.eval(&entry.expr)?,
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                map_system::evaluate(&mut self.model, *arity, &entries, expr.ty, *operation)
            }
            TypedForm::OddSystem(operation, arity, entries) => {
                let entries = entries
                    .iter()
                    .map(|entry| {
                        Ok(odd_massey::Entry {
                            from: entry.from,
                            to: entry.to,
                            ty: entry.expr.ty,
                            span: entry.span,
                            value: self.eval(&entry.expr)?,
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                odd_massey::evaluate(&mut self.model, *arity, &entries, expr.ty, *operation)
            }
            TypedForm::ClassBracket {
                family,
                target,
                inputs,
                query,
            } => {
                let mut eval_class = |term: &TypedExpr| -> Result<_> {
                    let Value::Class(coordinates) = self.eval(term)? else {
                        unreachable!()
                    };
                    Ok((term.ty, coordinates))
                };
                let inputs = inputs
                    .iter()
                    .map(&mut eval_class)
                    .collect::<Result<Vec<_>>>()?;
                let query = query
                    .as_deref()
                    .map(eval_class)
                    .transpose()?
                    .map(|(_, class)| class);
                match (family, query) {
                    (ClassBracketFamily::Map, Some(query)) => {
                        map_bracket::member(&mut self.model, *target, &inputs, &query)
                    }
                    (ClassBracketFamily::Map, None) => {
                        map_bracket::predicate(&mut self.model, *target, &inputs)
                    }
                    (ClassBracketFamily::Koszul, Some(query)) => rv_t_bracket::member(
                        &mut self.model,
                        *target,
                        inputs.as_slice().try_into().unwrap(),
                        &query,
                    ),
                    (ClassBracketFamily::Koszul, None) => rv_t_bracket::predicate(
                        &mut self.model,
                        *target,
                        inputs.as_slice().try_into().unwrap(),
                    ),
                }
            }
            TypedForm::DefiningSystem {
                arity,
                layout,
                cochains,
            } => {
                let Kind::CyclicH2(modulus) = *self.model.kind(expr.ty) else {
                    unreachable!()
                };
                let cochains = cochains
                    .iter()
                    .map(|cochain| {
                        let Value::Bar(support) = self.eval(cochain)? else {
                            unreachable!()
                        };
                        Ok(cyclic::Cochain {
                            support,
                            span: cochain.span,
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                cyclic::evaluate(&mut self.model, modulus, *arity, *layout, cochains)
            }
            TypedForm::Lambda(body) => {
                let Kind::Arrow(domain, range) = *self.model.kind(expr.ty) else {
                    unreachable!()
                };
                let inputs = self.model.elements(domain)?;
                self.model.allocate(1 + inputs.len())?;
                let mut rows = Vec::with_capacity(inputs.len());
                for input in inputs.iter() {
                    self.locals.push(input.clone());
                    let result = self.eval(body);
                    self.locals.pop();
                    rows.push(result?);
                }
                self.model.check_map(domain, range, &rows)?;
                Ok(Value::Map(rows.into()))
            }
            TypedForm::Call(function, argument) => {
                let function = self.eval(function)?;
                let value = self.eval(argument)?;
                let Value::Map(rows) = function else {
                    unreachable!()
                };
                let index = self.model.index(argument.ty, &value)?;
                Ok(rows[index].clone())
            }
            TypedForm::Pair(left, right) => {
                self.model.allocate(2)?;
                let left = self.eval(left)?;
                let right = self.eval(right)?;
                Ok(Value::Pair(Rc::new((left, right))))
            }
            TypedForm::Fst(pair) | TypedForm::Snd(pair) => {
                let Value::Pair(pair) = self.eval(pair)? else {
                    unreachable!()
                };
                if matches!(expr.kind, TypedForm::Fst(_)) {
                    Ok(pair.0.clone())
                } else {
                    Ok(pair.1.clone())
                }
            }
            TypedForm::Table(entries) => {
                let Kind::Arrow(domain, range) = *self.model.kind(expr.ty) else {
                    unreachable!()
                };
                self.model.allocate(1 + entries.len())?;
                let rows = entries
                    .iter()
                    .map(|entry| self.eval(entry))
                    .collect::<Result<Vec<_>>>()?;
                self.model.check_map(domain, range, &rows)?;
                Ok(Value::Map(rows.into()))
            }
            TypedForm::Restrict(base) => {
                let Kind::Arrow(star, _) = *self.model.kind(expr.ty) else {
                    unreachable!()
                };
                let Value::Map(base) = self.eval(base)? else {
                    unreachable!()
                };
                let points = self.model.parent_points(star);
                self.model.allocate(1 + points.len())?;
                self.model.spend(points.len())?;
                let rows = points.iter().map(|&point| base[point].clone()).collect();
                Ok(Value::Map(rows))
            }
            TypedForm::Sew(charts) => self.sew(expr.ty, charts),
        }
    }

    fn scalar_arguments(&mut self, terms: &[TypedExpr]) -> Result<Vec<transfer::Input>> {
        terms
            .iter()
            .map(|term| {
                let (_, s, t) = self.model.kind(term.ty).ext_grade().unwrap();
                let Value::Class(coordinates) = self.eval(term)? else {
                    unreachable!()
                };
                Ok(((s, t), coordinates))
            })
            .collect()
    }

    fn odd_arguments(&mut self, terms: &[TypedExpr]) -> Result<Vec<odd_transfer::Input>> {
        terms
            .iter()
            .map(|term| {
                let (_, _, s, t) = self.model.kind(term.ty).pext_grade().unwrap();
                let Value::OddClass(coordinates) = self.eval(term)? else {
                    unreachable!()
                };
                Ok(((s, t), coordinates))
            })
            .collect()
    }

    fn sew(&mut self, ty: TypeId, charts: &[(usize, TypedExpr)]) -> Result<Value> {
        let Kind::Arrow(domain, range) = *self.model.kind(ty) else {
            unreachable!()
        };
        let poset = self.model.poset(domain).unwrap();
        self.model.allocate(1 + poset.names.len())?;
        let mut sewn = vec![None::<(usize, Value)>; poset.names.len()];
        for (center, chart) in charts {
            let Value::Map(rows) = self.eval(chart)? else {
                unreachable!()
            };
            let Kind::Arrow(star, _) = *self.model.kind(chart.ty) else {
                unreachable!()
            };
            let points = self.model.parent_points(star);
            for (&point, image) in points.iter().zip(rows.iter()) {
                self.model.spend(1)?;
                if let Some((other_center, other_image)) = &sewn[point] {
                    if !self.model.equal(other_image, image)? {
                        let other_image = self.model.describe(range, other_image);
                        let image = self.model.describe(range, image);
                        return Err(Error::new(
                            chart.span,
                            format!(
                                "overlap conflict at {}.{}: chart '{}' gives {other_image}, chart '{}' gives {image}",
                                self.model.type_name(domain),
                                poset.names[point],
                                poset.names[*other_center],
                                poset.names[*center],
                            ),
                        ));
                    }
                } else {
                    sewn[point] = Some((*center, image.clone()));
                }
            }
        }

        // A star containing the lower endpoint of a comparison also contains
        // the upper endpoint, so the chart's monotonicity applies to both.
        let rows = sewn.into_iter().map(|row| row.unwrap().1).collect();
        Ok(Value::Map(rows))
    }
}

fn check_new_ideal(model: &Model, name: &Label) -> Result<()> {
    if model.rv_ideal_named(&name.text).is_some() {
        return Err(Error::new(
            name.span,
            format!("rank ideal '{}' is already defined", name.text),
        ));
    }
    Ok(())
}

struct Checker<'a> {
    model: &'a mut Model,
    bindings: &'a HashMap<String, Binding>,
    locals: Vec<(String, TypeId)>,
}

macro_rules! rank_arithmetic_expression {
    () => {
        ExprForm::HeckeBasis(..)
            | ExprForm::HeckeScalar(..)
            | ExprForm::HeckeNegate(..)
            | ExprForm::WzwBasis(..)
            | ExprForm::WzwSEntry(..)
            | ExprForm::WzwTEntry(..)
            | ExprForm::WzwUnary(..)
            | ExprForm::WzwScale(..)
            | ExprForm::LensTv(..)
            | ExprForm::WzwTv(..)
            | ExprForm::CycloRatio(..)
            | ExprForm::CycloRatioProduct(..)
            | ExprForm::CycloRatioGalois(..)
            | ExprForm::TlPath(..)
            | ExprForm::TlOperation(..)
            | ExprForm::TlCoordinate(..)
            | ExprForm::TlScale(..)
            | ExprForm::TlTrace(..)
            | ExprForm::TlJones(..)
            | ExprForm::TlColoredJones(..)
            | ExprForm::TlFramed(..)
            | ExprForm::TlSurgery(..)
            | ExprForm::WzwCoordinate(..)
            | ExprForm::CycloEmbed(..)
            | ExprForm::CycloGalois(..)
            | ExprForm::CycloRoot(..)
            | ExprForm::CycloProduct(..)
            | ExprForm::CycloScale(..)
            | ExprForm::CycloConjugate(..)
            | ExprForm::FusionHandle(..)
            | ExprForm::Verlinde(..)
            | ExprForm::FusionBasis(..)
            | ExprForm::FusionOne(..)
            | ExprForm::FusionProduct(..)
            | ExprForm::FusionScale(..)
            | ExprForm::FusionDual(..)
            | ExprForm::QuantumBasis(..)
            | ExprForm::ResSpace(..)
            | ExprForm::ResSpecialize(..)
            | ExprForm::ResSpecializeQuad(..)
            | ExprForm::ResFil(..)
            | ExprForm::ResOrthogonal(..)
            | ExprForm::ValMat(..)
            | ExprForm::ValMatQuad(..)
            | ExprForm::ValFil(..)
            | ExprForm::ValDual(..)
            | ExprForm::StableLine(..)
            | ExprForm::StableLineQuad(..)
            | ExprForm::StableProjection(..)
            | ExprForm::RVCTInitial(..)
            | ExprForm::RVCTQuery(..)
            | ExprForm::RVCTLiteral(..)
            | ExprForm::RVCTSection(..)
            | ExprForm::RVCTRead(..)
            | ExprForm::OreTerm(..)
            | ExprForm::OreProduct(..)
            | ExprForm::OreUnary(..)
            | ExprForm::OreApply(..)
            | ExprForm::QuantumReduce(..)
            | ExprForm::QuantumPSection(..)
            | ExprForm::QuantumCurvature(..)
            | ExprForm::QuantumNabla(..)
            | ExprForm::QuantumAnnihilator(..)
            | ExprForm::QuantumProduct(..)
            | ExprForm::QuantumUnary(..)
            | ExprForm::Hive(..)
            | ExprForm::SphericalBasis(..)
            | ExprForm::SphericalScalar(..)
            | ExprForm::SphericalProduct(..)
            | ExprForm::SatakeInverse(..)
            | ExprForm::SatakeCharacter(..)
            | ExprForm::SatakeTransform(..)
            | ExprForm::SatakeScalar(..)
            | ExprForm::SatakeProduct(..)
            | ExprForm::SphericalEmbedding(..)
            | ExprForm::HeckeProduct(..)
            | ExprForm::RVLaurentMonomial(..)
            | ExprForm::RVIwAdlv(..)
            | ExprForm::RVAdlv(..)
            | ExprForm::RVNewton(..)
            | ExprForm::RVScale(..)
            | ExprForm::RVField(..)
            | ExprForm::RVInverse(..)
            | ExprForm::RVHilbert90(..)
            | ExprForm::RVTwisted(..)
            | ExprForm::RVShintani(..)
            | ExprForm::RVSpecialize(..)
    };
}

macro_rules! bp_expression {
    () => {
        ExprForm::BpInteger(..)
            | ExprForm::BpVariable(..)
            | ExprForm::BpVPageGenerator(..)
            | ExprForm::BpVPageClass(..)
            | ExprForm::BpVPageOperation(..)
            | ExprForm::BpPageGenerator(..)
            | ExprForm::BpPageClass(..)
            | ExprForm::BpPageOperation(..)
            | ExprForm::BpModExtGen(..)
            | ExprForm::BpModExtClass(..)
            | ExprForm::BpModExtRep(..)
            | ExprForm::BpRaise(..)
            | ExprForm::BpTriple(..)
            | ExprForm::BpTripleBase(..)
            | ExprForm::BpHigher(..)
            | ExprForm::BpMassey(..)
            | ExprForm::BpHeightGenerator(..)
            | ExprForm::BpHeightOperation(..)
            | ExprForm::BpRelOperation(..)
            | ExprForm::BpRelGenerator(..)
            | ExprForm::BpExtGenerator(..)
            | ExprForm::BpExtOperation(..)
            | ExprForm::BpOperation(..)
    };
}

impl Checker<'_> {
    fn ty(&mut self, ty: ParsedType) -> Result<TypeId> {
        self.model.spend(1)?;
        match ty.kind {
            TypeForm::Name(name) => self
                .model
                .named(&name.text)
                .ok_or_else(|| Error::new(name.span, format!("unknown space '{}'", name.text))),
            TypeForm::SetNat(source, target) => {
                self.model.set_family(&source)?;
                self.model.set_family(&target)?;
                self.model
                    .set_nat_type(source.text.into(), target.text.into())
            }
            TypeForm::Star(base, point) => {
                let base = self.ty(*base)?;
                self.model.star(base, &point)
            }
            TypeForm::Cohomology(base, degree) => {
                let base = self.ty(*base)?;
                self.model.cohomology_type(base, degree)
            }
            TypeForm::Operations(degree) => self.model.operation_type(degree),
            TypeForm::Cobar(s, t) => self.model.cobar_type(s, t),
            TypeForm::MCobar(grade) => self.model.mcobar_type(grade),
            TypeForm::RMCobar(grade) => self.model.real_cobar_type(grade),
            TypeForm::RhoBock(grade) => self.model.rho_page_type(grade),
            TypeForm::RhoLoc(grade) => self.model.rho_loc_type(grade),
            TypeForm::RhoMin(grade) => self.model.rho_min_type(grade),
            TypeForm::MExt(grade) => self.model.mext_type(grade),
            TypeForm::SheafLocal(name, point, mode, degree) => {
                let data = self.model.sheaf_supports(&name)?;
                let point = data.point(self.model, &point)?;
                self.model.sheaf_local_type(
                    name.text.into(),
                    SupportGrade {
                        point,
                        mode,
                        degree,
                    },
                )
            }
            TypeForm::SheafPage(name, grade) => self.model.sheaf_page_type(&name, grade),
            TypeForm::SheafComplex(form, name, degree) => {
                self.model.sheaf_complex_type(form, &name, degree)
            }
            TypeForm::MixedHomology(name, b) => {
                let complex = self.model.mixed_id(&name)?;
                self.model.mixed_homology_type(mixed::Grade { complex, b })
            }
            TypeForm::MixedChain(name, b, total) => {
                let complex = self.model.mixed_id(&name)?;
                self.model
                    .mixed_chain_type(mixed::Grade { complex, b }, total)
            }
            TypeForm::MixedMinimal(name, b, total) => {
                let minimal = self.model.mixed_minimal_id(&name)?;
                self.model
                    .mixed_minimal_type(mixed_minimal::Grade { minimal, b }, total)
            }
            TypeForm::MixedPage(spec) => {
                let grade = self.mixed_grade(spec)?;
                self.model.mixed_page_type(grade)
            }
            TypeForm::Sheaf(form, name, degree) => self.model.sheaf_type(form, &name, degree),
            TypeForm::SheafDeformation(name, order) => {
                self.model.sheaf(&name)?;
                self.model
                    .known_sheaf_deformation_type(name.text.into(), order)
            }
            TypeForm::SheafHom(form, source, target, degree) => {
                self.model.sheaf_hom_type(form, &source, &target, degree)
            }
            TypeForm::RMExt(grade) => self.model.rmext_type(grade),
            TypeForm::RealRho(true, grade) => self.model.real_rho_ext_type(grade),
            TypeForm::RealRho(false, grade) => self.model.real_rho_cobar_type(grade),
            TypeForm::RMBar(grade) => self.model.rmbar_type(grade),
            TypeForm::RMTriple(grade) => self.model.rm_triple_type(grade),
            TypeForm::PExt(p, s, t) => {
                let prime = Prime::new(self.model, p)?;
                self.model.odd_ext_type(prime, s, t)
            }
            TypeForm::PExtAn(p, n, s, t) => {
                let prime = Prime::new(self.model, p)?;
                self.model.profile_odd_ext_type(prime, Profile::An(n), s, t)
            }
            TypeForm::PMay(p, r, s, t, f) => {
                let prime = Prime::new(self.model, p)?;
                self.model.pmay_type(prime, r, s, t, f)
            }
            TypeForm::PMayAn(p, n, r, s, t, f) => {
                let prime = Prime::new(self.model, p)?;
                self.model
                    .profile_pmay_type(prime, Profile::An(n), r, s, t, f)
            }
            TypeForm::PCeAn(p, n, r, s, t, f) => {
                let prime = Prime::new(self.model, p)?;
                self.model.pce_type(prime, n, r, s, t, f)
            }
            TypeForm::PCobar(p, s, t) => {
                let prime = Prime::new(self.model, p)?;
                self.model.pcobar_type(prime, s, t)
            }
            TypeForm::PCobarAn(p, n, s, t) => {
                let prime = Prime::new(self.model, p)?;
                self.model.profile_pcobar_type(prime, Profile::An(n), s, t)
            }
            TypeForm::CobarAn(n, s, t) => self.model.profile_cobar_type(Profile::An(n), s, t),
            TypeForm::Ext(s, t) => self.model.ext_type(s, t),
            TypeForm::ExtAn(n, s, t) => self.model.profile_ext_type(Profile::An(n), s, t),
            TypeForm::XExt(base, s, t) => {
                let base = self.ty(*base)?;
                self.model.xext_type(base, s, t)
            }
            TypeForm::RVIwHecke(rank) => self.model.iw_hecke_type(rank),
            TypeForm::RVSphHecke(rank) => self.model.sph_hecke_type(rank),
            TypeForm::RVSatake(rank) => self.model.satake_type(rank),
            TypeForm::RVWzw(rank, level) => {
                self.model.wzw_type(rv_wzw::Grade::checked(rank, level)?)
            }
            TypeForm::RvTl(level, strands, height) => self
                .model
                .tl_type(rv_tl::Grade::checked(level, strands, height)?),
            TypeForm::RVResSpace(p, n) => {
                let prime = crate::rv_phodge::residue_prime(self.model, p)?;
                let grade = ResGrade::new(prime, n)?;
                self.model.rv_res_space_type(grade)
            }
            TypeForm::StableTree(count) => self.model.rv_stable_tree_type(count),
            TypeForm::StableMarks(p, count) => {
                let prime = crate::rv_phodge::residue_prime(self.model, p)?;
                self.model.rv_stable_curve_type(prime.get(), count)
            }
            TypeForm::RVValMat(p, n, r) => {
                let prime = crate::rv_phodge::residue_prime(self.model, p)?;
                let grade = ValGrade::new(prime.get(), n, r)?;
                self.model.rv_valuated_matroid_type(grade)
            }
            TypeForm::RVQpQuad(name) => self.model.qp_quad_type(&name),
            TypeForm::RVCycloQ(order) => self.model.cyclotomic_rational_type(order),
            TypeForm::RVCyclo(order) => self.model.cyclotomic_type(order),
            TypeForm::RVFusion(rank, level) => self
                .model
                .fusion_type(rv_fusion::Grade::checked(rank, level)?),
            TypeForm::RVOre => self.model.ore_type(),
            TypeForm::RVCTState(name, p, s) => {
                let grade = self.model.rv_ct_state_grade(name, p, s)?;
                self.model.rv_ct_state_type(grade)
            }
            TypeForm::RVCTResidue(p, s) => {
                let level = PeriodLevel::new(self.model, p, s)?;
                self.model.rv_ct_residue_type(level)
            }
            TypeForm::RVQGrP(prime, rank, ambient) => {
                let grade = rv_pcurv::ModGrade::checked(self.model, prime, rank, ambient)?;
                self.model.quantum_mod_type(grade)
            }
            TypeForm::RVQGr(rank, ambient) => self
                .model
                .quantum_grassmannian_type(rv_qgr::Grade::checked(rank, ambient)?),
            TypeForm::RVHive(weights) => {
                let boundary = rv_hive::Boundary::checked(self.model, weights)?;
                self.model.hive_type(boundary)
            }
            TypeForm::RVLaurent(spec) => {
                let grade = self.rv_grade(spec)?;
                self.model.rv_laurent_type(grade)
            }
            TypeForm::RVHom(spec) => {
                let grade = self.rv_grade(spec)?;
                self.model.rv_hom_type(grade)
            }
            TypeForm::JetPage(spec) => {
                let grade = self.jet_page_grade(spec)?;
                self.model.jet_page_type(grade)
            }
            TypeForm::JetHom(spec) => {
                let grade = self.jet_grade(spec)?;
                self.model.jet_hom_type(grade)
            }
            TypeForm::QKoszul(spec) => {
                let grade = self.q_grade(spec)?;
                self.model.q_koszul_type(grade)
            }
            TypeForm::QExt(spec) => {
                let grade = self.q_grade(spec)?;
                self.model.q_ext_type(grade)
            }
            TypeForm::VKoszul(spec) => {
                let grade = self.q_grade(spec)?;
                self.model.v_koszul_type(grade)
            }
            TypeForm::VExt(spec) => {
                let grade = self.q_grade(spec)?;
                self.model.v_ext_type(grade)
            }
            TypeForm::VBock(spec) => self.vbock_type(spec),
            TypeForm::BPVPage(spec) => {
                let grade = self.bp_vpage_grade(spec)?;
                self.model.bp_vpage_type(grade)
            }
            TypeForm::BpPage(mode, spec) => {
                let grade = self.bp_page_grade(mode, spec)?;
                self.model.bp_page_type(grade)
            }
            TypeForm::BPHeight(spec) => {
                let grade = self.bp_height_grade(spec)?;
                self.model.bp_height_type(grade)
            }
            TypeForm::BPRelExt(spec, t) => {
                let shape = self.bp_relative_shape(spec)?;
                self.model.bp_rel_ext_type(bp_rel::Grade { shape, t })
            }
            TypeForm::BPRel(spec) => {
                let shape = self.bp_relative_shape(spec)?;
                self.model.bp_rel_type(shape)
            }
            TypeForm::BPTriple(crate::syntax::BpSpec(p, n, a), t) => {
                let shape = bp::Shape::new(self.model, p, n, a)?;
                self.model.bp_triple_type(bp_classes::Grade { shape, t })
            }
            TypeForm::BPExt(spec, t) => {
                let grade = self.bp_ext_grade(spec, t)?;
                self.model.bp_ext_type(grade)
            }
            TypeForm::BPModExt(spec) => {
                let grade = self.bp_mod_ext_grade(spec)?;
                self.model.bp_mod_ext_type(grade)
            }
            TypeForm::BPThick(crate::syntax::BpSpec(p, n, a), h, power) => {
                let shape = bp::Shape::new(self.model, p, n, a)?;
                if !(2..=shape.height + 1).contains(&h) {
                    return Err(Error::plain("BPThick requires 2 <= h <= n+1"));
                }
                self.model.bp_type(shape.quotient(h - 1)?.thick(power)?)
            }
            TypeForm::BP(spec, depth) => {
                let mut shape = bp::Shape::new(self.model, spec.0, spec.1, spec.2)?;
                if let Some(depth) = depth {
                    shape = shape.quotient(depth)?;
                }
                self.model.bp_type(shape)
            }
            TypeForm::RVUntruncatedCyclic(spec, mode, homology) => {
                let ordinary = self.rv_hh_grade(spec)?;
                self.model
                    .rv_cyclic_type(rv_cyclic_homology::Grade { ordinary, mode }, homology)
            }
            TypeForm::RVNegativeCyclic(spec, homology) => {
                let grade = self.rv_negative_grade(spec)?;
                self.model.rv_negative_type(grade, homology)
            }
            TypeForm::RVHochschild(spec, homology) => {
                let grade = self.rv_hh_grade(spec)?;
                if homology {
                    self.model.rv_hh_class_type(grade)
                } else {
                    self.model.rv_hh_type(grade)
                }
            }
            TypeForm::RVBarCobar(spec) => {
                let grade = self.rv_bar_grade(spec)?;
                self.model.rv_barcobar_type(grade)
            }
            TypeForm::RVBarDual(spec) => {
                let grade = self.rv_bar_grade(spec)?;
                self.model.rv_path_type(grade, false)
            }
            TypeForm::RVPathExt(spec) => {
                let grade = self.rv_bar_grade(spec)?;
                self.model.rv_path_type(grade, true)
            }
            TypeForm::RVBarE(spec) => {
                let grade = self.rv_bar_page_grade(spec)?;
                self.model.rv_bar_page_type(grade)
            }
            TypeForm::RVBar(spec) => {
                let grade = self.rv_bar_grade(spec)?;
                self.model.rv_bar_type(grade)
            }
            TypeForm::RVComplex(spec, homology) => {
                let grade = self.rv_koszul_grade(spec)?;
                self.model.rv_complex_type(grade, homology)
            }
            TypeForm::ASDeRham(cover, modulus) => self.model.rv_as_derham_type(&cover, modulus),
            TypeForm::EOType(genus) => self.model.eo_type(genus),
            TypeForm::BinaryScalar(modulus) => {
                let field = BinaryField::new(self.model, modulus)?;
                self.model.binary_scalar_type(field)
            }
            TypeForm::ASOddDiff(cover) => self.model.rv_as_odd_diff_type(&cover),
            TypeForm::ASTwoDiff(cover) => self.model.rv_as_two_diff_type(&cover),
            TypeForm::ASOddDeRham(cover) => self.model.rv_as_odd_derham_type(&cover),
            TypeForm::ASOddTate(cover, degree) => self.model.rv_as_odd_tate_type(&cover, degree),
            TypeForm::PrimeScalar(p) => {
                let prime = self.scalar_prime(p)?;
                self.model.prime_scalar_type(prime)
            }
            TypeForm::CurveDiff(family, cover, modulus) => {
                self.model.rv_curve_diff_type(family, &cover, modulus)
            }
            TypeForm::RVPoly(selected, degree) => {
                let selected = q_koszul::indices(self.model, selected)?;
                self.model.rv_poly_type(selected, degree)
            }
            TypeForm::MapMay(spec) => self.map_may_type(spec),
            TypeForm::MapExt(source, target, s, t) => {
                let source = self.ty(*source)?;
                let target = self.ty(*target)?;
                self.model.map_ext_type(source, target, s, t)
            }
            TypeForm::MapCobar(source, target, s, t) => {
                let source = self.ty(*source)?;
                let target = self.ty(*target)?;
                self.model.map_cobar_type(source, target, s, t)
            }
            TypeForm::May(r, s, t, f) => self.model.may_type(r, s, t, f),
            TypeForm::MayAn(n, r, s, t, f) => {
                self.model.profile_may_type(Profile::An(n), r, s, t, f)
            }
            TypeForm::CeAn(n, r, s, t, f) => self.model.ce_type(n, r, s, t, f),
            TypeForm::Bar(modulus) => self.model.cyclic_type(modulus, false),
            TypeForm::CyclicH2(modulus) => self.model.cyclic_type(modulus, true),
            TypeForm::Suspension(base) => {
                let base = self.ty(*base)?;
                self.model.suspension_type(base)
            }
            TypeForm::Stable(base, shift, degree) => {
                let base = self.ty(*base)?;
                self.model.stable_type(base, shift, degree)
            }
            TypeForm::Product(left, right) => {
                let left = self.ty(*left)?;
                let right = self.ty(*right)?;
                self.model.product(left, right)
            }
            TypeForm::Arrow(domain, range) => {
                let domain = self.ty(*domain)?;
                let range = self.ty(*range)?;
                self.model.arrow(domain, range)
            }
        }
        .map_err(|error| error.at(ty.span))
    }

    fn set_predicate(&self, ty: TypeId) -> Result<(Rc<str>, Rc<str>, Rc<set_logic::Omega>)> {
        let (source, range) = self.set_nat_ends(ty)?;
        let omega = self.model.set_omega(&range)?;
        Ok((source, range, omega))
    }

    fn set_logic(
        &mut self,
        operation: SetLogic,
        labels: Vec<Label>,
        inputs: Vec<ParsedExpr>,
    ) -> Result<(TypeId, TypedForm)> {
        let inputs = inputs
            .into_iter()
            .map(|input| self.expr(input))
            .collect::<Result<Vec<_>>>()?;
        let (domain, range, action) = match operation {
            SetLogic::Top | SetLogic::Bottom => {
                let domain = self.model.set_family(&labels[0])?;
                let omega = self.model.set_omega(&labels[1].text)?;
                if domain.base != omega.base {
                    return Err(Error::plain(
                        "truth classifier and family must use the same site",
                    ));
                }
                (
                    Rc::from(labels[0].text.as_str()),
                    Rc::from(labels[1].text.as_str()),
                    LogicAction::Constant(omega, domain, matches!(operation, SetLogic::Top)),
                )
            }
            SetLogic::And | SetLogic::Or | SetLogic::Implies | SetLogic::Not => {
                let (source, target, omega) = self.set_predicate(inputs[0].ty)?;
                if inputs
                    .get(1)
                    .is_some_and(|second| second.ty != inputs[0].ty)
                {
                    return Err(Error::plain("logical operands must be predicates on the same family with the same classifier"));
                }
                let domain = self
                    .model
                    .set_family_named(&source)
                    .expect("predicate domain");
                (
                    source,
                    target,
                    LogicAction::Connect(omega, domain, operation),
                )
            }
            SetLogic::Exists | SetLogic::Forall => {
                let (source, target) = self.set_nat_ends(inputs[0].ty)?;
                let (predicate, classifier, omega) = self.set_predicate(inputs[1].ty)?;
                if predicate != source {
                    return Err(Error::plain(
                        "quantified predicate must be on the natural map's source family",
                    ));
                }
                let a = self.model.set_family_named(&source).expect("map source");
                let b = self.model.set_family_named(&target).expect("map target");
                (
                    target,
                    classifier,
                    LogicAction::Quantify(omega, a, b, matches!(operation, SetLogic::Forall)),
                )
            }
            SetLogic::Classify => {
                let (source, target) = self.set_nat_ends(inputs[0].ty)?;
                let a = self.model.set_family_named(&source).expect("map source");
                let b = self.model.set_family_named(&target).expect("map target");
                let omega = self.model.set_omega(&labels[0].text)?;
                if omega.base != b.base {
                    return Err(Error::plain(
                        "truth classifier and family must use the same site",
                    ));
                }
                (
                    target,
                    Rc::from(labels[0].text.as_str()),
                    LogicAction::Classify(omega, a, b),
                )
            }
            SetLogic::Inclusion => {
                let name = Rc::from(labels[0].text.as_str());
                let subobject = self.model.set_subobject(&name)?;
                let ty = self
                    .model
                    .register_set_nat_type(name, subobject.family.clone(), true)?;
                return Ok((ty, TypedForm::Constant(subobject.inclusion())));
            }
        };
        let ty = self.model.register_set_nat_type(domain, range, true)?;
        Ok((ty, TypedForm::SetLogic(action, inputs)))
    }

    fn get_set_product(&self, name: &str) -> Result<Rc<dependent_sets::Product>> {
        match self.model.set_origin(name) {
            Some(dependent_sets::Origin::Product(origin)) => Ok(origin),
            _ => Err(Error::plain(format!(
                "'{name}' is not a dependent set product"
            ))),
        }
    }

    fn get_set_sum(&self, name: &str) -> Result<Rc<dependent_sets::Sum>> {
        match self.model.set_origin(name) {
            Some(dependent_sets::Origin::Sum(origin)) => Ok(origin),
            _ => Err(Error::plain(format!("'{name}' is not a dependent set sum"))),
        }
    }

    fn set_dependent(
        &mut self,
        operation: SetOperation,
        first: Label,
        second: Option<Label>,
        input: Option<Box<ParsedExpr>>,
    ) -> Result<(TypeId, TypedForm)> {
        let input = input
            .map(|expr| self.expr(*expr))
            .transpose()?
            .map(Box::new);
        let first: Rc<str> = first.text.into();
        let second = second.map(|label| Rc::<str>::from(label.text));
        let (source, target) = if let Some(operand) = &input {
            self.set_nat_ends(operand.ty)?
        } else {
            (Rc::from(""), Rc::from(""))
        };
        let (domain, range, action) = match operation {
            SetOperation::PullNatural => {
                let second = second.expect("parsed target pullback");
                let a = dependent_sets::pullback(self.model, &first)?;
                let b = dependent_sets::pullback(self.model, &second)?;
                if a.family != source || b.family != target || !a.same_map(&b) {
                    return Err(Error::plain(
                        "setPullNat requires pullbacks of the two families along the same map",
                    ));
                }
                let diagram = self
                    .model
                    .set_family_named(&source)
                    .expect("registered source");
                (first, second, DependentAction::PullNatural(a, diagram))
            }
            SetOperation::Evaluation => {
                let second = second.expect("parsed projection pullback");
                let product = self.get_set_product(&first)?;
                let lifted = dependent_sets::pullback(self.model, &second)?;
                product.projection.check(&lifted)?;
                if lifted.family != first {
                    return Err(Error::plain(
                        "setEval requires the pullback of the indicated product",
                    ));
                }
                (
                    second,
                    product.projection.body.clone(),
                    DependentAction::Evaluation(product),
                )
            }
            SetOperation::Transpose => {
                let product = self.get_set_product(&first)?;
                let lifted = dependent_sets::pullback(self.model, &source)?;
                product.projection.check(&lifted)?;
                if target != product.projection.body {
                    return Err(Error::plain(
                        "setTranspose requires a map to the product's body family",
                    ));
                }
                let a = self
                    .model
                    .set_family_named(&lifted.family)
                    .expect("registered source");
                let pulled = self
                    .model
                    .set_family_named(&source)
                    .expect("registered pullback");
                (
                    lifted.family.clone(),
                    first,
                    DependentAction::Transpose(product, a, pulled),
                )
            }
            SetOperation::Untranspose => {
                let product = self.get_set_product(&target)?;
                let lifted = dependent_sets::pullback(self.model, &first)?;
                product.projection.check(&lifted)?;
                if lifted.family != source {
                    return Err(Error::plain(
                        "setUntranspose requires the pullback of the map's source family",
                    ));
                }
                let a = self
                    .model
                    .set_family_named(&source)
                    .expect("registered source");
                let pulled = self
                    .model
                    .set_family_named(&first)
                    .expect("registered pullback");
                (
                    first,
                    product.projection.body.clone(),
                    DependentAction::Untranspose(product, a, pulled),
                )
            }
            SetOperation::SumLift => {
                let sum = self.get_set_sum(&first)?;
                let lifted = dependent_sets::pullback(self.model, &target)?;
                sum.projection.check(&lifted)?;
                if source != sum.projection.body {
                    return Err(Error::plain(
                        "setSumLift requires a map from the sum's body family",
                    ));
                }
                let body = self
                    .model
                    .set_family_named(&source)
                    .expect("registered body");
                (
                    first,
                    lifted.family.clone(),
                    DependentAction::SumLift(sum, body),
                )
            }
            SetOperation::SumRestrict => {
                let sum = self.get_set_sum(&source)?;
                let lifted = dependent_sets::pullback(self.model, &first)?;
                sum.projection.check(&lifted)?;
                if lifted.family != target {
                    return Err(Error::plain(
                        "setSumRestrict requires the pullback of the map's target family",
                    ));
                }
                let a = self
                    .model
                    .set_family_named(&source)
                    .expect("registered sum");
                let body = self
                    .model
                    .set_family_named(&sum.projection.body)
                    .expect("registered body");
                (
                    sum.projection.body.clone(),
                    first,
                    DependentAction::SumRestrict(sum, a, body),
                )
            }
        };
        let evaluation = matches!(operation, SetOperation::Evaluation);
        let ty = self
            .model
            .register_set_nat_type(domain, range, evaluation)?;
        Ok((ty, TypedForm::SetDependent(action, input)))
    }

    fn set_nat_ends(&self, ty: TypeId) -> Result<(Rc<str>, Rc<str>)> {
        if let Kind::SetNat(source, target) = self.model.kind(ty) {
            return Ok((source.clone(), target.clone()));
        }
        Err(Error::plain(format!(
            "expected a SetNat, found {}",
            self.model.type_name(ty)
        )))
    }

    fn expr(&mut self, expr: ParsedExpr) -> Result<TypedExpr> {
        let (ty, kind) = self
            .expr_inner(expr.kind, expr.span)
            .map_err(|error| error.at(expr.span))?;
        Ok(TypedExpr {
            span: expr.span,
            ty,
            kind,
        })
    }

    fn set_topology_expression(
        &mut self,
        operation: set_topology::Operation,
        name: Label,
        input: Option<Box<ParsedExpr>>,
    ) -> Result<(TypeId, TypedForm)> {
        use set_topology::Operation;
        if matches!(operation, Operation::Map) {
            let topology = self.model.set_topology(&name.text)?;
            let ty = self.model.register_set_nat_type(
                topology.omega_name.clone(),
                topology.omega_name.clone(),
                true,
            )?;
            return Ok((ty, TypedForm::Constant(topology.map())));
        }
        let reflection = self.model.set_reflection(&name.text)?;
        let source = reflection.source.clone();
        let reflected: Rc<str> = name.text.into();
        match operation {
            Operation::Map => unreachable!(),
            Operation::Unit => {
                let ty = self.model.register_set_nat_type(source, reflected, true)?;
                Ok((ty, TypedForm::Constant(reflection.unit())))
            }
            Operation::Inverse => {
                let family = self
                    .model
                    .set_family_named(&source)
                    .expect("reflection source");
                let inverse = reflection.inverse(self.model, &family)?;
                let ty = self.model.register_set_nat_type(reflected, source, true)?;
                Ok((ty, TypedForm::Constant(inverse)))
            }
            Operation::Extend => {
                let map = self.expr(*input.expect("map to extend"))?;
                let (from, target) = self.set_nat_ends(map.ty)?;
                if from != source {
                    return Err(Error::plain(
                        "setReflectExtend requires a map from the original family",
                    ));
                }
                let certificate = self
                    .model
                    .set_reflection_certificate(reflection.topology.clone(), target.clone())?
                    .ok_or_else(|| {
                        Error::plain("setReflectExtend target must be a sheaf for the topology")
                    })?;
                let a = self.model.set_family_named(&source).expect("map source");
                let b = self.model.set_family_named(&target).expect("map target");
                let ty = self.model.set_nat_type(reflected, target)?;
                Ok((
                    ty,
                    TypedForm::SetReflectExtend {
                        reflection,
                        source: a,
                        target: b,
                        certificate,
                        map: Box::new(map),
                    },
                ))
            }
        }
    }

    fn expr_inner(&mut self, expr: ExprForm, span: Span) -> Result<(TypeId, TypedForm)> {
        self.model.spend(1)?;
        // Keep the checker frames for distinct expression families separate.
        match expr {
            rank_arithmetic_expression!() => self.expr_rank_arithmetic(expr),
            bp_expression!() => self.expr_bp(expr, span),
            expr => self.expr_regular(expr, span),
        }
    }

    fn expr_regular(&mut self, expr: ExprForm, span: Span) -> Result<(TypeId, TypedForm)> {
        match expr {
            ExprForm::SetTopology(operation, name, input) => {
                self.set_topology_expression(operation, name, input)
            }
            ExprForm::SetJUnit(name, inverse) => {
                let reflection = self.model.set_sheafification(&name.text)?;
                let name: Rc<str> = name.text.into();
                let original = reflection.family.clone();
                let (source, target, value) = if inverse {
                    let certificate = self.model.set_sheaf_certificate(original.clone())?
                        .ok_or_else(|| Error::plain("setJInverse requires the original family to be a double-negation sheaf"))?;
                    (name, original, certificate.inverse())
                } else {
                    (original, name, reflection.unit())
                };
                let ty = self.model.register_set_nat_type(source, target, true)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SetJExtend(name, map) => {
                let reflection = self.model.set_sheafification(&name.text)?;
                let map = self.expr(*map)?;
                let (source, target) = self.set_nat_ends(map.ty)?;
                if source != reflection.family {
                    return Err(Error::plain(
                        "setJExtend requires a map from the original family",
                    ));
                }
                let certificate = self
                    .model
                    .set_sheaf_certificate(target.clone())?
                    .ok_or_else(|| {
                        Error::plain("setJExtend target must be a double-negation sheaf")
                    })?;
                let a = self.model.set_family_named(&source).expect("map source");
                let b = self.model.set_family_named(&target).expect("map target");
                let ty = self.model.set_nat_type(name.text.into(), target)?;
                Ok((
                    ty,
                    TypedForm::SetJExtend {
                        reflection,
                        source: a,
                        target: b,
                        certificate,
                        map: Box::new(map),
                    },
                ))
            }
            ExprForm::SetLogic(operation, labels, inputs) => {
                self.set_logic(operation, labels, inputs)
            }
            ExprForm::SetKan(operation, first, second, input) => {
                let input = input
                    .map(|expr| self.expr(*expr))
                    .transpose()?
                    .map(Box::new);
                let ends = input
                    .as_ref()
                    .map(|arg| self.set_nat_ends(arg.ty))
                    .transpose()?;
                let (ty, action) = operation.plan(
                    self.model,
                    first.text.into(),
                    second.map(|label| label.text.into()),
                    ends,
                )?;
                Ok((ty, TypedForm::SetKan(action, input)))
            }
            ExprForm::SetDependent(operation, first, second, input) => {
                self.set_dependent(operation, first, second, input)
            }
            ExprForm::SetProjection(name) => {
                let (domain, elements) = self.model.set_elements(&name)?;
                let diagram = self
                    .model
                    .set_family_named(&elements.family)
                    .expect("elements family");
                let ty = self.model.arrow(domain, diagram.base)?;
                Ok((ty, TypedForm::Constant(elements.projection())))
            }
            ExprForm::SetSection(name, index) => {
                let (range, elements) = self.model.set_elements(&name)?;
                let diagram = self
                    .model
                    .set_family_named(&elements.family)
                    .expect("elements family");
                let ty = self.model.arrow(diagram.base, range)?;
                Ok((ty, TypedForm::SetSection(elements, diagram, index)))
            }
            ExprForm::SetIdentity(name) => {
                let diagram = self.model.set_family(&name)?;
                let family: Rc<str> = name.text.into();
                let ty = self.model.set_nat_type(family.clone(), family)?;
                Ok((ty, TypedForm::SetIdentity(diagram)))
            }
            ExprForm::SetCompose(inputs) => {
                let [outer, inner] = *inputs;
                let outer = self.expr(outer)?;
                let inner = self.expr(inner)?;
                let (middle, target) = self.set_nat_ends(outer.ty)?;
                let (source, through) = self.set_nat_ends(inner.ty)?;
                if through != middle {
                    return Err(Error::plain("setComp requires matching middle families"));
                }
                let a = self
                    .model
                    .set_family_named(&source)
                    .expect("registered source");
                let b = self
                    .model
                    .set_family_named(&middle)
                    .expect("registered middle");
                let ty = self.model.set_nat_type(source, target)?;
                Ok((ty, TypedForm::SetCompose(a, b, Box::new([outer, inner]))))
            }
            ExprForm::SetTotal(source, target, input) => {
                let (domain, a) = self.model.set_elements(&source)?;
                let (range, b) = self.model.set_elements(&target)?;
                let input = self.expr(*input)?;
                let (source, target) = self.set_nat_ends(input.ty)?;
                if source != a.family || target != b.family {
                    return Err(Error::plain(
                        "setTotal requires the declared elements spaces of the map's two families",
                    ));
                }
                let target = self
                    .model
                    .set_family_named(&target)
                    .expect("registered target");
                let ty = self.model.arrow(domain, range)?;
                Ok((ty, TypedForm::SetTotal(a, target, Box::new(input))))
            }
            ExprForm::Name(name) => {
                if let Some((index, (_, ty))) = self
                    .locals
                    .iter()
                    .enumerate()
                    .rev()
                    .find(|(_, (local, _))| *local == name.text)
                {
                    return Ok((*ty, TypedForm::Local(index)));
                }
                let binding = self.bindings.get(&name.text).ok_or_else(|| {
                    Error::new(name.span, format!("unknown value '{}'", name.text))
                })?;
                Ok((binding.ty, TypedForm::Global(binding.index)))
            }
            ExprForm::Point(ty, point) => {
                let ty = self.ty(ty)?;
                let point = self.model.point(ty, &point)?;
                Ok((ty, TypedForm::Point(point)))
            }
            ExprForm::Class(base, degree, form) => {
                let base = self.ty(base)?;
                let ty = self.model.cohomology_type(base, degree)?;
                let value = match form {
                    ClassForm::Zero => cohomology::zero(),
                    ClassForm::Generator(index, at) => {
                        let rank = self.model.cohomology(ty)?.rank();
                        self.generator(ty, rank, index, at)?
                    }
                    ClassForm::Cocycle(chains) => {
                        cohomology::literal(self.model, ty, chains, span)?
                    }
                };
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::Zero(ty) => {
                let ty = self.ty(ty)?;
                let value = match self.model.kind(ty) {
                    Kind::Operations(_) => steenrod::zero(),
                    Kind::Cobar(..) | Kind::CobarAn(..) | Kind::MCobar(..) => cobar::zero(),
                    Kind::RMCobar(..) | Kind::RMRhoCobar(..) => real_motivic::zero(),
                    Kind::Sheaf(SheafForm::Cochain, ..) | Kind::SheafComplex(SheafForm::Cochain, ..) | Kind::SheafHom(SheafForm::Cochain, ..) =>
                        Value::SheafC([].into()),
                    Kind::SheafPage(..) | Kind::SheafLocal(..) => odd_ext::zero(),
                    Kind::MixedMinimum(..) => odd_ext::zero(),
                    Kind::MixedChain(..) | Kind::MixedPage(..) | Kind::MixedHomology(..) => odd_ext::zero(),
                    Kind::Sheaf(SheafForm::Class, ..) | Kind::SheafComplex(SheafForm::Class, ..) | Kind::SheafHom(SheafForm::Class, ..) =>
                        odd_ext::zero(),
                    Kind::RhoBock(..) | Kind::RMExt(..) | Kind::RMBar(..) | Kind::RMRhoExt(..) => cohomology::zero(),
                    Kind::RhoLoc(..) => rho_localized::zero(),
                    Kind::RhoMin(..) => rho_minimal::zero(),
                    Kind::PCobar(..) | Kind::PCobarAn(..) => odd_cobar::zero(),
                    Kind::PExt(..) | Kind::PExtAn(..) | Kind::PMay(..) | Kind::PMayAn(..) | Kind::PCeAn(..) | Kind::BPPage(..) | Kind::BPModExt(..) | Kind::BPVPage(..) => odd_ext::zero(),
                    Kind::MapCobar(..) => Value::MapCobar([].into()),
                    Kind::JetHom(..) => jet_hom::zero(),
                    Kind::JetPage(..) => jet_page::zero(),
                    Kind::QKoszul(..) | Kind::VKoszul(..) => Value::QKoszul([].into()),
                    Kind::RVBarE(..) => cohomology::zero(),
                    Kind::RvHochschild(..) => cohomology::zero(),
                    Kind::RvCyclicChain(..) => Value::RVNCChain([].into()),
                    Kind::RvCyclicHomology(..) => cohomology::zero(),
                    Kind::RVNCChain(..) => Value::RVNCChain([].into()),
                    Kind::RvNegativeCyclic(..) => cohomology::zero(),
                    Kind::RVHHChain(..) => Value::RVHHChain([].into()),
                    Kind::RVBarCobar(..) => Value::RVBarCobar([].into()),
                    Kind::RVPathExt(..) => cohomology::zero(),
                    Kind::RVBar(..) | Kind::RVBarDual(..) => Value::RVBar([].into()),
                    Kind::RVKoszul(..) => Value::RVKoszul([].into()),
                    Kind::RVTor(..) => cohomology::zero(),
                    Kind::RVPoly(..) => Value::RVPoly([].into()),
                    Kind::ASDeRham(..) => Value::ASDeRham([].into()),
                    Kind::BinaryScalar(..) => Value::BinaryScalar(0),
                    Kind::ASOddDiff(..) | Kind::ASTwoDiff(..) => odd_ext::zero(),
                    Kind::ASOddDeRham(..) | Kind::ASOddTate(..) => odd_ext::zero(),
                    Kind::PrimeScalar(..) => Value::PrimeScalar(0),
                    Kind::CurveDiff(..) => Value::CurveDiff([].into()),
                    Kind::BP(..) => bp::zero(),
                    Kind::BPRel(..) => bp_rel::zero(),
                    Kind::BPExt(..) | Kind::BPRelExt(..) | Kind::BPHeight(..) => bp_classes::zero(),
                    Kind::QExt(..) | Kind::VExt(..) | Kind::VBock(..) => cohomology::zero(),
                    Kind::RVHom(..) => rv_hom::zero(),
                    Kind::RVIwHecke(..) => rv_iwhecke::zero(),
                    Kind::RVSphHecke(..) => rv_sphhecke::zero(),
                    Kind::RVSatake(..) => rv_satake::zero(),
                    Kind::RVQGr(..) => rv_qgr::zero(),
                    Kind::RVQGrP(..) => rv_pcurv::zero(),
                    Kind::RVOre => rv_ore::zero(),
                    Kind::RVCTState(..) => Value::PeriodState([].into()),
                    Kind::RVCTResidue(..) => Value::PeriodResidue(0),
                    Kind::RVFusion(..) => rv_fusion::zero(),
                    Kind::RVCyclo(..) => rv_cyclotomic::zero(),
                    Kind::RVResSpace(grade) => ResSpace::zero(*grade).persist(self.model)?,
                    Kind::RVQpQuad(..) => Value::RVQpQuad(QuadElement::ZERO),
                    Kind::RVCycloQ(..) => rv_lens::zero(),
                    Kind::RvTl(..) => rv_tl::zero(),
                    Kind::RVWzw(..) => rv_wzw::zero(),
                    Kind::RVLaurent(..) => rv_laurent::zero(),
                    Kind::Cohomology(..) | Kind::Stable(..) | Kind::Ext(..) | Kind::ExtAn(..) | Kind::MExt(..) | Kind::XExt(..) | Kind::MapExt(..) | Kind::May(..) | Kind::MayAn(..) | Kind::MapMay(..) | Kind::CeAn(..) => {
                        cohomology::zero()
                    }
                    Kind::Bar(_) => Value::Bar([].into()),
                    Kind::CyclicH2(_) => Value::CyclicClass(false),
                    _ => {
                        return Err(Error::new(
                            span,
                            format!(
                                "zero requires a cohomology, Steenrod algebra, cobar, MapCobar, Ext, XExt, MapExt, May, MapMay, or bar cochain type, found {}",
                                self.model.type_name(ty),
                            ),
                        ))
                    }
                };
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::Square(degree) => {
                let ty = self.model.operation_type(degree)?;
                let square = steenrod::square(self.model, degree)?;
                Ok((ty, TypedForm::Constant(square)))
            }
            ExprForm::Milnor(exponents) => {
                let (degree, value) = milnor::literal(self.model, exponents)?;
                let ty = self.model.operation_type(degree)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::Coact(n, index, at, input) => {
                let input = self.expr(*input)?;
                let (base, q) = self.class_type(&input, "coact")?;
                let (s, t) = milnor::coact_bidegree(self.model, base, n, index, q)
                    .map_err(|error| error.at(at))?;
                let ty = self.model.cobar_type(s, t)?;
                Ok((ty, TypedForm::Coact(n, index, Box::new(input))))
            }
            ExprForm::CobarWord(factors) => {
                let ((s, t), value) = cobar::literal(self.model, factors)?;
                let ty = self.model.cobar_type(s, t)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::MCobarWord(coefficient, factors) => {
                let (grade, value) = motivic::literal(self.model, coefficient, factors)?;
                let ty = self.model.mcobar_type(grade)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RMCobarWord(rho, tau, factors) => {
                let (grade, value) = real_motivic::literal(self.model, rho, tau, factors)?;
                let ty = self.model.real_cobar_type(grade)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RhoLocWord(t, vertex, factors) => {
                let (grade, value) = rho_localized::literal(self.model, t, vertex, factors)?;
                let ty = self.model.rho_loc_type(grade)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RhoLocOperation(operation, input) => {
                let input = self.expr(*input)?;
                let ty = match operation {
                    RhoLocOp::Localize => {
                        let Kind::RMCobar(grade) = *self.model.kind(input.ty) else {
                            return Err(self.invalid_operand(
                                &input,
                                operation.name(),
                                "an RMCobar cochain",
                            ));
                        };
                        self.model.rho_loc_type(grade)?
                    }
                    RhoLocOp::Section(t) => {
                        let Kind::Cobar(s, k) = *self.model.kind(input.ty) else {
                            return Err(self.invalid_operand(
                                &input,
                                operation.name(),
                                "a classical Cobar cochain",
                            ));
                        };
                        self.model
                            .rho_loc_type(rho_localized::section_grade(s, t, k)?)?
                    }
                    RhoLocOp::Project => match *self.model.kind(input.ty) {
                        Kind::RhoLoc(grade) => self
                            .model
                            .cobar_type(grade.s, rho_localized::classical_degree(grade)?)?,
                        Kind::RMExt(grade) => self
                            .model
                            .ext_type(grade.s, rho_localized::classical_degree(grade)?)?,
                        _ => {
                            return Err(self.invalid_operand(
                                &input,
                                operation.name(),
                                "a RhoLoc cochain or RMExt class",
                            ))
                        }
                    },
                    RhoLocOp::Homotopy => {
                        let Kind::RhoLoc(grade) = *self.model.kind(input.ty) else {
                            return Err(self.invalid_operand(
                                &input,
                                operation.name(),
                                "a RhoLoc cochain",
                            ));
                        };
                        let previous = grade.previous().ok_or_else(|| {
                            Error::plain("rhoHomotopy requires positive cobar length")
                        })?;
                        self.model.rho_loc_type(previous)?
                    }
                };
                Ok((ty, TypedForm::RhoLocOperation(operation, Box::new(input))))
            }
            ExprForm::RhoGenerator(grade, index, at) => {
                let ty = self.model.rho_page_type(grade)?;
                let rank = self.model.rho_page(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RhoOperation(operation, input) => {
                let input = self.expr(*input)?;
                let ty = match operation {
                    RhoOp::Class(r, f) => {
                        let Kind::RMCobar(chain) = *self.model.kind(input.ty) else {
                            return Err(self.invalid_operand(
                                &input,
                                operation.name(),
                                "an RMCobar cochain",
                            ));
                        };
                        self.model.rho_page_type(rho_bock::Grade { chain, r, f })?
                    }
                    _ => {
                        let Kind::RhoBock(grade) = *self.model.kind(input.ty) else {
                            return Err(self.invalid_operand(
                                &input,
                                operation.name(),
                                "a RhoBock class",
                            ));
                        };
                        match operation {
                            RhoOp::Representative => self.model.real_cobar_type(grade.chain)?,
                            RhoOp::Differential => {
                                self.model.rho_page_type(grade.differential()?)?
                            }
                            RhoOp::Next => self.model.rho_page_type(grade.next()?)?,
                            RhoOp::Previous => self.model.rho_page_type(grade.previous()?)?,
                            RhoOp::Class(..) => unreachable!(),
                        }
                    }
                };
                Ok((ty, TypedForm::RhoOperation(operation, Box::new(input))))
            }
            ExprForm::RhoMinOperation(op, input) => {
                let input = self.expr(*input)?;
                let (grade, ty) = self.rho_min_result(op, &input)?;
                Ok((ty, TypedForm::RhoMinOperation(op, grade, Box::new(input))))
            }
            ExprForm::RhoShift(power, input) => {
                let input = self.expr(*input)?;
                if let Kind::RhoLoc(grade) = *self.model.kind(input.ty) {
                    let ty = self
                        .model
                        .rho_loc_type(rho_localized::shift(grade, power)?)?;
                    return Ok((ty, TypedForm::MotivicCast(Box::new(input))));
                }
                if let Kind::RhoMin(grade) = *self.model.kind(input.ty) {
                    let power = usize::try_from(power).map_err(|_| {
                        Error::plain("rhoShift on RhoMin requires an unsigned exponent")
                    })?;
                    let ty = self
                        .model
                        .rho_min_type(real_motivic::shifted(grade, power)?)?;
                    return Ok((ty, TypedForm::RhoShift(power, Box::new(input))));
                }
                let truncated = match *self.model.kind(input.ty) {
                    Kind::RMRhoExt(grade) => Some((true, grade)),
                    Kind::RMRhoCobar(grade) => Some((false, grade)),
                    _ => None,
                };
                if let Some((is_class, grade)) = truncated {
                    let power = usize::try_from(power).map_err(|_| {
                        Error::plain("rhoShift on a rho truncation requires an unsigned exponent")
                    })?;
                    let shifted = real_rho::Grade {
                        chain: real_motivic::shifted(grade.chain, power)?,
                        ..grade
                    };
                    let ty = if is_class {
                        self.model.real_rho_ext_type(shifted)?
                    } else {
                        self.model.real_rho_cobar_type(shifted)?
                    };
                    return Ok((ty, TypedForm::RhoShift(power, Box::new(input))));
                }
                let grade = self.real_grade(&input, "rhoShift")?;
                let power = usize::try_from(power).map_err(|_| {
                    let label = if matches!(self.model.kind(input.ty), Kind::RMExt(..)) {
                        "RMExt"
                    } else {
                        "RMCobar"
                    };
                    Error::plain(format!("rhoShift on {label} requires an unsigned exponent"))
                })?;
                let shifted = real_motivic::shifted(grade, power)?;
                let ty = if matches!(self.model.kind(input.ty), Kind::RMExt(..)) {
                    self.model.rmext_type(shifted)?
                } else {
                    self.model.real_cobar_type(shifted)?
                };
                Ok((ty, TypedForm::RhoShift(power, Box::new(input))))
            }
            ExprForm::RealComplex(input) => {
                let input = self.expr(*input)?;
                let ty = match *self.model.kind(input.ty) {
                    Kind::RMCobar(grade) => self.model.mcobar_type(grade)?,
                    Kind::RMExt(grade) => self.model.mext_type(grade)?,
                    Kind::RMRhoCobar(grade) => self.model.mcobar_type(grade.chain)?,
                    Kind::RMRhoExt(grade) => self.model.mext_type(grade.chain)?,
                    _ => {
                        return Err(self.invalid_operand(
                            &input,
                            "realComplex",
                            "an RMCobar cochain or RMExt class, or a rho truncation of either",
                        ))
                    }
                };
                Ok((ty, TypedForm::RealComplex(Box::new(input))))
            }
            ExprForm::RealRhoOperation(op, input) => {
                let input = self.expr(*input)?;
                let (ty, grade) = self.real_rho_result(op, &input)?;
                Ok((ty, TypedForm::RealRhoOperation(op, grade, Box::new(input))))
            }
            ExprForm::RhoNull(power, input) => {
                let input = self.expr(*input)?;
                let grade = self.real_grade(&input, "rhoNull")?;
                let previous = grade
                    .previous()
                    .ok_or_else(|| Error::plain("rhoNull requires positive cobar length"))?;
                let ty = self
                    .model
                    .real_cobar_type(real_motivic::shifted(previous, power)?)?;
                Ok((ty, TypedForm::RhoNull(power, Box::new(input))))
            }
            ExprForm::RealMassey(coset, inputs) => {
                let name = if coset { "rmTriple" } else { "rmMassey" };
                let [a, b, c] = *inputs;
                let inputs = [self.expr(a)?, self.expr(b)?, self.expr(c)?];
                let [a, b, c] = inputs.each_ref();
                let layout = real_massey::Layout::new([
                    self.real_grade(a, name)?,
                    self.real_grade(b, name)?,
                    self.real_grade(c, name)?,
                ])?;
                let ty = if coset {
                    self.model.rm_triple_type(layout.target)?
                } else {
                    self.model.real_cobar_type(layout.target)?
                };
                Ok((ty, TypedForm::RealMassey(coset, layout, Box::new(inputs))))
            }
            ExprForm::RealMasseyNull(inputs) => {
                let [left, right] = *inputs;
                let left = self.expr(left)?;
                let right = self.expr(right)?;
                let grade = self
                    .real_grade(&left, "rmMasseyNull")?
                    .product(self.real_grade(&right, "rmMasseyNull")?)?;
                let previous = grade.previous().ok_or_else(|| {
                    Error::new(span, "rmMasseyNull requires positive total cobar length")
                })?;
                let ty = self.model.real_cobar_type(previous)?;
                Ok((
                    ty,
                    TypedForm::RealMasseyNull(grade, Box::new([left, right])),
                ))
            }
            ExprForm::RealTripleClass(index, input) => {
                let input = self.expr(*input)?;
                let name = if index.is_some() {
                    "rmTripleDir"
                } else {
                    "rmTripleBase"
                };
                let grade = self.rm_triple_grade(&input, name)?;
                let ty = self.model.rmext_type(grade)?;
                Ok((ty, TypedForm::RealTripleClass(index, Box::new(input))))
            }
            ExprForm::TauShift(power, input) => {
                let input = self.expr(*input)?;
                let ty = match *self.model.kind(input.ty) {
                    Kind::MCobar(grade) => self.model.mcobar_type(grade.shift(power)?)?,
                    Kind::MExt(grade) => self.model.mext_type(grade.shift(power)?)?,
                    _ => {
                        return Err(self.invalid_operand(
                            &input,
                            "tauShift",
                            "an MCobar cochain or MExt class",
                        ));
                    }
                };
                Ok((ty, TypedForm::MotivicCast(Box::new(input))))
            }
            ExprForm::Realize(input) => {
                let input = self.expr(*input)?;
                let ty = match *self.model.kind(input.ty) {
                    Kind::MCobar(grade) => {
                        self.model.cobar_type(grade.s, grade.classical_degree()?)?
                    }
                    Kind::MExt(grade) => self.model.ext_type(grade.s, grade.classical_degree()?)?,
                    _ => {
                        return Err(self.invalid_operand(
                            &input,
                            "realize",
                            "an MCobar cochain or MExt class",
                        ));
                    }
                };
                Ok((ty, TypedForm::MotivicCast(Box::new(input))))
            }
            ExprForm::PCobarWord(p, factors) => {
                let prime = Prime::new(self.model, p)?;
                let ((s, t), value) = odd_cobar::literal(self.model, prime, factors)?;
                let ty = self.model.pcobar_type(prime, s, t)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RvHHToMixed(name, f, input) => {
                let input = self.expr(*input)?;
                let complex = self.model.mixed_id(&name)?;
                let origin = self.model.mixed_hh_origin(complex)?;
                let ty = match self.model.kind(input.ty) {
                    Kind::RVHHChain(grade) => {
                        origin.check_weight(grade)?;
                        if f.is_some() {
                            return Err(Error::plain(
                                "rvHHToMixed on a Hochschild chain takes only a mixed complex name",
                            ));
                        }
                        self.model.mixed_chain_type(
                            mixed::Grade {
                                complex,
                                b: grade.b,
                            },
                            false,
                        )?
                    }
                    Kind::RvHochschild(grade) => {
                        origin.check_weight(grade)?;
                        let f = f.unwrap_or(0);
                        let chain = mixed::Grade {
                            complex,
                            b: grade.b,
                        }
                        .shift(-2 * f as i128)?;
                        self.model
                            .mixed_page_type(mixed_page::Grade { chain, r: 1, f })?
                    }
                    _ => {
                        return Err(self.invalid_operand(
                            &input,
                            "rvHHToMixed",
                            "an RVHHChain or RVHH class",
                        ))
                    }
                };
                Ok((ty, TypedForm::RvHHToMixed(origin, Box::new(input))))
            }
            ExprForm::MixedUnit(name, source, target) => {
                let id = self.model.mixed_id(&name)?;
                let hom = self.model.mixed_hom(id)?;
                let x = self.model.mixed_complex(hom.source).index(&source)?;
                let y = self.model.mixed_complex(hom.target).index(&target)?;
                let index = hom.index(x, y);
                let b = self.model.mixed_complex(id).raw_degree(index);
                let ty = self
                    .model
                    .mixed_chain_type(mixed::Grade { complex: id, b }, false)?;
                let value = odd_ext::persist(self.model, vec![(index, 1)])?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::MixedIdentity(name) => {
                let id = self.model.mixed_id(&name)?;
                let entries = self.model.mixed_hom(id)?.identity(self.model)?;
                let ty = self
                    .model
                    .mixed_chain_type(mixed::Grade { complex: id, b: 0 }, false)?;
                let value = odd_ext::persist(self.model, entries)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::MixedConeMap(name, operation, input) => {
                self.check_mixed_cone_map(name, operation, *input)
            }
            ExprForm::MixedConeShear(name, homotopy, input) => {
                let homotopy = self.expr(*homotopy)?;
                let input = self.expr(*input)?;
                let (h, true) = self.mixed_chain(&homotopy, "mixConeShear")? else {
                    return Err(self.invalid_operand(
                        &homotopy,
                        "mixConeShear",
                        "a degree-one MixedU map",
                    ));
                };
                if h.b != 1 {
                    return Err(self.invalid_operand(
                        &homotopy,
                        "mixConeShear",
                        "a degree-one MixedU map",
                    ));
                }
                let (grade, true) = self.mixed_chain(&input, "mixConeShear")? else {
                    return Err(self.invalid_operand(
                        &input,
                        "mixConeShear",
                        "a MixedU cone chain",
                    ));
                };
                let source = self.model.mixed_cone(grade.complex)?;
                let id = self.model.mixed_id(&name)?;
                let target = self.model.mixed_cone(id)?;
                target.check_shear(&source, self.model, h.complex)?;
                let ty = self.model.mixed_chain_type(
                    mixed::Grade {
                        complex: id,
                        b: grade.b,
                    },
                    true,
                )?;
                Ok((
                    ty,
                    TypedForm::MixedConeShear(
                        target,
                        source,
                        h.complex,
                        Box::new(homotopy),
                        Box::new(input),
                    ),
                ))
            }
            ExprForm::MixedTensorProduct(name, left, right) => {
                let left = self.expr(*left)?;
                let right = self.expr(*right)?;
                let (a, total) = self.mixed_chain(&left, "mixTensor")?;
                let (b, right_total) = self.mixed_chain(&right, "mixTensor")?;
                if total != right_total {
                    return Err(Error::plain("mixTensor requires chains of the same form"));
                }
                let target = self.model.mixed_id(&name)?;
                let tensor = self.model.mixed_tensor(target)?;
                if tensor.left != a.complex || tensor.right != b.complex {
                    return Err(Error::plain(
                        "mixTensor requires the declared factor complexes in order",
                    ));
                }
                let grade = mixed::Grade {
                    complex: target,
                    b: a.b,
                }
                .shift(b.b.into())?;
                let ty = self.model.mixed_chain_type(grade, total)?;
                Ok((
                    ty,
                    TypedForm::MixedTensorProduct(tensor, Box::new(left), Box::new(right)),
                ))
            }
            ExprForm::MixedTensorSwap(name, input) => {
                let input = self.expr(*input)?;
                let (grade, total) = self.mixed_chain(&input, "mixSwap")?;
                let source = self.model.mixed_tensor(grade.complex)?;
                let target_id = self.model.mixed_id(&name)?;
                let target = self.model.mixed_tensor(target_id)?;
                if source.left != target.right || source.right != target.left {
                    return Err(Error::plain(
                        "mixSwap requires the declared tensor with factors reversed",
                    ));
                }
                let ty = self.model.mixed_chain_type(
                    mixed::Grade {
                        complex: target_id,
                        ..grade
                    },
                    total,
                )?;
                Ok((
                    ty,
                    TypedForm::MixedTensorSwap(source, target, Box::new(input)),
                ))
            }
            ExprForm::MixedTensorCertificate(name, i, j, kind) => {
                let complex = self.model.mixed_id(&name)?;
                let tensor = self.model.mixed_tensor(complex)?;
                let left = self.model.mixed_bars(tensor.left)?;
                let right = self.model.mixed_bars(tensor.right)?;
                let (b, support) =
                    tensor.certificate(self.model, kind, left.bar(i)?, right.bar(j)?)?;
                let ty = self
                    .model
                    .mixed_chain_type(mixed::Grade { complex, b }, true)?;
                let value = odd_ext::persist(self.model, support)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::MixedComposition(name, outer, inner) => {
                let outer = self.expr(*outer)?;
                let inner = self.expr(*inner)?;
                let (g, total) = self.mixed_chain(&outer, "mixComp")?;
                let (f, inner_total) = self.mixed_chain(&inner, "mixComp")?;
                if total != inner_total {
                    return Err(Error::plain("mixComp requires chains of the same form"));
                }
                let outer_hom = self.model.mixed_hom(g.complex)?;
                let inner_hom = self.model.mixed_hom(f.complex)?;
                let output = self.model.mixed_id(&name)?;
                let output_hom = self.model.mixed_hom(output)?;
                if inner_hom.target != outer_hom.source
                    || output_hom.source != inner_hom.source
                    || output_hom.target != outer_hom.target
                {
                    return Err(Error::plain("mixComp endpoint complexes do not match"));
                }
                let chain = mixed::Grade {
                    complex: output,
                    b: g.b,
                }
                .shift(f.b.into())?;
                let ty = self.model.mixed_chain_type(chain, total)?;
                Ok((
                    ty,
                    TypedForm::MixedComposition(
                        [output_hom, outer_hom, inner_hom],
                        Box::new(outer),
                        Box::new(inner),
                    ),
                ))
            }
            ExprForm::MixedLiteral(name, b, terms) => {
                let complex = self.model.mixed_id(&name)?;
                let terms =
                    self.model
                        .mixed_complex(complex)
                        .literal(self.model, b.into(), terms)?;
                let ty = self
                    .model
                    .mixed_chain_type(mixed::Grade { complex, b }, false)?;
                let value = odd_ext::persist(self.model, terms)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::MixedHGenerator(name, b, index, at) => {
                self.check_mixed_h_generator(name, b, index, at)
            }
            ExprForm::MixedBarTerm(name, index, null) => {
                let complex = self.model.mixed_id(&name)?;
                let bars = self.model.mixed_bars(complex)?;
                let bar = bars.bar(index)?;
                let (b, support) = if null {
                    let torsion = bar
                        .torsion
                        .as_ref()
                        .ok_or_else(|| Error::plain("mixBarNull requires a torsion bar"))?;
                    (torsion.null_degree, torsion.witness.as_ref())
                } else {
                    (bar.birth, bar.cycle.as_ref())
                };
                let ty = self
                    .model
                    .mixed_chain_type(mixed::Grade { complex, b }, true)?;
                let value = odd_ext::persist(self.model, support.to_vec())?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::MixedMinGenerator(name, b, index, at) => {
                let minimal = self.model.mixed_minimal_id(&name)?;
                let indices = self
                    .model
                    .mixed_minimal(minimal)
                    .basis(self.model, b, false)?;
                let &generator = indices
                    .get(index)
                    .ok_or_else(|| Error::new(at, "generator index out of range"))?;
                let ty = self
                    .model
                    .mixed_minimal_type(mixed_minimal::Grade { minimal, b }, false)?;
                let value = odd_ext::persist(self.model, vec![(generator, 1)])?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::MixedMinMap(name, operation, input) => {
                let input = self.expr(*input)?;
                let (grade, total) = self.mixed_chain(&input, operation.name())?;
                let minimal = self.model.mixed_minimal_id(&name)?;
                let data = self.model.mixed_minimal(minimal);
                if data.source != grade.complex {
                    return Err(Error::plain(
                        "minimal mixed model has a different source complex",
                    ));
                }
                let ty = match operation {
                    mixed_minimal::Operation::Projection => self.model.mixed_minimal_type(
                        mixed_minimal::Grade {
                            minimal,
                            b: grade.b,
                        },
                        total,
                    )?,
                    mixed_minimal::Operation::Homotopy => {
                        self.model.mixed_chain_type(grade.shift(1)?, total)?
                    }
                };
                Ok((
                    ty,
                    TypedForm::MixedMinMap(data, operation, total, Box::new(input)),
                ))
            }
            ExprForm::MixedDelta(order, input) => {
                let input = self.expr(*input)?;
                let degree = 2 * order as i128 - 1;
                match *self.model.kind(input.ty) {
                    Kind::MixedMinimum(grade, false) => {
                        let ty = self.model.mixed_minimal_type(grade.shift(degree)?, false)?;
                        let data = self.model.mixed_minimal(grade.minimal);
                        Ok((ty, TypedForm::MixedDelta(data, order, Box::new(input))))
                    }
                    Kind::MixedChain(grade, false) => {
                        let ty = self.model.mixed_chain_type(grade.shift(degree)?, false)?;
                        let data = self.model.mixed_complex(grade.complex);
                        Ok((ty, TypedForm::MixedRawDelta(data, order, Box::new(input))))
                    }
                    _ => Err(self.invalid_operand(
                        &input,
                        "mixDelta",
                        "a MixedMin class or a raw MixedChain",
                    )),
                }
            }
            ExprForm::MixedGenerator(spec, index, at) => {
                let grade = self.mixed_grade(spec)?;
                let ty = self.model.mixed_page_type(grade)?;
                let rank = self.model.mixed_page(ty)?.rank();
                let value = self.prime_generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::MixedEmbed(power, input) => {
                let input = self.expr(*input)?;
                let offset = -2 * power as i128;
                let ty = match *self.model.kind(input.ty) {
                    Kind::MixedHomology(grade) => {
                        let ty = self.model.mixed_homology_type(grade.shift(offset)?)?;
                        return Ok((ty, TypedForm::MixedHShift(Box::new(input))));
                    }
                    Kind::MixedChain(grade, _) => {
                        self.model.mixed_chain_type(grade.shift(offset)?, true)?
                    }
                    Kind::MixedMinimum(grade, _) => {
                        self.model.mixed_minimal_type(grade.shift(offset)?, true)?
                    }
                    _ => {
                        return Err(self.invalid_operand(
                            &input,
                            "mixU",
                            "a MixedChain, MixedU, MixedH, MixedMin, or MixedMinU value",
                        ))
                    }
                };
                Ok((ty, TypedForm::MixedEmbed(Box::new(input))))
            }
            ExprForm::MixedClass(r, f, input) => {
                let input = self.expr(*input)?;
                let Kind::MixedChain(chain, true) = *self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, "mixPageClass", "a MixedU chain"));
                };
                let ty = self
                    .model
                    .mixed_page_type(mixed_page::Grade { chain, r, f })?;
                Ok((ty, TypedForm::Cobar(CobarOp::MixClass, Box::new(input))))
            }
            ExprForm::SheafCellTerm(name, q, cell, index, coefficient) => {
                let SheafTotal::Cellular(data) = self.model.sheaf_total(&name)? else {
                    return Err(Error::new(
                        name.span,
                        "shCellTerm requires a declared cellular total complex",
                    ));
                };
                let (degree, value) = data.term(self.model, q, &cell, index, coefficient)?;
                let ty = self
                    .model
                    .sheaf_complex_type(SheafForm::Cochain, &name, degree)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SheafTerm(name, chain, index, coefficient) => {
                let data = self.model.sheaf(&name)?;
                let value = data.term(self.model, &chain, index, coefficient)?;
                let ty = self
                    .model
                    .sheaf_type(SheafForm::Cochain, &name, chain.len() - 1)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SheafGenerator(name, degree, index) => {
                let data = self.model.sheaf(&name)?;
                let value = data
                    .cohomology
                    .generator(self.model, data.as_ref(), degree, index)?;
                let ty = self.model.sheaf_type(SheafForm::Class, &name, degree)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SheafHyperGenerator(name, degree, index) => {
                let data = self.model.sheaf_total(&name)?;
                let value =
                    data.cohomology()
                        .generator(self.model, data.filtered(), degree, index)?;
                let ty = self
                    .model
                    .sheaf_complex_type(SheafForm::Class, &name, degree)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SheafLocalGenerator(name, point, mode, degree, index) => {
                let data = self.model.sheaf_supports(&name)?;
                let point = data.point(self.model, &point)?;
                let grade = SupportGrade {
                    point,
                    mode,
                    degree,
                };
                let rank = data.rank(self.model, grade)?;
                if index >= rank {
                    return Err(Error::new(
                        name.span,
                        format!("local sheaf generator {index} is out of range for rank {rank}"),
                    ));
                }
                let ty = self.model.sheaf_local_type(name.text.into(), grade)?;
                let value = odd_ext::persist(self.model, vec![(index, 1)])?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SheafLocalClass(name, point, mode, input) => {
                let data = self.model.sheaf_supports(&name)?;
                let point = data.point(self.model, &point)?;
                let input = self.expr(*input)?;
                let degree = match (data.input.form(), self.model.kind(input.ty)) {
                    (SheafInputForm::Sheaf, Kind::Sheaf(SheafForm::Cochain, source, degree))
                    | (
                        SheafInputForm::Complex,
                        Kind::SheafComplex(SheafForm::Cochain, source, degree),
                    ) if source.as_ref() == data.input.name() => *degree,
                    _ => {
                        return Err(Error::new(
                            input.span,
                            format!(
                                "shLocalClass requires an ambient cochain of '{}'",
                                data.input.name()
                            ),
                        ))
                    }
                };
                let grade = SupportGrade {
                    point,
                    mode,
                    degree,
                };
                let ty = self.model.sheaf_local_type(name.text.into(), grade)?;
                Ok((ty, TypedForm::SheafLocalClass(data, grade, Box::new(input))))
            }
            ExprForm::SheafLocalOperation(operation, input) => {
                let input = self.expr(*input)?;
                let Kind::SheafLocal(name, grade) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, operation.name(), "a ShLocal class"));
                };
                let name = name.clone();
                let grade = *grade;
                let data = self
                    .model
                    .sheaf_supports_named(&name)
                    .expect("registered support data");
                let ty = if let Some(output) = operation.output(grade)? {
                    self.model.sheaf_local_type(name, output)?
                } else {
                    data.cochain_type(self.model, grade.degree)?
                };
                Ok((
                    ty,
                    TypedForm::SheafLocalOperation(data, grade, operation, Box::new(input)),
                ))
            }
            ExprForm::SheafLerayIncl(name, input) => {
                let leray = self.model.leray_complex(&name)?;
                let input = self.expr(*input)?;
                let form = leray.source.form();
                let expected = match form {
                    SheafInputForm::Sheaf => "a ShC cochain",
                    SheafInputForm::Complex => "a ShTot cochain",
                };
                let (source, degree, description) = match (form, self.model.kind(input.ty)) {
                    (SheafInputForm::Sheaf, Kind::Sheaf(SheafForm::Cochain, name, degree)) => {
                        (name, *degree, "sheaf")
                    }
                    (
                        SheafInputForm::Complex,
                        Kind::SheafComplex(SheafForm::Cochain, name, degree),
                    ) => (name, *degree, "complex"),
                    _ => return Err(self.invalid_operand(&input, "shLerayIncl", expected)),
                };
                if source.as_ref() != leray.source.name() {
                    return Err(Error::new(
                        input.span,
                        format!(
                            "shLerayIncl[{}] requires source {description} '{}', found '{}'",
                            name.text,
                            leray.source.name(),
                            source,
                        ),
                    ));
                }
                let ty = self
                    .model
                    .sheaf_complex_type(SheafForm::Cochain, &name, degree)?;
                Ok((
                    ty,
                    TypedForm::SheafLerayIncl(leray, degree, Box::new(input)),
                ))
            }
            ExprForm::SheafHomTerm(source, target, chain, column, row, coefficient) => {
                let degree = chain
                    .len()
                    .checked_sub(1)
                    .ok_or_else(|| Error::plain("shHomTerm requires a nonempty strict chain"))?;
                let ty = self
                    .model
                    .sheaf_hom_type(SheafForm::Cochain, &source, &target, degree)?;
                let Kind::SheafHom(_, ends, _) = self.model.kind(ty) else {
                    unreachable!()
                };
                let data = self.model.sheaf_hom(ends);
                let value = data.term(self.model, &chain, column, row, coefficient)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SheafExtGenerator(source, target, degree, index) => {
                let ty = self
                    .model
                    .sheaf_hom_type(SheafForm::Class, &source, &target, degree)?;
                let Kind::SheafHom(_, ends, _) = self.model.kind(ty) else {
                    unreachable!()
                };
                let data = self.model.sheaf_hom(ends);
                let value = data
                    .cohomology
                    .generator(self.model, data.as_ref(), degree, index)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SheafHomIdentity(name) => {
                let ty = self
                    .model
                    .sheaf_hom_type(SheafForm::Cochain, &name, &name, 0)?;
                let data = self.model.sheaf(&name)?;
                let value = sheaf_ext::identity(self.model, &data)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SheafHomMap(name) => {
                let map = self.model.sheaf_map(&name)?;
                let ends = SheafEnds {
                    source: map.source.name.clone(),
                    target: map.target.name.clone(),
                };
                let ty = self
                    .model
                    .known_sheaf_hom_type(SheafForm::Cochain, ends, 0)?;
                let value = sheaf_ext::map_cochain(self.model, &map)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SheafHigher(operation, inputs, witnesses) => {
                self.sheaf_higher(operation, inputs, witnesses)
            }
            form @ (ExprForm::SheafDefStart(..)
            | ExprForm::SheafDefOperation(..)
            | ExprForm::SheafDefAppend(..)
            | ExprForm::SheafDefGauge(..)) => self.check_sheaf_deformation(form),
            ExprForm::SheafKuranishi(query) => self.sheaf_kuranishi_expression(*query),
            ExprForm::SheafCompose(arguments) => {
                let [after, before] = *arguments;
                let after = self.expr(after)?;
                let before = self.expr(before)?;
                let (left_form, a, m) = self.sheaf_hom_grade(&after)?;
                let (right_form, b, n) = self.sheaf_hom_grade(&before)?;
                if left_form != right_form || a.source != b.target {
                    return Err(Error::plain("shExtCompose requires matching middle sheaf and two cochains or two classes"));
                }
                let degree = n.checked_add(m).ok_or_else(|| {
                    Error::plain("derived sheaf Hom composition degree is out of range")
                })?;
                let ends = SheafEnds {
                    source: b.source.clone(),
                    target: a.target.clone(),
                };
                let ty = self.model.known_sheaf_hom_type(left_form, ends, degree)?;
                Ok((ty, TypedForm::SheafCompose(Box::new([after, before]))))
            }
            ExprForm::SheafApply(name, input) => {
                let map = self.model.sheaf_map(&name)?;
                let input = self.expr(*input)?;
                let Kind::Sheaf(form, source, degree) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(
                        &input,
                        "shApply",
                        "a ShC cochain or ShH class",
                    ));
                };
                if *source != map.source.name {
                    return Err(Error::new(
                        input.span,
                        format!(
                            "shApply[{}] requires source sheaf '{}', found '{}'",
                            name.text, map.source.name, source
                        ),
                    ));
                }
                let degree = *degree;
                let ty = self
                    .model
                    .known_sheaf_type(*form, map.target.name.clone(), degree)?;
                Ok((ty, TypedForm::SheafApply(map, degree, Box::new(input))))
            }
            ExprForm::SheafEmbed(name, q, input) => {
                let data = self.model.sheaf_complex(&name)?;
                let component = data.component(q)?;
                let input = self.expr(*input)?;
                let Kind::Sheaf(SheafForm::Cochain, source, p) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, "shTotEmbed", "a ShC cochain"));
                };
                if *source != component.name {
                    return Err(Error::new(
                        input.span,
                        format!(
                            "shTotEmbed[{}, {q}] requires sheaf '{}', found '{}'",
                            name.text, component.name, source
                        ),
                    ));
                }
                let degree = p
                    .checked_add(q)
                    .ok_or_else(|| Error::new(input.span, "sheaf total degree is out of range"))?;
                let ty = self
                    .model
                    .sheaf_complex_type(SheafForm::Cochain, &name, degree)?;
                Ok((ty, TypedForm::SheafEmbed(data, q, Box::new(input))))
            }
            ExprForm::SheafPart(q, input) => {
                let input = self.expr(*input)?;
                let Kind::SheafComplex(SheafForm::Cochain, name, degree) =
                    self.model.kind(input.ty)
                else {
                    return Err(self.invalid_operand(&input, "shTotPart", "a ShTot cochain"));
                };
                let SheafTotal::Roos(data) = self
                    .model
                    .sheaf_total_named(name)
                    .expect("registered sheaf complex")
                else {
                    return Err(Error::new(
                        input.span,
                        "shTotPart requires ordinary Roos cochains",
                    ));
                };
                let component = data.component(q)?;
                let p = degree.checked_sub(q).ok_or_else(|| {
                    Error::new(input.span, "shTotPart would have a negative Roos degree")
                })?;
                let ty =
                    self.model
                        .known_sheaf_type(SheafForm::Cochain, component.name.clone(), p)?;
                Ok((ty, TypedForm::SheafPart(data, q, Box::new(input))))
            }
            ExprForm::SheafPageGenerator(name, grade, index) => {
                let complex = self.model.sheaf_total(&name)?;
                let rank = complex.filtered().page(self.model, grade)?.rank();
                if index >= rank {
                    return Err(Error::new(
                        name.span,
                        format!("sheaf page generator {index} is out of range for rank {rank}"),
                    ));
                }
                let ty = self.model.sheaf_page_type(&name, grade)?;
                let value = odd_ext::persist(self.model, vec![(index, 1)])?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SheafPageLift(operation, r, p, input) => {
                let input = self.expr(*input)?;
                let Kind::SheafComplex(SheafForm::Cochain, name, total) = self.model.kind(input.ty)
                else {
                    return Err(self.invalid_operand(&input, operation.name(), "a ShTot cochain"));
                };
                let name = name.clone();
                let grade = SheafPageGrade::from_total(r, p, *total);
                let complex = self
                    .model
                    .sheaf_total_named(&name)
                    .expect("registered sheaf complex");
                let ty = match operation {
                    SheafPageLift::Class => self.model.known_sheaf_page_type(name, grade)?,
                    SheafPageLift::Correct => input.ty,
                };
                Ok((
                    ty,
                    TypedForm::SheafPageLift(operation, complex, grade, Box::new(input)),
                ))
            }
            ExprForm::SheafPageOperation(operation, input) => {
                let input = self.expr(*input)?;
                let Kind::SheafPage(name, grade) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, operation.name(), "a ShE class"));
                };
                let name = name.clone();
                let grade = *grade;
                let complex = self
                    .model
                    .sheaf_total_named(&name)
                    .expect("registered sheaf complex");
                let ty = match operation {
                    SheafPageOp::Differential => {
                        self.model.known_sheaf_page_type(name, grade.outgoing()?)?
                    }
                    SheafPageOp::Representative => self.model.known_sheaf_complex_type(
                        SheafForm::Cochain,
                        name,
                        grade.total()?,
                    )?,
                    SheafPageOp::Next => self.model.known_sheaf_page_type(name, grade.next()?)?,
                    SheafPageOp::Previous => {
                        self.model.known_sheaf_page_type(name, grade.previous()?)?
                    }
                };
                Ok((
                    ty,
                    TypedForm::SheafPageOperation(operation, complex, grade, Box::new(input)),
                ))
            }
            ExprForm::SheafOperation(operation, input) => {
                let input = self.expr(*input)?;
                if let Kind::SheafHom(form, ends, degree) = self.model.kind(input.ty).clone() {
                    let (output, grade) = match (operation, form) {
                        (SheafOp::Differential, SheafForm::Cochain) => (
                            form,
                            degree.checked_add(1).ok_or_else(|| {
                                Error::new(input.span, "sheaf differential degree is out of range")
                            })?,
                        ),
                        (SheafOp::Class, SheafForm::Cochain) => (SheafForm::Class, degree),
                        (SheafOp::Representative, SheafForm::Class) => (SheafForm::Cochain, degree),
                        (SheafOp::Scale(_), _) => (form, degree),
                        (SheafOp::Differential | SheafOp::Class, _) => {
                            return Err(self.invalid_operand(
                                &input,
                                operation.name(),
                                "a ShHomC cochain",
                            ));
                        }
                        (SheafOp::Representative, _) => {
                            return Err(self.invalid_operand(
                                &input,
                                operation.name(),
                                "a ShExt class",
                            ));
                        }
                    };
                    let data = SheafComputation::Hom(self.model.sheaf_hom(&ends));
                    let ty = self.model.known_sheaf_hom_type(output, ends, grade)?;
                    return Ok((
                        ty,
                        TypedForm::SheafOperation(operation, data, degree, Box::new(input)),
                    ));
                }
                if let (
                    SheafOp::Scale(scalar),
                    Kind::SheafPage(name, _) | Kind::SheafLocal(name, _),
                ) = (operation, self.model.kind(input.ty))
                {
                    let prime = self
                        .model
                        .sheaf_total_named(name)
                        .expect("registered sheaf complex")
                        .filtered()
                        .prime();
                    let scalar = sheaf::residue(prime, scalar) as usize;
                    return Ok((input.ty, TypedForm::PScale(prime, scalar, Box::new(input))));
                }
                let cochain_type = match self.model.kind(input.ty) {
                    Kind::SheafComplex(..) => "a ShTot cochain",
                    _ => "a ShC cochain or ShTot cochain",
                };
                let expected = match operation {
                    SheafOp::Differential | SheafOp::Class => cochain_type,
                    SheafOp::Representative => "a ShH class or ShHyper class",
                    SheafOp::Scale(_) => "a sheaf cochain or class",
                };
                let (form, name, degree, total) = match self.model.kind(input.ty) {
                    Kind::Sheaf(form, name, degree) => (*form, name.clone(), *degree, false),
                    Kind::SheafComplex(form, name, degree) => (*form, name.clone(), *degree, true),
                    _ => return Err(self.invalid_operand(&input, operation.name(), expected)),
                };
                let (output, grade) = match (operation, form) {
                    (SheafOp::Differential, SheafForm::Cochain) => (
                        form,
                        degree.checked_add(1).ok_or_else(|| {
                            Error::new(input.span, "sheaf differential degree is out of range")
                        })?,
                    ),
                    (SheafOp::Class, SheafForm::Cochain) => (SheafForm::Class, degree),
                    (SheafOp::Representative, SheafForm::Class) => (SheafForm::Cochain, degree),
                    (SheafOp::Scale(_), _) => (form, degree),
                    _ => return Err(self.invalid_operand(&input, operation.name(), expected)),
                };
                let (ty, data) = if total {
                    let data = self
                        .model
                        .sheaf_total_named(&name)
                        .expect("registered sheaf complex");
                    (
                        self.model.known_sheaf_complex_type(output, name, grade)?,
                        SheafComputation::Total(data),
                    )
                } else {
                    let data = self.model.sheaf_named(&name).expect("registered sheaf");
                    (
                        self.model.known_sheaf_type(output, name, grade)?,
                        SheafComputation::Ordinary(data),
                    )
                };
                Ok((
                    ty,
                    TypedForm::SheafOperation(operation, data, degree, Box::new(input)),
                ))
            }
            ExprForm::PScale(scalar, expression) => {
                let value = self.expr(*expression)?;
                let p = match *self.model.kind(value.ty) {
                    Kind::PCobar(p, ..)
                    | Kind::PCobarAn(p, ..)
                    | Kind::PExt(p, ..)
                    | Kind::PExtAn(p, ..)
                    | Kind::PMay(p, ..)
                    | Kind::PMayAn(p, ..)
                    | Kind::PCeAn(p, ..)
                    | Kind::ASOddDiff(_, p)
                    | Kind::ASTwoDiff(_, p)
                    | Kind::ASOddDeRham(_, p)
                    | Kind::ASOddTate(_, p, _) => p,
                    Kind::MixedChain(grade, _) | Kind::MixedHomology(grade) => self.model.mixed_complex(grade.complex).prime,
                    Kind::MixedMinimum(grade,_) => self.model.mixed_minimal(grade.minimal).prime(),
                    Kind::MixedPage(grade) => self.model.mixed_complex(grade.chain.complex).prime,
                    Kind::BPPage(grade) => Prime::validated(grade.shape.prime),
                    Kind::BPVPage(grade) => grade.prime(),
                    Kind::BPModExt(grade) => grade.prime(),
                    _ => return Err(self.invalid_operand(
                        &value,
                        "pscale",
                        "a PCobar cochain, PExt class, PMay class, PCEAn class, BP page class, mod-p BP Ext class, or mixed complex chain or page class",
                    )),
                };
                Ok((value.ty, TypedForm::PScale(p, scalar, Box::new(value))))
            }
            ExprForm::PExtGenerator(p, s, t, index, at) => {
                let prime = Prime::new(self.model, p)?;
                let ty = self.model.odd_ext_type(prime, s, t)?;
                let rank = self.model.odd_ext(ty)?.rank();
                let value = self.prime_generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::PExtAnGenerator(p, n, s, t, index, at) => {
                let prime = Prime::new(self.model, p)?;
                let ty = self
                    .model
                    .profile_odd_ext_type(prime, Profile::An(n), s, t)?;
                let rank = self.model.odd_ext(ty)?.rank();
                let value = self.prime_generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::PMayGenerator(p, r, s, t, f, index, at) => {
                let prime = Prime::new(self.model, p)?;
                let ty = self.model.pmay_type(prime, r, s, t, f)?;
                let rank = self.model.odd_page(ty)?.rank();
                let value = self.prime_generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::PMayAnGenerator(p, n, r, s, t, f, index, at) => {
                let prime = Prime::new(self.model, p)?;
                let ty = self
                    .model
                    .profile_pmay_type(prime, Profile::An(n), r, s, t, f)?;
                let rank = self.model.odd_page(ty)?.rank();
                let value = self.prime_generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::PCeAnGenerator(p, n, r, s, t, f, index, at) => {
                let prime = Prime::new(self.model, p)?;
                let ty = self.model.pce_type(prime, n, r, s, t, f)?;
                let rank = self.model.odd_page(ty)?.rank();
                let value = self.prime_generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RealBarTerm(grade, index, at) => {
                let group = real_extensions::bars(self.model, grade)?;
                group
                    .get(index)
                    .and_then(|bar| bar.exponent_at(grade.t))
                    .map_err(|error| error.at(at))?;
                let ty = self.model.rmbar_type(grade)?;
                self.model.allocate(1)?;
                Ok((ty, TypedForm::Constant(Value::Class([index].into()))))
            }
            ExprForm::RealExtension(op, input) => {
                let input = self.expr(*input)?;
                let grade = match op {
                    RealExtOp::Class | RealExtOp::Lift => {
                        let Kind::RMBar(grade) = *self.model.kind(input.ty) else {
                            return Err(self.invalid_operand(
                                &input,
                                op.name(),
                                "RMBar coordinates",
                            ));
                        };
                        grade
                    }
                    _ => self.rmext_grade(&input, op.name())?,
                };
                let ty = match op {
                    RealExtOp::Coordinates => self.model.rmbar_type(grade)?,
                    RealExtOp::Class => self.model.rmext_type(grade)?,
                    RealExtOp::Lift => self.model.real_cobar_type(grade)?,
                    RealExtOp::Divide(power) => self
                        .model
                        .rmext_type(real_extensions::quotient_grade(grade, power)?)?,
                    RealExtOp::Edge(r, f) => self
                        .model
                        .rho_page_type(real_extensions::edge_grade(grade, r, f)?)?,
                };
                Ok((ty, TypedForm::RealExtension(op, grade, Box::new(input))))
            }
            ExprForm::RhoProductNull(power, inputs) => {
                let [a, b] = *inputs;
                let a = self.expr(a)?;
                let b = self.expr(b)?;
                let grade = self
                    .rmext_grade(&a, "rhoProductNull")?
                    .product(self.rmext_grade(&b, "rhoProductNull")?)?;
                real_extensions::quotient_grade(grade, power)?;
                let previous = grade.previous().ok_or_else(|| {
                    Error::new(span, "rhoProductNull requires positive total cobar length")
                })?;
                let ty = self.model.real_cobar_type(previous)?;
                Ok((ty, TypedForm::RhoProductNull(power, Box::new([a, b]))))
            }
            ExprForm::RealBar(op, s, k, index, at) => {
                let group = self.model.real_bars(s, k)?;
                let bar = group.get(index).map_err(|error| error.at(at))?;
                let ty = match op {
                    RealBarOp::Generator => self.model.rmext_type(bar.grade(s, k)?)?,
                    RealBarOp::Representative => self.model.real_cobar_type(bar.grade(s, k)?)?,
                    RealBarOp::Witness => self.model.real_cobar_type(bar.witness_grade(s, k)?)?,
                };
                Ok((ty, TypedForm::RealBar(op, s, k, index)))
            }
            ExprForm::RMExtGenerator(grade, index, at) => {
                let ty = self.model.rmext_type(grade)?;
                let rank = self.model.real_ext(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RealRhoGenerator(grade, index, at) => {
                let ty = self.model.real_rho_ext_type(grade)?;
                let rank = self.model.real_ext(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::MExtGenerator(grade, index, at) => {
                let ty = self.model.mext_type(grade)?;
                let rank = self.model.adams(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::ExtGenerator(s, t, index, at) => {
                let ty = self.model.ext_type(s, t)?;
                let rank = self.model.adams(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::ExtAnGenerator(n, s, t, index, at) => {
                let ty = self.model.profile_ext_type(Profile::An(n), s, t)?;
                let rank = self.model.adams(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RestrictAlgebra(n, value) => {
                let value = self.expr(*value)?;
                let kind = self.model.kind(value.ty);
                let target = Profile::An(n);
                let ty = if let Some((source, r, s, t, f)) = kind.may_grade() {
                    require_restriction(source, target, span)?;
                    self.model.profile_may_type(target, r, s, t, f)?
                } else if let Some((prime, source, r, s, t, f)) = kind.pmay_grade() {
                    require_restriction(source, target, span)?;
                    self.model.profile_pmay_type(prime, target, r, s, t, f)?
                } else if let Some((source, s, t)) = kind.cobar_grade() {
                    require_restriction(source, target, span)?;
                    self.model.profile_cobar_type(target, s, t)?
                } else if let Some((source, s, t)) = kind.ext_grade() {
                    require_restriction(source, target, span)?;
                    self.model.profile_ext_type(target, s, t)?
                } else if let Some((p, source, s, t)) = kind.pcobar_grade() {
                    require_restriction(source, target, span)?;
                    self.model.profile_pcobar_type(p, target, s, t)?
                } else if let Some((p, source, s, t)) = kind.pext_grade() {
                    require_restriction(source, target, span)?;
                    self.model.profile_odd_ext_type(p, target, s, t)?
                } else {
                    return Err(self.invalid_operand(
                        &value,
                        "resA",
                        "a cobar cochain or Ext class, or a May class",
                    ));
                };
                Ok((ty, TypedForm::RestrictAlgebra(target, Box::new(value))))
            }
            ExprForm::CobarCup(i, left, right) => {
                let left = self.expr(*left)?;
                let right = self.expr(*right)?;
                let (profile, s, t) = self.profile_cobar_bidegree(&left, "cCup")?;
                let (other, u, v) = self.profile_cobar_bidegree(&right, "cCup")?;
                require_profile(profile, other, "cCup", right.span)?;
                let (length, degree) = cobar_squares::cup_bidegree(i, (s, t), (u, v))?;
                let ty = self.model.profile_cobar_type(profile, length, degree)?;
                Ok((
                    ty,
                    TypedForm::CobarCup(profile, i, Box::new(left), Box::new(right)),
                ))
            }
            ExprForm::ExtSquare(k, value) => {
                let value = self.expr(*value)?;
                let (profile, s, t) = self
                    .model
                    .kind(value.ty)
                    .ext_grade()
                    .ok_or_else(|| self.invalid_operand(&value, "extSq", "an Ext class"))?;
                let (length, degree) = cobar_squares::square_bidegree(k, s, t)?;
                let ty = self.model.profile_ext_type(profile, length, degree)?;
                Ok((ty, TypedForm::ExtSquare(k, Box::new(value))))
            }
            ExprForm::Nakamura(k, value) => {
                let value = self.expr(*value)?;
                let (profile, s, t) = self.profile_cobar_bidegree(&value, "cNakamura")?;
                let (s, t) = nakamura::bidegree(k, s, t)?;
                let ty = self.model.profile_cobar_type(profile, s, t)?;
                Ok((ty, TypedForm::Nakamura(k, Box::new(value))))
            }
            ExprForm::MayLiftTo(r, f, value) => {
                let value = self.expr(*value)?;
                let kind = self.model.kind(value.ty);
                if kind.cobar_grade().is_none() && !matches!(kind, Kind::MapCobar(..)) {
                    return Err(self.invalid_operand(
                        &value,
                        "mayLiftTo",
                        "a mod-two cobar cochain",
                    ));
                }
                Ok((value.ty, TypedForm::MayLiftTo(r, f, Box::new(value))))
            }
            ExprForm::XExtGenerator(base, s, t, index, at) => {
                let base = self.ty(base)?;
                let ty = self.model.xext_type(base, s, t)?;
                let rank = self.model.xext(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::MapExtGenerator(source, target, s, t, index, at) => {
                let source = self.ty(source)?;
                let target = self.ty(target)?;
                let ty = self.model.map_ext_type(source, target, s, t)?;
                let rank = self.model.map_ext(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::MapExtIdentity(base) => {
                let base = self.ty(base)?;
                let ty = self.model.map_ext_type(base, base, 0, 0)?;
                let value = map_ext::identity(self.model, ty)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::MapWord(source, target, input) => {
                let source_ty = self.ty(source.base)?;
                let target_ty = self.ty(target.base)?;
                let input = self.expr(*input)?;
                let (s, degree) = self.cobar_bidegree(&input, "mword")?;
                let t = map_cochain::internal_degree(source.degree, target.degree, degree)?;
                let source =
                    map_cochain::generator(self.model, source_ty, source.degree, source.index)
                        .map_err(|error| error.at(source.span))?;
                let target =
                    map_cochain::generator(self.model, target_ty, target.degree, target.index)
                        .map_err(|error| error.at(target.span))?;
                let ty = self.model.map_cobar_type(source_ty, target_ty, s, t)?;
                Ok((ty, TypedForm::MapWord(source, target, Box::new(input))))
            }
            ExprForm::JetTerm(literal) => {
                let source = self.q_endpoint(literal.source.base)?;
                let target = self.q_endpoint(literal.target.base)?;
                let path = self.model.jet_path_id(&literal.path)?;
                if literal.power >= self.model.jet_path(path).order {
                    return Err(Error::new(
                        span,
                        "jetterm exponent is outside the arc truncation",
                    ));
                }
                let term = jet_hom::Term {
                    power: literal.power,
                    source: q_module::generator(
                        self.model,
                        source,
                        literal.source.degree,
                        literal.source.index,
                    )
                    .map_err(|error| error.at(literal.source.span))?,
                    target: q_module::generator(
                        self.model,
                        target,
                        literal.target.degree,
                        literal.target.index,
                    )
                    .map_err(|error| error.at(literal.target.span))?,
                };
                let ty = self.model.jet_hom_type(jet_hom::Grade {
                    diagram: jet_hom::Diagram {
                        source,
                        target,
                        path,
                    },
                    parity: term.parity(),
                })?;
                let value = jet_hom::persist(self.model, vec![(term, 1)])?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::JetOperation(op, value) => {
                let value = self.expr(*value)?;
                let ty = self.jet_operation_type(op, &value)?;
                Ok((ty, TypedForm::JetOperation(op, Box::new(value))))
            }
            ExprForm::JetGenerator(spec, index, at) => {
                let grade = self.jet_page_grade(spec)?;
                let ty = self.model.jet_page_type(grade)?;
                let rank = self.model.jet_page(ty)?.rank();
                if index >= rank {
                    return Err(Error::new(
                        at,
                        format!("jetgen index {index} is out of range for rank {rank}"),
                    ));
                }
                let value = jet_page::persist(self.model, vec![(index, 1)])?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::QTerm(literal) => {
                let source = self.q_endpoint(literal.source.base)?;
                let target = self.q_endpoint(literal.target.base)?;
                let selected = q_koszul::indices(self.model, literal.selected)?;
                let (s, t) = q_module::term_bidegree(
                    self.model,
                    literal.source.degree,
                    literal.target.degree,
                    &selected,
                    &literal.powers,
                )?;
                let ty = self.model.q_koszul_type(Grade {
                    source,
                    target,
                    selected,
                    s,
                    t,
                })?;
                let source = q_module::generator(
                    self.model,
                    source,
                    literal.source.degree,
                    literal.source.index,
                )
                .map_err(|error| error.at(literal.source.span))?;
                let target = q_module::generator(
                    self.model,
                    target,
                    literal.target.degree,
                    literal.target.index,
                )
                .map_err(|error| error.at(literal.target.span))?;
                let powers = literal
                    .powers
                    .into_iter()
                    .map(|power| power as i128)
                    .collect();
                let terms = q_koszul::persist(
                    self.model,
                    vec![q_koszul::Term {
                        source,
                        target,
                        powers,
                    }],
                )?;
                Ok((ty, TypedForm::Constant(Value::QKoszul(terms))))
            }
            ExprForm::QIdentity(base, selected) => {
                let base = self.q_endpoint(base)?;
                let selected = q_koszul::indices(self.model, selected)?;
                let ty = self.model.q_ext_type(Grade {
                    source: base,
                    target: base,
                    selected,
                    s: 0,
                    t: 0,
                })?;
                let value = q_koszul::identity(self.model, ty)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RVGenerator(spec, index, at) => {
                let grade = self.rv_grade(spec)?;
                let ty = self.model.rv_hom_type(grade)?;
                let rank = self.model.rv_hom(ty)?.rank();
                self.check_generator(ty, rank, index, at)?;
                let value = rv_hom::persist(self.model, vec![(index, 1)])?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RVIdentity(base, selected, modulus, point) => {
                let base = self.q_endpoint(base)?;
                let selected = q_koszul::indices(self.model, selected)?;
                let grade = rv_hom::Grade::new(self.model, base, base, selected, modulus, point)?;
                let ty = self.model.rv_hom_type(grade)?;
                let value = rv_hom::identity(self.model, ty)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::QExtGenerator(spec, index, at) => {
                let grade = self.q_grade(spec)?;
                let ty = self.model.q_ext_type(grade)?;
                let rank = self.model.q_ext(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::VExtGenerator(spec, index, at) => {
                let grade = self.q_grade(spec)?;
                let ty = self.model.v_ext_type(grade)?;
                let rank = self.model.q_ext(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RvMonomial(selected, powers) => {
                let selected = q_koszul::indices(self.model, selected)?;
                if selected.len() != powers.len() {
                    return Err(Error::new(
                        span,
                        "rvmon requires one exponent per selected primitive",
                    ));
                }
                let ty = self
                    .model
                    .rv_poly_type(selected, rv_poly::degree(&powers)?)?;
                let polynomial = rv_poly::persist(self.model, rv_poly::monomial(powers))?;
                Ok((ty, TypedForm::Constant(Value::RVPoly(polynomial))))
            }
            ExprForm::RvMinor(spec, rows, columns) => {
                let base = self.q_endpoint(*spec.base)?;
                let selected = q_koszul::indices(self.model, spec.selected)?;
                let polynomial = rv_poly::minor(self.model, base, &selected, &rows, &columns)?;
                let ty = self.model.rv_poly_type(selected, rows.len())?;
                let polynomial = rv_poly::persist(self.model, polynomial)?;
                Ok((ty, TypedForm::Constant(Value::RVPoly(polynomial))))
            }
            ExprForm::RvHasse(orders, input) => {
                let input = self.expr(*input)?;
                let (selected, degree) = self.rv_poly_grade(&input, "rvHasse")?;
                if orders.len() != selected.len() {
                    return Err(Error::new(
                        span,
                        "rvHasse requires one order per selected primitive",
                    ));
                }
                let degree = degree
                    .checked_sub(rv_poly::degree(&orders)?)
                    .ok_or_else(|| {
                        Error::new(span, "rvHasse total order exceeds the polynomial degree")
                    })?;
                let ty = self.model.rv_poly_type(selected, degree)?;
                Ok((ty, TypedForm::RvHasse(orders, Box::new(input))))
            }
            ExprForm::RvVanish(modulus, point, input) => {
                let input = self.expr(*input)?;
                let (selected, _) = self.rv_poly_grade(&input, "rvVanish")?;
                let field = BinaryField::new(self.model, modulus)?;
                field.check_point(&point, selected.len(), "rvVanish")?;
                let ty = self.model.named("S").unwrap();
                Ok((ty, TypedForm::RvVanish(field, point, Box::new(input))))
            }
            ExprForm::RvCyclicGenerator(spec, mode, index, at) => {
                let ordinary = self.rv_hh_grade(spec)?;
                let ty = self
                    .model
                    .rv_cyclic_type(rv_cyclic_homology::Grade { ordinary, mode }, true)?;
                let rank = self.model.rv_cyclic_homology(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RvCyclicTerm(mode, power, input) => {
                let input = self.expr(*input)?;
                let Kind::RVHHChain(ordinary) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, mode.term_name(), "an RVHHChain"));
                };
                let grade = rv_cyclic_homology::Grade::term(ordinary.clone(), mode, power)?;
                let ty = self.model.rv_cyclic_type(grade, false)?;
                Ok((ty, TypedForm::RvNCTerm(power.into(), Box::new(input))))
            }
            ExprForm::RvCyclicPart(power, input) => {
                let input = self.expr(*input)?;
                let Kind::RvCyclicChain(grade) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, "rvCyclicPart", "a cyclic chain"));
                };
                let ty = self.model.rv_hh_type(grade.part(power)?)?;
                Ok((ty, TypedForm::RvNCPart(power.into(), Box::new(input))))
            }
            ExprForm::RvCyclicShift(shift, input) => {
                let input = self.expr(*input)?;
                let (grade, homology) = self.cyclic_grade(&input, shift.name())?;
                let ty = self.model.rv_cyclic_type(shift.target(&grade)?, homology)?;
                Ok((ty, TypedForm::RvCyclicShift(shift, Box::new(input))))
            }
            ExprForm::RvNCGenerator(spec, index, at) => {
                let grade = self.rv_negative_grade(spec)?;
                let ty = self.model.rv_negative_type(grade, true)?;
                let rank = self.model.rv_negative_homology(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RvNCTerm(length, power, input) => {
                let input = self.expr(*input)?;
                let Kind::RVHHChain(grade) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, "rvNCTerm", "an RVHHChain"));
                };
                let ty = self.model.rv_negative_type(
                    rv_negative::Grade::term(grade.clone(), length, power)?,
                    false,
                )?;
                Ok((ty, TypedForm::RvNCTerm(power as i128, Box::new(input))))
            }
            ExprForm::RvNCPart(power, input) => {
                let input = self.expr(*input)?;
                let Kind::RVNCChain(grade) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, "rvNCPart", "an RVNCChain"));
                };
                let ty = self.model.rv_hh_type(grade.part(power)?)?;
                Ok((ty, TypedForm::RvNCPart(power as i128, Box::new(input))))
            }
            ExprForm::RvNCTransfer(operation, input) => {
                let input = self.expr(*input)?;
                let (grade, homology) = match self.model.kind(input.ty) {
                    Kind::RVNCChain(grade) => (grade, false),
                    Kind::RvNegativeCyclic(grade) => (grade, true),
                    _ => {
                        return Err(self.invalid_operand(
                            &input,
                            operation.name(),
                            "an RVNCChain or RVNC class",
                        ))
                    }
                };
                let ty = self
                    .model
                    .rv_negative_type(operation.target(grade)?, homology)?;
                Ok((ty, TypedForm::RvNCTransfer(operation, Box::new(input))))
            }
            ExprForm::RvHHGenerator(spec, index, at) => {
                let grade = self.rv_hh_grade(spec)?;
                let ty = self.model.rv_hh_class_type(grade)?;
                let rank = self.model.rv_hh_homology(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RvPathGenerator(spec, index, at) => {
                let grade = self.rv_bar_grade(spec)?;
                let ty = self.model.rv_path_type(grade, true)?;
                let rank = self.model.rv_path_ext(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RvBarEGenerator(spec, index, at) => {
                let grade = self.rv_bar_page_grade(spec)?;
                let ty = self.model.rv_bar_page_type(grade)?;
                let rank = self.model.rv_bar_page(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RvBarClass(r, length, argument) => {
                let argument = self.expr(*argument)?;
                let Kind::RVBar(chain) = *self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(&argument, "rvBarClass", "an RVBar chain"));
                };
                let ty = self
                    .model
                    .rv_bar_page_type(rv_bar_page::Grade { chain, r, length })?;
                Ok((
                    ty,
                    TypedForm::Cobar(CobarOp::RvBarClass, Box::new(argument)),
                ))
            }
            ExprForm::RvHHWord(inputs) => {
                let inputs = inputs
                    .into_iter()
                    .map(|input| self.expr(input))
                    .collect::<Result<Vec<_>>>()?;
                let grade = rv_hh::word_type(self.model, inputs.iter().map(|input| input.ty))?;
                let ty = self.model.rv_hh_type(grade)?;
                Ok((ty, TypedForm::RvHHWord(inputs)))
            }
            ExprForm::RvHHLength(length, input) => {
                let input = self.expr(*input)?;
                if !matches!(self.model.kind(input.ty), Kind::RVHHChain(_)) {
                    return Err(self.invalid_operand(&input, "rvHHLength", "an RVHHChain"));
                }
                Ok((input.ty, TypedForm::RvHHLength(length, Box::new(input))))
            }
            ExprForm::RvBarWord(path, cuts, inputs) => {
                let path = self.model.rv_bar_path_id(&path)?;
                let inputs = inputs
                    .into_iter()
                    .map(|input| self.expr(input))
                    .collect::<Result<Vec<_>>>()?;
                let (grade, degrees) = rv_bar::word_type(
                    self.model,
                    path,
                    &cuts,
                    inputs.iter().map(|input| input.ty),
                )?;
                let ty = self.model.rv_bar_type(grade)?;
                Ok((
                    ty,
                    TypedForm::RvBarWord(cuts.into_iter().skip(1).collect(), degrees, inputs),
                ))
            }
            ExprForm::RvBarLength(length, input) => {
                let input = self.expr(*input)?;
                if !matches!(self.model.kind(input.ty), Kind::RVBar(_)) {
                    return Err(self.invalid_operand(&input, "rvBarLength", "an RVBar chain"));
                }
                Ok((input.ty, TypedForm::RvBarLength(length, Box::new(input))))
            }
            ExprForm::RvHFace(name, vertices, face) => {
                let (grade, terms) = rv_hochster::face(self.model, name, &vertices, &face)?;
                let ty = self.model.rv_complex_type(grade, false)?;
                Ok((ty, TypedForm::Constant(Value::RVKoszul(terms))))
            }
            ExprForm::RvMTerm(name, indices, input) => {
                let ExprForm::RvMonomial(selected, powers) = input.kind else {
                    return Err(Error::new(
                        input.span,
                        "rvmkterm requires an rvmon literal to infer its multidegree",
                    ));
                };
                let ideal = self.model.rv_ideal(&name)?;
                if selected.as_slice() != ideal.selected.as_ref() {
                    return Err(Error::new(
                        input.span,
                        "rvmkterm coefficient uses a different primitive selection",
                    ));
                }
                if selected.len() != powers.len() {
                    return Err(Error::new(
                        input.span,
                        "rvmon requires one exponent per selected primitive",
                    ));
                }
                let wedge = rv_koszul::mask(&ideal, &indices)?;
                let multidegree = rv_koszul::add_wedge(&powers, wedge)?;
                let grade = rv_koszul::Grade::local(
                    self.model,
                    &ideal,
                    name.text.into(),
                    indices.len(),
                    multidegree,
                )?;
                let ty = self.model.rv_complex_type(grade, false)?;
                let coefficient = rv_poly::monomial(powers);
                let terms = rv_koszul::term(self.model, &ideal, wedge, coefficient)?;
                Ok((ty, TypedForm::Constant(Value::RVKoszul(terms))))
            }
            ExprForm::RvMComponent(multidegree, input) => {
                let input = self.expr(*input)?;
                let (mut grade, homology) = self.rv_complex_family(&input, "rvMComponent")?;
                if grade.multidegree.is_some() {
                    return Err(Error::new(
                        input.span,
                        "rvMComponent requires a total Koszul chain or Tor class",
                    ));
                }
                let ideal = self.model.rv_ideal_named(&grade.ideal).unwrap();
                let internal = grade.t;
                grade =
                    rv_koszul::Grade::local(self.model, &ideal, grade.ideal, grade.s, multidegree)?;
                if grade.t != internal {
                    return Err(Error::new(
                        input.span,
                        "rvMComponent multidegree must have the same total internal degree",
                    ));
                }
                let ty = self.model.rv_complex_type(grade, homology)?;
                Ok((ty, TypedForm::RvMComponent(Box::new(input))))
            }
            ExprForm::RvKoszulTerm(name, indices, input) => {
                let ideal = self.model.rv_ideal(&name)?;
                let input = self.expr(*input)?;
                let (selected, degree) = self.rv_poly_grade(&input, "rvkterm")?;
                if selected != ideal.selected {
                    return Err(Error::new(
                        input.span,
                        "rvkterm coefficient uses a different primitive selection",
                    ));
                }
                let wedge = rv_koszul::mask(&ideal, &indices)?;
                let s = indices.len();
                let t = checked_degree(degree, s, "rvkterm", span)?;
                let grade = rv_koszul::Grade::total(name.text.into(), s, t);
                let ty = self.model.rv_complex_type(grade, false)?;
                Ok((ty, TypedForm::RvKoszulTerm(wedge, ideal, Box::new(input))))
            }
            ExprForm::RvTorGenerator(spec, index, at) => {
                let grade = self.rv_koszul_grade(spec)?;
                let ty = self.model.rv_complex_type(grade, true)?;
                let rank = self.model.rv_tor(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RvIdeal(op, name, input) => {
                let ideal = self.model.rv_ideal(&name)?;
                let input = self.expr(*input)?;
                let operation = match op {
                    RvIdealOp::Normal => "rvNormal",
                    RvIdealOp::Contains => "rvInIdeal",
                    RvIdealOp::Witness(_) => "rvWitness",
                };
                let (selected, degree) = self.rv_poly_grade(&input, operation)?;
                if selected != ideal.selected {
                    return Err(Error::new(
                        span,
                        format!("{operation} polynomial uses a different primitive selection"),
                    ));
                }
                let ty = match op {
                    RvIdealOp::Normal => input.ty,
                    RvIdealOp::Contains => self.model.named("S").unwrap(),
                    RvIdealOp::Witness(index) => {
                        let degree = ideal.witness_degree(index, degree)?;
                        self.model.rv_poly_type(selected, degree)?
                    }
                };
                Ok((ty, TypedForm::RvIdeal(op, ideal, Box::new(input))))
            }
            ExprForm::EOType(sequence) => {
                let ty = self.model.eo_type(sequence.len())?;
                let value = eo::persist(self.model, sequence).map_err(|error| error.at(span))?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::QpQuadElement(name, a, b, c) => {
                let ty = self.model.qp_quad_type(&name)?;
                let value = QuadElement::new(a, b, c).map_err(|error| error.at(span))?;
                Ok((
                    ty,
                    TypedForm::Constant(rv_qp_quad::persist(self.model, value)?),
                ))
            }
            ExprForm::QpQuadUnary(operation, input) => {
                let input = self.expr(*input)?;
                let field = self.qp_quad_field(&input, operation.name())?;
                Ok((
                    input.ty,
                    TypedForm::QpQuadUnary(field, operation, Box::new(input)),
                ))
            }
            ExprForm::EOOddCover(name) => {
                let cover = self.model.rv_as_odd_cover(&name)?;
                let sequence = cover.derham(self.model)?.final_type(self.model, &cover)?;
                let ty = self.model.eo_type(cover.genus())?;
                let value = eo::persist(self.model, sequence)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::CrystalEO(name) => {
                let crystal = self.model.rv_crystal(&name)?;
                let ty = self.model.eo_type(crystal.genus())?;
                let value = crystal
                    .eo_type(self.model)
                    .map_err(|error| error.at(span))?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::CrystalCongruent(left, right, level) => {
                let left = self.model.rv_crystal(&left)?;
                let right = self.model.rv_crystal(&right)?;
                let result = left
                    .congruent(self.model, &right, level)
                    .map_err(|error| error.at(span))?;
                let ty = self.model.named("S").unwrap();
                Ok((ty, TypedForm::Constant(Value::Point(usize::from(result)))))
            }
            ExprForm::EOCover(name) => {
                let cover = self.model.rv_as_cover(&name)?;
                let sequence = rv_as_derham::final_type(self.model, &cover)?;
                let ty = self.model.eo_type(cover.genus())?;
                let value = eo::persist(self.model, sequence)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::BinaryScalar(modulus, coefficient) => {
                let field = BinaryField::new(self.model, modulus)?;
                field.check_element(coefficient)?;
                let ty = self.model.binary_scalar_type(field)?;
                Ok((ty, TypedForm::Constant(Value::BinaryScalar(coefficient))))
            }
            ExprForm::AsDeRhamGen(name, modulus, index, at) => {
                let ty = self.model.rv_as_derham_type(&name, modulus)?;
                let cover = self.model.rv_as_cover(&name)?;
                let value = rv_as_derham::generator(self.model, cover.genus(), index)
                    .map_err(|error| error.at(at))?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::AsDeRhamOperation(operation, input) => {
                use rv_as_derham::Operation::{Frobenius, Include, Scale, Verschiebung};
                let input = self.expr(*input)?;
                let name = operation.name();
                let (curve, field) = match (operation, self.model.kind(input.ty)) {
                    (Include, Kind::CurveDiff(DiffFamily::ArtinSchreier, curve, field))
                    | (Frobenius | Verschiebung | Scale(_), Kind::ASDeRham(curve, field)) => {
                        (curve.clone(), *field)
                    }
                    (Include, _) => {
                        return Err(self.invalid_operand(&input, name, "an ASDiff differential"))
                    }
                    _ => return Err(self.invalid_operand(&input, name, "an ASDeRham class")),
                };
                if let Scale(coefficient) = operation {
                    field.check_element(coefficient)?;
                }
                let cover = self
                    .model
                    .rv_as_cover_named(&curve)
                    .expect("registered Artin–Schreier cover");
                let ty = if matches!(operation, Include) {
                    self.model
                        .rv_as_derham_field(curve, field)
                        .map_err(|error| error.at(input.span))?
                } else {
                    input.ty
                };
                Ok((
                    ty,
                    TypedForm::AsDeRhamOperation(operation, cover, field, Box::new(input)),
                ))
            }
            ExprForm::PrimeScalar(prime, coefficient) => {
                let prime = self.scalar_prime(prime)?;
                let ty = self.model.prime_scalar_type(prime)?;
                let value = coefficient.rem_euclid(i128::from(prime.get())) as u32;
                Ok((ty, TypedForm::Constant(Value::PrimeScalar(value))))
            }
            ExprForm::AsOddTateGen(name, degree, index) => {
                let cover = self.model.rv_as_odd_cover(&name)?;
                let ty =
                    self.model
                        .rv_as_odd_tate_named(name.text.as_str().into(), &cover, degree)?;
                let rank = cover
                    .tate(self.model)?
                    .group(self.model, &cover, degree)?
                    .rank();
                let value = self.prime_generator(ty, rank, index, span)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::AsOddTateOperation(operation, input) => {
                use rv_as_odd_tate::Operation::*;
                let input = self.expr(*input)?;
                let name = operation.name();
                let (cover, degree) = match (operation, self.model.kind(input.ty)) {
                    (Differential(degree) | Class(degree), Kind::ASOddDeRham(cover, _)) => {
                        (cover.clone(), degree)
                    }
                    (
                        Representative | Frobenius | Verschiebung,
                        Kind::ASOddTate(cover, _, degree),
                    ) => (cover.clone(), *degree),
                    (Differential(_) | Class(_), _) => {
                        return Err(self.invalid_operand(&input, name, "an ASOddDeRham class"));
                    }
                    _ => return Err(self.invalid_operand(&input, name, "an ASOddTate class")),
                };
                let data = self
                    .model
                    .rv_as_odd_cover_named(&cover)
                    .expect("registered odd Artin–Schreier cover");
                let ty = match operation {
                    Class(_) => self.model.rv_as_odd_tate_named(cover, &data, degree)?,
                    Representative => self.model.rv_as_odd_derham_named(cover, &data)?,
                    _ => input.ty,
                };
                Ok((
                    ty,
                    TypedForm::AsOddTateOperation(operation, data, degree, Box::new(input)),
                ))
            }
            ExprForm::AsOddDeRhamGen(name, index) => {
                let ty = self.model.rv_as_odd_derham_type(&name)?;
                let cover = self.model.rv_as_odd_cover(&name)?;
                let value = cover
                    .derham(self.model)?
                    .generator(self.model, index)
                    .map_err(|error| error.at(span))?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::AsOddDeRhamOperation(operation, input) => {
                use rv_as_odd_derham::Operation::*;
                let input = self.expr(*input)?;
                let expected = operation.name();
                let name = match (operation, self.model.kind(input.ty)) {
                    (Include, Kind::ASOddDiff(name, _))
                    | (Frobenius | Verschiebung | Translate(_), Kind::ASOddDeRham(name, _)) => {
                        name.clone()
                    }
                    (Include, _) => {
                        return Err(self.invalid_operand(
                            &input,
                            expected,
                            "an ASOddDiff differential",
                        ))
                    }
                    _ => {
                        return Err(self.invalid_operand(&input, expected, "an ASOddDeRham class"))
                    }
                };
                let cover = self
                    .model
                    .rv_as_odd_cover_named(&name)
                    .expect("registered odd Artin–Schreier cover");
                let ty = if matches!(operation, Include) {
                    self.model.rv_as_odd_derham_named(name, &cover)?
                } else {
                    input.ty
                };
                Ok((
                    ty,
                    TypedForm::AsOddDeRhamOperation(operation, cover, Box::new(input)),
                ))
            }
            ExprForm::AsTwoOmega(name, i, j) => {
                let cover = self.model.rv_as_two_cover(&name)?;
                let ty = self.model.rv_as_two_diff_type(&name)?;
                let value = cover
                    .cartier(self.model)?
                    .generator(self.model, i, j)
                    .map_err(|error| error.at(span))?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::AsTwoOrdinaryGen(name, index) => {
                let cover = self.model.rv_as_two_cover(&name)?;
                let ty = self.model.rv_as_two_diff_type(&name)?;
                let value = cover
                    .cartier(self.model)?
                    .ordinary_generator(self.model, index)
                    .map_err(|error| error.at(span))?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::AsTwoDiffOperation(operation, input) => {
                let input = self.expr(*input)?;
                let Kind::ASTwoDiff(name, _) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(
                        &input,
                        operation.name(),
                        "an ASTwoDiff differential",
                    ));
                };
                let cover = self
                    .model
                    .rv_as_two_cover_named(name)
                    .expect("registered two-pole cover");
                Ok((
                    input.ty,
                    TypedForm::AsTwoDiffOperation(operation, cover, Box::new(input)),
                ))
            }
            ExprForm::AsOddOmega(name, b, r) => {
                let cover = self.model.rv_as_odd_cover(&name)?;
                let ty = self.model.rv_as_odd_diff_type(&name)?;
                let value = cover.cartier(self.model)?.generator(self.model, b, r)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::AsOddDiffOperation(operation, input) => {
                let input = self.expr(*input)?;
                let Kind::ASOddDiff(name, _) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(
                        &input,
                        operation.name(),
                        "an ASOddDiff differential",
                    ));
                };
                let cover = self
                    .model
                    .rv_as_odd_cover_named(name)
                    .expect("registered odd Artin–Schreier cover");
                Ok((
                    input.ty,
                    TypedForm::AsOddDiffOperation(operation, cover, Box::new(input)),
                ))
            }
            ExprForm::AsTwoCoefficient(name, character, index) => {
                let cover = self.model.rv_as_two_cover(&name)?;
                let ty = self.model.cyclotomic_type(cover.order())?;
                let value = cover.coefficient(self.model, character, index)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::AsOddCoefficient(name, character, index) => {
                let cover = self.model.rv_as_odd_cover(&name)?;
                let ty = self.model.cyclotomic_type(cover.order())?;
                let value = cover.coefficient(self.model, character, index)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RvOmega(family, name, modulus, index) => {
                let ty = self.model.rv_curve_diff_type(family, &name, modulus)?;
                let value = match family {
                    DiffFamily::ArtinSchreier => {
                        let cover = self.model.rv_as_cover(&name)?;
                        rv_as_cartier::omega(self.model, &cover, index)?
                    }
                    DiffFamily::Witt => {
                        let cover = self.model.rv_witt_cover(&name)?;
                        cover.omega(self.model, index)?
                    }
                };
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RvDiffScale(family, scalar, input) => {
                let input = self.expr(*input)?;
                let field = match *self.model.kind(input.ty) {
                    Kind::CurveDiff(actual, _, field) if actual == family => field,
                    _ => {
                        return Err(self.invalid_operand(
                            &input,
                            family.scale_name(),
                            family.operand(),
                        ))
                    }
                };
                field.check_element(scalar)?;
                Ok((
                    input.ty,
                    TypedForm::DiffScale(field, scalar, Box::new(input)),
                ))
            }
            ExprForm::RvCartier(family, input) => {
                let input = self.expr(*input)?;
                let (name, field) = match self.model.kind(input.ty) {
                    Kind::CurveDiff(actual, name, field) if *actual == family => (name, *field),
                    _ => {
                        return Err(self.invalid_operand(
                            &input,
                            family.cartier_name(),
                            family.operand(),
                        ))
                    }
                };
                let ty = input.ty;
                let form = match family {
                    DiffFamily::ArtinSchreier => {
                        let cover = self
                            .model
                            .rv_as_cover_named(name)
                            .expect("registered Artin–Schreier cover");
                        TypedForm::AsCartier(cover, field, Box::new(input))
                    }
                    DiffFamily::Witt => {
                        let cover = self
                            .model
                            .rv_witt_cover_named(name)
                            .expect("registered Witt cover");
                        TypedForm::WittCartier(cover, field, Box::new(input))
                    }
                };
                Ok((ty, form))
            }
            ExprForm::RvFrobeniusWronskian(name, exponent) => {
                let curve = self.model.rv_curve(&name)?;
                let (degree, polynomial) = curve.frobenius_wronskian(self.model, exponent)?;
                let polynomial = rv_poly::persist(self.model, polynomial)?;
                let ty = self.model.rv_poly_type(curve.selected.clone(), degree)?;
                Ok((ty, TypedForm::Constant(Value::RVPoly(polynomial))))
            }
            ExprForm::RvOscMinor(name, k, columns) => {
                let curve = self.model.rv_curve(&name)?;
                let (degree, polynomial) = curve.minor(self.model, k, columns)?;
                let polynomial = rv_poly::persist(self.model, polynomial)?;
                let ty = self.model.rv_poly_type(curve.selected.clone(), degree)?;
                Ok((ty, TypedForm::Constant(Value::RVPoly(polynomial))))
            }
            ExprForm::RvGenerator(name, index, at) => {
                let ideal = self.model.rv_ideal(&name)?;
                let (degree, polynomial) = ideal.generator(index).ok_or_else(|| {
                    Error::new(
                        at,
                        format!("rvGenerator reduced basis index {index} is out of range"),
                    )
                })?;
                let ty = self.model.rv_poly_type(ideal.selected.clone(), degree)?;
                Ok((ty, TypedForm::Constant(Value::RVPoly(polynomial))))
            }
            ExprForm::RvInMinors(spec, size, input) => {
                let base = self.q_endpoint(*spec.base)?;
                let selected = q_koszul::indices(self.model, spec.selected)?;
                self.model.check_q_endpoints(base, base, &selected)?;
                let input = self.expr(*input)?;
                let expected = self.model.rv_poly_type(selected.clone(), size)?;
                self.require_type(expected, &input, "rvInMinors")?;
                let ty = self.model.named("S").unwrap();
                Ok((
                    ty,
                    TypedForm::RvInMinors(base, selected, size, Box::new(input)),
                ))
            }
            ExprForm::VBockGenerator(spec, index, at) => {
                let ty = self.vbock_type(spec)?;
                let rank = self.model.vbock_page(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::VBockClass(r, f, input) => {
                let input = self.expr(*input)?;
                let Kind::VKoszul(localized) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, "vbClass", "a VKoszul cochain"));
                };
                let grade = vbock::Grade {
                    localized: localized.clone(),
                    r,
                    f,
                };
                let ty = self.model.vbock_type(grade)?;
                Ok((ty, TypedForm::Cobar(CobarOp::VBockClass, Box::new(input))))
            }
            ExprForm::VShift(k, input) => {
                let input = self.expr(*input)?;
                let ty = if matches!(self.model.kind(input.ty), Kind::VBock(..)) {
                    vbock::shift_type(self.model, input.ty, k)?
                } else {
                    v_koszul::shift_type(self.model, input.ty, k)?
                };
                Ok((ty, TypedForm::VShift(k, Box::new(input))))
            }
            ExprForm::VRatio(k, input) => {
                let input = self.expr(*input)?;
                let ty = v_bars::ratio_type(self.model, input.ty, k)?;
                Ok((ty, TypedForm::VRatio(k, Box::new(input))))
            }
            ExprForm::VBarCertificate(operation, spec, degree, exponent, index) => {
                let family = self.vbar_family(spec)?;
                let ty = operation.ty(self.model, family, degree, exponent)?;
                let bars = self.model.v_bars(family)?;
                let value = bars.certificate(self.model, operation, degree, exponent, index)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::QRestrict(selected, input) => {
                let input = self.expr(*input)?;
                let selected = q_koszul::indices(self.model, selected)?;
                let ty = q_koszul::restriction_type(self.model, input.ty, selected)?;
                Ok((ty, TypedForm::QRestrict(Box::new(input))))
            }
            ExprForm::MapMayGenerator(spec, index, at) => {
                let ty = self.map_may_type(spec)?;
                let rank = self.model.map_may(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::MayAnGenerator(n, r, s, t, f, index, at) => {
                let ty = self.model.profile_may_type(Profile::An(n), r, s, t, f)?;
                let rank = self.model.may_page(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::CeGenerator(n, r, s, t, f, index, at) => {
                let ty = self.model.ce_type(n, r, s, t, f)?;
                let rank = self.model.ce_page(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::CeClass(r, f, value) => {
                let value = self.expr(*value)?;
                let ty = match *self.model.kind(value.ty) {
                    Kind::CobarAn(n, s, t) => self.model.ce_type(n, r, s, t, f)?,
                    Kind::PCobarAn(p, n, s, t) => self.model.pce_type(p, n, r, s, t, f)?,
                    _ => {
                        return Err(self.invalid_operand(
                            &value,
                            "ceClass",
                            "a CobarAn cochain or a PCobarAn cochain",
                        ))
                    }
                };
                Ok((ty, TypedForm::Cobar(CobarOp::CeClass, Box::new(value))))
            }
            ExprForm::MayGenerator(r, s, t, f, index, at) => {
                let ty = self.model.may_type(r, s, t, f)?;
                let rank = self.model.may_page(ty)?.rank();
                let value = self.generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::MayClass(r, f, value) => {
                let value = self.expr(*value)?;
                let ty = match *self.model.kind(value.ty) {
                    Kind::MapCobar(source, target, s, t) => {
                        self.model.map_may_type(source, target, r, s, t, f)?
                    }
                    Kind::PCobar(p, s, t) => self.model.pmay_type(p, r, s, t, f)?,
                    Kind::PCobarAn(p, n, s, t) => {
                        self.model
                            .profile_pmay_type(p, Profile::An(n), r, s, t, f)?
                    }
                    _ => {
                        let (profile, s, t) = self.profile_cobar_bidegree(&value, "mayClass")?;
                        self.model.profile_may_type(profile, r, s, t, f)?
                    }
                };
                Ok((ty, TypedForm::Cobar(CobarOp::MayClass, Box::new(value))))
            }
            ExprForm::Cyclic(modulus, at, form) => {
                let ty = self
                    .model
                    .cyclic_type(modulus, matches!(form, CyclicForm::Carry))
                    .map_err(|error| error.at(at))?;
                let value = match form {
                    CyclicForm::Support(residues) => {
                        cyclic::support(self.model, modulus, residues)?
                    }
                    CyclicForm::Binomial(index, index_at) => {
                        cyclic::binomial(self.model, modulus, index, index_at)?
                    }
                    CyclicForm::Carry => Value::CyclicClass(modulus & 1 == 0),
                };
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::Massey(inputs, witnesses) => {
                self.defining_system(inputs, Some(witnesses), span)
            }
            ExprForm::PowerDef(inputs) => self.defining_system(inputs, None, span),
            ExprForm::Bracket(inputs) => self.bracket(inputs, span),
            ExprForm::Transfer(operation, inputs) => self.transfer(operation, inputs),
            ExprForm::HigherRestriction(n, inputs) => self.higher_restriction(n, inputs, span),
            ExprForm::RvTBracket(operation, inputs) => self.rv_t_bracket(operation, inputs),
            ExprForm::MapBracket(operation, inputs) => self.map_bracket(operation, inputs),
            ExprForm::MapTransfer(operation, inputs) => self.map_transfer(operation, inputs),
            ExprForm::RvTTransfer(operation, inputs) => self.rv_t_transfer(operation, inputs),
            ExprForm::RvTSystem(operation, inputs, witnesses) => {
                self.rv_t_system(operation, inputs, witnesses, span)
            }
            ExprForm::MapSystem(operation, inputs, witnesses) => {
                self.map_system(operation, inputs, witnesses, span)
            }
            ExprForm::OddSystem(operation, inputs, witnesses) => {
                self.odd_system(operation, inputs, witnesses, span)
            }
            ExprForm::Suspension(base, form) => {
                let base = self.ty(base)?;
                let suspension = self.model.suspension_type(base)?;
                match form {
                    SuspensionForm::North => Ok((
                        suspension,
                        TypedForm::Constant(Value::Suspension(SuspensionPoint::North)),
                    )),
                    SuspensionForm::South => Ok((
                        suspension,
                        TypedForm::Constant(Value::Suspension(SuspensionPoint::South)),
                    )),
                    SuspensionForm::Equator => {
                        Ok((self.model.arrow(base, suspension)?, TypedForm::Equator))
                    }
                }
            }
            ExprForm::Shift(amount, value) => {
                let value = self.expr(*value)?;
                let (base, existing, degree) = match *self.model.kind(value.ty) {
                    Kind::Cohomology(base, degree) => (base, 0, degree),
                    Kind::Stable(base, shift, degree) => (base, shift, degree),
                    _ => return Err(self.invalid_operand(&value, "shift", "a cohomology class")),
                };
                let shift = checked_degree(existing, amount, "shift", span)?;
                let degree = checked_degree(degree, amount, "shift", span)?;
                let ty = self.model.stable_type(base, shift, degree)?;
                Ok((ty, TypedForm::Shift(Box::new(value))))
            }
            ExprForm::Binary(op, left, right) => {
                let left = self.expr(*left)?;
                let right = self.expr(*right)?;
                let ty = self.binary_result(op, &left, &right, span)?;
                Ok((ty, TypedForm::Binary(op, Box::new(left), Box::new(right))))
            }
            ExprForm::Lambda(parameter, domain, body) => {
                let domain = self.ty(domain)?;
                self.locals.push((parameter.text, domain));
                let body = self.expr(*body);
                self.locals.pop();
                let body = body?;
                let ty = self.model.arrow(domain, body.ty)?;
                Ok((ty, TypedForm::Lambda(Box::new(body))))
            }
            ExprForm::Call(function, argument) => {
                if let ExprForm::Name(name) = &function.kind {
                    if self.unbound(&name.text) {
                        if let Some(op) = CobarOp::from_name(&name.text) {
                            let argument = self.expr(*argument)?;
                            let ty = self.cobar_result(&op, &argument, &name.text, span)?;
                            return Ok((ty, TypedForm::Cobar(op, Box::new(argument))));
                        }
                        if let Some(operation) = transfer::Operation::from_name(&name.text) {
                            return self.transfer(operation, vec![*argument]);
                        }
                        if let Some(operation) = map_transfer::Operation::from_name(&name.text) {
                            return self.map_transfer(operation, vec![*argument]);
                        }
                        if let Some(operation) = rv_t_transfer::Operation::from_name(&name.text) {
                            return self.rv_t_transfer(operation, vec![*argument]);
                        }
                        if let Some(operation) = rv_t_system::Operation::from_name(&name.text) {
                            return self.rv_t_system(operation, vec![*argument], Vec::new(), span);
                        }
                        if let Some(operation) = map_system::Operation::from_name(&name.text) {
                            return self.map_system(operation, vec![*argument], Vec::new(), span);
                        }
                        if let Some(operation) = odd_massey::Operation::from_name(&name.text) {
                            return self.odd_system(operation, vec![*argument], Vec::new(), span);
                        }
                        if let Some(operation) = rv_t_bracket::Operation::from_name(&name.text) {
                            return self.rv_t_bracket(operation, vec![*argument]);
                        }
                        if let Some(operation) = map_bracket::Operation::from_name(&name.text) {
                            return self.map_bracket(operation, vec![*argument]);
                        }
                        if name.text == "emap" {
                            let argument = self.expr(*argument)?;
                            let Kind::Arrow(source, target) = *self.model.kind(argument.ty) else {
                                return Err(self.invalid_operand(&argument, "emap", "a function"));
                            };
                            let ty = self.model.map_ext_type(source, target, 0, 0)?;
                            return Ok((ty, TypedForm::MapExtInduced(Box::new(argument))));
                        }
                    }
                }
                if matches!(&function.kind, ExprForm::Name(name) if name.text == "susp")
                    && self.unbound("susp")
                {
                    let argument = self.expr(*argument)?;
                    let (base, degree) = self.class_type(&argument, "susp")?;
                    let degree = checked_degree(degree, 1, "susp", span)?;
                    let suspension = self.model.suspension_type(base)?;
                    let ty = self.model.cohomology_type(suspension, degree)?;
                    return Ok((ty, TypedForm::Suspend(Box::new(argument))));
                }
                if matches!(&function.kind, ExprForm::Name(name) if name.text == "powerdef")
                    && self.unbound("powerdef")
                {
                    return self.defining_system(vec![*argument], None, span);
                }
                if matches!(&function.kind, ExprForm::Name(name) if name.text == "bracket")
                    && self.unbound("bracket")
                {
                    return self.bracket(vec![*argument], span);
                }
                let function = self.expr(*function)?;
                let Kind::Arrow(domain, range) = *self.model.kind(function.ty) else {
                    return Err(Error::new(
                        function.span,
                        format!(
                            "application requires a function, found {}",
                            self.model.type_name(function.ty)
                        ),
                    ));
                };
                let argument = self.expr(*argument)?;
                self.require_type(domain, &argument, "function argument")?;
                Ok((
                    range,
                    TypedForm::Call(Box::new(function), Box::new(argument)),
                ))
            }
            ExprForm::Pair(left, right) => {
                let left = self.expr(*left)?;
                let right = self.expr(*right)?;
                let ty = self.model.product(left.ty, right.ty)?;
                Ok((ty, TypedForm::Pair(Box::new(left), Box::new(right))))
            }
            ExprForm::Fst(pair) => self.projection(*pair, true),
            ExprForm::Snd(pair) => self.projection(*pair, false),
            ExprForm::Table(domain, range, entries) => {
                let domain = self.ty(domain)?;
                let range = self.ty(range)?;
                let rows = self.table(domain, range, entries, span)?;
                let ty = self.model.arrow(domain, range)?;
                Ok((ty, TypedForm::Table(rows)))
            }
            ExprForm::Restrict(point, base) => {
                let base = self.expr(*base)?;
                let Kind::Arrow(domain, range) = *self.model.kind(base.ty) else {
                    return Err(Error::new(
                        base.span,
                        format!(
                            "res requires a function, found {}",
                            self.model.type_name(base.ty)
                        ),
                    ));
                };
                let star = self.model.star(domain, &point)?;
                let ty = self.model.arrow(star, range)?;
                Ok((ty, TypedForm::Restrict(Box::new(base))))
            }
            ExprForm::Sew(domain, range, entries) => {
                let domain = self.ty(domain)?;
                let range = self.ty(range)?;
                let charts = self.atlas(domain, range, entries, span)?;
                let ty = self.model.arrow(domain, range)?;
                Ok((ty, TypedForm::Sew(charts)))
            }
            rank_arithmetic_expression!() | bp_expression!() => unreachable!(),
        }
    }

    fn expr_rank_arithmetic(&mut self, expr: ExprForm) -> Result<(TypeId, TypedForm)> {
        match expr {
            ExprForm::HeckeBasis(canonical, window) => {
                let ty = self.model.iw_hecke_type(window.len())?;
                let value = if canonical {
                    rv_iwhecke::canonical(self.model, &window)?
                } else {
                    rv_iwhecke::standard(self.model, window)?
                };
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::HeckeScalar(rank, degree) => {
                let ty = self.model.iw_hecke_type(rank)?;
                let value = rv_iwhecke::scalar(self.model, rank, degree)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::HeckeNegate(input) => {
                let input = self.expr(*input)?;
                if matches!(self.model.kind(input.ty), Kind::RVQpQuad(..)) {
                    let field = self.qp_quad_field(&input, "neg")?;
                    return Ok((
                        input.ty,
                        TypedForm::QpQuadUnary(field, QuadUnary::Negate, Box::new(input)),
                    ));
                }
                if matches!(self.model.kind(input.ty), Kind::RVOre) {
                    return Ok((input.ty, TypedForm::HeckeNegate(Box::new(input))));
                }
                if let Kind::RVQGrP(grade) = *self.model.kind(input.ty) {
                    return Ok((
                        input.ty,
                        TypedForm::QuantumPUnary(grade, QuantumOp::Scale(-1), Box::new(input)),
                    ));
                }
                if matches!(self.model.kind(input.ty), Kind::RVQGr(..)) {
                    return Ok((
                        input.ty,
                        TypedForm::QuantumUnary(QuantumOp::Scale(-1), Box::new(input)),
                    ));
                }
                if let Kind::RvTl(grade) = *self.model.kind(input.ty) {
                    return Ok((
                        input.ty,
                        TypedForm::TlOperation(grade, TlOp::Negate, Box::new(input)),
                    ));
                }
                if matches!(self.model.kind(input.ty), Kind::RVCycloQ(..)) {
                    return Ok((input.ty, TypedForm::CycloRatioNegate(Box::new(input))));
                }
                if let Kind::RVWzw(grade) = *self.model.kind(input.ty) {
                    return Ok((
                        input.ty,
                        TypedForm::WzwUnary(grade, WzwOp::Negate, Box::new(input)),
                    ));
                }
                if matches!(self.model.kind(input.ty), Kind::RVCyclo(..)) {
                    return Ok((input.ty, TypedForm::CycloScale(-1, Box::new(input))));
                }
                if let Kind::RVFusion(grade) = *self.model.kind(input.ty) {
                    return Ok((input.ty, TypedForm::FusionScale(grade, -1, Box::new(input))));
                }
                if !matches!(
                    self.model.kind(input.ty),
                    Kind::RVIwHecke(..) | Kind::RVSphHecke(..) | Kind::RVSatake(..)
                ) {
                    return Err(self.invalid_operand(
                        &input,
                        "neg",
                        "an affine Hecke element, a spherical Hecke element, a Satake polynomial, a fusion element, or a cyclotomic integer",
                    ));
                }
                Ok((input.ty, TypedForm::HeckeNegate(Box::new(input))))
            }
            ExprForm::WzwBasis(rank, level, weight) => {
                let grade = rv_wzw::Grade::checked(rank, level)?;
                let ty = self.model.wzw_type(grade)?;
                let value = rv_wzw::basis(self.model, grade, weight)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::WzwSEntry(rank, level, left, right) => {
                let grade = rv_wzw::Grade::checked(rank, level)?;
                let ty = self.model.cyclotomic_type(grade.order)?;
                let value = rv_wzw::s_entry(self.model, grade, left, right)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::WzwTEntry(rank, level, weight, power) => {
                let grade = rv_wzw::Grade::checked(rank, level)?;
                let ty = self.model.cyclotomic_type(grade.order)?;
                let value = rv_wzw::t_entry(self.model, grade, weight, power)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::WzwUnary(operation, input) => {
                let input = self.expr(*input)?;
                let grade = self.wzw_grade(&input, operation.name())?;
                Ok((
                    input.ty,
                    TypedForm::WzwUnary(grade, operation, Box::new(input)),
                ))
            }
            ExprForm::WzwScale(arguments) => {
                let [scalar, input] = *arguments;
                let input = self.expr(input)?;
                let grade = self.wzw_grade(&input, "rvWzwScale")?;
                let scalar = self.expr(scalar)?;
                let ty = self.model.cyclotomic_type(grade.order)?;
                self.require_type(ty, &scalar, "rvWzwScale")?;
                Ok((
                    input.ty,
                    TypedForm::WzwScale(grade, Box::new([scalar, input])),
                ))
            }
            ExprForm::LensTv(rank, level, p, q) => {
                let grade = rv_wzw::Grade::checked(rank, level)?;
                let ty = self.model.cyclotomic_rational_type(grade.order)?;
                let value = rv_lens::lens(self.model, grade, p, q)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::WzwTv(rank, level, word) => {
                let grade = rv_wzw::Grade::checked(rank, level)?;
                let ty = self.model.cyclotomic_rational_type(grade.order)?;
                let value = rv_lens::word(self.model, grade, &word)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::TlPath(level, path) => {
                let Some(&height) = path.last() else {
                    return Err(Error::plain("Temperley–Lieb path must start at zero"));
                };
                let grade = rv_tl::Grade::checked(level, path.len() - 1, height)?;
                let ty = self.model.tl_type(grade)?;
                let value = rv_tl::basis(self.model, grade, path)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::TlTrace(level, strands, height, word) => {
                let grade = rv_tl::Grade::checked(level, strands, height)?;
                let ty = self.model.cyclotomic_rational_type(grade.order)?;
                let value = rv_tl::trace(self.model, grade, &word)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::TlJones(level, strands, word) => {
                let grade = rv_tl::Grade::checked(level, strands, 0)?;
                let ty = self.model.cyclotomic_rational_type(grade.order)?;
                let value = rv_tl::jones(self.model, grade, &word)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::TlColoredJones(level, strands, color, word) => {
                let grade = rv_tl::Grade::checked(level, strands, 0)?;
                let ty = self.model.cyclotomic_rational_type(grade.order)?;
                let value = rv_tl::colored_jones(self.model, grade, color, &word)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::TlFramed(level, strands, word, colors, framings) => {
                let grade = rv_tl::Grade::framed(level, strands)?;
                let ty = self.model.cyclotomic_rational_type(grade.order)?;
                let value = rv_tl::framed(self.model, grade, &word, &colors, &framings)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::TlSurgery(kind, level, strands, word, framings) => {
                let grade = rv_tl::Grade::framed(level, strands)?;
                let ty = self.model.cyclotomic_rational_type(grade.order)?;
                let value = rv_tl::surgery_value(self.model, grade, &word, &framings, kind)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::TlOperation(operation, input) => {
                let input = self.expr(*input)?;
                let grade = self.tl_grade(&input, "Temperley–Lieb action")?;
                Ok((
                    input.ty,
                    TypedForm::TlOperation(grade, operation, Box::new(input)),
                ))
            }
            ExprForm::TlCoordinate(path, input) => {
                let input = self.expr(*input)?;
                let grade = self.tl_grade(&input, "rvTLCoord")?;
                let ty = self.model.cyclotomic_rational_type(grade.order)?;
                Ok((ty, TypedForm::TlCoordinate(grade, path, Box::new(input))))
            }
            ExprForm::TlScale(arguments) => {
                let [scalar, input] = *arguments;
                let input = self.expr(input)?;
                let grade = self.tl_grade(&input, "rvTLScale")?;
                let scalar = self.expr(scalar)?;
                let ty = self.model.cyclotomic_rational_type(grade.order)?;
                self.require_type(ty, &scalar, "rvTLScale")?;
                Ok((
                    input.ty,
                    TypedForm::TlScale(grade, Box::new([scalar, input])),
                ))
            }
            ExprForm::CycloRatioProduct(arguments) => {
                let [left, right] = *arguments;
                let left = self.expr(left)?;
                let order = self.cyclotomic_ratio_order(&left, "rvCycloQMul")?;
                let right = self.expr(right)?;
                self.require_type(left.ty, &right, "rvCycloQMul")?;
                Ok((
                    left.ty,
                    TypedForm::CycloRatioProduct(order, Box::new([left, right])),
                ))
            }
            ExprForm::CycloRatioGalois(power, input) => {
                let input = self.expr(*input)?;
                let order = self.cyclotomic_ratio_order(&input, "rvCycloQGalois")?;
                Ok((
                    input.ty,
                    TypedForm::CycloRatioGalois(order, power, Box::new(input)),
                ))
            }
            ExprForm::CycloRatio(denominator, input) => {
                let input = self.expr(*input)?;
                let order = self.cyclotomic_order(&input, "rvCycloRat")?;
                let ty = self.model.cyclotomic_rational_type(order)?;
                Ok((ty, TypedForm::CycloRatio(denominator, Box::new(input))))
            }
            ExprForm::WzwCoordinate(weight, input) => {
                let input = self.expr(*input)?;
                let grade = self.wzw_grade(&input, "rvWzwCoord")?;
                let ty = self.model.cyclotomic_type(grade.order)?;
                Ok((ty, TypedForm::WzwCoordinate(grade, weight, Box::new(input))))
            }
            ExprForm::CycloGalois(power, input) => {
                let input = self.expr(*input)?;
                let order = self.cyclotomic_order(&input, "rvCycloGalois")?;
                Ok((
                    input.ty,
                    TypedForm::CycloGalois(order, power, Box::new(input)),
                ))
            }
            ExprForm::CycloEmbed(target, input) => {
                let input = self.expr(*input)?;
                let source = self.cyclotomic_order(&input, "rvCycloEmbed")?;
                let ty = self.model.cyclotomic_type(target)?;
                if !target.is_multiple_of(source) {
                    return Err(Error::new(
                        input.span,
                        "cyclotomic embedding requires the source order to divide the target order",
                    ));
                }
                Ok((ty, TypedForm::CycloEmbed(source, target, Box::new(input))))
            }
            ExprForm::CycloRoot(order, power) => {
                let ty = self.model.cyclotomic_type(order)?;
                let value = rv_cyclotomic::root(self.model, order, power)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::CycloProduct(arguments) => {
                let [left, right] = *arguments;
                let left = self.expr(left)?;
                let order = self.cyclotomic_order(&left, "rvCycloMul")?;
                let right = self.expr(right)?;
                self.require_type(left.ty, &right, "rvCycloMul")?;
                Ok((
                    left.ty,
                    TypedForm::CycloProduct(order, Box::new([left, right])),
                ))
            }
            ExprForm::CycloScale(factor, input) => {
                let input = self.expr(*input)?;
                self.cyclotomic_order(&input, "rvCycloScale")?;
                Ok((input.ty, TypedForm::CycloScale(factor, Box::new(input))))
            }
            ExprForm::CycloConjugate(input) => {
                let input = self.expr(*input)?;
                let order = self.cyclotomic_order(&input, "rvCycloConj")?;
                Ok((input.ty, TypedForm::CycloConjugate(order, Box::new(input))))
            }
            ExprForm::FusionHandle(rank, level) => {
                let grade = rv_fusion::Grade::checked(rank, level)?;
                let ty = self.model.fusion_type(grade)?;
                let value = rv_fusion::handle(self.model, grade)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::Verlinde(label, input) => {
                let input = self.expr(*input)?;
                let Kind::RVFusion(grade) = *self.model.kind(input.ty) else {
                    let name = if label.is_some() {
                        "rvVerlinde"
                    } else {
                        "rvQDim"
                    };
                    return Err(self.invalid_operand(&input, name, "a fusion element"));
                };
                let order = rv_verlinde::order(grade)?;
                let ty = self.model.cyclotomic_type(order)?;
                let point = rv_verlinde::torus(
                    self.model,
                    grade,
                    label.unwrap_or_else(|| vec![0; grade.rank]),
                )?;
                Ok((ty, TypedForm::Verlinde(order, point, Box::new(input))))
            }
            ExprForm::FusionBasis(level, weight) => {
                let grade = rv_fusion::Grade::checked(weight.len(), level)?;
                let ty = self.model.fusion_type(grade)?;
                let value = rv_fusion::basis(self.model, grade, weight)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::FusionOne(rank, level) => {
                let grade = rv_fusion::Grade::checked(rank, level)?;
                let ty = self.model.fusion_type(grade)?;
                let value = rv_fusion::basis(self.model, grade, vec![0; rank])?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::FusionProduct(arguments) => {
                let [left, right] = *arguments;
                let left = self.expr(left)?;
                let Kind::RVFusion(grade) = *self.model.kind(left.ty) else {
                    return Err(self.invalid_operand(&left, "rvFusionMul", "a fusion element"));
                };
                let right = self.expr(right)?;
                self.require_type(left.ty, &right, "rvFusionMul")?;
                Ok((
                    left.ty,
                    TypedForm::FusionProduct(grade, Box::new([left, right])),
                ))
            }
            ExprForm::FusionScale(factor, input) => {
                let input = self.expr(*input)?;
                let Kind::RVFusion(grade) = *self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, "rvFusionScale", "a fusion element"));
                };
                Ok((
                    input.ty,
                    TypedForm::FusionScale(grade, factor, Box::new(input)),
                ))
            }
            ExprForm::FusionDual(input) => {
                let input = self.expr(*input)?;
                let Kind::RVFusion(grade) = *self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, "rvFusionDual", "a fusion element"));
                };
                Ok((input.ty, TypedForm::FusionDual(grade, Box::new(input))))
            }
            ExprForm::ResSpace(p, n, rows) => {
                let prime = crate::rv_phodge::residue_prime(self.model, p)?;
                let grade = ResGrade::new(prime, n)?;
                let space = ResSpace::entries(self.model, grade, rows)?;
                let ty = self.model.rv_res_space_type(grade)?;
                Ok((ty, TypedForm::Constant(space.persist(self.model)?)))
            }
            ExprForm::ResSpecialize(p, n, weights, rows) => {
                let prime = crate::rv_phodge::residue_prime(self.model, p)?;
                let grade = ResGrade::new(prime, n)?;
                let space = rv_res_space::rational(self.model, grade, &weights, rows)?;
                let ty = self.model.rv_res_space_type(grade)?;
                Ok((ty, TypedForm::Constant(space.persist(self.model)?)))
            }
            ExprForm::ResSpecializeQuad(name, n, weights, rows) => {
                let field = self.model.qp_quad(&name)?;
                let field_type = self.model.qp_quad_type(&name)?;
                let grade = ResGrade::new(field.prime(), n)?;
                rv_res_space::validate_rows(grade, &rows)?;
                let rows = self.check_quadratic_rows(field_type, rows, "rvResSpecializeQuad")?;
                let ty = self.model.rv_res_space_type(grade)?;
                Ok((
                    ty,
                    TypedForm::ResSpecializeQuad(grade, field, weights, rows),
                ))
            }
            ExprForm::ResFil(name, cutoff, weights) => {
                let module = self.model.rv_phodge(&name)?;
                let space = module.residue_filtration(self.model, cutoff, &weights)?;
                let ty = self.model.rv_res_space_type(space.grade)?;
                Ok((ty, TypedForm::Constant(space.persist(self.model)?)))
            }
            ExprForm::ResOrthogonal(input) => {
                let input = self.expr(*input)?;
                self.require_res_space(&input, "rvResOrthogonal")?;
                Ok((input.ty, TypedForm::ResOrthogonal(Box::new(input))))
            }
            ExprForm::StableLine(p, marks) => {
                let prime = crate::rv_phodge::residue_prime(self.model, p)?;
                let ty = self.model.rv_stable_curve_type(prime.get(), marks.len())?;
                let value = StableCurve::rational(self.model, prime, marks)?;
                Ok((ty, TypedForm::Constant(Value::RVStableMarks(value))))
            }
            ExprForm::StableLineQuad(name, marks) => {
                let field = self.model.qp_quad(&name)?;
                let field_type = self.model.qp_quad_type(&name)?;
                let ty = self
                    .model
                    .rv_stable_curve_type(field.prime().get(), marks.len())?;
                let mut points = Vec::with_capacity(marks.len());
                for [first, second] in marks {
                    let a = self.expr(first)?;
                    let b = self.expr(second)?;
                    self.require_type(field_type, &a, "rvStableLineQuad")?;
                    self.require_type(field_type, &b, "rvStableLineQuad")?;
                    points.push([a, b]);
                }
                Ok((ty, TypedForm::StableLineQuad(field, points)))
            }
            ExprForm::StableProjection(input) => {
                let input = self.expr(*input)?;
                let count = match *self.model.kind(input.ty) {
                    Kind::RVValMat(grade) if grade.rank == 2 => grade.ground,
                    Kind::RVStableMarks(_, count) => count,
                    _ => {
                        return Err(self.invalid_operand(
                            &input,
                            "rvStableTree",
                            "a rank-two RVValMat or an RVStableMarks",
                        ))
                    }
                };
                let ty = self.model.rv_stable_tree_type(count)?;
                Ok((ty, TypedForm::StableProjection(Box::new(input))))
            }
            ExprForm::ValMat(p, n, rows) => {
                let prime = crate::rv_phodge::residue_prime(self.model, p)?;
                let grade = ValGrade::new(prime.get(), n, rows.len())?;
                let value = tropical_matroid::rational(self.model, grade, rows)?;
                let ty = self.model.rv_valuated_matroid_type(grade)?;
                Ok((ty, TypedForm::Constant(Value::RVValMat(value))))
            }
            ExprForm::ValMatQuad(name, n, rows) => {
                let field = self.model.qp_quad(&name)?;
                let field_type = self.model.qp_quad_type(&name)?;
                let grade = ValGrade::new(field.prime().get(), n, rows.len())?;
                tropical_matroid::check_rows(grade, &rows)?;
                let rows = self.check_quadratic_rows(field_type, rows, "rvValMatQuad")?;
                let ty = self.model.rv_valuated_matroid_type(grade)?;
                Ok((ty, TypedForm::ValMatQuad(grade, field, rows)))
            }
            ExprForm::ValFil(name, cutoff) => {
                let module = self.model.rv_phodge(&name)?;
                let value = module.valuated_filtration(self.model, cutoff)?;
                let ty = self.model.rv_valuated_matroid_type(value.grade)?;
                Ok((ty, TypedForm::Constant(Value::RVValMat(value))))
            }
            ExprForm::ValDual(input) => {
                let input = self.expr(*input)?;
                let Kind::RVValMat(grade) = *self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, "rvValDual", "an RVValMat"));
                };
                let ty = self.model.rv_valuated_matroid_type(grade.dual())?;
                Ok((ty, TypedForm::ValDual(Box::new(input))))
            }
            ExprForm::QuantumBasis(rank, ambient, shape) => {
                let grade = rv_qgr::Grade::checked(rank, ambient)?;
                let ty = self.model.quantum_grassmannian_type(grade)?;
                let value = rv_qgr::basis(self.model, grade, shape)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::RVCTInitial(name, p, s) => {
                let grade = self.model.rv_ct_state_grade(name, p, s)?;
                let ty = self.model.rv_ct_state_type(grade.clone())?;
                Ok((ty, TypedForm::PeriodInitial(grade)))
            }
            ExprForm::RVCTQuery(name, p, s, index) => {
                let grade = self.model.rv_ct_state_grade(name, p, s)?;
                let ty = self.model.rv_ct_residue_type(grade.level)?;
                Ok((ty, TypedForm::PeriodQuery(grade, index)))
            }
            ExprForm::RVCTLiteral(p, s, value) => {
                let level = PeriodLevel::new(self.model, p, s)?;
                let ty = self.model.rv_ct_residue_type(level)?;
                Ok((
                    ty,
                    TypedForm::Constant(Value::PeriodResidue(value % level.modulus)),
                ))
            }
            ExprForm::RVCTSection(digit, input) => {
                let input = self.expr(*input)?;
                let Kind::RVCTState(grade) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, "rvCTSection", "an RVCTState"));
                };
                let digit = grade.level.check_digit(digit)?;
                Ok((
                    input.ty,
                    TypedForm::PeriodSection(grade.clone(), digit, Box::new(input)),
                ))
            }
            ExprForm::RVCTRead(input) => {
                let input = self.expr(*input)?;
                let Kind::RVCTState(grade) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, "rvCTRead", "an RVCTState"));
                };
                let grade = grade.clone();
                let ty = self.model.rv_ct_residue_type(grade.level)?;
                Ok((ty, TypedForm::PeriodRead(grade, Box::new(input))))
            }
            ExprForm::OreTerm(q, theta, coefficient) => {
                let ty = self.model.ore_type()?;
                let value = rv_ore::term(self.model, q, theta, coefficient)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::OreUnary(operation, input) => {
                let input = self.expr(*input)?;
                if !matches!(self.model.kind(input.ty), Kind::RVOre) {
                    return Err(self.invalid_operand(&input, operation.name(), "an Ore operator"));
                }
                Ok((input.ty, TypedForm::OreUnary(operation, Box::new(input))))
            }
            ExprForm::OreProduct(arguments) => {
                let [a, b] = *arguments;
                let a = self.expr(a)?;
                if !matches!(self.model.kind(a.ty), Kind::RVOre) {
                    return Err(self.invalid_operand(&a, "rvOreMul", "an Ore operator"));
                }
                let b = self.expr(b)?;
                self.require_type(a.ty, &b, "rvOreMul")?;
                Ok((a.ty, TypedForm::OreProduct(Box::new([a, b]))))
            }
            ExprForm::OreApply(arguments) => {
                let [a, b] = *arguments;
                let a = self.expr(a)?;
                if !matches!(self.model.kind(a.ty), Kind::RVOre) {
                    return Err(self.invalid_operand(&a, "rvOreApply", "an Ore operator"));
                }
                let b = self.expr(b)?;
                let Kind::RVQGr(grade) = *self.model.kind(b.ty) else {
                    return Err(self.invalid_operand(
                        &b,
                        "rvOreApply",
                        "a quantum Grassmannian element",
                    ));
                };
                Ok((b.ty, TypedForm::OreApply(grade, Box::new([a, b]))))
            }
            ExprForm::QuantumReduce(p, input) => {
                let input = self.expr(*input)?;
                let Kind::RVQGr(grassmannian) = *self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(
                        &input,
                        "rvQGrMod",
                        "an integral quantum Grassmannian element",
                    ));
                };
                let prime = rv_pcurv::checked_prime(self.model, p)?;
                let grade = rv_pcurv::ModGrade {
                    prime,
                    grassmannian,
                };
                let ty = self.model.quantum_mod_type(grade)?;
                Ok((ty, TypedForm::QuantumReduce(grade, Box::new(input))))
            }
            ExprForm::QuantumPSection(p, k, n, shape, column) => {
                let grade = rv_pcurv::ModGrade::checked(self.model, p, k, n)?;
                let ty = self.model.quantum_mod_type(grade)?;
                let value = if column {
                    rv_pcurv::column(self.model, grade, shape)?
                } else {
                    rv_pcurv::basis(self.model, grade, shape)?
                };
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::QuantumCurvature(satake, input) => {
                let input = self.expr(*input)?;
                let Kind::RVQGrP(grade) = *self.model.kind(input.ty) else {
                    let op = if satake {
                        "rvQGrCurvatureSatake"
                    } else {
                        "rvQGrCurvature"
                    };
                    return Err(self.invalid_operand(
                        &input,
                        op,
                        "a finite-field Schubert section",
                    ));
                };
                Ok((
                    input.ty,
                    TypedForm::QuantumCurvature(grade, satake, Box::new(input)),
                ))
            }
            ExprForm::QuantumNabla(input) => {
                let input = self.expr(*input)?;
                if let Kind::RVQGrP(grade) = *self.model.kind(input.ty) {
                    return Ok((input.ty, TypedForm::QuantumPNabla(grade, Box::new(input))));
                }
                let Kind::RVQGr(grade) = *self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(
                        &input,
                        "rvQGrNabla",
                        "a quantum Grassmannian element",
                    ));
                };
                Ok((input.ty, TypedForm::QuantumNabla(grade, Box::new(input))))
            }
            ExprForm::QuantumAnnihilator(rank, ambient) => {
                let grade = rv_qgr::Grade::checked(rank, ambient)?;
                let ty = self.model.ore_type()?;
                let value = rv_qconnection::annihilator(self.model, grade)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::QuantumProduct(arguments) => {
                let [left, right] = *arguments;
                let left = self.expr(left)?;
                let Kind::RVQGr(grade) = *self.model.kind(left.ty) else {
                    return Err(self.invalid_operand(
                        &left,
                        "rvQGrMul",
                        "a quantum Grassmannian element",
                    ));
                };
                let right = self.expr(right)?;
                self.require_type(left.ty, &right, "rvQGrMul")?;
                Ok((
                    left.ty,
                    TypedForm::QuantumProduct(grade, Box::new([left, right])),
                ))
            }
            ExprForm::QuantumUnary(operation, input) => {
                let input = self.expr(*input)?;
                if let Kind::RVQGrP(grade) = *self.model.kind(input.ty) {
                    let ty = if matches!(operation, QuantumOp::Transpose) {
                        let grassmannian = grade.grassmannian.complement();
                        self.model.quantum_mod_type(rv_pcurv::ModGrade {
                            grassmannian,
                            ..grade
                        })?
                    } else {
                        input.ty
                    };
                    return Ok((
                        ty,
                        TypedForm::QuantumPUnary(grade, operation, Box::new(input)),
                    ));
                }
                let Kind::RVQGr(grade) = *self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(
                        &input,
                        operation.name(),
                        "a quantum Grassmannian element",
                    ));
                };
                let ty = if matches!(operation, QuantumOp::Transpose) {
                    self.model.quantum_grassmannian_type(grade.complement())?
                } else {
                    input.ty
                };
                Ok((ty, TypedForm::QuantumUnary(operation, Box::new(input))))
            }
            ExprForm::Hive(weights, interior) => {
                let boundary = rv_hive::Boundary::checked(self.model, weights)?;
                let ty = self.model.hive_type(boundary.clone())?;
                let value = if let Some(interior) = interior {
                    rv_hive::literal(self.model, &boundary, interior)?
                } else {
                    self.model.hive_witness(ty)?
                };
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SphericalBasis(coweight) => {
                let ty = self.model.sph_hecke_type(coweight.len())?;
                let value = rv_sphhecke::basis(self.model, coweight)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SphericalScalar(rank, degree) => {
                let ty = self.model.sph_hecke_type(rank)?;
                let value = rv_sphhecke::scalar(self.model, rank, degree)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SphericalProduct(arguments) => {
                let [left, right] = *arguments;
                let left = self.expr(left)?;
                if !matches!(self.model.kind(left.ty), Kind::RVSphHecke(..)) {
                    return Err(self.invalid_operand(
                        &left,
                        "rvSphMul",
                        "a spherical Hecke element",
                    ));
                }
                let right = self.expr(right)?;
                self.require_type(left.ty, &right, "rvSphMul")?;
                Ok((
                    left.ty,
                    TypedForm::SphericalProduct(Box::new([left, right])),
                ))
            }
            ExprForm::SatakeInverse(input) => {
                let input = self.expr(*input)?;
                let Kind::RVSatake(rank) = *self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(
                        &input,
                        "rvSatakeInverse",
                        "a Satake polynomial",
                    ));
                };
                let ty = self.model.sph_hecke_type(rank)?;
                Ok((ty, TypedForm::SatakeInverse(Box::new(input))))
            }
            ExprForm::SatakeCharacter(intersection, coweight) => {
                let ty = if intersection {
                    self.model.sph_hecke_type(coweight.len())?
                } else {
                    self.model.satake_type(coweight.len())?
                };
                let value = rv_satake::schur(self.model, &coweight, intersection)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SatakeTransform(input) => {
                let input = self.expr(*input)?;
                let Kind::RVSphHecke(rank) = *self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(
                        &input,
                        "rvSatake",
                        "a spherical Hecke element",
                    ));
                };
                let ty = self.model.satake_type(rank)?;
                Ok((ty, TypedForm::SatakeTransform(Box::new(input))))
            }
            ExprForm::SatakeScalar(rank, power) => {
                let ty = self.model.satake_type(rank)?;
                let value = rv_satake::scalar(self.model, rank, power)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::SatakeProduct(arguments) => {
                let [left, right] = *arguments;
                let left = self.expr(left)?;
                if !matches!(self.model.kind(left.ty), Kind::RVSatake(..)) {
                    return Err(self.invalid_operand(&left, "rvSatakeMul", "a Satake polynomial"));
                }
                let right = self.expr(right)?;
                self.require_type(left.ty, &right, "rvSatakeMul")?;
                Ok((left.ty, TypedForm::SatakeProduct(Box::new([left, right]))))
            }
            ExprForm::SphericalEmbedding(input) => {
                let input = self.expr(*input)?;
                let Kind::RVSphHecke(rank) = *self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(
                        &input,
                        "rvSphIwahori",
                        "a spherical Hecke element",
                    ));
                };
                let ty = self.model.iw_hecke_type(rank)?;
                Ok((ty, TypedForm::SphericalEmbedding(Box::new(input))))
            }
            ExprForm::HeckeProduct(arguments) => {
                let [left, right] = *arguments;
                let left = self.expr(left)?;
                if !matches!(self.model.kind(left.ty), Kind::RVIwHecke(..)) {
                    return Err(self.invalid_operand(
                        &left,
                        "rvHeckeMul",
                        "an affine Hecke element",
                    ));
                }
                let right = self.expr(right)?;
                self.require_type(left.ty, &right, "rvHeckeMul")?;
                Ok((left.ty, TypedForm::HeckeProduct(Box::new([left, right]))))
            }
            ExprForm::RVLaurentMonomial(power, input) => {
                let input = self.expr(*input)?;
                let Kind::RVHom(grade) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, "rvU", "an RVHom map"));
                };
                let ty = self.model.rv_laurent_type(grade.clone())?;
                Ok((ty, TypedForm::RVLaurentMonomial(power, Box::new(input))))
            }
            ExprForm::RVIwAdlv(boundary, modulus, root, window, arguments) => {
                let [a, x] = *arguments;
                let a = self.expr(a)?;
                let x = self.expr(x)?;
                let step = rv_hodge::relative_step(
                    self.model,
                    a.ty,
                    x.ty,
                    modulus,
                    root,
                    boundary.name(),
                )?;
                rv_iwahori::check_window(self.model, a.ty, &window)?;
                let ty = self.model.named("S").unwrap();
                Ok((
                    ty,
                    TypedForm::RVIwAdlv(boundary, step, window, Box::new([a, x])),
                ))
            }
            ExprForm::RVAdlv(boundary, modulus, root, mu, arguments) => {
                let [a, x] = *arguments;
                let a = self.expr(a)?;
                let x = self.expr(x)?;
                let step = rv_hodge::relative_step(
                    self.model,
                    a.ty,
                    x.ty,
                    modulus,
                    root,
                    boundary.name(),
                )?;
                rv_hodge::check_coweight(self.model, a.ty, &mu)?;
                let ty = self.model.named("S").unwrap();
                Ok((ty, TypedForm::RVAdlv(boundary, step, mu, Box::new([a, x]))))
            }
            ExprForm::RVNewton(operation, modulus, root, input) => {
                let input = self.expr(*input)?;
                let (ty, embedding) =
                    rv_newton::plan(self.model, input.ty, modulus, root, operation)?;
                Ok((
                    ty,
                    TypedForm::RVNewton(operation, embedding, Box::new(input)),
                ))
            }
            ExprForm::RVScale(scalar, input) => {
                let input = self.expr(*input)?;
                let field = match self.model.kind(input.ty) {
                    Kind::RVHom(grade) | Kind::RVLaurent(grade) => grade.field,
                    Kind::JetHom(grade) => self.model.jet_path(grade.diagram.path).field,
                    Kind::JetPage(grade) => self.model.jet_path(grade.chain.diagram.path).field,
                    _ => {
                        return Err(self.invalid_operand(
                            &input,
                            "rvScale",
                            "an RVHom or RVLaurent map, JetHom cochain, or JetPage class",
                        ))
                    }
                };
                field.check_element(scalar)?;
                Ok((input.ty, TypedForm::RVScale(field, scalar, Box::new(input))))
            }
            ExprForm::RVField(operation, input) => {
                let input = self.expr(*input)?;
                if let (Kind::RVLaurent(grade), rv_field::Operation::Frobenius(iterations)) =
                    (self.model.kind(input.ty), operation)
                {
                    let field = grade.field;
                    let ty = rv_laurent::frobenius_type(self.model, grade.clone(), iterations)?;
                    return Ok((
                        ty,
                        TypedForm::RVLaurentFrobenius(field, iterations, Box::new(input)),
                    ));
                }
                let Kind::RVHom(grade) = self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, operation.name(), "an RVHom map"));
                };
                let (ty, map) = rv_field::plan(self.model, grade.clone(), operation)?;
                Ok((ty, TypedForm::RVField(map, Box::new(input))))
            }
            ExprForm::RVInverse(input) => {
                let input = self.expr(*input)?;
                let ty = rv_descent::inverse_type(self.model, input.ty)?;
                Ok((ty, TypedForm::RVInverse(Box::new(input))))
            }
            ExprForm::RVHilbert90(modulus, root, input) => {
                let input = self.expr(*input)?;
                let embedding = rv_descent::hilbert_embedding(self.model, input.ty, modulus, root)?;
                Ok((input.ty, TypedForm::RVHilbert90(embedding, Box::new(input))))
            }
            ExprForm::RVTwisted(operation, modulus, root, arguments) => {
                let [left, input, right] = *arguments;
                let left = self.expr(left)?;
                let input = self.expr(input)?;
                let right = self.expr(right)?;
                let types = [left.ty, input.ty, right.ty];
                let (ty, embedding) =
                    rv_descent::twist_plan(self.model, operation, types, modulus, root)?;
                Ok((
                    ty,
                    TypedForm::RVTwisted(operation, embedding, Box::new([left, input, right])),
                ))
            }
            ExprForm::RVShintani(operation, modulus, root, input) => {
                let input = self.expr(*input)?;
                let (ty, embedding) =
                    rv_shintani::plan(self.model, input.ty, modulus, root, operation)?;
                Ok((
                    ty,
                    TypedForm::RVShintani(operation, embedding, Box::new(input)),
                ))
            }
            ExprForm::RVSpecialize(modulus, point, input) => {
                let input = self.expr(*input)?;
                let ty = rv_hom::specialize_type(self.model, input.ty, modulus, point)?;
                Ok((ty, TypedForm::RVSpecialize(Box::new(input))))
            }
            _ => unreachable!(),
        }
    }

    fn expr_bp(&mut self, expr: ExprForm, span: Span) -> Result<(TypeId, TypedForm)> {
        match expr {
            ExprForm::BpInteger(spec, integer) => {
                let shape = bp::Shape::new(self.model, spec.0, spec.1, spec.2)?;
                let ty = self.model.bp_type(shape)?;
                let value = bp::persist(self.model, BpPolynomial::integer(integer))?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::BpVariable(spec, slot, index) => {
                let shape = bp::Shape::new(self.model, spec.0, spec.1, spec.2)?;
                if slot == Some(0) {
                    return Err(Error::new(span, "BP arrow slots start at 1"));
                }
                let ty = self.model.bp_type(shape)?;
                let polynomial = shape.variable(slot.unwrap_or(0), index)?;
                let value = bp::persist(self.model, polynomial)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::BpVPageGenerator(spec, index, at) => {
                let grade = self.bp_vpage_grade(spec)?;
                let ty = self.model.bp_vpage_type(grade)?;
                let rank = self.model.bp_vpage(ty)?.rank();
                let value = self.prime_generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::BpVPageClass((page, t, filtration), input) => {
                let input = self.expr(*input)?;
                let shape = self.bp_shape(&input, "bpVClass")?;
                if shape.depth == 0 || shape.thickness.is_some() || shape.depth > shape.height {
                    return Err(self.invalid_operand(
                        &input,
                        "bpVClass",
                        "a BPI cochain with an invariant generator v_h",
                    ));
                }
                let grade = bp_vpage::Grade::new(shape, page, t, filtration)?;
                let ty = self.model.bp_vpage_type(grade)?;
                Ok((ty, TypedForm::BpVPageClass(Box::new(input))))
            }
            ExprForm::BpVPageOperation(name, operation, input) => {
                let input = self.expr(*input)?;
                let Kind::BPVPage(grade) = *self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, &name.text, "a BPVPage class"));
                };
                let ty = if matches!(operation, bp_vpage::Operation::Rep) {
                    self.model.bp_type(grade.shape)?
                } else {
                    self.model.bp_vpage_type(operation.target(grade)?)?
                };
                Ok((ty, TypedForm::BpVPageOperation(operation, Box::new(input))))
            }
            ExprForm::BpPageGenerator(mode, spec, index, at) => {
                let grade = self.bp_page_grade(mode, spec)?;
                let ty = self.model.bp_page_type(grade)?;
                let rank = self.model.bp_page(ty)?.rank();
                let value = self.prime_generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::BpPageClass(mode, (page, t, filtration), input) => {
                let input = self.expr(*input)?;
                let name = match mode {
                    bp_page::Mode::Bockstein => "bpBockClass",
                    bp_page::Mode::Novikov => "bpNovClass",
                };
                let shape = self.bp_shape(&input, name)?;
                if shape.depth != 0 {
                    return Err(self.invalid_operand(&input, name, "an integral BP cochain"));
                }
                let grade = bp_page::Grade {
                    shape,
                    mode,
                    page,
                    t,
                    filtration,
                };
                let ty = self.model.bp_page_type(grade)?;
                Ok((ty, TypedForm::BpPageClass(Box::new(input))))
            }
            ExprForm::BpPageOperation(name, mode, operation, input) => {
                let input = self.expr(*input)?;
                let Kind::BPPage(grade) = *self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, &name.text, mode.type_name()));
                };
                if grade.mode != mode {
                    return Err(self.invalid_operand(&input, &name.text, mode.type_name()));
                }
                let ty = if matches!(operation, bp_page::Operation::Representative) {
                    self.model.bp_type(grade.shape)?
                } else {
                    self.model.bp_page_type(operation.target(grade)?)?
                };
                Ok((ty, TypedForm::BpPageOperation(operation, Box::new(input))))
            }
            ExprForm::BpModExtGen(spec, index, at) => {
                let grade = self.bp_mod_ext_grade(spec)?;
                let ty = self.model.bp_mod_ext_type(grade)?;
                let rank = self.model.bp_mod_ext(ty)?.rank();
                let value = self.prime_generator(ty, rank, index, at)?;
                Ok((ty, TypedForm::Constant(value)))
            }
            ExprForm::BpModExtClass(thick, t, input) => {
                let (name, expected) = if thick {
                    ("bpThickExtClass", "a BPThick cochain")
                } else {
                    ("bpIExtClass", "a BPI cochain")
                };
                let input = self.expr(*input)?;
                let shape = self.bp_shape(&input, name)?;
                if shape.depth == 0 || shape.thickness.is_some() != thick {
                    return Err(self.invalid_operand(&input, name, expected));
                }
                let ty = self.model.bp_mod_ext_type(bp_mod_ext::Grade { shape, t })?;
                Ok((ty, TypedForm::BpModExtClass(Box::new(input))))
            }
            ExprForm::BpModExtRep(thick, input) => {
                let (name, expected) = if thick {
                    ("bpThickExtRep", "a BPThickExt class")
                } else {
                    ("bpIExtRep", "a BPIExt class")
                };
                let input = self.expr(*input)?;
                let Kind::BPModExt(grade) = *self.model.kind(input.ty) else {
                    return Err(self.invalid_operand(&input, name, expected));
                };
                if grade.shape.thickness.is_some() != thick {
                    return Err(self.invalid_operand(&input, name, expected));
                }
                let ty = self.model.bp_type(grade.shape)?;
                Ok((ty, TypedForm::BpModExtRep(Box::new(input))))
            }
            ExprForm::BpRaise(height, input) => {
                let input = self.expr(*input)?;
                let ty = match *self.model.kind(input.ty) {
                    Kind::BP(shape) => self.model.bp_type(shape.raise(height)?)?,
                    Kind::BPPage(mut grade) => {
                        grade.shape = grade.shape.raise(height)?;
                        self.model.bp_page_type(grade)?
                    }
                    Kind::BPExt(mut grade) => {
                        grade.shape = grade.shape.raise(height)?;
                        self.model.bp_ext_type(grade)?
                    }
                    _ => {
                        return Err(self.invalid_operand(
                            &input,
                            "bpRaise",
                            "a BP cochain, page class, or BPExt class",
                        ))
                    }
                };
                Ok((ty, TypedForm::BpRaise(Box::new(input))))
            }
            ExprForm::BpTriple(inputs) => {
                let (layout, inputs) = self.bp_triple_inputs(*inputs, "bpTriple")?;
                let ty = self.model.bp_triple_type(layout.target)?;
                Ok((ty, TypedForm::BpTriple(layout, inputs)))
            }
            ExprForm::BpTripleBase(input) => {
                let input = self.expr(*input)?;
                let grade = self.bp_triple_grade(&input, "bpTripleBase")?;
                let ty = self.model.bp_ext_type(grade)?;
                Ok((ty, TypedForm::BpTripleBase(Box::new(input))))
            }
            ExprForm::BpHigher(arity, mut inputs) => {
                let (operation, witnesses) = if let Some(arity) = arity {
                    ("bpHigherSystem", Some(inputs.split_off(arity)))
                } else {
                    ("bpChosenHigher", None)
                };
                let inputs = inputs
                    .into_iter()
                    .map(|input| self.expr(input))
                    .collect::<Result<Vec<_>>>()?;
                let grades = inputs
                    .iter()
                    .map(|input| self.bp_ext_operand(input, operation))
                    .collect::<Result<Vec<_>>>()?;
                let layout = bp_higher::Layout::new(self.model, grades)?;
                let witnesses = witnesses
                    .map(|witnesses| {
                        witnesses
                            .into_iter()
                            .zip(layout.witness_grades())
                            .map(|(witness, grade)| {
                                let witness = self.expr(witness)?;
                                let ty = self.model.bp_type(grade.shape)?;
                                self.require_type(ty, &witness, operation)?;
                                Ok(witness)
                            })
                            .collect::<Result<Vec<_>>>()
                    })
                    .transpose()?;
                let ty = self.model.bp_ext_type(layout.target)?;
                Ok((ty, TypedForm::BpHigher(layout, inputs, witnesses)))
            }
            ExprForm::BpMassey(inputs) => {
                let inputs = inputs
                    .into_iter()
                    .map(|input| self.expr(input))
                    .collect::<Result<Vec<_>>>()?;
                let shapes = inputs
                    .iter()
                    .map(|input| self.bp_shape(input, "bpMassey"))
                    .collect::<Result<Vec<_>>>()?;
                let layout = bp::Massey::new(shapes.try_into().unwrap_or_else(|_| unreachable!()))?;
                Ok((
                    self.model.bp_type(layout.target)?,
                    TypedForm::BpMassey(layout, inputs),
                ))
            }
            ExprForm::BpHeightGenerator(spec, index, at) => {
                let grade = self.bp_height_grade(spec)?;
                let ty = self.model.bp_height_type(grade)?;
                let rank = self.model.bp_ext(ty)?.rank();
                if index >= rank {
                    return Err(Error::new(
                        at,
                        format!("bpHeightgen index {index} is out of range for rank {rank}"),
                    ));
                }
                Ok((
                    ty,
                    TypedForm::Constant(bp_classes::generator(self.model, index)?),
                ))
            }
            ExprForm::BpHeightOperation(operation, input) => {
                let input = self.expr(*input)?;
                let ty = if let BpHeightOp::Class(r, h, t) = operation {
                    let top = self.bp_shape(&input, operation.name())?;
                    top.require_integral(operation.name())?;
                    self.model
                        .bp_height_type(bp_height::ClassGrade { top, r, h, t })?
                } else {
                    let Kind::BPHeight(grade) = *self.model.kind(input.ty) else {
                        return Err(self.invalid_operand(
                            &input,
                            operation.name(),
                            "a BPHeight class",
                        ));
                    };
                    match operation {
                        BpHeightOp::Representative => self.model.bp_type(grade.top)?,
                        BpHeightOp::Differential => self.model.bp_height_type(grade.target()?)?,
                        BpHeightOp::Next => self.model.bp_height_type(grade.next_page()?)?,
                        BpHeightOp::Previous => {
                            self.model.bp_height_type(grade.previous_page()?)?
                        }
                        BpHeightOp::Class(..) => unreachable!(),
                    }
                };
                Ok((ty, TypedForm::BpHeightOperation(operation, Box::new(input))))
            }
            ExprForm::BpRelOperation(operation, input) => {
                let input = self.expr(*input)?;
                let (ty, shape) = match (operation, self.model.kind(input.ty).clone()) {
                    (BpRelOp::Quotient(lower), Kind::BPExt(grade)) => {
                        let shape = bp_rel::Shape::with_lower(grade.shape, lower)?;
                        (
                            self.model
                                .bp_rel_ext_type(bp_rel::Grade { shape, t: grade.t })?,
                            shape,
                        )
                    }
                    (BpRelOp::Quotient(lower), _) => {
                        let upper = self.bp_shape(&input, "bpRel")?;
                        let shape = bp_rel::Shape::with_lower(upper, lower)?;
                        (self.model.bp_rel_type(shape)?, shape)
                    }
                    (BpRelOp::Lift, Kind::BPRel(shape)) => {
                        (self.model.bp_type(shape.upper)?, shape)
                    }
                    (BpRelOp::Class(t), Kind::BPRel(shape)) => (
                        self.model.bp_rel_ext_type(bp_rel::Grade { shape, t })?,
                        shape,
                    ),
                    (BpRelOp::Representative, Kind::BPRelExt(grade)) => {
                        (self.model.bp_rel_type(grade.shape)?, grade.shape)
                    }
                    (BpRelOp::HeightDelta, Kind::BPRel(shape)) => {
                        (self.model.bp_type(shape.next()?.lower_shape())?, shape)
                    }
                    (BpRelOp::HeightDelta, Kind::BPRelExt(grade)) => {
                        let shape = grade.shape.next()?.lower_shape();
                        let ty = self
                            .model
                            .bp_ext_type(bp_classes::Grade { shape, t: grade.t })?;
                        (ty, grade.shape)
                    }
                    (BpRelOp::Lift, _) => {
                        return Err(self.invalid_operand(&input, "bpRelLift", "a BPRel cochain"))
                    }
                    (BpRelOp::Class(_), _) => {
                        return Err(self.invalid_operand(&input, "bpRelClass", "a BPRel cochain"))
                    }
                    (BpRelOp::Representative, _) => {
                        return Err(self.invalid_operand(&input, "bpRelRep", "a BPRelExt class"))
                    }
                    (BpRelOp::HeightDelta, _) => {
                        return Err(self.invalid_operand(
                            &input,
                            "bpHeightDelta",
                            "a BPRel cochain or BPRelExt class",
                        ))
                    }
                    _ => unreachable!(),
                };
                Ok((
                    ty,
                    TypedForm::BpRelOperation(operation, shape, Box::new(input)),
                ))
            }
            ExprForm::BpRelGenerator(spec, t, index, at) => {
                let shape = self.bp_relative_shape(spec)?;
                let ty = self.model.bp_rel_ext_type(bp_rel::Grade { shape, t })?;
                let rank = self.model.bp_ext(ty)?.rank();
                if index >= rank {
                    return Err(Error::new(
                        at,
                        format!("bprelgen index {index} is out of range for rank {rank}"),
                    ));
                }
                Ok((
                    ty,
                    TypedForm::Constant(bp_classes::generator(self.model, index)?),
                ))
            }
            ExprForm::BpExtGenerator(spec, t, index, at) => {
                let grade = self.bp_ext_grade(spec, t)?;
                let ty = self.model.bp_ext_type(grade)?;
                let rank = self.model.bp_ext(ty)?.rank();
                if index >= rank {
                    return Err(Error::new(
                        at,
                        format!("bpextgen index {index} is out of range for rank {rank}"),
                    ));
                }
                Ok((
                    ty,
                    TypedForm::Constant(bp_classes::generator(self.model, index)?),
                ))
            }
            ExprForm::BpExtOperation(operation, input) => {
                let input = self.expr(*input)?;
                let ty = match operation {
                    BpExtOp::Class(t) => {
                        let shape = self.bp_shape(&input, "bpExtClass")?;
                        shape.require_integral("bpExtClass")?;
                        self.model.bp_ext_type(bp_classes::Grade { shape, t })?
                    }
                    BpExtOp::Representative => {
                        let grade = self.bp_ext_operand(&input, "bpExtRep")?;
                        self.model.bp_type(grade.shape)?
                    }
                    BpExtOp::Scale(_) => {
                        if !matches!(
                            self.model.kind(input.ty),
                            Kind::BPExt(..)
                                | Kind::BPRel(..)
                                | Kind::BPRelExt(..)
                                | Kind::BPHeight(..)
                        ) {
                            return Err(self.invalid_operand(
                                &input,
                                "bpScale",
                                "a BPExt, BPRelExt, or BPHeight class or a BPRel cochain",
                            ));
                        }
                        input.ty
                    }
                };
                Ok((ty, TypedForm::BpExtOperation(operation, Box::new(input))))
            }
            ExprForm::BpOperation(operation, input) => {
                let input = self.expr(*input)?;
                let kind = self.model.kind(input.ty);
                if let Kind::BPModExt(grade) = *kind {
                    let grade = grade.map(operation)?;
                    let ty = if grade.shape.depth == 0 {
                        self.model.bp_ext_type(bp_classes::Grade {
                            shape: grade.shape,
                            t: grade.t,
                        })?
                    } else {
                        self.model.bp_mod_ext_type(grade)?
                    };
                    return Ok((ty, TypedForm::BpModExtOperation(operation, Box::new(input))));
                }
                if matches!(operation, bp::Operation::Negate)
                    && matches!(
                        kind,
                        Kind::BPExt(..) | Kind::BPRel(..) | Kind::BPRelExt(..) | Kind::BPHeight(..)
                    )
                {
                    return Ok((
                        input.ty,
                        TypedForm::BpExtOperation(BpExtOp::Scale(-1), Box::new(input)),
                    ));
                }
                if let Kind::BPRel(shape) = *kind {
                    let (target, operation) = match operation {
                        bp::Operation::Differential => (shape.next()?, BpRelOp::Differential),
                        bp::Operation::Null(t) => (shape.previous()?, BpRelOp::Null(t)),
                        _ => {
                            return Err(self.invalid_operand(
                                &input,
                                operation.name(),
                                "a BP polynomial",
                            ))
                        }
                    };
                    let ty = self.model.bp_rel_type(target)?;
                    return Ok((
                        ty,
                        TypedForm::BpRelOperation(operation, shape, Box::new(input)),
                    ));
                }
                let shape = self.bp_shape(&input, operation.name())?;
                let ty = self.model.bp_type(operation.target(shape)?)?;
                Ok((
                    ty,
                    TypedForm::BpOperation(operation, shape, Box::new(input)),
                ))
            }
            _ => unreachable!(),
        }
    }

    fn class_type(&self, value: &TypedExpr, op: &str) -> Result<(TypeId, usize)> {
        match *self.model.kind(value.ty) {
            Kind::Cohomology(base, degree) => Ok((base, degree)),
            _ => Err(Error::new(
                value.span,
                format!(
                    "{op} requires a cohomology class, found {}",
                    self.model.type_name(value.ty),
                ),
            )),
        }
    }

    fn unbound(&self, name: &str) -> bool {
        !self.bindings.contains_key(name) && !self.locals.iter().any(|(local, _)| local == name)
    }

    fn map_bracket(
        &mut self,
        operation: map_bracket::Operation,
        mut inputs: Vec<ParsedExpr>,
    ) -> Result<(TypeId, TypedForm)> {
        operation.require_arity(inputs.len())?;
        let query = match operation {
            map_bracket::Operation::Predicate => None,
            map_bracket::Operation::Membership => Some(inputs.pop().unwrap()),
        };
        let inputs = inputs
            .into_iter()
            .map(|input| self.expr(input))
            .collect::<Result<Vec<_>>>()?;
        let target = map_bracket::result_type(
            self.model,
            inputs.iter().map(|input| input.ty),
            operation.name(),
        )?;
        self.class_bracket(
            ClassBracketFamily::Map,
            target,
            inputs,
            query,
            operation.name(),
        )
    }

    fn rv_t_bracket(
        &mut self,
        operation: rv_t_bracket::Operation,
        mut inputs: Vec<ParsedExpr>,
    ) -> Result<(TypeId, TypedForm)> {
        operation.require_arity(inputs.len())?;
        let query = match operation {
            rv_t_bracket::Operation::Predicate => None,
            rv_t_bracket::Operation::Membership => inputs.pop(),
        };
        let inputs = inputs
            .into_iter()
            .map(|input| self.expr(input))
            .collect::<Result<Vec<_>>>()?;
        let [a, b, c] = inputs.as_slice() else {
            unreachable!();
        };
        let target = rv_t_bracket::result_type(self.model, [a.ty, b.ty, c.ty], operation.name())?;
        self.class_bracket(
            ClassBracketFamily::Koszul,
            target,
            inputs,
            query,
            operation.name(),
        )
    }

    fn class_bracket(
        &mut self,
        family: ClassBracketFamily,
        target: TypeId,
        inputs: Vec<TypedExpr>,
        query: Option<ParsedExpr>,
        name: &str,
    ) -> Result<(TypeId, TypedForm)> {
        let range = self.model.named("S").unwrap();
        let (ty, query) = if let Some(query) = query {
            let query = self.expr(query)?;
            self.require_type(target, &query, name)?;
            (range, Some(Box::new(query)))
        } else {
            (self.model.arrow(target, range)?, None)
        };
        Ok((
            ty,
            TypedForm::ClassBracket {
                family,
                target,
                inputs,
                query,
            },
        ))
    }

    fn transfer(
        &mut self,
        operation: transfer::Operation,
        inputs: Vec<ParsedExpr>,
    ) -> Result<(TypeId, TypedForm)> {
        let arity = inputs.len();
        let terms = inputs
            .into_iter()
            .map(|input| self.expr(input))
            .collect::<Result<Vec<_>>>()?;
        if terms
            .iter()
            .any(|term| matches!(self.model.kind(term.ty), Kind::RMExt(..)))
        {
            let ty = real_transfer::result_type(
                self.model,
                operation,
                terms.iter().map(|term| term.ty),
            )?;
            return Ok((ty, TypedForm::RealTransfer(operation, terms)));
        }
        operation.require_arity(arity)?;
        let (field, degrees) = self.scalar_degrees(operation.name(), &terms)?;
        let (s, t) = operation.bidegree(degrees, arity)?;
        let ty = match (field, operation) {
            (ScalarField::ModTwo(profile), transfer::Operation::Product) => {
                self.model.profile_ext_type(profile, s, t)?
            }
            (ScalarField::ModTwo(profile), transfer::Operation::Lift) => {
                self.model.profile_cobar_type(profile, s, t)?
            }
            (ScalarField::Odd(prime, profile), transfer::Operation::Product) => {
                self.model.profile_odd_ext_type(prime, profile, s, t)?
            }
            (ScalarField::Odd(prime, profile), transfer::Operation::Lift) => {
                self.model.profile_pcobar_type(prime, profile, s, t)?
            }
        };
        Ok((ty, TypedForm::Transfer(field, operation, terms)))
    }

    fn higher_restriction(
        &mut self,
        n: usize,
        inputs: Vec<ParsedExpr>,
        span: Span,
    ) -> Result<(TypeId, TypedForm)> {
        let arity = inputs.len();
        let terms = inputs
            .into_iter()
            .map(|input| self.expr(input))
            .collect::<Result<Vec<_>>>()?;
        let (field, degrees) = self.scalar_degrees("resHigher", &terms)?;
        let source = match field {
            ScalarField::ModTwo(profile) | ScalarField::Odd(_, profile) => profile,
        };
        let target = Profile::An(n);
        if !source.can_restrict_to(target) {
            return Err(Error::new(
                span,
                "resHigher cannot restrict from a smaller Steenrod algebra",
            ));
        }
        let (s, t) = restriction_transfer::bidegree(degrees, arity)?;
        let ty = match field {
            ScalarField::ModTwo(_) => self.model.profile_ext_type(target, s, t)?,
            ScalarField::Odd(prime, _) => self.model.profile_odd_ext_type(prime, target, s, t)?,
        };
        Ok((ty, TypedForm::HigherRestriction(field, target, terms)))
    }

    fn scalar_degrees(
        &self,
        operation: &str,
        terms: &[TypedExpr],
    ) -> Result<(ScalarField, Vec<transfer::Bidegree>)> {
        let mut degrees = Vec::with_capacity(terms.len());
        let mut source = None;
        for term in terms {
            let kind = self.model.kind(term.ty);
            let (field, s, t) = if let Some((prime, profile, s, t)) = kind.pext_grade() {
                (ScalarField::Odd(prime, profile), s, t)
            } else {
                let expected = if matches!(source, Some(ScalarField::Odd(..))) {
                    "a PExt class"
                } else {
                    "an Ext class"
                };
                let (profile, s, t) = kind
                    .ext_grade()
                    .ok_or_else(|| self.invalid_operand(term, operation, expected))?;
                (ScalarField::ModTwo(profile), s, t)
            };
            if let Some(previous) = source {
                match (previous, field) {
                    (ScalarField::ModTwo(p), ScalarField::ModTwo(q)) => {
                        require_profile(p, q, operation, term.span)?;
                    }
                    (ScalarField::Odd(p, profile), ScalarField::Odd(q, other)) => {
                        if p != q {
                            return Err(Error::new(
                                term.span,
                                format!("{operation} requires PExt classes at the same prime"),
                            ));
                        }
                        require_profile(profile, other, operation, term.span)?;
                    }
                    (ScalarField::Odd(..), ScalarField::ModTwo(_)) => {
                        return Err(self.invalid_operand(term, operation, "a PExt class"));
                    }
                    (ScalarField::ModTwo(_), ScalarField::Odd(..)) => {
                        return Err(self.invalid_operand(term, operation, "an Ext class"));
                    }
                }
            }
            source = Some(field);
            degrees.push((s, t));
        }
        let field = source
            .ok_or_else(|| Error::plain(format!("{operation} requires at least one Ext class")))?;
        Ok((field, degrees))
    }

    fn map_transfer(
        &mut self,
        operation: map_transfer::Operation,
        inputs: Vec<ParsedExpr>,
    ) -> Result<(TypeId, TypedForm)> {
        let arity = inputs.len();
        operation.require_arity(arity)?;
        let terms = inputs
            .into_iter()
            .map(|input| self.expr(input))
            .collect::<Result<Vec<_>>>()?;
        let ty = operation.result_type(self.model, terms.iter().map(|term| term.ty), arity)?;
        Ok((ty, TypedForm::MapTransfer(operation, terms)))
    }

    fn rv_t_transfer(
        &mut self,
        operation: rv_t_transfer::Operation,
        inputs: Vec<ParsedExpr>,
    ) -> Result<(TypeId, TypedForm)> {
        operation.require_arity(inputs.len())?;
        let inputs = inputs
            .into_iter()
            .map(|input| self.expr(input))
            .collect::<Result<Vec<_>>>()?;
        let ty =
            rv_t_transfer::result_type(self.model, operation, inputs.iter().map(|term| term.ty))?;
        Ok((ty, TypedForm::RvTTransfer(operation, inputs)))
    }

    fn rv_t_system(
        &mut self,
        operation: rv_t_system::Operation,
        inputs: Vec<ParsedExpr>,
        witnesses: Vec<Witness>,
        span: Span,
    ) -> Result<(TypeId, TypedForm)> {
        let arity = inputs.len();
        let inputs = inputs
            .into_iter()
            .map(|input| self.expr(input))
            .collect::<Result<Vec<_>>>()?;
        let signature = rv_t_system::Signature::new(
            self.model,
            operation,
            inputs.iter().map(|input| (input.ty, input.span)),
            span,
        )?;
        let ty = signature.result_type(self.model, span)?;
        let mut entries = Vec::with_capacity(arity + witnesses.len());
        for (from, expr) in inputs.into_iter().enumerate() {
            entries.push(DefiningEntry {
                from,
                to: from + 1,
                span: expr.span,
                expr,
            });
        }
        self.check_witnesses(witnesses, &mut entries, |model, from, to, at| {
            signature.witness_type(model, from, to, at)
        })?;
        Ok((ty, TypedForm::RvTSystem(operation, arity, entries)))
    }

    fn map_system(
        &mut self,
        operation: map_system::Operation,
        inputs: Vec<ParsedExpr>,
        witnesses: Vec<Witness>,
        span: Span,
    ) -> Result<(TypeId, TypedForm)> {
        let arity = inputs.len();
        let inputs = inputs
            .into_iter()
            .map(|input| self.expr(input))
            .collect::<Result<Vec<_>>>()?;
        let signature = map_system::Signature::new(
            self.model,
            operation,
            inputs.iter().map(|input| (input.ty, input.span)),
            span,
        )?;
        let ty = signature.result_type(self.model, span)?;
        let mut entries = Vec::with_capacity(arity + witnesses.len());
        for (index, expr) in inputs.into_iter().enumerate() {
            let to = arity - index;
            entries.push(DefiningEntry {
                from: to - 1,
                to,
                span: expr.span,
                expr,
            });
        }
        self.check_witnesses(witnesses, &mut entries, |model, from, to, at| {
            signature.witness_type(model, from, to, at)
        })?;
        Ok((ty, TypedForm::MapSystem(operation, arity, entries)))
    }

    fn odd_system(
        &mut self,
        operation: odd_massey::Operation,
        inputs: Vec<ParsedExpr>,
        witnesses: Vec<Witness>,
        span: Span,
    ) -> Result<(TypeId, TypedForm)> {
        let inputs = inputs
            .into_iter()
            .map(|input| self.expr(input))
            .collect::<Result<Vec<_>>>()?;
        let signature = odd_massey::Signature::new(
            self.model,
            operation,
            inputs.iter().map(|input| (input.ty, input.span)),
            span,
        )?;
        let ty = signature.result_type(self.model, span)?;
        let arity = signature.arity();
        let mut entries = Vec::with_capacity(inputs.len() + witnesses.len());
        for (index, expr) in inputs.into_iter().enumerate() {
            let to = index + 1;
            let from = if operation.explicit() { index } else { 0 };
            if !operation.explicit() && index > 0 {
                let expected = signature.witness_type(self.model, from, to, expr.span)?;
                self.require_type(expected, &expr, &format!("pPowerDef witness for gap {to}"))?;
            }
            entries.push(DefiningEntry {
                from,
                to,
                span: expr.span,
                expr,
            });
        }
        self.check_witnesses(witnesses, &mut entries, |model, from, to, at| {
            signature.witness_type(model, from, to, at)
        })?;
        Ok((ty, TypedForm::OddSystem(operation, arity, entries)))
    }

    fn sheaf_higher(
        &mut self,
        operation: SheafHigherOp,
        mut inputs: Vec<ParsedExpr>,
        witnesses: Vec<Witness>,
    ) -> Result<(TypeId, TypedForm)> {
        let target = if operation.targeted() {
            Some(inputs.pop().ok_or_else(|| {
                Error::plain(format!(
                    "{} requires a final ShExt outcome",
                    operation.name()
                ))
            })?)
        } else {
            None
        };
        let args = inputs
            .into_iter()
            .map(|value| self.expr(value))
            .collect::<Result<Vec<_>>>()?;
        let positions: Vec<_> = args.iter().map(|value| (value.ty, value.span)).collect();
        let signature = sheaf_higher::Signature::new(self.model, &positions)?;
        let target = if let Some(target) = target {
            let value = self.expr(target)?;
            let ty = signature.top_type(self.model, SheafForm::Class)?;
            self.require_type(ty, &value, operation.name())?;
            Some(Box::new(value))
        } else {
            None
        };
        let ty = match operation {
            SheafHigherOp::Curvature => signature.top_type(self.model, SheafForm::Cochain)?,
            SheafHigherOp::Contains => self.model.named("S").unwrap(),
            SheafHigherOp::Witness { from, to, .. } => {
                signature.witness_type(self.model, from, to)?
            }
            _ => signature.top_type(self.model, SheafForm::Class)?,
        };
        let mut entries = Vec::new();
        self.check_witnesses(witnesses, &mut entries, |model, from, to, span| {
            signature
                .witness_type(model, from, to)
                .map_err(|error| error.at(span))
        })?;
        Ok((
            ty,
            TypedForm::SheafHigher(operation, signature, args, entries, target),
        ))
    }

    fn check_witnesses(
        &mut self,
        witnesses: Vec<Witness>,
        entries: &mut Vec<DefiningEntry>,
        mut witness_type: impl FnMut(&mut Model, usize, usize, Span) -> Result<TypeId>,
    ) -> Result<()> {
        let mut supplied = HashSet::new();
        for Witness {
            from,
            to,
            span,
            value,
        } in witnesses
        {
            let expected = witness_type(self.model, from, to, span)?;
            if !supplied.insert((from, to)) {
                return Err(Error::new(
                    span,
                    format!("duplicate witness for interval [{from},{to}]"),
                ));
            }
            let expr = self.expr(value)?;
            self.require_type(
                expected,
                &expr,
                &format!("witness for interval [{from},{to}]"),
            )?;
            entries.push(DefiningEntry {
                from,
                to,
                span,
                expr,
            });
        }
        Ok(())
    }

    fn vbar_family(&mut self, spec: BarFamilySpec) -> Result<v_bars::Family> {
        let source = self.q_endpoint(spec.source)?;
        let target = self.q_endpoint(spec.target)?;
        let selected = q_koszul::indices(self.model, spec.indices)?;
        v_bars::Family::new(source, target, &selected)
    }

    fn bracket(&mut self, inputs: Vec<ParsedExpr>, span: Span) -> Result<(TypeId, TypedForm)> {
        if inputs.len() < 2 {
            return Err(Error::new(
                span,
                "bracket requires at least two cohomology classes",
            ));
        }
        let mut terms = Vec::with_capacity(inputs.len());
        let mut base = None;
        let mut degree = 2usize;
        for expression in inputs {
            let term = self.expr(expression)?;
            let (space, input_degree) = self.class_type(&term, "bracket")?;
            if let Some(expected) = base {
                if expected != space {
                    return Err(Error::new(
                        term.span,
                        format!(
                            "bracket requires classes on the same space, found {} and {}",
                            self.model.type_name(expected),
                            self.model.type_name(space),
                        ),
                    ));
                }
            }
            let reduced = input_degree.checked_sub(1).ok_or_else(|| {
                Error::new(
                    term.span,
                    "bracket requires positive-degree cohomology classes",
                )
            })?;
            degree = degree
                .checked_add(reduced)
                .ok_or_else(|| Error::new(term.span, "bracket result degree overflow"))?;
            base = Some(space);
            terms.push(term);
        }
        let domain = self.model.cohomology_type(base.unwrap(), degree)?;
        let range = self.model.named("S").unwrap();
        let ty = self.model.arrow(domain, range)?;
        Ok((ty, TypedForm::Bracket(terms)))
    }

    fn defining_system(
        &mut self,
        inputs: Vec<ParsedExpr>,
        witnesses: Option<Vec<Witness>>,
        span: Span,
    ) -> Result<(TypeId, TypedForm)> {
        let (arity, layout) = if witnesses.is_some() {
            if inputs.len() < 2 {
                return Err(Error::new(
                    span,
                    "massey requires at least two adjacent cochains",
                ));
            }
            (inputs.len(), Layout::Explicit)
        } else {
            if inputs.is_empty() {
                return Err(Error::new(
                    span,
                    "powerdef requires at least one bar cochain",
                ));
            }
            (inputs.len() + 1, Layout::Repeated)
        };
        let mut input_iter = inputs.into_iter();
        let first = self.expr(input_iter.next().unwrap())?;
        let op = if witnesses.is_some() {
            "massey"
        } else {
            "powerdef"
        };
        let Kind::Bar(modulus) = *self.model.kind(first.ty) else {
            return Err(self.invalid_operand(&first, op, "a bar cochain"));
        };
        let expected = first.ty;
        let mut cochains = vec![first];
        for argument in input_iter {
            let value = self.expr(argument)?;
            self.require_type(expected, &value, op)?;
            cochains.push(value);
        }

        if let Some(witnesses) = witnesses {
            let mut entries = HashMap::new();
            for Witness {
                from,
                to,
                span,
                value,
            } in witnesses
            {
                if from >= to || to > arity {
                    return Err(Error::new(
                        span,
                        format!(
                            "invalid interval [{from},{to}]: expected 0 <= start < end <= {arity}",
                        ),
                    ));
                }
                let gap = to - from;
                if gap == 1 {
                    return Err(Error::new(
                        span,
                        format!("adjacent interval [{from},{to}] already has an input"),
                    ));
                }
                if gap == arity {
                    return Err(Error::new(
                        span,
                        format!("top interval [{from},{to}] cannot have a witness"),
                    ));
                }
                let std::collections::hash_map::Entry::Vacant(entry) = entries.entry((gap, from))
                else {
                    return Err(Error::new(
                        span,
                        format!("duplicate witness for interval [{from},{to}]"),
                    ));
                };
                let value = self.expr(value)?;
                self.require_type(
                    expected,
                    &value,
                    &format!("witness for interval [{from},{to}]"),
                )?;
                entry.insert(value);
            }
            for gap in 2..arity {
                for from in 0..=arity - gap {
                    let to = from + gap;
                    let value = entries.remove(&(gap, from)).ok_or_else(|| {
                        Error::new(span, format!("missing witness for interval [{from},{to}]"))
                    })?;
                    cochains.push(value);
                }
            }
        }
        let ty = self.model.cyclic_type(modulus, true)?;
        Ok((
            ty,
            TypedForm::DefiningSystem {
                arity,
                layout,
                cochains,
            },
        ))
    }

    fn operation_degree(&self, value: &TypedExpr, op: &str) -> Result<usize> {
        match *self.model.kind(value.ty) {
            Kind::Operations(degree) => Ok(degree),
            _ => Err(self.invalid_operand(value, op, "a Steenrod operation")),
        }
    }

    fn cobar_bidegree(&self, value: &TypedExpr, op: &str) -> Result<(usize, usize)> {
        match *self.model.kind(value.ty) {
            Kind::Cobar(s, t) => Ok((s, t)),
            _ => Err(self.invalid_operand(value, op, "a cobar cochain")),
        }
    }

    fn ext_bidegree(&self, value: &TypedExpr, op: &str) -> Result<(usize, usize)> {
        match *self.model.kind(value.ty) {
            Kind::Ext(s, t) => Ok((s, t)),
            _ => Err(self.invalid_operand(value, op, "an Ext class")),
        }
    }

    fn bp_vpage_grade(
        &mut self,
        (p, n, page, a, h, t, filtration): crate::syntax::BpVPageSpec,
    ) -> Result<bp_vpage::Grade> {
        let shape = bp::Shape::new(self.model, p, n, a)?.quotient(h)?;
        bp_vpage::Grade::new(shape, page, t, filtration)
    }

    fn bp_page_grade(
        &mut self,
        mode: bp_page::Mode,
        (p, n, page, a, t, filtration): crate::syntax::BpPageSpec,
    ) -> Result<bp_page::Grade> {
        let shape = bp::Shape::new(self.model, p, n, a)?;
        Ok(bp_page::Grade {
            mode,
            shape,
            page,
            t,
            filtration,
        })
    }

    fn bp_height_grade(
        &mut self,
        BpHeightSpec(spec, r, h, t): BpHeightSpec,
    ) -> Result<bp_height::ClassGrade> {
        let top = self.bp_ext_grade(spec, t)?.shape;
        Ok(bp_height::ClassGrade { top, r, h, t })
    }

    fn bp_relative_shape(&mut self, spec: crate::syntax::BpRelSpec) -> Result<bp_rel::Shape> {
        bp_rel::Shape::new(self.model, spec.0, spec.1, spec.2, spec.3)
    }

    fn bp_mod_ext_grade(&mut self, spec: crate::syntax::BpModSpec) -> Result<bp_mod_ext::Grade> {
        let crate::syntax::BpSpec(p, n, a) = spec.shape;
        let shape = bp::Shape::new(self.model, p, n, a)?;
        let shape = if let Some(j) = spec.exponent {
            if !(2..=shape.height + 1).contains(&spec.h) {
                return Err(Error::plain("BPThickExt requires 2 <= h <= n+1"));
            }
            shape.quotient(spec.h - 1)?.thick(j)?
        } else {
            shape.quotient(spec.h)?
        };
        Ok(bp_mod_ext::Grade { shape, t: spec.t })
    }

    fn bp_ext_grade(&mut self, spec: crate::syntax::BpSpec, t: usize) -> Result<bp_classes::Grade> {
        let shape = bp::Shape::new(self.model, spec.0, spec.1, spec.2)?;
        Ok(bp_classes::Grade { shape, t })
    }

    fn bp_ext_operand(&self, value: &TypedExpr, name: &str) -> Result<bp_classes::Grade> {
        let Kind::BPExt(grade) = *self.model.kind(value.ty) else {
            return Err(self.invalid_operand(value, name, "a BPExt class"));
        };
        Ok(grade)
    }

    fn bp_triple_grade(&self, value: &TypedExpr, name: &str) -> Result<bp_classes::Grade> {
        match *self.model.kind(value.ty) {
            Kind::BPTriple(grade) => Ok(grade),
            _ => Err(self.invalid_operand(value, name, "a BPTriple coset")),
        }
    }

    fn bp_triple_inputs(
        &mut self,
        [a, b, c]: [ParsedExpr; 3],
        name: &str,
    ) -> Result<(bp_triple::Layout, Box<[TypedExpr; 3]>)> {
        let inputs = Box::new([self.expr(a)?, self.expr(b)?, self.expr(c)?]);
        let [a, b, c] = &*inputs;
        let layout = bp_triple::Layout::new([
            self.bp_ext_operand(a, name)?,
            self.bp_ext_operand(b, name)?,
            self.bp_ext_operand(c, name)?,
        ])?;
        Ok((layout, inputs))
    }

    fn bp_shape(&self, value: &TypedExpr, name: &str) -> Result<bp::Shape> {
        match self.model.kind(value.ty) {
            Kind::BP(shape) => Ok(*shape),
            _ => Err(self.invalid_operand(value, name, "a BP polynomial")),
        }
    }

    fn rv_bar_page_grade(&self, spec: RvBarPageSpec) -> Result<rv_bar_page::Grade> {
        let chain = self.rv_bar_grade(spec.chain)?;
        Ok(rv_bar_page::Grade {
            chain,
            r: spec.r,
            length: spec.length,
        })
    }

    fn cyclic_grade(
        &self,
        value: &TypedExpr,
        name: &str,
    ) -> Result<(rv_cyclic_homology::Grade, bool)> {
        match self.model.kind(value.ty) {
            Kind::RvCyclicChain(grade) => Ok((grade.clone(), false)),
            Kind::RvCyclicHomology(grade) => Ok((grade.clone(), true)),
            _ => Err(self.invalid_operand(value, name, "a cyclic chain or class")),
        }
    }

    fn rv_negative_grade(
        &mut self,
        (ideal, length, b, weight): RvNCSpec,
    ) -> Result<rv_negative::Grade> {
        rv_negative::Grade::new(self.rv_hh_grade((ideal, b, weight))?, length)
    }

    fn check_mixed_cone_map(
        &mut self,
        name: Label,
        operation: mixed_cone::Operation,
        input: ParsedExpr,
    ) -> Result<(TypeId, TypedForm)> {
        use mixed_cone::Operation;
        let input = self.expr(input)?;
        let id = self.model.mixed_id(&name)?;
        let cone = self.model.mixed_cone(id)?;
        let ty = match *self.model.kind(input.ty) {
            Kind::MixedHomology(grade) if !matches!(operation, Operation::Lift) => self
                .model
                .mixed_homology_type(cone.map_grade(id, grade, operation)?)?,
            Kind::MixedChain(grade, total) if !matches!(operation, Operation::Connect) => self
                .model
                .mixed_chain_type(cone.map_grade(id, grade, operation)?, total)?,
            _ => {
                let (name, expected) = match operation {
                    Operation::Connect => ("mixConeConnect", "a MixedH class"),
                    Operation::Lift => ("mixConeLift", "a MixedChain or MixedU chain"),
                    _ => (
                        "mixed mapping cone operation",
                        "a MixedChain, MixedU, or MixedH value",
                    ),
                };
                return Err(self.invalid_operand(&input, name, expected));
            }
        };
        Ok((
            ty,
            TypedForm::MixedConeMap(cone, operation, Box::new(input)),
        ))
    }

    fn check_mixed_h_generator(
        &mut self,
        name: Label,
        b: i64,
        index: usize,
        at: Span,
    ) -> Result<(TypeId, TypedForm)> {
        let complex = self.model.mixed_id(&name)?;
        let grade = mixed::Grade { complex, b };
        let ty = self.model.mixed_homology_type(grade)?;
        let rank = self.model.mixed_homology(grade)?.rank();
        let value = self.prime_generator(ty, rank, index, at)?;
        Ok((ty, TypedForm::Constant(value)))
    }

    fn check_mixed_h_unary(
        &mut self,
        op: &CobarOp,
        argument: &TypedExpr,
        name: &str,
    ) -> Result<TypeId> {
        match op {
            CobarOp::MixHClass | CobarOp::MixHPrimitive => {
                let Kind::MixedChain(grade, true) = *self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "a MixedU chain"));
                };
                if matches!(op, CobarOp::MixHClass) {
                    self.model.mixed_homology_type(grade)
                } else {
                    self.model.mixed_chain_type(grade.shift(1)?, true)
                }
            }
            CobarOp::MixHRep => {
                let Kind::MixedHomology(grade) = *self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "a MixedH class"));
                };
                self.model.mixed_chain_type(grade, true)
            }
            _ => unreachable!(),
        }
    }

    fn check_mixed_induce(&mut self, left: &TypedExpr, right: &TypedExpr) -> Result<TypeId> {
        let Kind::MixedChain(map, true) = *self.model.kind(left.ty) else {
            return Err(self.invalid_operand(left, "mixInduce", "a total Hom map"));
        };
        let Kind::MixedHomology(source) = *self.model.kind(right.ty) else {
            return Err(self.invalid_operand(right, "mixInduce", "a MixedH class"));
        };
        let hom = self.model.mixed_hom(map.complex)?;
        if source.complex != hom.source {
            return Err(Error::plain(
                "mixInduce class belongs to the wrong source complex",
            ));
        }
        let grade = mixed::Grade {
            complex: hom.target,
            b: source.b,
        }
        .shift(map.b.into())?;
        self.model.mixed_homology_type(grade)
    }

    fn mixed_chain(&self, value: &TypedExpr, name: &str) -> Result<(mixed::Grade, bool)> {
        match *self.model.kind(value.ty) {
            Kind::MixedChain(grade, total) => Ok((grade, total)),
            _ => Err(self.invalid_operand(value, name, "a MixedChain or MixedU chain")),
        }
    }

    fn mixed_grade(&self, (name, r, b, f): MixedSpec) -> Result<mixed_page::Grade> {
        let complex = self.model.mixed_id(&name)?;
        Ok(mixed_page::Grade {
            chain: mixed::Grade { complex, b },
            r,
            f,
        })
    }

    fn rv_hh_grade(&mut self, (ideal, b, weight): RvHHSpec) -> Result<rv_hh::Grade> {
        let basis = self.model.rv_ideal(&ideal)?;
        let koszul = rv_koszul::Grade::local(self.model, &basis, ideal.text.into(), 0, weight)?;
        Ok(rv_hh::Grade::from_tor(koszul, b))
    }

    fn rv_bar_grade(&self, spec: RvBarSpec) -> Result<rv_bar::Grade> {
        let path = self.model.rv_bar_path_id(&spec.path)?;
        self.model
            .rv_bar_path(path)
            .check_interval(spec.start, spec.end)?;
        Ok(rv_bar::Grade {
            path,
            start: spec.start,
            end: spec.end,
            b: spec.b,
        })
    }

    fn rv_koszul_grade(&mut self, spec: RvGradeSpec) -> Result<rv_koszul::Grade> {
        match spec.degree {
            RvDegree::Total(t) => self.model.rv_grade(spec.ideal, spec.s, t),
            RvDegree::Multi(multidegree) => {
                let ideal = self.model.rv_ideal(&spec.ideal)?;
                rv_koszul::Grade::local(
                    self.model,
                    &ideal,
                    spec.ideal.text.into(),
                    spec.s,
                    multidegree,
                )
            }
        }
    }

    fn rv_complex_family(&self, value: &TypedExpr, name: &str) -> Result<(rv_koszul::Grade, bool)> {
        match self.model.kind(value.ty) {
            Kind::RVKoszul(grade) => Ok((grade.clone(), false)),
            Kind::RVTor(grade) => Ok((grade.clone(), true)),
            _ => Err(self.invalid_operand(value, name, "a rank Koszul chain or Tor class")),
        }
    }

    fn rv_complex_grade(
        &self,
        value: &TypedExpr,
        homology: bool,
        name: &str,
    ) -> Result<rv_koszul::Grade> {
        match (self.model.kind(value.ty), homology) {
            (Kind::RVKoszul(grade), false) | (Kind::RVTor(grade), true) => Ok(grade.clone()),
            _ => {
                let expected = if homology {
                    "an RVTor class"
                } else {
                    "an RVKoszul chain"
                };
                Err(self.invalid_operand(value, name, expected))
            }
        }
    }

    fn rv_poly_grade(&self, value: &TypedExpr, name: &str) -> Result<(q_koszul::Indices, usize)> {
        match self.model.kind(value.ty) {
            Kind::RVPoly(selected, degree) => Ok((selected.clone(), *degree)),
            _ => Err(self.invalid_operand(value, name, "an RVPoly polynomial")),
        }
    }

    fn rv_grade(&mut self, spec: crate::syntax::RankVarietySpec) -> Result<rv_hom::Grade> {
        let source = self.q_endpoint(*spec.source)?;
        let target = self.q_endpoint(*spec.target)?;
        let selected = q_koszul::indices(self.model, spec.selected)?;
        rv_hom::Grade::new(
            self.model,
            source,
            target,
            selected,
            spec.modulus,
            spec.point.unwrap(),
        )
    }

    fn jet_page_grade(&mut self, spec: JetPageSpec) -> Result<jet_page::Grade> {
        let chain = self.jet_grade(spec.chain)?;
        Ok(jet_page::Grade {
            chain,
            r: spec.r,
            f: spec.f,
        })
    }

    fn jet_operation_type(&mut self, op: JetOp, value: &TypedExpr) -> Result<TypeId> {
        let kind = self.model.kind(value.ty);
        match (op, kind) {
            (JetOp::Class(r, f), &Kind::JetHom(chain)) => {
                self.model.jet_page_type(jet_page::Grade { chain, r, f })
            }
            (JetOp::Differential, &Kind::JetHom(grade)) => {
                self.model.jet_hom_type(grade.differential())
            }
            (JetOp::Differential, &Kind::JetPage(grade)) => {
                self.model.jet_page_type(grade.differential()?)
            }
            (JetOp::Representative, &Kind::JetPage(grade)) => self.model.jet_hom_type(grade.chain),
            (JetOp::Next, &Kind::JetPage(grade)) => self.model.jet_page_type(grade.next()?),
            (JetOp::Prev, &Kind::JetPage(grade)) => self.model.jet_page_type(grade.prev()?),
            (JetOp::Class(..), _) => {
                Err(self.invalid_operand(value, "jetClass", "a JetHom cochain"))
            }
            (JetOp::Differential, _) => {
                Err(self.invalid_operand(value, "jetD", "a JetHom cochain or JetPage class"))
            }
            (JetOp::Representative, _) => {
                Err(self.invalid_operand(value, "jetRep", "a JetPage class"))
            }
            (JetOp::Next, _) => Err(self.invalid_operand(value, "jetNext", "a JetPage class")),
            (JetOp::Prev, _) => Err(self.invalid_operand(value, "jetPrev", "a JetPage class")),
        }
    }

    fn jet_grade(&mut self, spec: JetHomSpec) -> Result<jet_hom::Grade> {
        let source = self.q_endpoint(*spec.source)?;
        let target = self.q_endpoint(*spec.target)?;
        let path = self.model.jet_path_id(&spec.path)?;
        Ok(jet_hom::Grade {
            diagram: jet_hom::Diagram {
                source,
                target,
                path,
            },
            parity: spec.parity,
        })
    }

    fn q_endpoint(&mut self, ty: ParsedType) -> Result<TypeId> {
        if let TypeForm::Name(name) = &ty.kind {
            if let Some(base) = self.model.q_endpoint_named(&name.text) {
                return Ok(base);
            }
        }
        self.ty(ty)
    }

    fn q_grade<S>(&mut self, spec: QSpec<S>) -> Result<Grade<S>> {
        let source = self.q_endpoint(*spec.source)?;
        let target = self.q_endpoint(*spec.target)?;
        let selected = q_koszul::indices(self.model, spec.selected)?;
        Ok(Grade {
            source,
            target,
            selected,
            s: spec.s,
            t: spec.t,
        })
    }

    fn vbock_type(&mut self, spec: VBockSpec) -> Result<TypeId> {
        let localized = self.q_grade(spec.localized)?;
        self.model.vbock_type(vbock::Grade {
            localized,
            r: spec.r,
            f: spec.f,
        })
    }

    fn map_may_type(&mut self, spec: MapMaySpec) -> Result<TypeId> {
        let source = self.ty(*spec.source)?;
        let target = self.ty(*spec.target)?;
        self.model
            .map_may_type(source, target, spec.r, spec.s, spec.t, spec.f)
    }

    fn may_degree(
        &self,
        value: &TypedExpr,
        op: &str,
    ) -> Result<(Profile, usize, usize, usize, usize)> {
        self.model
            .kind(value.ty)
            .may_grade()
            .ok_or_else(|| self.invalid_operand(value, op, "a May page class"))
    }

    fn ce_degree(
        &self,
        value: &TypedExpr,
        op: &str,
    ) -> Result<(usize, usize, usize, usize, usize)> {
        match *self.model.kind(value.ty) {
            Kind::CeAn(n, r, s, t, f) => Ok((n, r, s, t, f)),
            _ => Err(self.invalid_operand(value, op, "a CEAn class")),
        }
    }

    fn rmext_grade(&self, input: &TypedExpr, operation: &str) -> Result<motivic::Grade> {
        if let Kind::RMExt(grade) = *self.model.kind(input.ty) {
            Ok(grade)
        } else {
            Err(self.invalid_operand(input, operation, "an RMExt class"))
        }
    }

    fn real_grade(&self, input: &TypedExpr, operation: &str) -> Result<motivic::Grade> {
        self.model.kind(input.ty).real_grade().ok_or_else(|| {
            let expected = if operation == "rhoShift" {
                "an RMCobar cochain, an RMExt class, or a RhoLoc cochain or RhoMin polynomial"
            } else {
                "an RMCobar cochain or RMExt class"
            };
            self.invalid_operand(input, operation, expected)
        })
    }

    fn real_rho_result(
        &mut self,
        op: RealRhoOp,
        input: &TypedExpr,
    ) -> Result<(TypeId, motivic::Grade)> {
        let kind = self.model.kind(input.ty);
        if let RealRhoOp::Reduce(exponent) = op {
            let (chain,is_class) = match *kind {
                Kind::RMCobar(grade) => (grade,false),
                Kind::RMExt(grade) => (grade,true),
                Kind::RMRhoCobar(grade) | Kind::RMRhoExt(grade) => {
                    if exponent > grade.exponent {
                        return Err(Error::plain("rhoReduce cannot increase the truncation exponent"));
                    }
                    (grade.chain,matches!(kind,Kind::RMRhoExt(..)))
                }
                Kind::MCobar(grade) if exponent == 1 => (grade,false),
                Kind::MExt(grade) if exponent == 1 => (grade,true),
                _ => return Err(self.invalid_operand(input,op.name(),"a real-motivic cochain or class, a rho truncation, or a complex-motivic value when reducing to rho^1")),
            };
            let grade = real_rho::Grade::new(exponent, chain)?;
            let ty = if is_class {
                self.model.real_rho_ext_type(grade)?
            } else {
                self.model.real_rho_cobar_type(grade)?
            };
            return Ok((ty, chain));
        }
        let grade = match *kind {
            Kind::RMRhoExt(grade) => grade,
            Kind::RMRhoCobar(grade) if matches!(op, RealRhoOp::Lift) => grade,
            _ => {
                return Err(self.invalid_operand(
                    input,
                    op.name(),
                    if matches!(op, RealRhoOp::Lift) {
                        "an RMRhoCobar cochain or RMRhoExt class"
                    } else {
                        "an RMRhoExt class"
                    },
                ))
            }
        };
        let ty = match op {
            RealRhoOp::Lift => self.model.real_cobar_type(grade.chain)?,
            RealRhoOp::Connect => self.model.rmext_type(grade.connecting()?)?,
            RealRhoOp::ConnectRep => self.model.real_cobar_type(grade.connecting()?)?,
            RealRhoOp::Beta => self.model.real_rho_ext_type(real_rho::Grade {
                chain: grade.connecting()?,
                ..grade
            })?,
            RealRhoOp::BetaRep | RealRhoOp::BetaNull => {
                let mut chain = grade.connecting()?;
                if matches!(op, RealRhoOp::BetaNull) {
                    chain = chain
                        .previous()
                        .expect("connecting differential increases length");
                }
                self.model
                    .real_rho_cobar_type(real_rho::Grade { chain, ..grade })?
            }
            RealRhoOp::Reduce(..) => unreachable!(),
        };
        Ok((ty, grade.chain))
    }

    fn rho_min_result(
        &mut self,
        op: RhoMinOp,
        input: &TypedExpr,
    ) -> Result<(motivic::Grade, TypeId)> {
        let kind = self.model.kind(input.ty);
        let grade = match (kind, op) {
            (Kind::MExt(grade), RhoMinOp::Inject(_))
            | (Kind::RMCobar(grade), RhoMinOp::Project | RhoMinOp::Homotopy)
            | (
                Kind::RhoMin(grade),
                RhoMinOp::Coefficient(_) | RhoMinOp::Differential | RhoMinOp::Lift,
            ) => *grade,
            _ => {
                let expected = match op {
                    RhoMinOp::Inject(_) => "an MExt class",
                    RhoMinOp::Project | RhoMinOp::Homotopy => "an RMCobar cochain",
                    _ => "a RhoMin polynomial",
                };
                return Err(self.invalid_operand(input, op.name(), expected));
            }
        };
        let ty = match op {
            RhoMinOp::Inject(f) => self.model.rho_min_type(real_motivic::shifted(grade, f)?)?,
            RhoMinOp::Coefficient(f) => self
                .model
                .mext_type(rho_minimal::coefficient_grade(grade, f)?)?,
            RhoMinOp::Differential => self.model.rho_min_type(grade.differential()?)?,
            RhoMinOp::Lift => self.model.real_cobar_type(grade)?,
            RhoMinOp::Project => self.model.rho_min_type(grade)?,
            RhoMinOp::Homotopy => {
                let previous = grade
                    .previous()
                    .ok_or_else(|| Error::plain("rhoMinH requires positive cobar length"))?;
                self.model.real_cobar_type(previous)?
            }
        };
        Ok((grade, ty))
    }

    fn rm_triple_grade(&self, input: &TypedExpr, operation: &str) -> Result<motivic::Grade> {
        match *self.model.kind(input.ty) {
            Kind::RMTriple(grade) => Ok(grade),
            _ => Err(self.invalid_operand(input, operation, "an RMTriple coset")),
        }
    }

    fn cobar_result(
        &mut self,
        op: &CobarOp,
        argument: &TypedExpr,
        name: &str,
        span: Span,
    ) -> Result<TypeId> {
        if let Kind::RhoLoc(grade) = *self.model.kind(argument.ty) {
            if matches!(op, CobarOp::Differential) {
                return self.model.rho_loc_type(grade.differential()?);
            }
        }
        if let Kind::RMRhoCobar(grade) = *self.model.kind(argument.ty) {
            match op {
                CobarOp::Differential => {
                    return self.model.real_rho_cobar_type(grade.differential()?)
                }
                CobarOp::Class => return self.model.real_rho_ext_type(grade),
                CobarOp::Primitive => {
                    let chain = grade.chain.previous().ok_or_else(|| {
                        Error::new(span, "extPrimitive requires positive cobar length")
                    })?;
                    return self
                        .model
                        .real_rho_cobar_type(real_rho::Grade { chain, ..grade });
                }
                _ => {}
            }
        }
        if let Kind::RMRhoExt(grade) = *self.model.kind(argument.ty) {
            if matches!(op, CobarOp::Representative) {
                return self.model.real_rho_cobar_type(grade);
            }
        }
        if let Kind::RMCobar(grade) = *self.model.kind(argument.ty) {
            match op {
                CobarOp::Differential => return self.model.real_cobar_type(grade.differential()?),
                CobarOp::Class => return self.model.rmext_type(grade),
                CobarOp::Primitive => {
                    let previous = grade.previous().ok_or_else(|| {
                        Error::new(span, "extPrimitive requires positive cobar length")
                    })?;
                    return self.model.real_cobar_type(previous);
                }
                _ => {}
            }
        }
        if let Kind::RMExt(grade) = *self.model.kind(argument.ty) {
            if matches!(op, CobarOp::Representative) {
                return self.model.real_cobar_type(grade);
            }
        }
        if let Kind::MCobar(grade) = *self.model.kind(argument.ty) {
            match op {
                CobarOp::Differential => return self.model.mcobar_type(grade.differential()?),
                CobarOp::Class => return self.model.mext_type(grade),
                CobarOp::Primitive => {
                    let previous = grade.previous().ok_or_else(|| {
                        Error::new(span, "extPrimitive requires positive cobar length")
                    })?;
                    return self.model.mcobar_type(previous);
                }
                _ => {}
            }
        }
        if let Kind::MExt(grade) = *self.model.kind(argument.ty) {
            if matches!(op, CobarOp::Representative) {
                return self.model.mcobar_type(grade);
            }
        }
        if let Some((p, profile, s, t)) = self.model.kind(argument.ty).pcobar_grade() {
            match op {
                CobarOp::Differential => {
                    let s = checked_degree(s, 1, name, span)?;
                    return self.model.profile_pcobar_type(p, profile, s, t);
                }
                CobarOp::Class => return self.model.profile_odd_ext_type(p, profile, s, t),
                CobarOp::Primitive => {
                    let s = s.checked_sub(1).ok_or_else(|| {
                        Error::new(span, "extPrimitive requires positive cobar length")
                    })?;
                    return self.model.profile_pcobar_type(p, profile, s, t);
                }
                _ => {}
            }
        }
        match op {
            CobarOp::VBockDifferential | CobarOp::VBockRepresentative => {
                let Kind::VBock(grade) = self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "a VBock class"));
                };
                if matches!(op, CobarOp::VBockRepresentative) {
                    return self.model.v_koszul_type(grade.localized.clone());
                }
                return vbock::differential_type(self.model, argument.ty);
            }
            CobarOp::PageTransition(family, direction) => {
                return direction
                    .result_type(self.model, argument.ty, *family)
                    .map_err(|error| error.at(argument.span));
            }
            CobarOp::CeRepresentative => {
                if let Some((p, n, _, s, t, _)) = self.model.kind(argument.ty).pce_grade() {
                    return self.model.profile_pcobar_type(p, Profile::An(n), s, t);
                }
                let (n, _, s, t, _) = self.ce_degree(argument, name)?;
                return self.model.profile_cobar_type(Profile::An(n), s, t);
            }
            CobarOp::CeDifferential => {
                if let Some((p, n, r, s, t, f)) = self.model.kind(argument.ty).pce_grade() {
                    let s = checked_degree(s, 1, name, span)?;
                    let f = checked_degree(f, r, name, span)?;
                    return self.model.pce_type(p, n, r, s, t, f);
                }
                let (n, r, s, t, f) = self.ce_degree(argument, name)?;
                let s = checked_degree(s, 1, name, span)?;
                let f = checked_degree(f, r, name, span)?;
                return self.model.ce_type(n, r, s, t, f);
            }
            CobarOp::RvBarNext | CobarOp::RvBarPrev | CobarOp::RvBarRep => {
                let Kind::RVBarE(grade) = *self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "an RVBarE class"));
                };
                return match op {
                    CobarOp::RvBarNext => self.model.rv_bar_page_type(grade.next()?),
                    CobarOp::RvBarPrev => self.model.rv_bar_page_type(grade.previous()?),
                    _ => self.model.rv_bar_type(grade.chain),
                };
            }
            CobarOp::RvCyclicD | CobarOp::RvCyclicClass => {
                let Kind::RvCyclicChain(grade) = self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "a cyclic chain"));
                };
                if matches!(op, CobarOp::RvCyclicD) {
                    return self.model.rv_cyclic_type(grade.differential()?, false);
                }
                return self.model.rv_cyclic_type(grade.clone(), true);
            }
            CobarOp::RvCyclicRep => {
                let Kind::RvCyclicHomology(grade) = self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "a cyclic homology class"));
                };
                return self.model.rv_cyclic_type(grade.clone(), false);
            }
            CobarOp::RvCyclicI => {
                let (ordinary, homology) = match self.model.kind(argument.ty) {
                    Kind::RVHHChain(grade) => (grade.clone(), false),
                    Kind::RvHochschild(grade) => (grade.clone(), true),
                    _ => {
                        return Err(self.invalid_operand(
                            argument,
                            name,
                            "an RVHHChain or RVHH class",
                        ))
                    }
                };
                let mode = rv_cyclic_homology::Mode::Ordinary;
                return self
                    .model
                    .rv_cyclic_type(rv_cyclic_homology::Grade { ordinary, mode }, homology);
            }
            CobarOp::RvCyclicB => {
                let (grade, homology) = self.cyclic_grade(argument, name)?;
                let result = grade.connes_output()?;
                return if homology {
                    self.model.rv_hh_class_type(result)
                } else {
                    self.model.rv_hh_type(result)
                };
            }
            CobarOp::RvNCd | CobarOp::RvNCClass => {
                let Kind::RVNCChain(grade) = self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "an RVNCChain"));
                };
                if matches!(op, CobarOp::RvNCd) {
                    return self.model.rv_negative_type(grade.differential()?, false);
                }
                return self.model.rv_negative_type(grade.clone(), true);
            }
            CobarOp::RvNCRep => {
                let Kind::RvNegativeCyclic(grade) = self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "an RVNC class"));
                };
                return self.model.rv_negative_type(grade.clone(), false);
            }
            CobarOp::RvMixedToHH => {
                return match *self.model.kind(argument.ty) {
                    Kind::MixedChain(grade, false) => {
                        let origin = self.model.mixed_hh_origin(grade.complex)?;
                        self.model
                            .rv_hh_type(origin.original_grade(grade.b.into())?)
                    }
                    Kind::MixedPage(grade) if grade.r == 1 => {
                        let origin = self.model.mixed_hh_origin(grade.chain.complex)?;
                        let b = i128::from(grade.chain.b) + 2 * grade.f as i128;
                        self.model.rv_hh_class_type(origin.original_grade(b)?)
                    }
                    _ => Err(self.invalid_operand(
                        argument,
                        name,
                        "an imported MixedChain or page-1 MixedE class",
                    )),
                };
            }
            CobarOp::MixHClass | CobarOp::MixHPrimitive | CobarOp::MixHRep => {
                return self.check_mixed_h_unary(op, argument, name);
            }
            CobarOp::MixMinI => {
                let Kind::MixedMinimum(grade, total) = *self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(
                        argument,
                        name,
                        "a MixedMin or MixedMinU value",
                    ));
                };
                let data = self.model.mixed_minimal(grade.minimal);
                return self.model.mixed_chain_type(
                    mixed::Grade {
                        complex: data.source,
                        b: grade.b,
                    },
                    total,
                );
            }
            CobarOp::MixDown | CobarOp::MixUp => {
                let Kind::MixedChain(grade, false) = *self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "a MixedChain"));
                };
                if !self.model.mixed_complex(grade.complex).is_strict() {
                    return Err(Error::plain("mixb and mixB require a strict mixed complex. Use mixDelta on a higher complex."));
                }
                let offset = if matches!(op, CobarOp::MixUp) { 1 } else { -1 };
                return self.model.mixed_chain_type(grade.shift(offset)?, false);
            }
            CobarOp::MixD => {
                return match *self.model.kind(argument.ty) {
                    Kind::MixedChain(grade, true) => {
                        self.model.mixed_chain_type(grade.shift(-1)?, true)
                    }
                    Kind::MixedMinimum(grade, true) => {
                        self.model.mixed_minimal_type(grade.shift(-1)?, true)
                    }
                    Kind::MixedPage(grade) => self.model.mixed_page_type(grade.differential()?),
                    _ => Err(self.invalid_operand(
                        argument,
                        name,
                        "a MixedU, MixedE, or MixedMinU value",
                    )),
                };
            }
            CobarOp::MixRep | CobarOp::MixNext | CobarOp::MixPrev => {
                let Kind::MixedPage(grade) = *self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "a MixedE class"));
                };
                return match op {
                    CobarOp::MixRep => self.model.mixed_chain_type(grade.chain, true),
                    CobarOp::MixNext => self.model.mixed_page_type(grade.next()?),
                    _ => self.model.mixed_page_type(grade.previous()?),
                };
            }
            CobarOp::RvHHd | CobarOp::RvConnesB | CobarOp::RvHHClass | CobarOp::RvHHRep => {
                return match (op, self.model.kind(argument.ty)) {
                    (CobarOp::RvHHd, Kind::RVHHChain(grade)) => {
                        self.model.rv_hh_type(grade.differential()?)
                    }
                    (CobarOp::RvConnesB, Kind::RVHHChain(grade)) => {
                        self.model.rv_hh_type(grade.connes()?)
                    }
                    (CobarOp::RvConnesB, Kind::RvHochschild(grade)) => {
                        self.model.rv_hh_class_type(grade.connes()?)
                    }
                    (CobarOp::RvHHClass, Kind::RVHHChain(grade)) => {
                        self.model.rv_hh_class_type(grade.clone())
                    }
                    (CobarOp::RvHHRep, Kind::RvHochschild(grade)) => {
                        self.model.rv_hh_type(grade.clone())
                    }
                    (CobarOp::RvConnesB, _) => {
                        Err(self.invalid_operand(argument, name, "an RVHHChain or RVHH class"))
                    }
                    (CobarOp::RvHHRep, _) => {
                        Err(self.invalid_operand(argument, name, "an RVHH class"))
                    }
                    _ => Err(self.invalid_operand(argument, name, "an RVHHChain")),
                };
            }
            CobarOp::RvBCDesuspend => {
                let Kind::RVBar(grade) = *self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "an RVBar chain"));
                };
                return self
                    .model
                    .rv_barcobar_type(rv_barcobar::suspend_down(grade, name)?);
            }
            CobarOp::RvBCD | CobarOp::RvBCProject | CobarOp::RvBCHomotopy => {
                let Kind::RVBarCobar(grade) = *self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "an RVBarCobar chain"));
                };
                return match op {
                    CobarOp::RvBCD => self
                        .model
                        .rv_barcobar_type(rv_barcobar::suspend_down(grade, name)?),
                    CobarOp::RvBCHomotopy => self
                        .model
                        .rv_barcobar_type(rv_barcobar::homotopy_grade(grade)?),
                    _ => {
                        let tor = self.model.rv_bar_path(grade.path).tor_grade(
                            grade.start,
                            grade.end,
                            grade.b,
                        )?;
                        self.model.rv_complex_type(tor, true)
                    }
                };
            }
            CobarOp::RvBarDualize => {
                let Kind::RVBar(grade) = *self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "an RVBar chain"));
                };
                return self.model.rv_path_type(grade, false);
            }
            CobarOp::RvBarDualD | CobarOp::RvPathClass => {
                let Kind::RVBarDual(grade) = *self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "an RVBarDual cochain"));
                };
                if matches!(op, CobarOp::RvBarDualD) {
                    return self
                        .model
                        .rv_path_type(rv_bar_dual::differential_grade(grade)?, false);
                }
                return self.model.rv_path_type(grade, true);
            }
            CobarOp::RvPathRep => {
                let Kind::RVPathExt(grade) = *self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "an RVPathExt class"));
                };
                return self.model.rv_path_type(grade, false);
            }
            CobarOp::RvBarD => {
                return match *self.model.kind(argument.ty) {
                    Kind::RVBar(grade) => self.model.rv_bar_type(grade.differential()?),
                    Kind::RVBarE(grade) => self.model.rv_bar_page_type(grade.differential()?),
                    _ => {
                        Err(self.invalid_operand(argument, name, "an RVBar chain or RVBarE class"))
                    }
                };
            }

            CobarOp::RvMForget => {
                let (mut grade, homology) = self.rv_complex_family(argument, name)?;
                if grade.multidegree.take().is_none() {
                    return Err(Error::new(
                        argument.span,
                        "rvMForget requires a multigraded Koszul chain or Tor class",
                    ));
                }
                return self.model.rv_complex_type(grade, homology);
            }
            CobarOp::RvKoszulDifferential
            | CobarOp::RvKoszulPrimitive
            | CobarOp::RvTorClass
            | CobarOp::RvTorRepresentative => {
                let underlying_tor = matches!(op, CobarOp::RvTorRepresentative);
                let grade = self.rv_complex_grade(argument, underlying_tor, name)?;
                let grade = match op {
                    CobarOp::RvKoszulDifferential => grade.differential()?,
                    CobarOp::RvKoszulPrimitive => grade.primitive()?,
                    _ => grade,
                };
                return self
                    .model
                    .rv_complex_type(grade, matches!(op, CobarOp::RvTorClass));
            }
            CobarOp::QDifferential | CobarOp::QClass => {
                if let Kind::VKoszul(grade) = self.model.kind(argument.ty) {
                    let mut grade = grade.clone();
                    if matches!(op, CobarOp::QClass) {
                        return self.model.v_ext_type(grade);
                    }
                    grade.s = grade
                        .s
                        .checked_add(1)
                        .ok_or_else(|| Error::new(span, "qD signed filtration is out of range"))?;
                    return self.model.v_koszul_type(grade);
                }
                let Kind::QKoszul(grade) = self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(
                        argument,
                        name,
                        "a QKoszul or VKoszul cochain",
                    ));
                };
                let mut grade = grade.clone();
                if matches!(op, CobarOp::QClass) {
                    return self.model.q_ext_type(grade);
                }
                grade.s = checked_degree(grade.s, 1, name, span)?;
                return self.model.q_koszul_type(grade);
            }
            CobarOp::QRepresentative => {
                if let Kind::VExt(grade) = self.model.kind(argument.ty) {
                    return self.model.v_koszul_type(grade.clone());
                }
                let Kind::QExt(grade) = self.model.kind(argument.ty) else {
                    return Err(self.invalid_operand(argument, name, "a QExt or VExt class"));
                };
                return self.model.q_koszul_type(grade.clone());
            }
            CobarOp::VLocal => return v_koszul::localization_type(self.model, argument.ty),
            CobarOp::MapDifferential => {
                let (source, target, s, t) = self.map_cobar_bidegree(argument, name)?;
                let s = checked_degree(s, 1, name, span)?;
                return self.model.map_cobar_type(source, target, s, t);
            }
            CobarOp::MapClass => {
                let (source, target, s, t) = self.map_cobar_bidegree(argument, name)?;
                return self.model.map_ext_type(source, target, s, t);
            }
            CobarOp::MapRepresentative => {
                let (source, target, s, t) = self.map_ext_bidegree(argument, name)?;
                return self.model.map_cobar_type(source, target, s, t);
            }
            CobarOp::Representative => {
                if let Some((p, profile, s, t)) = self.model.kind(argument.ty).pext_grade() {
                    return self.model.profile_pcobar_type(p, profile, s, t);
                }
                let (profile, s, t) = self
                    .model
                    .kind(argument.ty)
                    .ext_grade()
                    .ok_or_else(|| self.invalid_operand(argument, name, "an Ext class"))?;
                return self.model.profile_cobar_type(profile, s, t);
            }
            CobarOp::MayRepresentative => {
                if let Some((p, profile, _, s, t, _)) = self.model.kind(argument.ty).pmay_grade() {
                    return self.model.profile_pcobar_type(p, profile, s, t);
                }
                if let Kind::MapMay(source, target, _, s, t, _) = *self.model.kind(argument.ty) {
                    return self.model.map_cobar_type(source, target, s, t);
                }
                let (profile, _, s, t, _) = self.may_degree(argument, name)?;
                return self.model.profile_cobar_type(profile, s, t);
            }
            CobarOp::MayDifferential => {
                if let Some((p, profile, r, s, t, f)) = self.model.kind(argument.ty).pmay_grade() {
                    let s = checked_degree(s, 1, name, span)?;
                    let f = checked_degree(f, r, name, span)?;
                    return self.model.profile_pmay_type(p, profile, r, s, t, f);
                }
                if let Kind::MapMay(source, target, r, s, t, f) = *self.model.kind(argument.ty) {
                    let s = checked_degree(s, 1, name, span)?;
                    let f = checked_degree(f, r, name, span)?;
                    return self.model.map_may_type(source, target, r, s, t, f);
                }
                let (profile, r, s, t, f) = self.may_degree(argument, name)?;
                let s = checked_degree(s, 1, name, span)?;
                let f = checked_degree(f, r, name, span)?;
                return self.model.profile_may_type(profile, r, s, t, f);
            }
            _ => {}
        }
        let (profile, s, t) = self.profile_cobar_bidegree(argument, name)?;
        match op {
            CobarOp::Differential => {
                self.model
                    .profile_cobar_type(profile, checked_degree(s, 1, name, span)?, t)
            }
            CobarOp::Class => self.model.profile_ext_type(profile, s, t),
            CobarOp::Primitive => {
                let previous = s.checked_sub(1).ok_or_else(|| {
                    Error::new(span, "extPrimitive requires positive cobar length")
                })?;
                self.model.profile_cobar_type(profile, previous, t)
            }
            _ => unreachable!(),
        }
    }

    fn profile_cobar_bidegree(
        &self,
        value: &TypedExpr,
        op: &str,
    ) -> Result<(Profile, usize, usize)> {
        self.model
            .kind(value.ty)
            .cobar_grade()
            .ok_or_else(|| self.invalid_operand(value, op, "a cobar cochain"))
    }

    fn prime_generator(
        &mut self,
        ty: TypeId,
        rank: usize,
        index: usize,
        at: Span,
    ) -> Result<Value> {
        self.check_generator(ty, rank, index, at)?;
        odd_ext::persist(self.model, vec![(index, 1)])
    }

    fn generator(&mut self, ty: TypeId, rank: usize, index: usize, at: Span) -> Result<Value> {
        self.check_generator(ty, rank, index, at)?;
        self.model.allocate(1)?;
        Ok(Value::Class([index].into()))
    }

    fn check_generator(&self, ty: TypeId, rank: usize, index: usize, at: Span) -> Result<()> {
        if index >= rank {
            return Err(Error::new(
                at,
                format!(
                    "generator index {index} is out of range for {} of rank {rank}",
                    self.model.type_name(ty),
                ),
            ));
        }
        Ok(())
    }

    fn map_ext_bidegree(
        &self,
        value: &TypedExpr,
        op: &str,
    ) -> Result<(TypeId, TypeId, usize, i64)> {
        match *self.model.kind(value.ty) {
            Kind::MapExt(source, target, s, t) => Ok((source, target, s, t)),
            _ => Err(self.invalid_operand(value, op, "a MapExt class")),
        }
    }

    fn map_cobar_bidegree(
        &self,
        value: &TypedExpr,
        op: &str,
    ) -> Result<(TypeId, TypeId, usize, i64)> {
        match *self.model.kind(value.ty) {
            Kind::MapCobar(source, target, s, t) => Ok((source, target, s, t)),
            _ => Err(self.invalid_operand(value, op, "a MapCobar cochain")),
        }
    }

    fn wzw_grade(&self, input: &TypedExpr, operation: &str) -> Result<rv_wzw::Grade> {
        match *self.model.kind(input.ty) {
            Kind::RVWzw(grade) => Ok(grade),
            _ => Err(self.invalid_operand(input, operation, "a WZW torus state")),
        }
    }

    fn tl_grade(&self, input: &TypedExpr, operation: &str) -> Result<rv_tl::Grade> {
        match *self.model.kind(input.ty) {
            Kind::RvTl(grade) => Ok(grade),
            _ => Err(self.invalid_operand(input, operation, "a Temperley–Lieb path vector")),
        }
    }

    fn cyclotomic_ratio_order(&self, input: &TypedExpr, operation: &str) -> Result<usize> {
        match *self.model.kind(input.ty) {
            Kind::RVCycloQ(order) => Ok(order),
            _ => Err(self.invalid_operand(input, operation, "a cyclotomic rational")),
        }
    }

    fn cyclotomic_order(&self, input: &TypedExpr, operation: &str) -> Result<usize> {
        match *self.model.kind(input.ty) {
            Kind::RVCyclo(order) => Ok(order),
            _ => Err(self.invalid_operand(input, operation, "a cyclotomic integer")),
        }
    }

    fn sheaf_kuranishi_expression(
        &mut self,
        query: sheaf_kuranishi::Query,
    ) -> Result<(TypeId, TypedForm)> {
        use sheaf_kuranishi::{Coefficient, Query};
        let (ty, value) = match query {
            Query::Coefficient(name, powers, coefficient) => {
                let chart = self.model.sheaf_kuranishi(&name)?;
                let (form, degree) = match coefficient {
                    Coefficient::Lift => (SheafForm::Cochain, 1),
                    Coefficient::Residual => (SheafForm::Cochain, 2),
                    Coefficient::Obstruction => (SheafForm::Class, 2),
                };
                let ends = SheafEnds::endomorphisms(chart.sheaf.clone());
                let ty = self.model.known_sheaf_hom_type(form, ends, degree)?;
                let value = chart.coefficient(self.model, &powers, coefficient)?;
                (ty, value)
            }
            Query::Arc(name, amplitudes, weights, order, class) => {
                let chart = self.model.sheaf_kuranishi(&name)?;
                let ends = SheafEnds::endomorphisms(chart.sheaf.clone());
                let hom = self.model.sheaf_hom(&ends);
                let ty = if class {
                    self.model.known_sheaf_hom_type(SheafForm::Class, ends, 2)?
                } else {
                    self.model
                        .known_sheaf_deformation_type(chart.sheaf.clone(), order)?
                };
                let value =
                    chart.specialized(self.model, &hom, &amplitudes, &weights, order, class)?;
                (ty, value)
            }
        };
        Ok((ty, TypedForm::Constant(value)))
    }

    fn check_sheaf_deformation(&mut self, form: ExprForm) -> Result<(TypeId, TypedForm)> {
        match form {
            ExprForm::SheafDefStart(name) => {
                self.model.sheaf(&name)?;
                let ty = self
                    .model
                    .known_sheaf_deformation_type(name.text.into(), 0)?;
                let jet = sheaf_deform::Jet::new(self.model, Vec::new())?;
                Ok((ty, TypedForm::Constant(Value::SheafDeformation(jet))))
            }
            ExprForm::SheafDefOperation(operation, input) => {
                use sheaf_deform::Operation;
                let input = self.expr(*input)?;
                let (name, order) = self.sheaf_deformation(&input, operation.name())?;
                let ends = SheafEnds::endomorphisms(name.clone());
                let hom = self.model.sheaf_hom(&ends);
                let ty = match operation {
                    Operation::Term(k) => {
                        if k > order {
                            return Err(Error::plain(
                                "shDefTerm coefficient is beyond the known truncation",
                            ));
                        }
                        self.model
                            .known_sheaf_hom_type(SheafForm::Cochain, ends, 1)?
                    }
                    Operation::Curvature => {
                        self.model
                            .known_sheaf_hom_type(SheafForm::Cochain, ends, 2)?
                    }
                    Operation::Obstruction => {
                        self.model.known_sheaf_hom_type(SheafForm::Class, ends, 2)?
                    }
                    Operation::Lift => self
                        .model
                        .known_sheaf_deformation_type(name, next_def_order(order)?)?,
                };
                Ok((
                    ty,
                    TypedForm::SheafDefOperation(hom, operation, Box::new(input)),
                ))
            }
            ExprForm::SheafDefAppend(args) => {
                let [jet, term] = *args;
                let jet = self.expr(jet)?;
                let term = self.expr(term)?;
                let (name, order) = self.sheaf_deformation(&jet, "shDefAppend")?;
                let ends = SheafEnds::endomorphisms(name.clone());
                self.require_sheaf_def_cochain(&term, &ends, 1, "shDefAppend")?;
                let hom = self.model.sheaf_hom(&ends);
                let ty = self
                    .model
                    .known_sheaf_deformation_type(name, next_def_order(order)?)?;
                Ok((ty, TypedForm::SheafDefAppend(hom, Box::new([jet, term]))))
            }
            ExprForm::SheafDefGauge(exponent, witness, args) => {
                let [jet, term] = *args;
                let jet = self.expr(jet)?;
                let term = self.expr(term)?;
                let name = if witness {
                    "shDefGaugeWitness"
                } else {
                    "shDefGauge"
                };
                let (sheaf, order) = self.sheaf_deformation(&jet, name)?;
                if exponent == 0 || exponent > order {
                    return Err(Error::plain("sheaf deformation gauge exponent must lie between 1 and the truncation order"));
                }
                let ends = SheafEnds::endomorphisms(sheaf);
                self.require_sheaf_def_cochain(&term, &ends, 0, name)?;
                let hom = self.model.sheaf_hom(&ends);
                let ty = if witness {
                    self.model
                        .known_sheaf_hom_type(SheafForm::Cochain, ends, 1)?
                } else {
                    jet.ty
                };
                Ok((
                    ty,
                    TypedForm::SheafDefGauge(hom, exponent, witness, Box::new([jet, term])),
                ))
            }
            _ => unreachable!(),
        }
    }

    fn sheaf_deformation(&self, jet: &TypedExpr, op: &str) -> Result<(Rc<str>, usize)> {
        match self.model.kind(jet.ty) {
            Kind::SheafDeformation(name, order) => Ok((name.clone(), *order)),
            _ => Err(self.invalid_operand(jet, op, "a ShDef jet")),
        }
    }

    fn require_sheaf_def_cochain(
        &self,
        term: &TypedExpr,
        ends: &SheafEnds,
        degree: usize,
        op: &str,
    ) -> Result<()> {
        if matches!(self.model.kind(term.ty), Kind::SheafHom(SheafForm::Cochain, found, d) if found == ends && *d == degree)
        {
            Ok(())
        } else {
            let expected = format!(
                "a ShHomC({}, {}, {degree}) cochain",
                ends.source, ends.target
            );
            Err(self.invalid_operand(term, op, &expected))
        }
    }

    fn sheaf_hom_grade(&self, value: &TypedExpr) -> Result<(SheafForm, &SheafEnds, usize)> {
        match self.model.kind(value.ty) {
            Kind::SheafHom(form, ends, degree) => Ok((*form, ends, *degree)),
            _ => {
                Err(self.invalid_operand(value, "shExtCompose", "a ShHomC cochain or ShExt class"))
            }
        }
    }

    fn invalid_operand(&self, value: &TypedExpr, op: &str, expected: &str) -> Error {
        Error::new(
            value.span,
            format!(
                "{op} requires {expected}, found {}",
                self.model.type_name(value.ty),
            ),
        )
    }

    fn binary_result(
        &mut self,
        op: BinaryOp,
        left: &TypedExpr,
        right: &TypedExpr,
        span: Span,
    ) -> Result<TypeId> {
        match op {
            BinaryOp::Add => {
                if !matches!(
                    self.model.kind(left.ty),
                    Kind::Cohomology(..)
                        | Kind::Stable(..)
                        | Kind::Operations(..)
                        | Kind::Cobar(..)
                        | Kind::MCobar(..)
                        | Kind::RMCobar(..)
                        | Kind::RMRhoCobar(..)
                        | Kind::RhoBock(..)
                        | Kind::RhoLoc(..)
                        | Kind::RhoMin(..)
                        | Kind::MExt(..)
                        | Kind::MixedChain(..)
                        | Kind::MixedHomology(..)
                        | Kind::MixedMinimum(..)
                        | Kind::MixedPage(..)
                        | Kind::Sheaf(..)
                        | Kind::SheafComplex(..)
                        | Kind::SheafHom(..)
                        | Kind::SheafPage(..)
                        | Kind::RMExt(..)
                        | Kind::RMRhoExt(..)
                        | Kind::RMBar(..)
                        | Kind::SheafLocal(..)
                        | Kind::PCobar(..)
                        | Kind::PCobarAn(..)
                        | Kind::PExt(..)
                        | Kind::PExtAn(..)
                        | Kind::PMay(..)
                        | Kind::PMayAn(..)
                        | Kind::PCeAn(..)
                        | Kind::BPPage(..)
                        | Kind::BPVPage(..)
                        | Kind::BPModExt(..)
                        | Kind::CobarAn(..)
                        | Kind::MapCobar(..)
                        | Kind::QKoszul(..)
                        | Kind::JetHom(..)
                        | Kind::RVBarE(..)
                        | Kind::RvCyclicChain(..)
                        | Kind::RvCyclicHomology(..)
                        | Kind::RVNCChain(..)
                        | Kind::RvNegativeCyclic(..)
                        | Kind::RvHochschild(..)
                        | Kind::RVHHChain(..)
                        | Kind::RVBarCobar(..)
                        | Kind::RVBarDual(..)
                        | Kind::RVPathExt(..)
                        | Kind::RVBar(..)
                        | Kind::RVKoszul(..)
                        | Kind::RVTor(..)
                        | Kind::JetPage(..)
                        | Kind::RVPoly(..)
                        | Kind::CurveDiff(..)
                        | Kind::ASOddDiff(..)
                        | Kind::ASTwoDiff(..)
                        | Kind::ASOddDeRham(..)
                        | Kind::ASOddTate(..)
                        | Kind::PrimeScalar(..)
                        | Kind::ASDeRham(..)
                        | Kind::BinaryScalar(..)
                        | Kind::BP(..)
                        | Kind::BPRel(..)
                        | Kind::BPExt(..)
                        | Kind::BPRelExt(..)
                        | Kind::BPHeight(..)
                        | Kind::QExt(..)
                        | Kind::RVHom(..)
                        | Kind::RVLaurent(..)
                        | Kind::RVIwHecke(..)
                        | Kind::RVSphHecke(..)
                        | Kind::RVSatake(..)
                        | Kind::RVQGr(..)
                        | Kind::RVQGrP(..)
                        | Kind::RVOre
                        | Kind::RVCTState(..)
                        | Kind::RVCTResidue(..)
                        | Kind::RVFusion(..)
                        | Kind::RVCyclo(..)
                        | Kind::RVCycloQ(..)
                        | Kind::RVQpQuad(..)
                        | Kind::RvTl(..)
                        | Kind::RVWzw(..)
                        | Kind::VKoszul(..)
                        | Kind::VExt(..)
                        | Kind::VBock(..)
                        | Kind::Ext(..)
                        | Kind::ExtAn(..)
                        | Kind::XExt(..)
                        | Kind::MapExt(..)
                        | Kind::MapMay(..)
                        | Kind::May(..)
                        | Kind::MayAn(..)
                        | Kind::CeAn(..)
                        | Kind::Bar(..)
                        | Kind::CyclicH2(..)
                ) {
                    return Err(self.invalid_operand(
                        left,
                        "add",
                        "a cohomology class, Steenrod operation, cobar cochain, MapCobar cochain, Ext class, XExt class, MapExt class, May class, MapMay class, or bar cochain",
                    ));
                }
                self.require_type(left.ty, right, "add")?;
                Ok(left.ty)
            }
            BinaryOp::Mul => {
                let left = self.operation_degree(left, "mul")?;
                let right = self.operation_degree(right, "mul")?;
                let degree = checked_degree(left, right, "mul", span)?;
                self.model.operation_type(degree)
            }
            BinaryOp::MilnorPair => {
                let degree = self.operation_degree(left, "milnorPair")?;
                match *self.model.kind(right.ty) {
                    Kind::Cobar(1, t) if t == degree => {}
                    Kind::Cobar(0, 0) if degree == 0 => {}
                    _ => {
                        let expected = if degree == 0 {
                            "a cochain in Cobar(1, 0) or Cobar(0, 0)".into()
                        } else {
                            format!("a cochain in Cobar(1, {degree})")
                        };
                        return Err(self.invalid_operand(right, "milnorPair", &expected));
                    }
                }
                Ok(self.model.named("S").unwrap())
            }
            BinaryOp::Concat => {
                if let Kind::RMRhoCobar(a) = *self.model.kind(left.ty) {
                    let Kind::RMRhoCobar(b) = *self.model.kind(right.ty) else {
                        return Err(self.invalid_operand(right, "concat", "an RMRhoCobar cochain"));
                    };
                    return self.model.real_rho_cobar_type(a.product(b)?);
                }
                if matches!(self.model.kind(right.ty), Kind::RMRhoCobar(..)) {
                    return Err(self.invalid_operand(left, "concat", "an RMRhoCobar cochain"));
                }
                if let Kind::RhoLoc(a) = *self.model.kind(left.ty) {
                    let Kind::RhoLoc(b) = *self.model.kind(right.ty) else {
                        return Err(self.invalid_operand(right, "concat", "a RhoLoc cochain"));
                    };
                    return self.model.rho_loc_type(a.product(b)?);
                }
                if matches!(self.model.kind(right.ty), Kind::RhoLoc(..)) {
                    return Err(self.invalid_operand(left, "concat", "a RhoLoc cochain"));
                }
                if let Kind::RMCobar(a) = *self.model.kind(left.ty) {
                    let Kind::RMCobar(b) = *self.model.kind(right.ty) else {
                        return Err(self.invalid_operand(right, "concat", "an RMCobar cochain"));
                    };
                    return self.model.real_cobar_type(a.product(b)?);
                }
                if matches!(self.model.kind(right.ty), Kind::RMCobar(..)) {
                    return Err(self.invalid_operand(left, "concat", "an RMCobar cochain"));
                }
                if let Kind::MCobar(a) = *self.model.kind(left.ty) {
                    let Kind::MCobar(b) = *self.model.kind(right.ty) else {
                        return Err(self.invalid_operand(right, "concat", "an MCobar cochain"));
                    };
                    return self.model.mcobar_type(a.product(b)?);
                }
                if matches!(self.model.kind(right.ty), Kind::MCobar(..)) {
                    return Err(self.invalid_operand(left, "concat", "an MCobar cochain"));
                }
                if let Some((p, profile, s, t)) = self.model.kind(left.ty).pcobar_grade() {
                    let Some((q, other, u, v)) = self.model.kind(right.ty).pcobar_grade() else {
                        return Err(self.invalid_operand(right, "concat", "a PCobar cochain"));
                    };
                    if p != q {
                        return Err(Error::new(right.span, "concat requires the same odd prime"));
                    }
                    require_profile(profile, other, "concat", right.span)?;
                    let s = checked_degree(s, u, "concat", span)?;
                    let t = checked_degree(t, v, "concat", span)?;
                    return self.model.profile_pcobar_type(p, profile, s, t);
                }
                if self.model.kind(right.ty).pcobar_grade().is_some() {
                    return Err(self.invalid_operand(left, "concat", "a PCobar cochain"));
                }
                let (profile, s, t) = self.profile_cobar_bidegree(left, "concat")?;
                let (other, u, v) = self.profile_cobar_bidegree(right, "concat")?;
                require_profile(profile, other, "concat", right.span)?;
                let length = checked_degree(s, u, "concat", span)?;
                let degree = checked_degree(t, v, "concat", span)?;
                self.model.profile_cobar_type(profile, length, degree)
            }
            BinaryOp::Yoneda => {
                if let Kind::RMRhoExt(a) = *self.model.kind(left.ty) {
                    let Kind::RMRhoExt(b) = *self.model.kind(right.ty) else {
                        return Err(self.invalid_operand(right, "yoneda", "an RMRhoExt class"));
                    };
                    return self.model.real_rho_ext_type(a.product(b)?);
                }
                if matches!(self.model.kind(right.ty), Kind::RMRhoExt(..)) {
                    return Err(self.invalid_operand(left, "yoneda", "an RMRhoExt class"));
                }
                if let Kind::RMExt(a) = *self.model.kind(left.ty) {
                    let Kind::RMExt(b) = *self.model.kind(right.ty) else {
                        return Err(self.invalid_operand(right, "yoneda", "an RMExt class"));
                    };
                    return self.model.rmext_type(a.product(b)?);
                }
                if matches!(self.model.kind(right.ty), Kind::RMExt(..)) {
                    return Err(self.invalid_operand(left, "yoneda", "an RMExt class"));
                }
                if let Kind::MExt(a) = *self.model.kind(left.ty) {
                    let Kind::MExt(b) = *self.model.kind(right.ty) else {
                        return Err(self.invalid_operand(right, "yoneda", "an MExt class"));
                    };
                    return self.model.mext_type(a.product(b)?);
                }
                if matches!(self.model.kind(right.ty), Kind::MExt(..)) {
                    return Err(self.invalid_operand(left, "yoneda", "an MExt class"));
                }
                if let Some((p, profile, s, t)) = self.model.kind(left.ty).pext_grade() {
                    let Some((q, other, u, v)) = self.model.kind(right.ty).pext_grade() else {
                        return Err(self.invalid_operand(right, "yoneda", "a PExt class"));
                    };
                    if p != q {
                        return Err(Error::new(right.span, "yoneda requires the same odd prime"));
                    }
                    require_profile(profile, other, "yoneda", right.span)?;
                    let length = checked_degree(s, u, "yoneda", span)?;
                    let degree = checked_degree(t, v, "yoneda", span)?;
                    return self.model.profile_odd_ext_type(p, profile, length, degree);
                }
                let (profile, s, t) = self
                    .model
                    .kind(left.ty)
                    .ext_grade()
                    .ok_or_else(|| self.invalid_operand(left, "yoneda", "an Ext class"))?;
                let (other, u, v) = self
                    .model
                    .kind(right.ty)
                    .ext_grade()
                    .ok_or_else(|| self.invalid_operand(right, "yoneda", "an Ext class"))?;
                require_profile(profile, other, "yoneda", right.span)?;
                let length = checked_degree(s, u, "yoneda", span)?;
                let degree = checked_degree(t, v, "yoneda", span)?;
                self.model.profile_ext_type(profile, length, degree)
            }
            BinaryOp::XAction => {
                let (s, t) = self.ext_bidegree(left, "xaction")?;
                let Kind::XExt(base, u, v) = *self.model.kind(right.ty) else {
                    return Err(self.invalid_operand(right, "xaction", "an XExt class"));
                };
                let length = checked_degree(s, u, "xaction", span)?;
                let degree = checked_degree(t, v, "xaction", span)?;
                self.model.xext_type(base, length, degree)
            }
            BinaryOp::EOBruhat | BinaryOp::EOClosure => {
                let operation = if matches!(op, BinaryOp::EOBruhat) {
                    "eoBruhat"
                } else {
                    "eoClosure"
                };
                self.require_eo(left, operation)?;
                self.require_type(left.ty, right, operation)?;
                Ok(self.model.named("S").unwrap())
            }
            BinaryOp::AsDeRhamPair => {
                let op = "rvASDeRhamPair";
                let Kind::ASDeRham(_, field) = *self.model.kind(left.ty) else {
                    return Err(self.invalid_operand(left, op, "an ASDeRham class"));
                };
                self.require_type(left.ty, right, op)?;
                self.model.binary_scalar_type(field)
            }
            BinaryOp::AsOddTatePair => {
                let op = "rvASOddTatePair";
                let Kind::ASOddTate(first, prime, a) = self.model.kind(left.ty) else {
                    return Err(self.invalid_operand(left, op, "an ASOddTate class"));
                };
                let Kind::ASOddTate(second, _, b) = self.model.kind(right.ty) else {
                    return Err(self.invalid_operand(right, op, "an ASOddTate class"));
                };
                if first != second {
                    return Err(Error::new(
                        span,
                        "rvASOddTatePair requires classes of the same cover",
                    ));
                }
                if a.rem_euclid(2) == b.rem_euclid(2) {
                    return Err(Error::new(
                        span,
                        "rvASOddTatePair requires opposite Tate parities",
                    ));
                }
                self.model.prime_scalar_type(*prime)
            }
            BinaryOp::AsOddDeRhamPair => {
                let op = "rvASOddDeRhamPair";
                let Kind::ASOddDeRham(_, prime) = *self.model.kind(left.ty) else {
                    return Err(self.invalid_operand(left, op, "an ASOddDeRham class"));
                };
                self.require_type(left.ty, right, op)?;
                self.model.prime_scalar_type(prime)
            }
            operation @ (BinaryOp::ResSum | BinaryOp::ResIntersection) => {
                let op = match operation {
                    BinaryOp::ResSum => "rvResSum",
                    _ => "rvResIntersection",
                };
                self.require_res_space(left, op)?;
                self.require_type(left.ty, right, op)?;
                Ok(left.ty)
            }
            BinaryOp::QpQuadMultiply => {
                self.qp_quad_field(left, "rvQpQuadMul")?;
                self.require_type(left.ty, right, "rvQpQuadMul")?;
                Ok(left.ty)
            }
            BinaryOp::BpCup => {
                if matches!(self.model.kind(left.ty), Kind::BPExt(..))
                    || matches!(self.model.kind(right.ty), Kind::BPExt(..))
                {
                    let left = self.bp_ext_operand(left, "bpCup")?;
                    let right = self.bp_ext_operand(right, "bpCup")?;
                    return self.model.bp_ext_type(left.cup(right)?);
                }
                if let (&Kind::BPModExt(left), &Kind::BPModExt(right)) =
                    (self.model.kind(left.ty), self.model.kind(right.ty))
                {
                    return self.model.bp_mod_ext_type(left.cup(right)?);
                }
                if let (&Kind::BPVPage(left), &Kind::BPVPage(right)) =
                    (self.model.kind(left.ty), self.model.kind(right.ty))
                {
                    return self.model.bp_vpage_type(left.cup(right)?);
                }
                if let (&Kind::BPPage(left), &Kind::BPPage(right)) =
                    (self.model.kind(left.ty), self.model.kind(right.ty))
                {
                    return self.model.bp_page_type(left.cup(right)?);
                }
                let left = self.bp_shape(left, "bpCup")?;
                let right = self.bp_shape(right, "bpCup")?;
                self.model.bp_type(left.cup(right)?)
            }
            BinaryOp::BpMultiply => {
                self.bp_shape(left, "bpMul")?;
                self.bp_shape(right, "bpMul")?;
                self.require_type(left.ty, right, "bpMul")?;
                Ok(left.ty)
            }
            BinaryOp::RvBCMul => {
                let Kind::RVBarCobar(a) = *self.model.kind(left.ty) else {
                    return Err(self.invalid_operand(left, "rvBCMul", "an RVBarCobar chain"));
                };
                let Kind::RVBarCobar(b) = *self.model.kind(right.ty) else {
                    return Err(self.invalid_operand(right, "rvBCMul", "an RVBarCobar chain"));
                };
                self.model
                    .rv_barcobar_type(rv_barcobar::product_grade(a, b)?)
            }
            BinaryOp::MixInduce => self.check_mixed_induce(left, right),
            BinaryOp::MixEval => {
                let (map, total) = self.mixed_chain(left, "mixEval")?;
                let (argument, argument_total) = self.mixed_chain(right, "mixEval")?;
                let hom = self.model.mixed_hom(map.complex)?;
                if argument.complex != hom.source {
                    return Err(Error::plain(
                        "mixEval argument belongs to the wrong source complex",
                    ));
                }
                if total != argument_total {
                    return Err(Error::plain("mixEval requires chains of the same form"));
                }
                let grade = mixed::Grade {
                    complex: hom.target,
                    b: map.b,
                }
                .shift(argument.b.into())?;
                self.model.mixed_chain_type(grade, total)
            }
            BinaryOp::RvBarPair => {
                let Kind::RVBarDual(grade) = *self.model.kind(left.ty) else {
                    return Err(self.invalid_operand(left, "rvBarPair", "an RVBarDual cochain"));
                };
                let raw = self.model.rv_bar_type(grade)?;
                self.require_type(raw, right, "rvBarPair")?;
                Ok(self.model.named("S").unwrap())
            }
            BinaryOp::RvPathCup => {
                let (a, b, homology) =
                    match (self.model.kind(left.ty), self.model.kind(right.ty)) {
                        (Kind::RVBarDual(a), Kind::RVBarDual(b)) => (*a, *b, false),
                        (Kind::RVPathExt(a), Kind::RVPathExt(b)) => (*a, *b, true),
                        _ => return Err(Error::new(
                            span,
                            "rvPathCup requires two RVBarDual cochains or two RVPathExt classes",
                        )),
                    };
                self.model
                    .rv_path_type(rv_bar_dual::cup_grade(a, b)?, homology)
            }
            BinaryOp::RvKoszulProduct | BinaryOp::RvTorProduct => {
                let homology = matches!(op, BinaryOp::RvTorProduct);
                let name = if homology { "rvTorProd" } else { "rvkProd" };
                let left = self.rv_complex_grade(left, homology, name)?;
                let right = self.rv_complex_grade(right, homology, name)?;
                self.model.rv_complex_type(left.product(&right)?, homology)
            }
            BinaryOp::RvMultiply => {
                let (selected, s) = self.rv_poly_grade(left, "rvMul")?;
                let (other, t) = self.rv_poly_grade(right, "rvMul")?;
                if selected != other {
                    return Err(Error::new(
                        span,
                        "rvMul requires the same selected primitives",
                    ));
                }
                let degree = checked_degree(s, t, "rvMul", span)?;
                self.model.rv_poly_type(selected, degree)
            }
            BinaryOp::XMap => {
                let Kind::Arrow(domain, range) = *self.model.kind(left.ty) else {
                    return Err(self.invalid_operand(left, "xmap", "a function"));
                };
                let Kind::XExt(base, s, t) = *self.model.kind(right.ty) else {
                    return Err(self.invalid_operand(right, "xmap", "an XExt class"));
                };
                if base != domain {
                    let expected = self.model.xext_type(domain, s, t)?;
                    self.require_type(expected, right, "xmap")?;
                }
                self.model.xext_type(range, s, t)
            }
            BinaryOp::MapExtCompose => {
                if matches!(self.model.kind(left.ty), Kind::RVLaurent(..))
                    || matches!(self.model.kind(right.ty), Kind::RVLaurent(..))
                {
                    return rv_laurent::composition_type(self.model, left.ty, right.ty)
                        .map_err(|error| error.at(right.span));
                }
                if matches!(self.model.kind(left.ty), Kind::RVHom(..))
                    || matches!(self.model.kind(right.ty), Kind::RVHom(..))
                {
                    return rv_hom::composition_type(self.model, left.ty, right.ty)
                        .map_err(|error| error.at(right.span));
                }
                let left_kind = self.model.kind(left.ty);
                let right_kind = self.model.kind(right.ty);
                if matches!(left_kind, Kind::JetPage(..)) || matches!(right_kind, Kind::JetPage(..))
                {
                    let Kind::JetPage(outer) = *left_kind else {
                        return Err(self.invalid_operand(left, "ecomp", "a JetPage class"));
                    };
                    let Kind::JetPage(inner) = *right_kind else {
                        return Err(self.invalid_operand(right, "ecomp", "a JetPage class"));
                    };
                    return self
                        .model
                        .jet_page_type(jet_page::Grade::compose(outer, inner)?);
                }
                if matches!(left_kind, Kind::JetHom(..)) || matches!(right_kind, Kind::JetHom(..)) {
                    let Kind::JetHom(outer) = *left_kind else {
                        return Err(self.invalid_operand(left, "ecomp", "a JetHom cochain"));
                    };
                    let Kind::JetHom(inner) = *right_kind else {
                        return Err(self.invalid_operand(right, "ecomp", "a JetHom cochain"));
                    };
                    return self
                        .model
                        .jet_hom_type(jet_hom::Grade::compose(outer, inner)?);
                }

                if matches!(self.model.kind(left.ty), Kind::VBock(..))
                    || matches!(self.model.kind(right.ty), Kind::VBock(..))
                {
                    return vbock::composition_type(self.model, left.ty, right.ty);
                }
                if matches!(
                    self.model.kind(left.ty),
                    Kind::QKoszul(..) | Kind::QExt(..) | Kind::VKoszul(..) | Kind::VExt(..)
                ) || matches!(
                    self.model.kind(right.ty),
                    Kind::QKoszul(..) | Kind::QExt(..) | Kind::VKoszul(..) | Kind::VExt(..)
                ) {
                    return q_koszul::composition_type(self.model, left.ty, right.ty)
                        .map_err(|error| error.at(right.span));
                }
                if let Kind::MapMay(middle, target, r, gs, gt, gf) = *self.model.kind(left.ty) {
                    let Kind::MapMay(source, inner_target, page, fs, ft, ff) =
                        *self.model.kind(right.ty)
                    else {
                        return Err(self.invalid_operand(right, "ecomp", "a MapMay class"));
                    };
                    if r != page {
                        return Err(Error::new(
                            right.span,
                            "ecomp requires MapMay classes on the same page",
                        ));
                    }
                    if inner_target != middle {
                        return Err(Error::new(
                            right.span,
                            "ecomp requires composable MapMay endpoints",
                        ));
                    }
                    let s = checked_degree(fs, gs, "ecomp", span)?;
                    let t = ft.checked_add(gt).ok_or_else(|| {
                        Error::new(span, "ecomp signed internal degree is out of range")
                    })?;
                    let f = checked_degree(ff, gf, "ecomp", span)?;
                    return self.model.map_may_type(source, target, r, s, t, f);
                }
                let raw = matches!(self.model.kind(left.ty), Kind::MapCobar(..));
                let bidegree = if raw {
                    Self::map_cobar_bidegree
                } else {
                    Self::map_ext_bidegree
                };
                let (middle, target, gs, gt) = bidegree(self, left, "ecomp")?;
                let (source, inner_target, fs, ft) = bidegree(self, right, "ecomp")?;
                if inner_target != middle {
                    return Err(Error::new(
                        right.span,
                        format!(
                            "ecomp requires the inner target {} to match the outer source {}",
                            self.model.type_name(inner_target),
                            self.model.type_name(middle),
                        ),
                    ));
                }
                let s = checked_degree(fs, gs, "ecomp", span)?;
                let t = ft.checked_add(gt).ok_or_else(|| {
                    Error::new(span, "ecomp signed internal degree is out of range")
                })?;
                if raw {
                    self.model.map_cobar_type(source, target, s, t)
                } else {
                    self.model.map_ext_type(source, target, s, t)
                }
            }
            BinaryOp::CeMul => {
                if let Some((p, n, r, s, t, f)) = self.model.kind(left.ty).pce_grade() {
                    let Some((q, m, page, u, v, g)) = self.model.kind(right.ty).pce_grade() else {
                        return Err(self.invalid_operand(right, "ceMul", "a PCEAn class"));
                    };
                    if p != q {
                        return Err(Error::new(right.span, "ceMul requires the same prime"));
                    }
                    require_profile(Profile::An(n), Profile::An(m), "ceMul", right.span)?;
                    if r != page {
                        return Err(Error::new(
                            right.span,
                            "ceMul requires classes on the same page",
                        ));
                    }
                    let s = checked_degree(s, u, "ceMul", span)?;
                    let t = checked_degree(t, v, "ceMul", span)?;
                    let f = checked_degree(f, g, "ceMul", span)?;
                    return self.model.pce_type(p, n, r, s, t, f);
                }
                let (n, r, s, t, f) = self.ce_degree(left, "ceMul")?;
                let (m, page, u, v, g) = self.ce_degree(right, "ceMul")?;
                require_profile(Profile::An(n), Profile::An(m), "ceMul", right.span)?;
                if r != page {
                    return Err(Error::new(
                        right.span,
                        "ceMul requires classes on the same page",
                    ));
                }
                let s = checked_degree(s, u, "ceMul", span)?;
                let t = checked_degree(t, v, "ceMul", span)?;
                let f = checked_degree(f, g, "ceMul", span)?;
                self.model.ce_type(n, r, s, t, f)
            }
            BinaryOp::RhoMul => {
                let Kind::RhoBock(a) = *self.model.kind(left.ty) else {
                    return Err(self.invalid_operand(left, "rhoMul", "a RhoBock class"));
                };
                let Kind::RhoBock(b) = *self.model.kind(right.ty) else {
                    return Err(self.invalid_operand(right, "rhoMul", "a RhoBock class"));
                };
                self.model.rho_page_type(a.product(b)?)
            }
            BinaryOp::MayMul => {
                if let Some((p, profile, r, s, t, f)) = self.model.kind(left.ty).pmay_grade() {
                    let Some((q, other, page, u, v, g)) = self.model.kind(right.ty).pmay_grade()
                    else {
                        return Err(self.invalid_operand(right, "mayMul", "a PMay class"));
                    };
                    if p != q {
                        return Err(Error::new(right.span, "mayMul requires the same prime"));
                    }
                    require_profile(profile, other, "mayMul", right.span)?;
                    if r != page {
                        return Err(Error::new(
                            right.span,
                            "mayMul requires classes on the same page",
                        ));
                    }
                    let s = checked_degree(s, u, "mayMul", span)?;
                    let t = checked_degree(t, v, "mayMul", span)?;
                    let f = checked_degree(f, g, "mayMul", span)?;
                    return self.model.profile_pmay_type(p, profile, r, s, t, f);
                }
                let (profile, r, s, t, f) = self.may_degree(left, "mayMul")?;
                let (other, page, u, v, g) = self.may_degree(right, "mayMul")?;
                require_profile(profile, other, "mayMul", right.span)?;
                if r != page {
                    return Err(Error::new(
                        right.span,
                        "mayMul requires classes on the same page",
                    ));
                }
                let s = checked_degree(s, u, "mayMul", span)?;
                let t = checked_degree(t, v, "mayMul", span)?;
                let f = checked_degree(f, g, "mayMul", span)?;
                self.model.profile_may_type(profile, r, s, t, f)
            }
            BinaryOp::Act => {
                let operation = self.operation_degree(left, "act")?;
                match *self.model.kind(right.ty) {
                    Kind::Cohomology(base, degree) => {
                        let degree = checked_degree(degree, operation, "act", span)?;
                        self.model.cohomology_type(base, degree)
                    }
                    Kind::Stable(base, shift, degree) => {
                        let degree = checked_degree(degree, operation, "act", span)?;
                        self.model.stable_type(base, shift, degree)
                    }
                    _ => Err(self.invalid_operand(right, "act", "a cohomology class")),
                }
            }
            BinaryOp::Cup => {
                let (base, p) = self.class_type(left, "cup")?;
                let (other, q) = self.class_type(right, "cup")?;
                if base != other {
                    return Err(Error::new(
                        right.span,
                        format!(
                            "cup requires classes on the same space, found {} and {}",
                            self.model.type_name(base),
                            self.model.type_name(other),
                        ),
                    ));
                }
                let degree = p
                    .checked_add(q)
                    .ok_or_else(|| Error::new(span, "cup degree is out of range"))?;
                self.model.cohomology_type(base, degree)
            }
            BinaryOp::Pull => {
                let Kind::Arrow(domain, range) = *self.model.kind(left.ty) else {
                    return Err(Error::new(
                        left.span,
                        format!(
                            "pull requires a function, found {}",
                            self.model.type_name(left.ty),
                        ),
                    ));
                };
                let (_, degree) = self.class_type(right, "pull")?;
                let expected = self.model.cohomology_type(range, degree)?;
                self.require_type(expected, right, "pull class")?;
                self.model.cohomology_type(domain, degree)
            }
        }
    }

    fn projection(&mut self, pair: ParsedExpr, first: bool) -> Result<(TypeId, TypedForm)> {
        let pair = self.expr(pair)?;
        let Kind::Product(left, right) = *self.model.kind(pair.ty) else {
            let projection = if first { "fst" } else { "snd" };
            return Err(Error::new(
                pair.span,
                format!(
                    "{projection} requires a pair, found {}",
                    self.model.type_name(pair.ty)
                ),
            ));
        };
        if first {
            Ok((left, TypedForm::Fst(Box::new(pair))))
        } else {
            Ok((right, TypedForm::Snd(Box::new(pair))))
        }
    }

    fn require_eo(&self, input: &TypedExpr, operation: &str) -> Result<()> {
        if !matches!(self.model.kind(input.ty), Kind::EOType(..)) {
            return Err(self.invalid_operand(input, operation, "an EO final type"));
        }
        Ok(())
    }

    fn scalar_prime(&mut self, value: usize) -> Result<Prime> {
        if value == 2 {
            Ok(Prime::validated(2))
        } else {
            Prime::new(self.model, value)
        }
    }

    fn check_quadratic_rows(
        &mut self,
        field_type: TypeId,
        rows: Vec<Vec<ParsedExpr>>,
        op: &str,
    ) -> Result<Vec<Vec<TypedExpr>>> {
        rows.into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|entry| {
                        let entry = self.expr(entry)?;
                        self.require_type(field_type, &entry, op)?;
                        Ok(entry)
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()
    }

    fn require_res_space(&self, input: &TypedExpr, op: &str) -> Result<()> {
        if !matches!(self.model.kind(input.ty), Kind::RVResSpace(..)) {
            return Err(self.invalid_operand(input, op, "an RVResSpace"));
        }
        Ok(())
    }

    fn qp_quad_field(&self, input: &TypedExpr, operation: &str) -> Result<QuadField> {
        let Kind::RVQpQuad(name) = self.model.kind(input.ty) else {
            return Err(self.invalid_operand(input, operation, "an RVQpQuad element"));
        };
        Ok(self
            .model
            .qp_quad_named(name)
            .expect("registered quadratic field"))
    }

    fn require_type(&self, expected: TypeId, actual: &TypedExpr, context: &str) -> Result<()> {
        if expected != actual.ty {
            return Err(Error::new(
                actual.span,
                format!(
                    "{context}: expected type {}, found {}",
                    self.model.type_name(expected),
                    self.model.type_name(actual.ty),
                ),
            ));
        }
        Ok(())
    }

    fn labeled_domain(&self, domain: TypeId, form: &str, span: Span) -> Result<Rc<Poset>> {
        self.model.poset(domain).ok_or_else(|| {
            Error::new(
                span,
                format!(
                    "{form} requires a named or Star domain, found {}",
                    self.model.type_name(domain)
                ),
            )
        })
    }

    fn table(
        &mut self,
        domain: TypeId,
        range: TypeId,
        entries: Entries,
        span: Span,
    ) -> Result<Vec<TypedExpr>> {
        let poset = self.labeled_domain(domain, "tab", span)?;
        let mut rows: Vec<Option<TypedExpr>> = (0..poset.names.len()).map(|_| None).collect();
        for (label, expr) in entries {
            let index = self.model.point(domain, &label)?;
            if rows[index].is_some() {
                return Err(Error::new(
                    label.span,
                    format!("duplicate table entry for '{}'", label.text),
                ));
            }
            let expr = self.expr(expr)?;
            self.require_type(range, &expr, "table entry")?;
            rows[index] = Some(expr);
        }
        rows.into_iter()
            .enumerate()
            .map(|(index, row)| {
                row.ok_or_else(|| {
                    Error::new(
                        span,
                        format!(
                            "table is missing point {}.{}",
                            self.model.type_name(domain),
                            poset.names[index],
                        ),
                    )
                })
            })
            .collect()
    }

    fn atlas(
        &mut self,
        domain: TypeId,
        range: TypeId,
        entries: Entries,
        span: Span,
    ) -> Result<Vec<(usize, TypedExpr)>> {
        let poset = self.labeled_domain(domain, "sew", span)?;
        let mut seen = vec![false; poset.names.len()];
        let mut covered = vec![false; poset.names.len()];
        let mut charts = Vec::new();
        for (label, expr) in entries {
            let center = self.model.point(domain, &label)?;
            if seen[center] {
                return Err(Error::new(
                    label.span,
                    format!("duplicate chart center '{}'", label.text),
                ));
            }
            seen[center] = true;
            let star = self.model.star(domain, &label)?;
            let expected = self.model.arrow(star, range)?;
            let chart = self.expr(expr)?;
            self.require_type(expected, &chart, "chart")?;
            for &point in self.model.parent_points(star).iter() {
                covered[point] = true;
            }
            charts.push((center, chart));
        }
        if let Some(index) = covered.iter().position(|covered| !covered) {
            return Err(Error::new(
                span,
                format!(
                    "atlas does not cover point {}.{}",
                    self.model.type_name(domain),
                    poset.names[index],
                ),
            ));
        }
        Ok(charts)
    }
}

fn checked_degree(left: usize, right: usize, op: &str, span: Span) -> Result<usize> {
    left.checked_add(right)
        .ok_or_else(|| Error::new(span, format!("{op} degree is out of range")))
}

fn require_restriction(source: Profile, target: Profile, span: Span) -> Result<()> {
    if source.can_restrict_to(target) {
        return Ok(());
    }
    Err(Error::new(
        span,
        "resA cannot restrict from a smaller Steenrod algebra",
    ))
}

fn require_profile(left: Profile, right: Profile, op: &str, span: Span) -> Result<()> {
    if left == right {
        return Ok(());
    }
    Err(Error::new(
        span,
        format!("{op} requires the same Steenrod algebra profile"),
    ))
}

fn next_def_order(order: usize) -> Result<usize> {
    order
        .checked_add(1)
        .ok_or_else(|| Error::plain("sheaf deformation order is out of range"))
}
