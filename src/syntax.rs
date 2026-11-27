use crate::bp::Operation as BpOperation;
use crate::bp_page::{Mode as BpPageMode, Operation as BpPageOp};
use crate::bp_vpage::Operation as BpVPageOp;
use crate::comodule::Generator;
use crate::eo::Report as EoReport;
use crate::error::{Error, Result, Span};
use crate::map_bracket::Operation as MapBracketOp;
use crate::map_system::Operation as MapSystemOp;
use crate::map_transfer::Operation as MapTransferOp;
use crate::model::DiffFamily;
use crate::motivic::Grade as MotivicGrade;
use crate::odd_cobar::FactorLiteral;
use crate::odd_massey::Operation as OddSystemOp;
use crate::rho_bock;
use crate::rho_localized;
use crate::rho_minimal::Operation as RhoMinOp;
use crate::rv_arithmetic::Series as RankSeries;
use crate::rv_as_cover::Query as AsQuery;
use crate::rv_as_derham::{Operator as AsOperator, WordQuery as AsWordQuery};
use crate::rv_as_odd::Query as AsOddQuery;
use crate::rv_as_odd_cartier::Operation as AsOddDiffOp;
use crate::rv_as_odd_derham::Operation as AsOddDeRhamOp;
use crate::rv_as_odd_tate::Operation as AsOddTateOp;
use crate::rv_gauss::Operation as GaussOp;
use crate::rv_rees::Operation as ReesOp;
use crate::rv_t_bracket::Operation as RvTBracketOp;
use crate::rv_t_system::Operation as RvTSystemOp;
use crate::rv_t_transfer::Operation as RvTTransferOp;
use crate::rv_witt_cover::Query as WittQuery;
use crate::rv_witt_tower::Query as TowerQuery;
use crate::sheaf::{self, Form as SheafForm, Operation as SheafOp};
use crate::sheaf_image::Operation as SheafImageOp;
use crate::sheaf_leray::InputForm as SheafInputForm;
use crate::sheaf_page::{Grade as SheafPageGrade, Lift as SheafPageLift, Operation as SheafPageOp};
use crate::sheaf_support::{Mode as SupportMode, Operation as SupportOp};
use crate::transfer::Operation as TransferOp;
use crate::v_bars::Certificate as BarCertificate;

const MAX_NESTING: usize = 64;
pub(crate) const MAX_IDENTIFIER_BYTES: usize = 128;

#[derive(Clone)]
pub(crate) struct Label {
    pub(crate) text: String,
    pub(crate) span: Span,
}

pub(crate) struct ParsedType {
    pub(crate) span: Span,
    pub(crate) kind: TypeForm,
    depth: usize,
}

pub(crate) enum TypeForm {
    Name(Label),
    Star(Box<ParsedType>, Label),
    Cohomology(Box<ParsedType>, usize),
    MixedChain(Label, i64, bool),
    MixedMinimal(Label, i64, bool),
    MixedPage(MixedSpec),
    Sheaf(SheafForm, Label, usize),
    SheafHom(SheafForm, Label, Label, usize),
    SheafComplex(SheafForm, Label, usize),
    SheafPage(Label, SheafPageGrade),
    SheafLocal(Label, Label, SupportMode, usize),
    Operations(usize),
    Cobar(usize, usize),
    MCobar(MotivicGrade),
    RMCobar(MotivicGrade),
    RhoBock(rho_bock::Grade),
    RhoLoc(MotivicGrade),
    RhoMin(MotivicGrade),
    MExt(MotivicGrade),
    RMExt(MotivicGrade),
    RMTriple(MotivicGrade),
    PCobar(usize, usize, usize),
    PCobarAn(usize, usize, usize, usize),
    PExt(usize, usize, usize),
    PExtAn(usize, usize, usize, usize),
    PMay(usize, usize, usize, usize, usize),
    PMayAn(usize, usize, usize, usize, usize, usize),
    PCeAn(usize, usize, usize, usize, usize, usize),
    CobarAn(usize, usize, usize),
    Ext(usize, usize),
    ExtAn(usize, usize, usize),
    XExt(Box<ParsedType>, usize, usize),
    MapCobar(Box<ParsedType>, Box<ParsedType>, usize, i64),
    MapExt(Box<ParsedType>, Box<ParsedType>, usize, i64),
    MapMay(MapMaySpec),
    QKoszul(QSpec),
    QExt(QSpec),
    RVHom(RankVarietySpec),
    RVLaurent(RankVarietySpec),
    RVIwHecke(usize),
    RVSphHecke(usize),
    RVSatake(usize),
    RVHive(HiveSpec),
    RVQGr(usize, usize),
    RVQGrP(usize, usize, usize),
    RVOre,
    RVFusion(usize, usize),
    RVCyclo(usize),
    RVCycloQ(usize),
    RVWzw(usize, usize),
    JetHom(JetHomSpec),
    JetPage(JetPageSpec),
    VKoszul(QSpec<i64>),
    VExt(QSpec<i64>),
    VBock(VBockSpec),
    RVPoly(Vec<usize>, usize),
    RVBarE(RvBarPageSpec),
    RVUntruncatedCyclic(RvHHSpec, crate::rv_cyclic_homology::Mode, bool),
    RVNegativeCyclic(RvNCSpec, bool),
    RVHochschild(RvHHSpec, bool),
    RVBarCobar(RvBarSpec),
    RVBarDual(RvBarSpec),
    RVPathExt(RvBarSpec),
    RVBar(RvBarSpec),
    CurveDiff(DiffFamily, Label, usize),
    ASOddDiff(Label),
    ASOddDeRham(Label),
    ASOddTate(Label, i64),
    PrimeScalar(usize),
    ASDeRham(Label, usize),
    BinaryScalar(usize),
    EOType(usize),
    RVComplex(RvGradeSpec, bool),
    BP(BpSpec, Option<usize>),
    BPThick(BpSpec, usize, usize),
    BPModExt(BpModSpec),
    BPRel(BpRelSpec),
    BPRelExt(BpRelSpec, usize),
    BPHeight(BpHeightSpec),
    BpPage(BpPageMode, BpPageSpec),
    BPVPage(BpVPageSpec),
    BPExt(BpSpec, usize),
    BPTriple(BpSpec, usize),
    May(usize, usize, usize, usize),
    MayAn(usize, usize, usize, usize, usize),
    CeAn(usize, usize, usize, usize, usize),
    Suspension(Box<ParsedType>),
    Stable(Box<ParsedType>, usize, usize),
    Bar(usize),
    CyclicH2(usize),
    Product(Box<ParsedType>, Box<ParsedType>),
    Arrow(Box<ParsedType>, Box<ParsedType>),
}

#[derive(Clone, Copy)]
pub(crate) enum RvIdealOp {
    Normal,
    Contains,
    Witness(usize),
}

pub(crate) type RvNCSpec = (Label, usize, i64, Vec<usize>);
pub(crate) type HiveSpec = [Vec<i128>; 3];

#[derive(Clone, Copy)]
pub(crate) enum WzwOp {
    S(bool),
    T(i128),
    Charge,
    Negate,
}

impl WzwOp {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::S(false) => "rvWzwS",
            Self::S(true) => "rvWzwSAdjoint",
            Self::T(_) => "rvWzwT",
            Self::Charge => "rvWzwCharge",
            Self::Negate => "neg",
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum QuantumOp {
    Scale(i128),
    Shift(usize),
    Transpose,
    Classical,
}

impl QuantumOp {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Scale(_) => "rvQGrScale",
            Self::Shift(_) => "rvQGrShift",
            Self::Transpose => "rvQGrTranspose",
            Self::Classical => "rvQGrClassical",
        }
    }
}

pub(crate) type RvHHSpec = (Label, i64, Vec<usize>);

pub(crate) type MixedSpec = (Label, usize, i64, usize);

pub(crate) struct BpSpec(pub(crate) usize, pub(crate) usize, pub(crate) usize);

pub(crate) struct BpHeightSpec(
    pub(crate) BpSpec,
    pub(crate) usize,
    pub(crate) i128,
    pub(crate) usize,
);

#[derive(Clone, Copy)]
pub(crate) enum BpHeightOp {
    Class(usize, i128, usize),
    Representative,
    Differential,
    Next,
    Previous,
}

impl BpHeightOp {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Class(..) => "bpHeightClass",
            Self::Representative => "bpHeightRep",
            Self::Differential => "bpHeightD",
            Self::Next => "bpHeightNext",
            Self::Previous => "bpHeightPrev",
        }
    }
}

pub(crate) struct BpRelSpec(
    pub(crate) usize,
    pub(crate) usize,
    pub(crate) usize,
    pub(crate) usize,
);

#[derive(Clone, Copy)]
pub(crate) enum BpRelOp {
    Quotient(usize),
    Lift,
    Differential,
    Null(usize),
    Class(usize),
    Representative,
    HeightDelta,
}

pub(crate) struct BpModSpec {
    pub(crate) shape: BpSpec,
    pub(crate) h: usize,
    pub(crate) exponent: Option<usize>,
    pub(crate) t: usize,
}

pub(crate) type BpPageSpec = (usize, usize, usize, usize, usize, usize);
pub(crate) type BpVPageSpec = (usize, usize, usize, usize, usize, usize, usize);

#[derive(Clone, Copy)]
pub(crate) enum BpExtOp {
    Class(usize),
    Representative,
    Scale(i128),
}

pub(crate) struct RvBarSpec {
    pub(crate) path: Label,
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) b: usize,
}

pub(crate) struct RvBarPageSpec {
    pub(crate) chain: RvBarSpec,
    pub(crate) r: usize,
    pub(crate) length: i64,
}

pub(crate) struct RvGradeSpec {
    pub(crate) ideal: Label,
    pub(crate) s: usize,
    pub(crate) degree: RvDegree,
}

pub(crate) enum RvDegree {
    Total(usize),
    Multi(Vec<usize>),
}

#[derive(Clone, Copy)]
pub(crate) enum RvSchemeOp {
    Colon,
    Saturate,
}

pub(crate) struct RvMatrixSpec {
    pub(crate) base: Box<ParsedType>,
    pub(crate) selected: Vec<usize>,
}

pub(crate) struct RankVarietySpec {
    pub(crate) source: Box<ParsedType>,
    pub(crate) target: Box<ParsedType>,
    pub(crate) selected: Vec<usize>,
    pub(crate) modulus: usize,
    pub(crate) point: Option<Vec<usize>>,
}

pub(crate) struct RankJetSpec {
    pub(crate) source: Box<ParsedType>,
    pub(crate) target: Box<ParsedType>,
    pub(crate) path: JetPathArg,
}

pub(crate) struct JetPathSpec {
    pub(crate) selected: Vec<usize>,
    pub(crate) modulus: usize,
    pub(crate) order: usize,
    pub(crate) coefficients: Vec<Vec<usize>>,
}

pub(crate) struct RankArithmeticSpec {
    pub(crate) source: Box<ParsedType>,
    pub(crate) target: Box<ParsedType>,
    pub(crate) selected: Vec<usize>,
    pub(crate) max_degree: usize,
}

pub(crate) enum JetPathArg {
    Inline(JetPathSpec),
    Named(Label),
}

pub(crate) struct JetHomSpec {
    pub(crate) source: Box<ParsedType>,
    pub(crate) target: Box<ParsedType>,
    pub(crate) path: Label,
    pub(crate) parity: usize,
}

pub(crate) struct JetPageSpec {
    pub(crate) chain: JetHomSpec,
    pub(crate) r: usize,
    pub(crate) f: usize,
}

#[derive(Clone, Copy)]
pub(crate) enum RealBarOp {
    Generator,
    Representative,
    Witness,
}

#[derive(Clone, Copy)]
pub(crate) enum RhoLocOp {
    Localize,
    Project,
    Section(i128),
    Homotopy,
}

impl RhoLocOp {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Localize => "rhoLocalize",
            Self::Project => "rhoProject",
            Self::Section(_) => "rhoSection",
            Self::Homotopy => "rhoHomotopy",
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum RhoOp {
    Class(usize, usize),
    Representative,
    Differential,
    Next,
    Previous,
}

impl RhoOp {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Class(..) => "rhoClass",
            Self::Representative => "rhoRep",
            Self::Differential => "rhoD",
            Self::Next => "rhoNext",
            Self::Previous => "rhoPrev",
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum JetOp {
    Class(usize, usize),
    Representative,
    Differential,
    Next,
    Prev,
}

pub(crate) struct JetTerm {
    pub(crate) source: MapGenerator,
    pub(crate) target: MapGenerator,
    pub(crate) path: Label,
    pub(crate) power: usize,
}

pub(crate) struct QSpec<S = usize> {
    pub(crate) source: Box<ParsedType>,
    pub(crate) target: Box<ParsedType>,
    pub(crate) selected: Vec<usize>,
    pub(crate) s: S,
    pub(crate) t: i64,
}

pub(crate) struct VBockSpec {
    pub(crate) localized: QSpec<i64>,
    pub(crate) r: usize,
    pub(crate) f: usize,
}

pub(crate) struct BarFamilySpec {
    pub(crate) source: ParsedType,
    pub(crate) target: ParsedType,
    pub(crate) indices: Vec<usize>,
}

pub(crate) struct MapMaySpec {
    pub(crate) source: Box<ParsedType>,
    pub(crate) target: Box<ParsedType>,
    pub(crate) r: usize,
    pub(crate) s: usize,
    pub(crate) t: i64,
    pub(crate) f: usize,
}

pub(crate) struct ParsedExpr {
    pub(crate) span: Span,
    pub(crate) kind: ExprForm,
    depth: usize,
}

pub(crate) type Entries = Vec<(Label, ParsedExpr)>;

#[derive(Clone, Copy)]
pub(crate) enum BinaryOp {
    Add,
    Cup,
    Pull,
    Mul,
    Act,
    Concat,
    Yoneda,
    XAction,
    MapExtCompose,
    XMap,
    MayMul,
    RhoMul,
    CeMul,
    MilnorPair,
    RvMultiply,
    RvBCMul,
    MixEval,
    RvBarPair,
    AsOddDeRhamPair,
    AsOddTatePair,
    AsDeRhamPair,
    EOBruhat,
    EOClosure,
    RvPathCup,
    RvKoszulProduct,
    RvTorProduct,
    BpMultiply,
    BpCup,
}

impl BinaryOp {
    fn from_name(name: &str) -> Option<Self> {
        match name {
            "add" => Some(Self::Add),
            "cup" => Some(Self::Cup),
            "pull" => Some(Self::Pull),
            "mul" => Some(Self::Mul),
            "act" => Some(Self::Act),
            "concat" => Some(Self::Concat),
            "yoneda" => Some(Self::Yoneda),
            "xaction" => Some(Self::XAction),
            "ecomp" => Some(Self::MapExtCompose),
            "xmap" => Some(Self::XMap),
            "mayMul" => Some(Self::MayMul),
            "rhoMul" => Some(Self::RhoMul),
            "ceMul" => Some(Self::CeMul),
            "milnorPair" => Some(Self::MilnorPair),
            "rvMul" => Some(Self::RvMultiply),
            "rvBCMul" => Some(Self::RvBCMul),
            "mixEval" => Some(Self::MixEval),
            "rvBarPair" => Some(Self::RvBarPair),
            "rvASOddDeRhamPair" => Some(Self::AsOddDeRhamPair),
            "rvASOddTatePair" => Some(Self::AsOddTatePair),
            "rvASDeRhamPair" => Some(Self::AsDeRhamPair),
            "eoBruhat" => Some(Self::EOBruhat),
            "eoClosure" => Some(Self::EOClosure),
            "rvPathCup" => Some(Self::RvPathCup),
            "rvkProd" => Some(Self::RvKoszulProduct),
            "rvTorProd" => Some(Self::RvTorProduct),
            "bpMul" => Some(Self::BpMultiply),
            "bpCup" => Some(Self::BpCup),
            _ => None,
        }
    }
}

pub(crate) enum ClassForm {
    Zero,
    Generator(usize, Span),
    Cocycle(Vec<Vec<Label>>),
}

pub(crate) enum SuspensionForm {
    North,
    South,
    Equator,
}

pub(crate) enum CyclicForm {
    Support(Vec<(usize, Span)>),
    Binomial(usize, Span),
    Carry,
}

pub(crate) struct Witness {
    pub(crate) from: usize,
    pub(crate) to: usize,
    pub(crate) span: Span,
    pub(crate) value: ParsedExpr,
}

pub(crate) fn validate_witness_interval(
    arity: usize,
    from: usize,
    to: usize,
    span: Span,
) -> Result<()> {
    if from >= to || to > arity {
        return Err(Error::new(span, format!("invalid interval [{from},{to}]")));
    }
    if to - from == 1 {
        return Err(Error::new(
            span,
            format!("adjacent interval [{from},{to}] already has an input"),
        ));
    }
    if to - from == arity {
        return Err(Error::new(
            span,
            format!("top interval [{from},{to}] cannot have a witness"),
        ));
    }
    Ok(())
}

pub(crate) struct QTerm {
    pub(crate) source: MapGenerator,
    pub(crate) target: MapGenerator,
    pub(crate) selected: Vec<usize>,
    pub(crate) powers: Vec<usize>,
}

pub(crate) struct MapGenerator {
    pub(crate) base: ParsedType,
    pub(crate) degree: usize,
    pub(crate) index: usize,
    pub(crate) span: Span,
}

pub(crate) enum ExprForm {
    Name(Label),
    Point(ParsedType, Label),
    Class(ParsedType, usize, ClassForm),
    Zero(ParsedType),
    Square(usize),
    Milnor(Vec<usize>),
    Coact(usize, usize, Span, Box<ParsedExpr>),
    CobarWord(Vec<(Vec<usize>, Span)>),
    MCobarWord(usize, Vec<FactorLiteral>),
    RMCobarWord(usize, usize, Vec<FactorLiteral>),
    RhoGenerator(rho_bock::Grade, usize, Span),
    RhoOperation(RhoOp, Box<ParsedExpr>),
    RhoLocWord(i128, usize, Vec<rho_localized::FactorLiteral>),
    RhoLocOperation(RhoLocOp, Box<ParsedExpr>),
    RhoMinOperation(RhoMinOp, Box<ParsedExpr>),
    TauShift(usize, Box<ParsedExpr>),
    RhoShift(i128, Box<ParsedExpr>),
    RealComplex(Box<ParsedExpr>),
    RhoNull(usize, Box<ParsedExpr>),
    RealMassey(bool, Box<[ParsedExpr; 3]>),
    RealMasseyNull(Box<[ParsedExpr; 2]>),
    RealTripleClass(Option<usize>, Box<ParsedExpr>),
    Realize(Box<ParsedExpr>),
    PCobarWord(usize, Vec<FactorLiteral>),
    PExtGenerator(usize, usize, usize, usize, Span),
    PExtAnGenerator(usize, usize, usize, usize, usize, Span),
    PMayGenerator(usize, usize, usize, usize, usize, usize, Span),
    PMayAnGenerator(usize, usize, usize, usize, usize, usize, usize, Span),
    PCeAnGenerator(usize, usize, usize, usize, usize, usize, usize, Span),
    MixedLiteral(Label, i64, crate::mixed::Terms),
    MixedUnit(Label, Label, Label),
    MixedIdentity(Label),
    MixedConeMap(Label, crate::mixed_cone::Operation, Box<ParsedExpr>),
    MixedConeShear(Label, Box<ParsedExpr>, Box<ParsedExpr>),
    MixedTensorProduct(Label, Box<ParsedExpr>, Box<ParsedExpr>),
    MixedTensorSwap(Label, Box<ParsedExpr>),
    MixedTensorCertificate(Label, usize, usize, crate::mixed_tensor::Certificate),
    MixedComposition(Label, Box<ParsedExpr>, Box<ParsedExpr>),
    MixedEmbed(usize, Box<ParsedExpr>),
    RvHHToMixed(Label, Option<usize>, Box<ParsedExpr>),
    MixedBarTerm(Label, usize, bool),
    MixedMinGenerator(Label, i64, usize, Span),
    MixedMinMap(Label, crate::mixed_minimal::Operation, Box<ParsedExpr>),
    MixedDelta(usize, Box<ParsedExpr>),
    MixedGenerator(MixedSpec, usize, Span),
    MixedClass(usize, usize, Box<ParsedExpr>),
    SheafTerm(Label, Vec<Label>, usize, i64),
    SheafCellTerm(Label, usize, Label, usize, i64),
    SheafGenerator(Label, usize, usize),
    SheafLerayIncl(Label, Box<ParsedExpr>),
    SheafLocalGenerator(Label, Label, SupportMode, usize, usize),
    SheafLocalClass(Label, Label, SupportMode, Box<ParsedExpr>),
    SheafLocalOperation(SupportOp, Box<ParsedExpr>),
    SheafApply(Label, Box<ParsedExpr>),
    SheafEmbed(Label, usize, Box<ParsedExpr>),
    SheafPart(usize, Box<ParsedExpr>),
    SheafHyperGenerator(Label, usize, usize),
    SheafPageGenerator(Label, SheafPageGrade, usize),
    SheafPageOperation(SheafPageOp, Box<ParsedExpr>),
    SheafPageLift(SheafPageLift, usize, usize, Box<ParsedExpr>),
    SheafOperation(SheafOp, Box<ParsedExpr>),
    SheafHomTerm(Label, Label, Vec<Label>, usize, usize, i64),
    SheafExtGenerator(Label, Label, usize, usize),
    SheafHomMap(Label),
    SheafHomIdentity(Label),
    SheafCompose(Box<[ParsedExpr; 2]>),
    PScale(usize, Box<ParsedExpr>),
    RVScale(usize, Box<ParsedExpr>),
    RVLaurentMonomial(i64, Box<ParsedExpr>),
    HeckeBasis(bool, Vec<i128>),
    HeckeScalar(usize, usize),
    HeckeNegate(Box<ParsedExpr>),
    HeckeProduct(Box<[ParsedExpr; 2]>),
    SphericalBasis(Vec<i128>),
    SphericalScalar(usize, usize),
    SphericalProduct(Box<[ParsedExpr; 2]>),
    SphericalEmbedding(Box<ParsedExpr>),
    SatakeTransform(Box<ParsedExpr>),
    SatakeInverse(Box<ParsedExpr>),
    SatakeCharacter(bool, Vec<i128>),
    SatakeScalar(usize, i128),
    Hive(HiveSpec, Option<Vec<i128>>),
    QuantumBasis(usize, usize, Vec<usize>),
    QuantumProduct(Box<[ParsedExpr; 2]>),
    QuantumNabla(Box<ParsedExpr>),
    QuantumReduce(usize, Box<ParsedExpr>),
    QuantumPSection(usize, usize, usize, Vec<usize>, bool),
    QuantumCurvature(bool, Box<ParsedExpr>),
    QuantumAnnihilator(usize, usize),
    OreTerm(usize, usize, i128),
    OreProduct(Box<[ParsedExpr; 2]>),
    OreApply(Box<[ParsedExpr; 2]>),
    QuantumUnary(QuantumOp, Box<ParsedExpr>),
    FusionBasis(usize, Vec<i128>),
    FusionOne(usize, usize),
    FusionHandle(usize, usize),
    Verlinde(Option<Vec<i128>>, Box<ParsedExpr>),
    CycloRoot(usize, i128),
    CycloProduct(Box<[ParsedExpr; 2]>),
    CycloScale(i128, Box<ParsedExpr>),
    CycloConjugate(Box<ParsedExpr>),
    CycloEmbed(usize, Box<ParsedExpr>),
    WzwBasis(usize, usize, Vec<i128>),
    WzwSEntry(usize, usize, Vec<i128>, Vec<i128>),
    WzwTEntry(usize, usize, Vec<i128>, i128),
    WzwUnary(WzwOp, Box<ParsedExpr>),
    WzwScale(Box<[ParsedExpr; 2]>),
    WzwCoordinate(Vec<i128>, Box<ParsedExpr>),
    LensTv(usize, usize, usize, i128),
    WzwTv(usize, usize, Vec<crate::rv_lens::Step>),
    CycloRatio(u128, Box<ParsedExpr>),
    FusionProduct(Box<[ParsedExpr; 2]>),
    FusionScale(i128, Box<ParsedExpr>),
    FusionDual(Box<ParsedExpr>),
    SatakeProduct(Box<[ParsedExpr; 2]>),
    RVIwAdlv(
        crate::rv_iwahori::Boundary,
        usize,
        usize,
        Vec<i128>,
        Box<[ParsedExpr; 2]>,
    ),
    RVNewton(crate::rv_newton::Operation, usize, usize, Box<ParsedExpr>),
    RVAdlv(
        crate::rv_hodge::Boundary,
        usize,
        usize,
        Vec<i64>,
        Box<[ParsedExpr; 2]>,
    ),
    JetOperation(JetOp, Box<ParsedExpr>),
    JetGenerator(JetPageSpec, usize, Span),
    JetTerm(JetTerm),
    ExtGenerator(usize, usize, usize, Span),
    MExtGenerator(MotivicGrade, usize, Span),
    RMExtGenerator(MotivicGrade, usize, Span),
    RealBar(RealBarOp, usize, i128, usize, Span),
    ExtAnGenerator(usize, usize, usize, usize, Span),
    RestrictAlgebra(usize, Box<ParsedExpr>),
    CobarCup(usize, Box<ParsedExpr>, Box<ParsedExpr>),
    ExtSquare(usize, Box<ParsedExpr>),
    Nakamura(usize, Box<ParsedExpr>),
    MayLiftTo(usize, usize, Box<ParsedExpr>),
    HigherRestriction(usize, Vec<ParsedExpr>),
    XExtGenerator(ParsedType, usize, usize, usize, Span),
    MapWord(MapGenerator, MapGenerator, Box<ParsedExpr>),
    MapExtGenerator(ParsedType, ParsedType, usize, i64, usize, Span),
    MapExtIdentity(ParsedType),
    MapMayGenerator(MapMaySpec, usize, Span),
    QExtGenerator(QSpec, usize, Span),
    QTerm(QTerm),
    QIdentity(ParsedType, Vec<usize>),
    RVGenerator(RankVarietySpec, usize, Span),
    RVIdentity(ParsedType, Vec<usize>, usize, Vec<usize>),
    RVSpecialize(usize, Vec<usize>, Box<ParsedExpr>),
    RVField(crate::rv_field::Operation, Box<ParsedExpr>),
    RVInverse(Box<ParsedExpr>),
    RVHilbert90(usize, usize, Box<ParsedExpr>),
    RVShintani(crate::rv_shintani::Operation, usize, usize, Box<ParsedExpr>),
    RVTwisted(crate::rv_descent::Twist, usize, usize, Box<[ParsedExpr; 3]>),
    VExtGenerator(QSpec<i64>, usize, Span),
    VBockGenerator(VBockSpec, usize, Span),
    VBockClass(usize, usize, Box<ParsedExpr>),
    RvCyclicGenerator(RvHHSpec, crate::rv_cyclic_homology::Mode, usize, Span),
    RvCyclicTerm(crate::rv_cyclic_homology::Mode, i64, Box<ParsedExpr>),
    RvCyclicPart(i64, Box<ParsedExpr>),
    RvCyclicShift(crate::rv_cyclic_homology::Shift, Box<ParsedExpr>),
    RvNCGenerator(RvNCSpec, usize, Span),
    RvNCTerm(usize, usize, Box<ParsedExpr>),
    RvNCPart(usize, Box<ParsedExpr>),
    RvNCTransfer(crate::rv_negative::Transfer, Box<ParsedExpr>),
    RvHHGenerator(RvHHSpec, usize, Span),
    RvPathGenerator(RvBarSpec, usize, Span),
    RvBarEGenerator(RvBarPageSpec, usize, Span),
    RvBarClass(usize, i64, Box<ParsedExpr>),
    RvHHWord(Vec<ParsedExpr>),
    RvHHLength(usize, Box<ParsedExpr>),
    RvBarWord(Label, Vec<usize>, Vec<ParsedExpr>),
    RvBarLength(usize, Box<ParsedExpr>),
    RvHFace(Label, Vec<usize>, Vec<usize>),
    RvMTerm(Label, Vec<usize>, Box<ParsedExpr>),
    RvMComponent(Vec<usize>, Box<ParsedExpr>),
    RvKoszulTerm(Label, Vec<usize>, Box<ParsedExpr>),
    RvTorGenerator(RvGradeSpec, usize, Span),
    RvMonomial(Vec<usize>, Vec<usize>),
    BpInteger(BpSpec, i128),
    BpVariable(BpSpec, Option<usize>, usize),
    BpOperation(BpOperation, Box<ParsedExpr>),
    BpRelOperation(BpRelOp, Box<ParsedExpr>),
    BpRelGenerator(BpRelSpec, usize, usize, Span),
    BpHeightGenerator(BpHeightSpec, usize, Span),
    BpHeightOperation(BpHeightOp, Box<ParsedExpr>),
    BpExtOperation(BpExtOp, Box<ParsedExpr>),
    BpExtGenerator(BpSpec, usize, usize, Span),
    BpMassey(Vec<ParsedExpr>),
    BpHigher(Option<usize>, Vec<ParsedExpr>),
    BpTriple(Box<[ParsedExpr; 3]>),
    BpTripleBase(Box<ParsedExpr>),
    BpRaise(usize, Box<ParsedExpr>),
    BpPageGenerator(BpPageMode, BpPageSpec, usize, Span),
    BpPageClass(BpPageMode, (usize, usize, usize), Box<ParsedExpr>),
    BpModExtGen(BpModSpec, usize, Span),
    BpModExtClass(bool, usize, Box<ParsedExpr>),
    BpModExtRep(bool, Box<ParsedExpr>),
    BpPageOperation(Label, BpPageMode, BpPageOp, Box<ParsedExpr>),
    BpVPageGenerator(BpVPageSpec, usize, Span),
    BpVPageClass((usize, usize, usize), Box<ParsedExpr>),
    BpVPageOperation(Label, BpVPageOp, Box<ParsedExpr>),
    RvMinor(RvMatrixSpec, Vec<Generator>, Vec<Generator>),
    RvHasse(Vec<usize>, Box<ParsedExpr>),
    RvVanish(usize, Vec<usize>, Box<ParsedExpr>),
    RvInMinors(RvMatrixSpec, usize, Box<ParsedExpr>),
    RvIdeal(RvIdealOp, Label, Box<ParsedExpr>),
    RvGenerator(Label, usize, Span),
    RvOscMinor(Label, usize, Vec<usize>),
    RvFrobeniusWronskian(Label, usize),
    AsDeRhamGen(Label, usize, usize, Span),
    BinaryScalar(usize, usize),
    EOType(Vec<usize>),
    EOCover(Label),
    EOOddCover(Label),
    AsDeRhamOperation(crate::rv_as_derham::Operation, Box<ParsedExpr>),
    AsOddCoefficient(Label, usize, usize),
    AsOddOmega(Label, usize, usize),
    AsOddDiffOperation(AsOddDiffOp, Box<ParsedExpr>),
    AsOddDeRhamGen(Label, usize),
    AsOddTateGen(Label, i64, usize),
    AsOddTateOperation(AsOddTateOp, Box<ParsedExpr>),
    AsOddDeRhamOperation(AsOddDeRhamOp, Box<ParsedExpr>),
    PrimeScalar(usize, i128),
    RvOmega(DiffFamily, Label, usize, usize),
    RvDiffScale(DiffFamily, usize, Box<ParsedExpr>),
    RvCartier(DiffFamily, Box<ParsedExpr>),
    VShift(i64, Box<ParsedExpr>),
    VRatio(usize, Box<ParsedExpr>),
    VBarCertificate(BarCertificate, BarFamilySpec, i64, usize, usize),
    QRestrict(Vec<usize>, Box<ParsedExpr>),
    MayClass(usize, usize, Box<ParsedExpr>),
    MayGenerator(usize, usize, usize, usize, usize, Span),
    MayAnGenerator(usize, usize, usize, usize, usize, usize, Span),
    CeGenerator(usize, usize, usize, usize, usize, usize, Span),
    CeClass(usize, usize, Box<ParsedExpr>),
    Suspension(ParsedType, SuspensionForm),
    Shift(usize, Box<ParsedExpr>),
    Cyclic(usize, Span, CyclicForm),
    Massey(Vec<ParsedExpr>, Vec<Witness>),
    PowerDef(Vec<ParsedExpr>),
    Bracket(Vec<ParsedExpr>),
    Transfer(TransferOp, Vec<ParsedExpr>),
    RvTBracket(RvTBracketOp, Vec<ParsedExpr>),
    MapBracket(MapBracketOp, Vec<ParsedExpr>),
    RvTTransfer(RvTTransferOp, Vec<ParsedExpr>),
    RvTSystem(RvTSystemOp, Vec<ParsedExpr>, Vec<Witness>),
    MapSystem(MapSystemOp, Vec<ParsedExpr>, Vec<Witness>),
    OddSystem(OddSystemOp, Vec<ParsedExpr>, Vec<Witness>),
    MapTransfer(MapTransferOp, Vec<ParsedExpr>),
    Binary(BinaryOp, Box<ParsedExpr>, Box<ParsedExpr>),
    Lambda(Label, ParsedType, Box<ParsedExpr>),
    Call(Box<ParsedExpr>, Box<ParsedExpr>),
    Pair(Box<ParsedExpr>, Box<ParsedExpr>),
    Fst(Box<ParsedExpr>),
    Snd(Box<ParsedExpr>),
    Table(ParsedType, ParsedType, Entries),
    Restrict(Label, Box<ParsedExpr>),
    Sew(ParsedType, ParsedType, Entries),
}

pub(crate) struct Statement {
    pub(crate) span: Span,
    pub(crate) kind: StatementForm,
}

pub(crate) enum StatementForm {
    Space(Label, Vec<Vec<Label>>),
    MixedComplex(
        Label,
        usize,
        Vec<(Label, i64)>,
        crate::mixed::Rules,
        crate::mixed::Rules,
    ),
    RvHHMixed(Label, Label, Vec<usize>),
    MixedBars(Label),
    MixedBarRank(Label, i64, Option<i64>),
    MixedPeriodicRank(Label, i64),
    MixedHigher(Label, usize, Vec<(Label, i64)>, crate::mixed::HigherRules),
    MixedCone(Label, ParsedExpr),
    MixedTensor(Label, Label, Label),
    MixedHom(Label, Label, Label),
    MixedMinimal(Label, Label),
    MixedMinDimension(Label, i64),
    MixedDimension(MixedSpec),
    Sheaf(Label, Label, usize, sheaf::Presentation),
    SheafMap(Label, Label, Label, Vec<(Label, Vec<i64>)>),
    SheafExtDimension(Label, Label, usize),
    SheafComplex(Label, Vec<Label>),
    SheafResolution(Label, Label, usize, crate::sheaf_cellular::Presentation),
    SheafCellular(Label, Label, Label),
    SheafHyperDimension(Label, usize),
    SheafPageDimension(Label, SheafPageGrade),
    SheafExtra(SheafInputForm, Label, Label, Label),
    SheafLocalDimension(Label, Label, SupportMode, usize),
    SheafDimension(Label, usize),
    SheafSerre(Label, Label),
    SheafSerreShift(Label),
    SheafStalkCohomology(Label, Label, usize),
    SheafLeray(SheafInputForm, Label, Label, Label),
    SheafTruncate(Label, Label, usize),
    SheafDeligne(SheafInputForm, Label, Label, Label, usize),
    SheafImage(SheafImageOp, Label, Label, Label),
    QModule(Label, Vec<usize>, crate::q_module::Presentation),
    Let(Label, ParsedExpr),
    Show(ParsedExpr),
    RVInvariantFactors(usize, usize, ParsedExpr),
    RVNewton(usize, usize, ParsedExpr),
    RVHodge(ParsedExpr),
    RVIw(bool, ParsedExpr),
    RVIwKL(Vec<i128>, Vec<i128>),
    HeckeCoefficient(bool, Option<i128>, Vec<i128>, ParsedExpr),
    SatakeQuery(Option<i128>, Vec<i128>, ParsedExpr),
    HiveReport(HiveSpec, bool),
    HiveVertex(usize, usize, ParsedExpr),
    QuantumCoefficient(usize, Vec<usize>, ParsedExpr),
    GromovWitten(usize, usize, usize, [Vec<usize>; 3]),
    OreCoefficient(usize, usize, ParsedExpr),
    QuantumPeriod(usize, usize, usize),
    QuantumCurvature(crate::rv_pcurv::Report, usize, usize, usize),
    FusionCoefficient(Vec<i128>, ParsedExpr),
    CycloCoefficient(usize, ParsedExpr),
    ConformalDimension(usize, usize, usize, Vec<Vec<i128>>),
    RVIwSchubertKL(Vec<i128>, ParsedExpr),
    RVIwPosition(Box<[ParsedExpr; 2]>),
    RVIwAt(usize, usize, Box<[ParsedExpr; 2]>),
    RVTree(crate::rv_tree::Report, usize, usize, usize, ParsedExpr),
    RVPosition(Box<[ParsedExpr; 2]>),
    RVMuAt(usize, usize, Box<[ParsedExpr; 2]>),
    Card(ParsedExpr),
    Assert(ParsedExpr, ParsedExpr),
    Betti(ParsedType),
    BpExt(BpSpec, usize),
    BpRelExt(BpRelSpec, usize),
    BpOrder(usize, ParsedExpr),
    BpTripleDefined(Box<[ParsedExpr; 3]>),
    BpTripleContains(Box<[ParsedExpr; 2]>),
    BpTripleSize(ParsedExpr),
    ExtDimension(usize, usize),
    MExtDimension(MotivicGrade),
    RMExtDimension(MotivicGrade),
    RhoOrder(Box<ParsedExpr>),
    RealBars(usize, i128),
    RealBarRank(usize, i128, i128, usize),
    RealTripleDim(Box<ParsedExpr>),
    RealTripleContains(Box<[ParsedExpr; 2]>),
    RhoDimension(rho_bock::Grade),
    PExtDimension(usize, usize, usize),
    PExtAnDimension(usize, usize, usize, usize),
    BpPageDimension(BpPageMode, BpPageSpec),
    BpVPageDimension(BpVPageSpec),
    BpModExtDimension(BpModSpec),
    BpHeightPage(bool, BpPageSpec),
    PMayDimension(usize, usize, usize, usize, usize),
    PMayAnDimension(usize, usize, usize, usize, usize, usize),
    PCeAnDimension(usize, usize, usize, usize, usize, usize),
    ExtAnDimension(usize, usize, usize),
    XExtDimension(ParsedType, usize, usize),
    MapExtDimension(ParsedType, ParsedType, usize, i64),
    MapMayDimension(MapMaySpec),
    QExtDimension(QSpec),
    RankVariety(RankVarietySpec),
    RankOrbit(usize, Vec<usize>),
    RankArithmetic(RankSeries, RankArithmeticSpec),
    RvBarPath(Label, Label, Vec<Vec<usize>>),
    RvCyclicDimension(RvHHSpec, crate::rv_cyclic_homology::Mode),
    RvNCDimension(RvNCSpec),
    RvHHDimension(RvHHSpec),
    RvPathDimension(RvBarSpec),
    RvBarPageDimension(RvBarPageSpec),
    RvBarDimension(RvBarSpec, Option<usize>),
    RvHochster(Label, Vec<usize>, i64),
    RvBetti(RvGradeSpec),
    RvMinors(RvMatrixSpec, usize),
    RvIdeal(Label, Vec<usize>, Vec<ParsedExpr>),
    RvRees(ReesOp, Label, Label, Vec<usize>, Vec<ParsedExpr>),
    RvGauss(GaussOp, Label, Vec<usize>, ParsedExpr),
    RvCurve(Label, Vec<usize>, Vec<ParsedExpr>),
    RvAsCover(Label, Vec<usize>, Vec<ParsedExpr>),
    RvAsQuery(AsQuery, Label),
    RvAsOdd(Label, usize, Vec<(usize, i128)>),
    RvAsOddQuery(AsOddQuery, Label),
    EOReport(EoReport, ParsedExpr),
    EOWitness(ParsedExpr, ParsedExpr),
    RvAsWordQuery(AsWordQuery, Label, Vec<Vec<AsOperator>>),
    RvWittCover(Label, Vec<usize>, Vec<ParsedExpr>),
    RvWittQuery(WittQuery, Label),
    RvWittTower(Label, Vec<usize>, Vec<ParsedExpr>),
    RvWittTowerQuery(TowerQuery, Label),
    RvOrders(Label),
    RvFrobeniusOrders(Label, usize),
    RvOrderAt(usize, Vec<usize>, ParsedExpr),
    RvOsculating(Label, Label, usize, Vec<usize>),
    RvFitting(Label, RvMatrixSpec, usize),
    RvGroebner(Label),
    RvIntersect(Label, Label, Label),
    RvScheme(RvSchemeOp, Label, Label, ParsedExpr),
    RvProjective(Label, Label),
    RvHilbert(Label, usize),
    RankJet(RankJetSpec),
    JetPath(Label, JetPathSpec),
    JetDimension(JetPageSpec),
    VExtDimension(QSpec<i64>),
    VBars(BarFamilySpec),
    VBarDimension(QSpec<i64>),
    VBockDimension(VBockSpec),
    MayDimension(usize, usize, usize, usize),
    MayAnDimension(usize, usize, usize, usize, usize),
    CeDimension(usize, usize, usize, usize, usize),
    Homotopy(ParsedExpr, ParsedExpr),
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Token<'src> {
    Word(&'src str),
    Number(&'src str),
    LBrace,
    RBrace,
    LParen,
    RParen,
    LBracket,
    RBracket,
    Colon,
    Semi,
    Comma,
    Dot,
    Less,
    Star,
    Slash,
    Equal,
    EqualEqual,
    Minus,
    Arrow,
    FatArrow,
    Tilde,
    Eof,
}

impl Token<'_> {
    fn description(self) -> String {
        match self {
            Self::Word(word) => format!("'{word}'"),
            Self::Number(number) => format!("'{number}'"),
            Self::Eof => "end of input".into(),
            symbol => {
                let spelling = match symbol {
                    Self::LBrace => "{",
                    Self::RBrace => "}",
                    Self::LParen => "(",
                    Self::RParen => ")",
                    Self::LBracket => "[",
                    Self::RBracket => "]",
                    Self::Colon => ":",
                    Self::Semi => ";",
                    Self::Comma => ",",
                    Self::Dot => ".",
                    Self::Less => "<",
                    Self::Star => "*",
                    Self::Slash => "\\",
                    Self::Equal => "=",
                    Self::EqualEqual => "==",
                    Self::Minus => "-",
                    Self::Arrow => "->",
                    Self::FatArrow => "=>",
                    Self::Tilde => "~",
                    Self::Word(_) | Self::Number(_) | Self::Eof => unreachable!(),
                };
                format!("'{spelling}'")
            }
        }
    }
}

struct Lexer<'src> {
    source: &'src str,
    offset: usize,
    span: Span,
}

impl<'src> Lexer<'src> {
    fn new(source: &'src str) -> Self {
        Self {
            source,
            offset: 0,
            span: Span { line: 1, column: 1 },
        }
    }

    fn peek(&self) -> Option<u8> {
        self.source.as_bytes().get(self.offset).copied()
    }

    fn step(&mut self) {
        if self.peek() == Some(b'\n') {
            self.span.line += 1;
            self.span.column = 1;
        } else {
            self.span.column += 1;
        }
        self.offset += 1;
    }

    fn eat(&mut self, byte: u8) -> bool {
        if self.peek() == Some(byte) {
            self.step();
            return true;
        }
        false
    }

    fn next(&mut self) -> Result<(Token<'src>, Span)> {
        loop {
            match self.peek() {
                Some(byte) if byte.is_ascii_whitespace() => self.step(),
                Some(b'/') if self.source.as_bytes().get(self.offset + 1) == Some(&b'/') => {
                    while self.peek().is_some_and(|byte| byte != b'\n') {
                        self.step();
                    }
                }
                _ => break,
            }
        }

        let span = self.span;
        let Some(byte) = self.peek() else {
            return Ok((Token::Eof, span));
        };
        let start = self.offset;
        self.step();
        let token = match byte {
            b'a'..=b'z' | b'A'..=b'Z' | b'_' => {
                while self
                    .peek()
                    .is_some_and(|byte| byte.is_ascii_alphanumeric() || byte == b'_')
                {
                    self.step();
                }
                if self.offset - start > MAX_IDENTIFIER_BYTES {
                    return Err(Error::new(span, "identifier exceeds 128 bytes"));
                }
                Token::Word(&self.source[start..self.offset])
            }
            b'0'..=b'9' => {
                while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    self.step();
                }
                Token::Number(&self.source[start..self.offset])
            }
            b'{' => Token::LBrace,
            b'}' => Token::RBrace,
            b'(' => Token::LParen,
            b')' => Token::RParen,
            b'[' => Token::LBracket,
            b']' => Token::RBracket,
            b':' => Token::Colon,
            b';' => Token::Semi,
            b',' => Token::Comma,
            b'.' => Token::Dot,
            b'<' => Token::Less,
            b'*' => Token::Star,
            b'\\' => Token::Slash,
            b'=' if self.eat(b'=') => Token::EqualEqual,
            b'=' if self.eat(b'>') => Token::FatArrow,
            b'=' => Token::Equal,
            b'-' if self.eat(b'>') => Token::Arrow,
            b'-' => Token::Minus,
            b'~' => Token::Tilde,
            _ => {
                let character = self.source[start..].chars().next().unwrap();
                return Err(Error::new(
                    span,
                    format!("unexpected character {character:?}"),
                ));
            }
        };
        Ok((token, span))
    }
}

pub(crate) struct Parser<'src> {
    lexer: Lexer<'src>,
    lookahead: (Token<'src>, Span),
}

impl<'src> Parser<'src> {
    pub(crate) fn new(source: &'src str) -> Result<Self> {
        let mut lexer = Lexer::new(source);
        let lookahead = lexer.next()?;
        Ok(Self { lexer, lookahead })
    }

    fn token(&self) -> Token<'src> {
        self.lookahead.0
    }

    fn span(&self) -> Span {
        self.lookahead.1
    }

    fn advance(&mut self) -> Result<()> {
        self.lookahead = self.lexer.next()?;
        Ok(())
    }

    fn eat(&mut self, token: Token<'src>) -> Result<bool> {
        if self.token() == token {
            self.advance()?;
            return Ok(true);
        }
        Ok(false)
    }

    fn expect(&mut self, token: Token<'src>) -> Result<()> {
        if !self.eat(token)? {
            return Err(self.expected(&token.description()));
        }
        Ok(())
    }

    fn expected(&self, expected: &str) -> Error {
        Error::new(
            self.span(),
            format!("expected {expected}, found {}", self.token().description()),
        )
    }

    fn label(&mut self) -> Result<Label> {
        let Token::Word(text) = self.token() else {
            return Err(self.expected("an identifier"));
        };
        if matches!(
            text,
            "space"
                | "let"
                | "show"
                | "assert"
                | "betti"
                | "homotopy"
                | "Star"
                | "fst"
                | "snd"
                | "tab"
                | "res"
                | "sew"
        ) {
            return Err(self.expected("an identifier"));
        }
        let label = Label {
            text: text.into(),
            span: self.span(),
        };
        self.advance()?;
        Ok(label)
    }

    fn number(&mut self) -> Result<(usize, Span)> {
        let span = self.span();
        let Token::Number(text) = self.token() else {
            return Err(self.expected("a nonnegative decimal integer"));
        };
        let number = text
            .parse()
            .map_err(|_| Error::new(span, "decimal integer is out of range"))?;
        self.advance()?;
        Ok((number, span))
    }

    fn rho_grade(&mut self) -> Result<rho_bock::Grade> {
        let (r, _) = self.number()?;
        self.expect(Token::Comma)?;
        let chain = self.motivic_grade()?;
        self.expect(Token::Comma)?;
        let (f, _) = self.number()?;
        Ok(rho_bock::Grade { chain, r, f })
    }

    fn motivic_grade(&mut self) -> Result<MotivicGrade> {
        let (s, _) = self.number()?;
        self.expect(Token::Comma)?;
        let t = self.signed_wide("a signed topological degree", "motivic topological degree")?;
        self.expect(Token::Comma)?;
        let w = self.signed_number("motivic weight")?;
        Ok(MotivicGrade { s, t, w })
    }

    fn rv_nc_spec(&mut self) -> Result<RvNCSpec> {
        let ideal = self.label()?;
        self.expect(Token::Comma)?;
        let length = self.number()?.0;
        self.expect(Token::Comma)?;
        let b = self.signed_degree()?;
        self.expect(Token::Comma)?;
        let weight = self.exponents()?;
        Ok((ideal, length, b, weight))
    }

    fn rv_hh_spec(&mut self) -> Result<RvHHSpec> {
        let ideal = self.label()?;
        self.expect(Token::Comma)?;
        let b = self.signed_degree()?;
        self.expect(Token::Comma)?;
        let weight = self.exponents()?;
        Ok((ideal, b, weight))
    }

    fn signed_degree(&mut self) -> Result<i64> {
        self.signed_number("signed internal degree")
    }

    fn signed_number(&mut self, name: &str) -> Result<i64> {
        let span = self.span();
        let negative = self.eat(Token::Minus)?;
        let Token::Number(text) = self.token() else {
            return Err(self.expected(&format!("a {name}")));
        };
        let out_of_range = || Error::new(span, format!("{name} is out of range"));
        let magnitude = i128::from(text.parse::<u64>().map_err(|_| out_of_range())?);
        let value = i64::try_from(if negative { -magnitude } else { magnitude })
            .map_err(|_| out_of_range())?;
        self.advance()?;
        Ok(value)
    }

    fn bp_relative_shape(&mut self) -> Result<BpRelSpec> {
        let (p, n, m) = self.prime_bidegree()?;
        self.expect(Token::Comma)?;
        let (a, _) = self.number()?;
        Ok(BpRelSpec(p, n, m, a))
    }

    fn bp_mod_ext_grade(&mut self, thick: bool) -> Result<BpModSpec> {
        let (p, n, a) = self.prime_bidegree()?;
        self.expect(Token::Comma)?;
        let (h, _) = self.number()?;
        self.expect(Token::Comma)?;
        let exponent = if thick {
            let (power, _) = self.number()?;
            self.expect(Token::Comma)?;
            Some(power)
        } else {
            None
        };
        let (t, _) = self.number()?;
        Ok(BpModSpec {
            shape: BpSpec(p, n, a),
            h,
            exponent,
            t,
        })
    }

    fn bp_ext_grade(&mut self) -> Result<(BpSpec, usize)> {
        let (prime, height, arrows) = self.prime_bidegree()?;
        self.expect(Token::Comma)?;
        let (t, _) = self.number()?;
        Ok((BpSpec(prime, height, arrows), t))
    }

    fn bp_height_spec(&mut self) -> Result<BpHeightSpec> {
        let (prime, height, r) = self.prime_bidegree()?;
        self.expect(Token::Comma)?;
        let h = self.signed_wide("a signed BP height", "BP height filtration")?;
        self.expect(Token::Comma)?;
        let (a, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (t, _) = self.number()?;
        Ok(BpHeightSpec(BpSpec(prime, height, a), r, h, t))
    }

    fn bp_coefficient(&mut self) -> Result<i128> {
        self.signed_wide("a signed BP integer", "BP coefficient")
    }

    fn signed_wide(&mut self, expected: &str, label: &str) -> Result<i128> {
        let span = self.span();
        let negative = self.eat(Token::Minus)?;
        let Token::Number(text) = self.token() else {
            return Err(self.expected(expected));
        };
        let overflow = || Error::new(span, format!("{label} exceeds the signed i128 range"));
        let magnitude = text.parse::<u128>().map_err(|_| overflow())?;
        let value = if negative && magnitude == i128::MIN.unsigned_abs() {
            i128::MIN
        } else {
            let value = i128::try_from(magnitude).map_err(|_| overflow())?;
            if negative {
                -value
            } else {
                value
            }
        };
        self.advance()?;
        Ok(value)
    }

    fn exponents(&mut self) -> Result<Vec<usize>> {
        self.expect(Token::LBracket)?;
        self.unsigned_list(Token::RBracket)
    }

    fn unsigned_list(&mut self, closing: Token<'src>) -> Result<Vec<usize>> {
        let mut exponents = Vec::new();
        while self.token() != closing {
            exponents.push(self.number()?.0);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(closing)?;
        Ok(exponents)
    }

    fn q_module_body(&mut self) -> Result<crate::q_module::Presentation> {
        self.expect(Token::LBrace)?;
        let mut generators = Vec::new();
        while self.token() != Token::RBrace {
            let name = self.label()?;
            self.expect(Token::Colon)?;
            generators.push((name, self.number()?.0));
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        self.expect(Token::LBrace)?;
        let mut rules = Vec::new();
        while self.token() != Token::RBrace {
            let (primitive, span) = self.number()?;
            self.expect(Token::LParen)?;
            let source = self.label()?;
            self.expect(Token::RParen)?;
            self.expect(Token::FatArrow)?;
            self.expect(Token::LBracket)?;
            let mut targets = Vec::new();
            while self.token() != Token::RBracket {
                targets.push(self.label()?);
                if !self.eat(Token::Comma)? {
                    break;
                }
            }
            self.expect(Token::RBracket)?;
            rules.push(crate::q_module::Rule {
                primitive,
                span,
                source,
                targets,
            });
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(crate::q_module::Presentation { generators, rules })
    }

    fn sheaf_incidence(&mut self, end: Token<'src>) -> Result<Vec<(Label, i64)>> {
        let mut terms = Vec::new();
        while self.token() != end {
            let label = self.label()?;
            self.expect(Token::Colon)?;
            terms.push((label, self.signed_number("signed incidence coefficient")?));
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(end)?;
        Ok(terms)
    }

    fn sheaf_resolution(&mut self) -> Result<crate::sheaf_cellular::Presentation> {
        self.expect(Token::LBrace)?;
        let mut cells = Vec::new();
        while self.token() != Token::RBrace {
            let label = self.label()?;
            self.expect(Token::Colon)?;
            self.expect(Token::LParen)?;
            let degree = self.number()?.0;
            self.expect(Token::Comma)?;
            let point = self.label()?;
            self.expect(Token::RParen)?;
            cells.push((label, degree, point));
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        self.expect(Token::LBrace)?;
        let mut boundaries = Vec::new();
        while self.token() != Token::RBrace {
            let label = self.label()?;
            self.expect(Token::Colon)?;
            self.expect(Token::LBracket)?;
            boundaries.push((label, self.sheaf_incidence(Token::RBracket)?));
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        self.expect(Token::LBrace)?;
        let augmentation = self.sheaf_incidence(Token::RBrace)?;
        Ok(crate::sheaf_cellular::Presentation {
            cells,
            boundaries,
            augmentation,
        })
    }

    fn sheaf_local_site(&mut self) -> Result<(Label, Label, SupportMode)> {
        let complex = self.label()?;
        self.expect(Token::Comma)?;
        let point = self.label()?;
        self.expect(Token::Comma)?;
        let mode = SupportMode::from_label(&self.label()?)?;
        Ok((complex, point, mode))
    }

    fn sheaf_page_grade(&mut self) -> Result<SheafPageGrade> {
        let r = self.number()?.0;
        self.expect(Token::Comma)?;
        let p = self.number()?.0;
        self.expect(Token::Comma)?;
        let q = self.signed_wide("a signed stalk degree", "sheaf page stalk degree")?;
        Ok(SheafPageGrade { r, p, q })
    }

    fn sheaf_hom_grade(&mut self) -> Result<(Label, Label, usize)> {
        let source = self.label()?;
        self.expect(Token::Comma)?;
        let target = self.label()?;
        self.expect(Token::Comma)?;
        let (degree, _) = self.number()?;
        Ok((source, target, degree))
    }

    fn sheaf_chain(&mut self) -> Result<Vec<Label>> {
        self.expect(Token::LBracket)?;
        let mut chain = Vec::new();
        while self.token() != Token::RBracket {
            chain.push(self.label()?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBracket)?;
        Ok(chain)
    }

    fn sheaf_matrix(&mut self) -> Result<Vec<i64>> {
        self.expect(Token::LBracket)?;
        let mut entries = Vec::new();
        while self.token() != Token::RBracket {
            entries.push(self.signed_number("signed sheaf matrix entry")?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBracket)?;
        Ok(entries)
    }

    fn mixed_terms(&mut self) -> Result<crate::mixed::Terms> {
        self.expect(Token::LBracket)?;
        let mut terms = Vec::new();
        while self.token() != Token::RBracket {
            self.expect(Token::LParen)?;
            let name = self.label()?;
            self.expect(Token::Comma)?;
            let coefficient = self.signed_number("mixed complex coefficient")?;
            self.expect(Token::RParen)?;
            terms.push((name, coefficient));
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBracket)?;
        Ok(terms)
    }

    fn mixed_rules(&mut self) -> Result<crate::mixed::Rules> {
        self.expect(Token::LBrace)?;
        let mut rules = Vec::new();
        while self.token() != Token::RBrace {
            let generator = self.label()?;
            self.expect(Token::FatArrow)?;
            let terms = self.mixed_terms()?;
            rules.push((generator, terms));
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(rules)
    }

    fn mixed_higher_rules(&mut self) -> Result<crate::mixed::HigherRules> {
        self.expect(Token::LBrace)?;
        let mut rules = Vec::new();
        while self.token() != Token::RBrace {
            let generator = self.label()?;
            self.expect(Token::FatArrow)?;
            self.expect(Token::LBracket)?;
            let mut terms = Vec::new();
            while self.token() != Token::RBracket {
                self.expect(Token::LParen)?;
                let power = self.number()?.0;
                self.expect(Token::Comma)?;
                let target = self.label()?;
                self.expect(Token::Comma)?;
                let coefficient = self.signed_number("higher mixed coefficient")?;
                self.expect(Token::RParen)?;
                terms.push((power, target, coefficient));
                if !self.eat(Token::Comma)? {
                    break;
                }
            }
            self.expect(Token::RBracket)?;
            rules.push((generator, terms));
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(rules)
    }

    fn mixed_spec(&mut self) -> Result<MixedSpec> {
        let complex = self.label()?;
        self.expect(Token::Comma)?;
        let r = self.number()?.0;
        self.expect(Token::Comma)?;
        let degree = self.signed_degree()?;
        self.expect(Token::Comma)?;
        let filtration = self.number()?.0;
        Ok((complex, r, degree, filtration))
    }

    fn sheaf_body(&mut self) -> Result<sheaf::Presentation> {
        self.expect(Token::LBrace)?;
        let mut stalks = Vec::new();
        while self.token() != Token::RBrace {
            let point = self.label()?;
            self.expect(Token::Colon)?;
            stalks.push((point, self.number()?.0));
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        self.expect(Token::LBrace)?;
        let mut maps = Vec::new();
        while self.token() != Token::RBrace {
            let source = self.label()?;
            self.expect(Token::Less)?;
            let target = self.label()?;
            self.expect(Token::Colon)?;
            let entries = self.sheaf_matrix()?;
            maps.push(sheaf::Cover {
                source,
                target,
                entries,
            });
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(sheaf::Presentation { stalks, maps })
    }

    fn chains(&mut self) -> Result<Vec<Vec<Label>>> {
        self.expect(Token::LBrace)?;
        let mut chains = Vec::new();
        while self.token() != Token::RBrace {
            let mut chain = vec![self.label()?];
            while self.eat(Token::Less)? {
                chain.push(self.label()?);
            }
            chains.push(chain);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(chains)
    }

    pub(crate) fn statement(&mut self) -> Result<Option<Statement>> {
        let span = self.span();
        let kind = match self.token() {
            Token::Eof => return Ok(None),
            Token::Word("qmodule") => {
                self.advance()?;
                let name = self.label()?;
                let selected = self.exponents()?;
                let presentation = self.q_module_body()?;
                StatementForm::QModule(name, selected, presentation)
            }
            Token::Word("rvHHMixed") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let ideal = self.label()?;
                self.expect(Token::Comma)?;
                let weight = self.exponents()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvHHMixed(name, ideal, weight)
            }
            Token::Word("mixcone") => {
                self.advance()?;
                let name = self.label()?;
                let map = self.paren_expr(0)?;
                StatementForm::MixedCone(name, map)
            }
            Token::Word("mixtensor") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let left = self.label()?;
                self.expect(Token::Comma)?;
                let right = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::MixedTensor(name, left, right)
            }
            Token::Word("mixhom") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let source = self.label()?;
                self.expect(Token::Comma)?;
                let target = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::MixedHom(name, source, target)
            }
            Token::Word("mixminimal") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let source = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::MixedMinimal(name, source)
            }
            Token::Word("mixmindim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::Comma)?;
                let b = self.signed_degree()?;
                self.expect(Token::RBracket)?;
                StatementForm::MixedMinDimension(name, b)
            }
            Token::Word(command @ ("mixedcomplex" | "mixedhigher")) => {
                self.advance()?;
                let complex = self.label()?;
                self.expect(Token::LBracket)?;
                let prime = self.number()?.0;
                self.expect(Token::RBracket)?;
                self.expect(Token::LBrace)?;
                let mut generators = Vec::new();
                while self.token() != Token::RBrace {
                    let label = self.label()?;
                    self.expect(Token::Colon)?;
                    generators.push((label, self.signed_degree()?));
                    if !self.eat(Token::Comma)? {
                        break;
                    }
                }
                self.expect(Token::RBrace)?;
                if command == "mixedhigher" {
                    let rules = self.mixed_higher_rules()?;
                    StatementForm::MixedHigher(complex, prime, generators, rules)
                } else {
                    let b = self.mixed_rules()?;
                    let big_b = self.mixed_rules()?;
                    StatementForm::MixedComplex(complex, prime, generators, b, big_b)
                }
            }
            Token::Word("mixbars") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::MixedBars(name)
            }
            Token::Word(command @ ("mixbardim" | "mixbarrank" | "mixbarperiodicdim")) => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::Comma)?;
                let source = self.signed_degree()?;
                let target = if command == "mixbarrank" {
                    self.expect(Token::Comma)?;
                    Some(self.signed_degree()?)
                } else {
                    None
                };
                self.expect(Token::RBracket)?;
                if command == "mixbarperiodicdim" {
                    StatementForm::MixedPeriodicRank(name, source)
                } else {
                    StatementForm::MixedBarRank(name, source, target)
                }
            }
            Token::Word("mixdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.mixed_spec()?;
                self.expect(Token::RBracket)?;
                StatementForm::MixedDimension(spec)
            }
            Token::Word("sheaf") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let base = self.label()?;
                self.expect(Token::Comma)?;
                let prime = self.number()?.0;
                self.expect(Token::RBracket)?;
                StatementForm::Sheaf(name, base, prime, self.sheaf_body()?)
            }
            Token::Word("shtruncate" | "shdeligne" | "shdelignecomplex") => {
                let source = match self.token() {
                    Token::Word("shdeligne") => Some(SheafInputForm::Sheaf),
                    Token::Word("shdelignecomplex") => Some(SheafInputForm::Complex),
                    _ => None,
                };
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let input = self.label()?;
                self.expect(Token::Comma)?;
                let statement = if let Some(form) = source {
                    let sheaf = self.label()?;
                    self.expect(Token::Comma)?;
                    StatementForm::SheafDeligne(form, name, input, sheaf, self.number()?.0)
                } else {
                    StatementForm::SheafTruncate(name, input, self.number()?.0)
                };
                self.expect(Token::RBracket)?;
                statement
            }
            Token::Word("shextra" | "shextracomplex") => {
                let form = if self.token() == Token::Word("shextra") {
                    SheafInputForm::Sheaf
                } else {
                    SheafInputForm::Complex
                };
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let map = self.label()?;
                self.expect(Token::Comma)?;
                let source = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SheafExtra(form, name, map, source)
            }
            Token::Word("shlocaldim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (complex, point, mode) = self.sheaf_local_site()?;
                self.expect(Token::Comma)?;
                let degree = self.number()?.0;
                self.expect(Token::RBracket)?;
                StatementForm::SheafLocalDimension(complex, point, mode, degree)
            }
            Token::Word("shleray" | "shleraycomplex") => {
                let form = if self.token() == Token::Word("shleray") {
                    SheafInputForm::Sheaf
                } else {
                    SheafInputForm::Complex
                };
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let map = self.label()?;
                self.expect(Token::Comma)?;
                let source = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SheafLeray(form, name, map, source)
            }
            Token::Word("shserre") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let sheaf = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SheafSerre(name, sheaf)
            }
            Token::Word("shserreshift") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SheafSerreShift(name)
            }
            Token::Word("shcohomology") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let complex = self.label()?;
                self.expect(Token::Comma)?;
                let degree = self.number()?.0;
                self.expect(Token::RBracket)?;
                StatementForm::SheafStalkCohomology(name, complex, degree)
            }
            Token::Word("shpull" | "shpush") => {
                let pull = self.token() == Token::Word("shpull");
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let function = self.label()?;
                self.expect(Token::Comma)?;
                let sheaf = self.label()?;
                let operation = if pull {
                    SheafImageOp::Pull
                } else {
                    self.expect(Token::Comma)?;
                    SheafImageOp::Push(self.number()?.0)
                };
                self.expect(Token::RBracket)?;
                StatementForm::SheafImage(operation, name, function, sheaf)
            }
            Token::Word("shmap") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let source = self.label()?;
                self.expect(Token::Comma)?;
                let target = self.label()?;
                self.expect(Token::RBracket)?;
                self.expect(Token::LBrace)?;
                let mut matrices = Vec::new();
                while self.token() != Token::RBrace {
                    let point = self.label()?;
                    self.expect(Token::Colon)?;
                    matrices.push((point, self.sheaf_matrix()?));
                    if !self.eat(Token::Comma)? {
                        break;
                    }
                }
                self.expect(Token::RBrace)?;
                StatementForm::SheafMap(name, source, target, matrices)
            }
            Token::Word("shresolution") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let base = self.label()?;
                self.expect(Token::Comma)?;
                let prime = self.number()?.0;
                self.expect(Token::RBracket)?;
                StatementForm::SheafResolution(name, base, prime, self.sheaf_resolution()?)
            }
            Token::Word("shcellular") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let resolution = self.label()?;
                self.expect(Token::Comma)?;
                let complex = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SheafCellular(name, resolution, complex)
            }
            Token::Word("shcomplex") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBrace)?;
                let mut entries = Vec::new();
                while self.token() != Token::RBrace {
                    entries.push(self.label()?);
                    if !self.eat(Token::Comma)? {
                        break;
                    }
                }
                self.expect(Token::RBrace)?;
                StatementForm::SheafComplex(name, entries)
            }
            Token::Word("shedim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::Comma)?;
                let grade = self.sheaf_page_grade()?;
                self.expect(Token::RBracket)?;
                StatementForm::SheafPageDimension(name, grade)
            }
            Token::Word("shhyperdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::Comma)?;
                let degree = self.number()?.0;
                self.expect(Token::RBracket)?;
                StatementForm::SheafHyperDimension(name, degree)
            }
            Token::Word("shExtDim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (source, target, degree) = self.sheaf_hom_grade()?;
                self.expect(Token::RBracket)?;
                StatementForm::SheafExtDimension(source, target, degree)
            }
            Token::Word("shdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::Comma)?;
                let degree = self.number()?.0;
                self.expect(Token::RBracket)?;
                StatementForm::SheafDimension(name, degree)
            }
            Token::Word("space") => {
                self.advance()?;
                let name = self.label()?;
                StatementForm::Space(name, self.chains()?)
            }
            Token::Word("let") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::Equal)?;
                StatementForm::Let(name, self.expr(0)?)
            }
            Token::Word(
                name @ ("rvTreeShell" | "rvTreeBall" | "rvTreeCount" | "rvTreeClosureCount"),
            ) => {
                use crate::rv_hodge::Boundary;
                use crate::rv_tree::Report;
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (modulus, _) = self.number()?;
                self.expect(Token::Comma)?;
                let (root, _) = self.number()?;
                self.expect(Token::Comma)?;
                let (radius, _) = self.number()?;
                let report = match name {
                    "rvTreeShell" => Report::Shell,
                    "rvTreeBall" => Report::Ball,
                    _ => {
                        self.expect(Token::Comma)?;
                        let mu = self.coweight()?;
                        let boundary = if name == "rvTreeCount" {
                            Boundary::Exact
                        } else {
                            Boundary::Closure
                        };
                        Report::Count(boundary, mu)
                    }
                };
                self.expect(Token::RBracket)?;
                StatementForm::RVTree(report, modulus, root, radius, self.paren_expr(0)?)
            }
            Token::Word("rvCycloCoeff") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (power, _) = self.number()?;
                self.expect(Token::RBracket)?;
                StatementForm::CycloCoefficient(power, self.paren_expr(0)?)
            }
            Token::Word("rvConformalDim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (rank, _) = self.number()?;
                self.expect(Token::Comma)?;
                let (level, _) = self.number()?;
                self.expect(Token::Comma)?;
                let (genus, _) = self.number()?;
                self.expect(Token::Comma)?;
                self.expect(Token::LBracket)?;
                let mut punctures = Vec::new();
                while self.token() != Token::RBracket {
                    punctures.push(self.signed_tuple("a signed fusion weight", "fusion weight")?);
                    if !self.eat(Token::Comma)? {
                        break;
                    }
                }
                self.expect(Token::RBracket)?;
                self.expect(Token::RBracket)?;
                StatementForm::ConformalDimension(rank, level, genus, punctures)
            }
            Token::Word("rvFusionCoeff") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let weight = self.signed_tuple("a signed fusion weight", "fusion weight")?;
                self.expect(Token::RBracket)?;
                StatementForm::FusionCoefficient(weight, self.paren_expr(0)?)
            }
            Token::Word("rvOreCoeff" | "rvQGrPeriod") => {
                let period = self.token() == Token::Word("rvQGrPeriod");
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (a, b) = self.quantum_grade()?;
                if period {
                    self.expect(Token::Comma)?;
                    let (degree, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    StatementForm::QuantumPeriod(a, b, degree)
                } else {
                    self.expect(Token::RBracket)?;
                    StatementForm::OreCoefficient(a, b, self.paren_expr(0)?)
                }
            }
            Token::Word(
                "rvQGrCurvatureChar" | "rvQGrCurvatureCoeff" | "rvQGrCurvatureHorizontal",
            ) => {
                use crate::rv_pcurv::Report;
                let kind = self.token();
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (prime, rank, ambient) = self.quantum_mod_grade()?;
                let report = match kind {
                    Token::Word("rvQGrCurvatureCoeff") => {
                        self.expect(Token::Comma)?;
                        Report::Coefficient(self.number()?.0)
                    }
                    Token::Word("rvQGrCurvatureHorizontal") => Report::Horizontal,
                    _ => Report::Characteristic,
                };
                self.expect(Token::RBracket)?;
                StatementForm::QuantumCurvature(report, prime, rank, ambient)
            }
            Token::Word("rvQGrCoeff") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (degree, _) = self.number()?;
                self.expect(Token::Comma)?;
                let shape = self.exponents()?;
                self.expect(Token::RBracket)?;
                StatementForm::QuantumCoefficient(degree, shape, self.paren_expr(0)?)
            }
            Token::Word("rvQGrGW") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (rank, ambient) = self.quantum_grade()?;
                self.expect(Token::Comma)?;
                let (degree, _) = self.number()?;
                self.expect(Token::Comma)?;
                let a = self.exponents()?;
                self.expect(Token::Comma)?;
                let b = self.exponents()?;
                self.expect(Token::Comma)?;
                let c = self.exponents()?;
                self.expect(Token::RBracket)?;
                StatementForm::GromovWitten(rank, ambient, degree, [a, b, c])
            }
            Token::Word("rvHiveCount" | "rvHorn") => {
                let horn = self.token() == Token::Word("rvHorn");
                self.advance()?;
                self.expect(Token::LBracket)?;
                let weights = self.hive_spec()?;
                self.expect(Token::RBracket)?;
                StatementForm::HiveReport(weights, horn)
            }
            Token::Word("rvHiveAt") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (i, _) = self.number()?;
                self.expect(Token::Comma)?;
                let (j, _) = self.number()?;
                self.expect(Token::RBracket)?;
                StatementForm::HiveVertex(i, j, self.paren_expr(0)?)
            }
            Token::Word("rvSatakeCoeff" | "rvSatakeAt") => {
                let numeric = self.token() == Token::Word("rvSatakeAt");
                self.advance()?;
                self.expect(Token::LBracket)?;
                let v = if numeric {
                    let value =
                        self.signed_wide("a signed v specialization", "v specialization")?;
                    self.expect(Token::Comma)?;
                    Some(value)
                } else {
                    None
                };
                let coordinates =
                    self.signed_tuple("a signed Satake coordinate", "Satake coordinate")?;
                self.expect(Token::RBracket)?;
                StatementForm::SatakeQuery(v, coordinates, self.paren_expr(0)?)
            }
            Token::Word("rvHeckeCoeff" | "rvHeckeAt" | "rvSphCoeff" | "rvSphAt") => {
                let evaluate = matches!(self.token(), Token::Word("rvHeckeAt" | "rvSphAt"));
                let spherical = matches!(self.token(), Token::Word("rvSphCoeff" | "rvSphAt"));
                self.advance()?;
                self.expect(Token::LBracket)?;
                let parameter = if evaluate {
                    let parameter =
                        self.signed_wide("a signed Hecke parameter", "Hecke parameter")?;
                    self.expect(Token::Comma)?;
                    Some(parameter)
                } else {
                    None
                };
                let window = self.affine_window()?;
                self.expect(Token::RBracket)?;
                StatementForm::HeckeCoefficient(spherical, parameter, window, self.paren_expr(0)?)
            }
            Token::Word("rvIwKL") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let lower = self.affine_window()?;
                self.expect(Token::Comma)?;
                let upper = self.affine_window()?;
                self.expect(Token::RBracket)?;
                StatementForm::RVIwKL(lower, upper)
            }
            Token::Word("rvIwSchubertKL") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let upper = self.affine_window()?;
                self.expect(Token::RBracket)?;
                StatementForm::RVIwSchubertKL(upper, self.paren_expr(0)?)
            }
            Token::Word("rvIw" | "rvIwLength") => {
                let length = self.token() == Token::Word("rvIwLength");
                self.advance()?;
                StatementForm::RVIw(length, self.paren_expr(0)?)
            }
            Token::Word("rvIwPosition") => {
                self.advance()?;
                StatementForm::RVIwPosition(Box::new(self.paren_pair(0)?))
            }
            Token::Word("rvIwAt") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (modulus, root) = self.embedding_parameters()?;
                StatementForm::RVIwAt(modulus, root, Box::new(self.paren_pair(0)?))
            }
            Token::Word("rvHodge") => {
                self.advance()?;
                StatementForm::RVHodge(self.paren_expr(0)?)
            }
            Token::Word("rvPosition") => {
                self.advance()?;
                StatementForm::RVPosition(Box::new(self.paren_pair(0)?))
            }
            Token::Word("rvMuAt") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (modulus, root) = self.embedding_parameters()?;
                StatementForm::RVMuAt(modulus, root, Box::new(self.paren_pair(0)?))
            }
            Token::Word("rvNewton") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (modulus, root) = self.embedding_parameters()?;
                StatementForm::RVNewton(modulus, root, self.paren_expr(0)?)
            }
            Token::Word("rvInvariantFactors") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (modulus, root) = self.embedding_parameters()?;
                StatementForm::RVInvariantFactors(modulus, root, self.paren_expr(0)?)
            }
            Token::Word("show") => {
                self.advance()?;
                StatementForm::Show(self.expr(0)?)
            }
            Token::Word("card") => {
                self.advance()?;
                StatementForm::Card(self.expr(0)?)
            }
            Token::Word("assert") => {
                self.advance()?;
                let left = self.expr(0)?;
                self.expect(Token::EqualEqual)?;
                StatementForm::Assert(left, self.expr(0)?)
            }
            Token::Word("betti") => {
                self.advance()?;
                StatementForm::Betti(self.ty(0)?)
            }
            Token::Word("bprelext") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let shape = self.bp_relative_shape()?;
                self.expect(Token::Comma)?;
                let (t, _) = self.number()?;
                self.expect(Token::RBracket)?;
                StatementForm::BpRelExt(shape, t)
            }
            Token::Word("bpext") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (spec, t) = self.bp_ext_grade()?;
                self.expect(Token::RBracket)?;
                StatementForm::BpExt(spec, t)
            }
            Token::Word("bpTripleDefined") => {
                self.advance()?;
                self.expect(Token::LParen)?;
                StatementForm::BpTripleDefined(self.bp_triple_arguments(
                    0,
                    span,
                    "bpTripleDefined",
                )?)
            }
            Token::Word("bpTripleContains") => {
                self.advance()?;
                self.expect(Token::LParen)?;
                let pair: [ParsedExpr; 2] = self.arguments(0)?.try_into().map_err(|_| {
                    Error::new(span, "bpTripleContains requires a bracket and a class")
                })?;
                StatementForm::BpTripleContains(Box::new(pair))
            }
            Token::Word("bpTripleSize") => {
                self.advance()?;
                StatementForm::BpTripleSize(self.paren_expr(0)?)
            }
            Token::Word("bporder") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (t, _) = self.number()?;
                self.expect(Token::RBracket)?;
                StatementForm::BpOrder(t, self.paren_expr(0)?)
            }
            Token::Word("pextdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (p, s, t) = self.prime_bidegree()?;
                self.expect(Token::RBracket)?;
                StatementForm::PExtDimension(p, s, t)
            }
            Token::Word("bpHeightPage" | "bpHeightImage") => {
                let image = self.token() == Token::Word("bpHeightImage");
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.prime_profile_page_degree()?;
                self.expect(Token::RBracket)?;
                StatementForm::BpHeightPage(image, spec)
            }
            Token::Word("bpIExtdim") | Token::Word("bpThickExtdim") => {
                let thick = self.token() == Token::Word("bpThickExtdim");
                self.advance()?;
                self.expect(Token::LBracket)?;
                let grade = self.bp_mod_ext_grade(thick)?;
                self.expect(Token::RBracket)?;
                StatementForm::BpModExtDimension(grade)
            }
            Token::Word("bpVdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.bp_vpage_spec()?;
                self.expect(Token::RBracket)?;
                StatementForm::BpVPageDimension(spec)
            }
            Token::Word("bpBockdim") | Token::Word("bpNovdim") => {
                let mode = if self.token() == Token::Word("bpBockdim") {
                    BpPageMode::Bockstein
                } else {
                    BpPageMode::Novikov
                };
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.prime_profile_page_degree()?;
                self.expect(Token::RBracket)?;
                StatementForm::BpPageDimension(mode, spec)
            }
            Token::Word("pextandim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (p, n, s, t) = self.prime_profile_bidegree()?;
                self.expect(Token::RBracket)?;
                StatementForm::PExtAnDimension(p, n, s, t)
            }
            Token::Word("rhobockdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let grade = self.rho_grade()?;
                self.expect(Token::RBracket)?;
                StatementForm::RhoDimension(grade)
            }
            Token::Word("rmbars" | "rmbarRank") => {
                let is_rank = self.token() == Token::Word("rmbarRank");
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (s, _) = self.number()?;
                self.expect(Token::Comma)?;
                let k = self.signed_wide("a signed degree t-w", "rho bar degree")?;
                let form = if is_rank {
                    self.expect(Token::Comma)?;
                    let t = self
                        .signed_wide("a signed topological degree", "rho bar topological degree")?;
                    self.expect(Token::Comma)?;
                    let (power, _) = self.number()?;
                    StatementForm::RealBarRank(s, k, t, power)
                } else {
                    StatementForm::RealBars(s, k)
                };
                self.expect(Token::RBracket)?;
                form
            }
            Token::Word("rhoorder") => {
                self.advance()?;
                let input = self.paren_expr(0)?;
                StatementForm::RhoOrder(Box::new(input))
            }
            Token::Word("rmTripleDim") => {
                self.advance()?;
                StatementForm::RealTripleDim(Box::new(self.paren_expr(0)?))
            }
            Token::Word("rmTripleContains") => {
                let span = self.span();
                self.advance()?;
                self.expect(Token::LParen)?;
                StatementForm::RealTripleContains(self.exact_arguments(
                    0,
                    span,
                    "rmTripleContains",
                    "a bracket and a candidate",
                )?)
            }
            Token::Word("rmextdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let grade = self.motivic_grade()?;
                self.expect(Token::RBracket)?;
                StatementForm::RMExtDimension(grade)
            }
            Token::Word("mextdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let grade = self.motivic_grade()?;
                self.expect(Token::RBracket)?;
                StatementForm::MExtDimension(grade)
            }
            Token::Word("extdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (s, t) = self.bidegree()?;
                self.expect(Token::RBracket)?;
                StatementForm::ExtDimension(s, t)
            }
            Token::Word("extandim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (n, s, t) = self.profile_bidegree()?;
                self.expect(Token::RBracket)?;
                StatementForm::ExtAnDimension(n, s, t)
            }
            Token::Word("xextdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (base, s, t) = self.space_bidegree(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::XExtDimension(base, s, t)
            }
            Token::Word("mapextdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (source, target, s, t) = self.map_bidegree(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::MapExtDimension(source, target, s, t)
            }
            Token::Word("rvideal") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let selected = self.exponents()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvIdeal(name, selected, self.rv_generators()?)
            }
            Token::Word("rvwittcover" | "rvwitttower") => {
                let tower = self.token() == Token::Word("rvwitttower");
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let variables = self.exponents()?;
                self.expect(Token::RBracket)?;
                let coordinates = self.rv_generators()?;
                if tower {
                    StatementForm::RvWittTower(name, variables, coordinates)
                } else {
                    StatementForm::RvWittCover(name, variables, coordinates)
                }
            }
            Token::Word(
                "rvWittTowerBreaks"
                | "rvWittTowerGenus"
                | "rvWittTowerPoints"
                | "rvWittTowerCharacterL"
                | "rvWittTowerCharacterNewton"
                | "rvWittTowerCharacterBaseL"
                | "rvWittTowerFrobPoints"
                | "rvWittTowerLPolynomial"
                | "rvWittTowerNewton",
            ) => {
                let operation = self.token();
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                let query = match operation {
                    Token::Word("rvWittTowerBreaks") => TowerQuery::Breaks,
                    Token::Word("rvWittTowerGenus") => TowerQuery::Genus,
                    Token::Word("rvWittTowerPoints") => {
                        self.expect(Token::Comma)?;
                        TowerQuery::Points(self.number()?.0)
                    }
                    Token::Word("rvWittTowerCharacterL") => {
                        self.expect(Token::Comma)?;
                        TowerQuery::CharacterL(self.number()?.0)
                    }
                    Token::Word("rvWittTowerCharacterNewton") => {
                        self.expect(Token::Comma)?;
                        TowerQuery::CharacterNewton(self.number()?.0)
                    }
                    Token::Word("rvWittTowerCharacterBaseL") => {
                        self.expect(Token::Comma)?;
                        let height = self.number()?.0;
                        self.expect(Token::Comma)?;
                        TowerQuery::CharacterBaseL(height, self.number()?.0)
                    }
                    Token::Word("rvWittTowerFrobPoints") => {
                        self.expect(Token::Comma)?;
                        TowerQuery::FrobPoints(self.number()?.0)
                    }
                    Token::Word("rvWittTowerLPolynomial") => TowerQuery::LPolynomial,
                    Token::Word("rvWittTowerNewton") => TowerQuery::Newton,
                    _ => unreachable!(),
                };
                self.expect(Token::RBracket)?;
                StatementForm::RvWittTowerQuery(query, name)
            }
            Token::Word(
                "rvWittBreaks" | "rvWittGenus" | "rvWittPoints" | "rvWittLPolynomial"
                | "rvWittNewton" | "rvWittCharacterL" | "rvWittGaps" | "rvWittBasis"
                | "rvWittCartierRank" | "rvWittaNumber",
            ) => {
                let operation = self.token();
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                let query = match operation {
                    Token::Word("rvWittBreaks") => WittQuery::Breaks,
                    Token::Word("rvWittGenus") => WittQuery::Genus,
                    Token::Word("rvWittPoints") => {
                        self.expect(Token::Comma)?;
                        WittQuery::Points(self.number()?.0)
                    }
                    Token::Word("rvWittLPolynomial") => WittQuery::LPolynomial,
                    Token::Word("rvWittNewton") => WittQuery::Newton,
                    Token::Word("rvWittCharacterL") => WittQuery::CharacterL,
                    Token::Word("rvWittGaps") => WittQuery::Gaps,
                    Token::Word("rvWittBasis") => WittQuery::Basis,
                    Token::Word("rvWittCartierRank") => {
                        self.expect(Token::Comma)?;
                        WittQuery::CartierRank(self.number()?.0)
                    }
                    Token::Word("rvWittaNumber") => WittQuery::ANumber,
                    _ => unreachable!(),
                };
                self.expect(Token::RBracket)?;
                StatementForm::RvWittQuery(query, name)
            }
            Token::Word("rvasodd") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let prime = self.number()?.0;
                self.expect(Token::RBracket)?;
                self.expect(Token::LBrace)?;
                let mut terms = Vec::new();
                while self.token() != Token::RBrace {
                    self.expect(Token::LParen)?;
                    let exponent = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let coefficient =
                        self.signed_wide("a signed coefficient", "odd Artin–Schreier coefficient")?;
                    self.expect(Token::RParen)?;
                    terms.push((exponent, coefficient));
                    if !self.eat(Token::Comma)? {
                        break;
                    }
                }
                self.expect(Token::RBrace)?;
                StatementForm::RvAsOdd(name, prime, terms)
            }
            Token::Word(
                "rvASOddReduced"
                | "rvASOddBasis"
                | "rvASOddCartierRank"
                | "rvASOddaNumber"
                | "rvASOddEO"
                | "rvASOddGram"
                | "rvASOddTatedim"
                | "rvASOddTateJordan"
                | "rvASOddGenus"
                | "rvASOddPrank"
                | "rvASOddCharacterL"
                | "rvASOddCharacterNewton"
                | "rvASOddLPolynomial"
                | "rvASOddNewton"
                | "rvASOddPoints",
            ) => {
                let operation = self.token();
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                let query = match operation {
                    Token::Word("rvASOddReduced") => AsOddQuery::Reduced,
                    Token::Word("rvASOddBasis") => AsOddQuery::Basis,
                    Token::Word("rvASOddaNumber") => AsOddQuery::ANumber,
                    Token::Word("rvASOddEO") => AsOddQuery::FinalType,
                    Token::Word("rvASOddGram") => AsOddQuery::Gram,
                    Token::Word("rvASOddTateJordan") => AsOddQuery::TateJordan,
                    Token::Word("rvASOddTatedim") => {
                        self.expect(Token::Comma)?;
                        AsOddQuery::TateDim(self.signed_number("Tate degree")?)
                    }
                    Token::Word("rvASOddCartierRank") => {
                        self.expect(Token::Comma)?;
                        AsOddQuery::CartierRank(self.number()?.0)
                    }
                    Token::Word("rvASOddGenus") => AsOddQuery::Genus,
                    Token::Word("rvASOddPrank") => AsOddQuery::PRank,
                    Token::Word("rvASOddLPolynomial") => AsOddQuery::LPolynomial,
                    Token::Word("rvASOddNewton") => AsOddQuery::Newton,
                    Token::Word(
                        "rvASOddCharacterL" | "rvASOddCharacterNewton" | "rvASOddPoints",
                    ) => {
                        self.expect(Token::Comma)?;
                        let index = self.number()?.0;
                        match operation {
                            Token::Word("rvASOddCharacterL") => AsOddQuery::Character(index),
                            Token::Word("rvASOddCharacterNewton") => {
                                AsOddQuery::CharacterNewton(index)
                            }
                            _ => AsOddQuery::Points(index),
                        }
                    }
                    _ => unreachable!(),
                };
                self.expect(Token::RBracket)?;
                StatementForm::RvAsOddQuery(query, name)
            }
            Token::Word("rvascover") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let variables = self.exponents()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvAsCover(name, variables, self.rv_generators()?)
            }
            Token::Word("eoDim" | "eoWeyl" | "eoFirstSlope" | "eoSlopeCycles") => {
                let report = match self.token() {
                    Token::Word("eoDim") => EoReport::Dimension,
                    Token::Word("eoWeyl") => EoReport::Weyl,
                    Token::Word("eoFirstSlope") => EoReport::FirstSlope,
                    Token::Word("eoSlopeCycles") => EoReport::SlopeCycles,
                    _ => unreachable!(),
                };
                self.advance()?;
                StatementForm::EOReport(report, self.paren_expr(0)?)
            }
            Token::Word("eoClosureWitness") => {
                self.advance()?;
                self.expect(Token::LParen)?;
                let lower = self.expr(0)?;
                self.expect(Token::Comma)?;
                let upper = self.expr(0)?;
                self.expect(Token::RParen)?;
                StatementForm::EOWitness(lower, upper)
            }
            Token::Word("rvASWordRank" | "rvASCommonKer" | "rvASSumImages") => {
                let query = match self.token() {
                    Token::Word("rvASWordRank") => AsWordQuery::Rank,
                    Token::Word("rvASCommonKer") => AsWordQuery::CommonKernel,
                    Token::Word("rvASSumImages") => AsWordQuery::SumImages,
                    _ => unreachable!(),
                };
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::Comma)?;
                let words = if matches!(query, AsWordQuery::Rank) {
                    vec![self.as_operator_word()?]
                } else {
                    self.as_operator_words()?
                };
                self.expect(Token::RBracket)?;
                StatementForm::RvAsWordQuery(query, name, words)
            }
            Token::Word(
                "rvASGenus" | "rvASPrank" | "rvASPoints" | "rvASLPolynomial" | "rvASNewton"
                | "rvASCartierRank" | "rvASaNumber" | "rvASEO" | "rvASPoleOrders",
            ) => {
                let operation = self.token();
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                let query = match operation {
                    Token::Word("rvASGenus") => AsQuery::Genus,
                    Token::Word("rvASPrank") => AsQuery::PRank,
                    Token::Word("rvASPoints") => {
                        self.expect(Token::Comma)?;
                        AsQuery::Points(self.number()?.0)
                    }
                    Token::Word("rvASLPolynomial") => AsQuery::LPolynomial,
                    Token::Word("rvASNewton") => AsQuery::Newton,
                    Token::Word("rvASCartierRank") => {
                        self.expect(Token::Comma)?;
                        AsQuery::CartierRank(self.number()?.0)
                    }
                    Token::Word("rvASaNumber") => AsQuery::ANumber,
                    Token::Word("rvASEO") => AsQuery::FinalType,
                    Token::Word("rvASPoleOrders") => AsQuery::PoleOrders,
                    _ => unreachable!(),
                };
                self.expect(Token::RBracket)?;
                StatementForm::RvAsQuery(query, name)
            }
            Token::Word("rvcurve") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let variables = self.exponents()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvCurve(name, variables, self.rv_generators()?)
            }
            Token::Word("rvFrobOrders") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::Comma)?;
                let (exponent, _) = self.number()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvFrobeniusOrders(name, exponent)
            }
            Token::Word("rvOrdAt") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (modulus, _) = self.number()?;
                self.expect(Token::Comma)?;
                let point = self.exponents()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvOrderAt(modulus, point, self.paren_expr(0)?)
            }
            Token::Word("rvorders") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvOrders(name)
            }
            Token::Word("rvosculating") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let curve = self.label()?;
                self.expect(Token::Comma)?;
                let (k, _) = self.number()?;
                self.expect(Token::Comma)?;
                let variables = self.exponents()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvOsculating(name, curve, k, variables)
            }
            Token::Word("rvgauss" | "rvdual") => {
                let operation = if self.token() == Token::Word("rvgauss") {
                    GaussOp::Graph
                } else {
                    GaussOp::Dual
                };
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let variables = self.exponents()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvGauss(operation, name, variables, self.paren_expr(0)?)
            }
            Token::Word("rvrees" | "rvnormalcone") => {
                let operation = if self.token() == Token::Word("rvrees") {
                    ReesOp::Rees
                } else {
                    ReesOp::NormalCone
                };
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let base = self.label()?;
                self.expect(Token::Comma)?;
                let variables = self.exponents()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvRees(operation, name, base, variables, self.rv_generators()?)
            }
            Token::Word("rvfitting") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let spec = self.rv_matrix(0)?;
                let (size, _) = self.number()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvFitting(name, spec, size)
            }
            Token::Word("rvintersect" | "rvcolon" | "rvsaturate" | "rvprojective") => {
                let operation = self.token();
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let left = self.label()?;
                if operation == Token::Word("rvintersect") {
                    self.expect(Token::Comma)?;
                    let right = self.label()?;
                    self.expect(Token::RBracket)?;
                    StatementForm::RvIntersect(name, left, right)
                } else {
                    self.expect(Token::RBracket)?;
                    match operation {
                        Token::Word("rvprojective") => StatementForm::RvProjective(name, left),
                        Token::Word("rvcolon") => StatementForm::RvScheme(
                            RvSchemeOp::Colon,
                            name,
                            left,
                            self.paren_expr(0)?,
                        ),
                        Token::Word("rvsaturate") => StatementForm::RvScheme(
                            RvSchemeOp::Saturate,
                            name,
                            left,
                            self.paren_expr(0)?,
                        ),
                        _ => unreachable!(),
                    }
                }
            }
            Token::Word("rvgroebner" | "rvhilbert") => {
                let hilbert = self.token() == Token::Word("rvhilbert");
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                let form = if hilbert {
                    self.expect(Token::Comma)?;
                    StatementForm::RvHilbert(name, self.number()?.0)
                } else {
                    StatementForm::RvGroebner(name)
                };
                self.expect(Token::RBracket)?;
                form
            }
            Token::Word("rvbarpath") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let ideal = self.label()?;
                self.expect(Token::RBracket)?;
                self.expect(Token::LBrace)?;
                let mut steps = Vec::new();
                while self.token() != Token::RBrace {
                    steps.push(self.exponents()?);
                    if !self.eat(Token::Comma)? {
                        break;
                    }
                }
                self.expect(Token::RBrace)?;
                StatementForm::RvBarPath(name, ideal, steps)
            }
            Token::Word(name @ ("rvHNdim" | "rvHPdim" | "rvHCdim")) => {
                use crate::rv_cyclic_homology::Mode;
                let mode = match name {
                    "rvHNdim" => Mode::Negative,
                    "rvHPdim" => Mode::Periodic,
                    _ => Mode::Ordinary,
                };
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.rv_hh_spec()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvCyclicDimension(spec, mode)
            }
            Token::Word("rvNCdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.rv_nc_spec()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvNCDimension(spec)
            }
            Token::Word("rvHHdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.rv_hh_spec()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvHHDimension(spec)
            }
            Token::Word("rvpathextdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.rv_bar_spec()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvPathDimension(spec)
            }
            Token::Word("rvbaredim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.rv_bar_page_spec()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvBarPageDimension(spec)
            }
            Token::Word("rvbardim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let grade = self.rv_bar_spec()?;
                let length = if self.eat(Token::Comma)? {
                    Some(self.number()?.0)
                } else {
                    None
                };
                self.expect(Token::RBracket)?;
                StatementForm::RvBarDimension(grade, length)
            }
            Token::Word("rvhochster") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let ideal = self.label()?;
                self.expect(Token::Comma)?;
                let vertices = self.exponents()?;
                self.expect(Token::Comma)?;
                let degree = self.signed_degree()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvHochster(ideal, vertices, degree)
            }
            Token::Word(name @ ("rvbetti" | "rvmbetti")) => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let grade = self.rv_grade(name == "rvmbetti")?;
                self.expect(Token::RBracket)?;
                StatementForm::RvBetti(grade)
            }
            Token::Word("rvminors") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.rv_matrix(0)?;
                let (size, _) = self.number()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvMinors(spec, size)
            }
            Token::Word("jetdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.jet_page(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::JetDimension(spec)
            }
            Token::Word("rvpath") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let path = self.jet_path()?;
                self.expect(Token::RBracket)?;
                StatementForm::JetPath(name, path)
            }
            Token::Word("rvjet") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (source, target) = self.map_endpoints(0)?;
                let path = if self.token() == Token::LBracket {
                    JetPathArg::Inline(self.jet_path()?)
                } else {
                    JetPathArg::Named(self.label()?)
                };
                self.expect(Token::RBracket)?;
                StatementForm::RankJet(RankJetSpec {
                    source: Box::new(source),
                    target: Box::new(target),
                    path,
                })
            }
            Token::Word("rvorbit") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (modulus, _) = self.number()?;
                self.expect(Token::Comma)?;
                let point = self.exponents()?;
                self.expect(Token::RBracket)?;
                StatementForm::RankOrbit(modulus, point)
            }
            Token::Word(name @ ("rvclosed" | "rvzeta")) => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (source, target) = self.map_endpoints(0)?;
                let selected = self.exponents()?;
                self.expect(Token::Comma)?;
                let (max_degree, _) = self.number()?;
                self.expect(Token::RBracket)?;
                let series = if name == "rvclosed" {
                    RankSeries::Closed
                } else {
                    RankSeries::Zeta
                };
                StatementForm::RankArithmetic(
                    series,
                    RankArithmeticSpec {
                        source: Box::new(source),
                        target: Box::new(target),
                        selected,
                        max_degree,
                    },
                )
            }
            Token::Word("rvdim" | "rvspectrum") => {
                let affine = self.token() == Token::Word("rvdim");
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.rank_variety_spec(0, affine)?;
                self.expect(Token::RBracket)?;
                StatementForm::RankVariety(spec)
            }
            Token::Word("qextdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.q_spec(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::QExtDimension(spec)
            }
            Token::Word("vextdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.v_spec(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::VExtDimension(spec)
            }
            Token::Word("vbockdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.vbock_spec(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::VBockDimension(spec)
            }
            Token::Word("vbars") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.bar_family(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::VBars(spec)
            }
            Token::Word("vbardim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.v_spec(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::VBarDimension(spec)
            }
            Token::Word("mapmaydim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let spec = self.map_may_spec(0)?;
                self.expect(Token::RBracket)?;
                StatementForm::MapMayDimension(spec)
            }
            Token::Word(name @ ("mayandim" | "ceandim")) => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (n, r, s, t, f) = self.profile_page_degree()?;
                self.expect(Token::RBracket)?;
                if name == "ceandim" {
                    StatementForm::CeDimension(n, r, s, t, f)
                } else {
                    StatementForm::MayAnDimension(n, r, s, t, f)
                }
            }
            Token::Word("maydim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (r, s, t, f) = self.may_degree()?;
                self.expect(Token::RBracket)?;
                StatementForm::MayDimension(r, s, t, f)
            }
            Token::Word("pmaydim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (p, r, s, t, f) = self.pmay_degree()?;
                self.expect(Token::RBracket)?;
                StatementForm::PMayDimension(p, r, s, t, f)
            }
            Token::Word(name @ ("pmayandim" | "pceandim")) => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (p, n, r, s, t, f) = self.prime_profile_page_degree()?;
                self.expect(Token::RBracket)?;
                if name == "pceandim" {
                    StatementForm::PCeAnDimension(p, n, r, s, t, f)
                } else {
                    StatementForm::PMayAnDimension(p, n, r, s, t, f)
                }
            }
            Token::Word("homotopy") => {
                self.advance()?;
                let left = self.expr(0)?;
                self.expect(Token::Tilde)?;
                StatementForm::Homotopy(left, self.expr(0)?)
            }
            _ => return Err(self.expected("a statement")),
        };
        self.expect(Token::Semi)?;
        Ok(Some(Statement { span, kind }))
    }

    fn prime_bidegree(&mut self) -> Result<(usize, usize, usize)> {
        let (p, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (s, t) = self.bidegree()?;
        Ok((p, s, t))
    }

    fn prime_profile_bidegree(&mut self) -> Result<(usize, usize, usize, usize)> {
        let (p, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (n, s, t) = self.profile_bidegree()?;
        Ok((p, n, s, t))
    }

    fn bidegree(&mut self) -> Result<(usize, usize)> {
        let (s, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (t, _) = self.number()?;
        Ok((s, t))
    }

    fn profile_bidegree(&mut self) -> Result<(usize, usize, usize)> {
        let (n, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (s, t) = self.bidegree()?;
        Ok((n, s, t))
    }

    fn space_bidegree(&mut self, nest: usize) -> Result<(ParsedType, usize, usize)> {
        let base = self.ty(nest + 1)?;
        self.expect(Token::Comma)?;
        let (s, t) = self.bidegree()?;
        Ok((base, s, t))
    }

    fn rv_generators(&mut self) -> Result<Vec<ParsedExpr>> {
        self.expect(Token::LBrace)?;
        let mut generators = Vec::new();
        while self.token() != Token::RBrace {
            generators.push(self.expr(0)?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(generators)
    }

    fn rv_matrix(&mut self, nest: usize) -> Result<RvMatrixSpec> {
        let base = self.ty(nest + 1)?;
        self.expect(Token::Comma)?;
        let selected = self.exponents()?;
        self.expect(Token::Comma)?;
        Ok(RvMatrixSpec {
            base: Box::new(base),
            selected,
        })
    }

    fn rv_coordinates(&mut self) -> Result<Vec<Generator>> {
        self.expect(Token::LBracket)?;
        let mut result = Vec::new();
        while self.token() != Token::RBracket {
            self.expect(Token::LParen)?;
            let (degree, index) = self.bidegree()?;
            self.expect(Token::RParen)?;
            result.push(Generator { degree, index });
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBracket)?;
        Ok(result)
    }

    fn rank_variety_spec(&mut self, nest: usize, affine: bool) -> Result<RankVarietySpec> {
        let (source, target) = self.map_endpoints(nest)?;
        let selected = self.exponents()?;
        self.expect(Token::Comma)?;
        let (modulus, _) = self.number()?;
        let point = if affine {
            self.expect(Token::Comma)?;
            Some(self.exponents()?)
        } else {
            None
        };
        Ok(RankVarietySpec {
            source: Box::new(source),
            target: Box::new(target),
            selected,
            modulus,
            point,
        })
    }

    fn jet_page(&mut self, nest: usize) -> Result<JetPageSpec> {
        let (source, target) = self.map_endpoints(nest)?;
        let path = self.label()?;
        self.expect(Token::Comma)?;
        let (r, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (parity, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (f, _) = self.number()?;
        Ok(JetPageSpec {
            chain: JetHomSpec {
                source: Box::new(source),
                target: Box::new(target),
                path,
                parity,
            },
            r,
            f,
        })
    }

    fn jet_path(&mut self) -> Result<JetPathSpec> {
        let selected = self.exponents()?;
        self.expect(Token::Comma)?;
        let (modulus, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (order, _) = self.number()?;
        self.expect(Token::Comma)?;
        self.expect(Token::LBracket)?;
        let mut coefficients = Vec::new();
        while self.token() != Token::RBracket {
            coefficients.push(self.exponents()?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBracket)?;
        Ok(JetPathSpec {
            selected,
            modulus,
            order,
            coefficients,
        })
    }

    fn map_endpoints(&mut self, nest: usize) -> Result<(ParsedType, ParsedType)> {
        let source = self.ty(nest + 1)?;
        self.expect(Token::Comma)?;
        let target = self.ty(nest + 1)?;
        self.expect(Token::Comma)?;
        Ok((source, target))
    }

    fn map_bidegree(&mut self, nest: usize) -> Result<(ParsedType, ParsedType, usize, i64)> {
        let (source, target) = self.map_endpoints(nest)?;
        let (s, _) = self.number()?;
        self.expect(Token::Comma)?;
        let t = self.signed_degree()?;
        Ok((source, target, s, t))
    }

    fn q_spec(&mut self, nest: usize) -> Result<QSpec> {
        self.koszul_spec(nest, |parser| parser.number().map(|(s, _)| s))
    }

    fn v_spec(&mut self, nest: usize) -> Result<QSpec<i64>> {
        self.koszul_spec(nest, Self::signed_degree)
    }

    fn vbock_spec(&mut self, nest: usize) -> Result<VBockSpec> {
        let QSpec {
            source,
            target,
            selected,
            s: (r, s),
            t,
        } = self.koszul_spec(nest, |parser| {
            let (r, _) = parser.number()?;
            parser.expect(Token::Comma)?;
            Ok((r, parser.signed_degree()?))
        })?;
        self.expect(Token::Comma)?;
        let (f, _) = self.number()?;
        Ok(VBockSpec {
            localized: QSpec {
                source,
                target,
                selected,
                s,
                t,
            },
            r,
            f,
        })
    }

    fn bar_family(&mut self, nest: usize) -> Result<BarFamilySpec> {
        let (source, target) = self.map_endpoints(nest)?;
        let indices = self.exponents()?;
        Ok(BarFamilySpec {
            source,
            target,
            indices,
        })
    }

    fn koszul_spec<S>(
        &mut self,
        nest: usize,
        filtration: impl FnOnce(&mut Self) -> Result<S>,
    ) -> Result<QSpec<S>> {
        let (source, target) = self.map_endpoints(nest)?;
        let selected = self.exponents()?;
        self.expect(Token::Comma)?;
        let s = filtration(self)?;
        self.expect(Token::Comma)?;
        let t = self.signed_degree()?;
        Ok(QSpec {
            source: Box::new(source),
            target: Box::new(target),
            selected,
            s,
            t,
        })
    }

    fn map_may_spec(&mut self, nest: usize) -> Result<MapMaySpec> {
        let (source, target) = self.map_endpoints(nest)?;
        let (r, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (s, _) = self.number()?;
        self.expect(Token::Comma)?;
        let t = self.signed_degree()?;
        self.expect(Token::Comma)?;
        let (f, _) = self.number()?;
        Ok(MapMaySpec {
            source: Box::new(source),
            target: Box::new(target),
            r,
            s,
            t,
            f,
        })
    }

    fn map_generator(&mut self, nest: usize) -> Result<MapGenerator> {
        let base = self.ty(nest + 1)?;
        self.expect(Token::Comma)?;
        let (degree, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (index, span) = self.number()?;
        Ok(MapGenerator {
            base,
            degree,
            index,
            span,
        })
    }

    fn profile_page_degree(&mut self) -> Result<(usize, usize, usize, usize, usize)> {
        let (n, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (r, s, t, f) = self.may_degree()?;
        Ok((n, r, s, t, f))
    }

    fn may_degree(&mut self) -> Result<(usize, usize, usize, usize)> {
        let (r, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (s, t) = self.bidegree()?;
        self.expect(Token::Comma)?;
        let (f, _) = self.number()?;
        Ok((r, s, t, f))
    }

    fn pmay_degree(&mut self) -> Result<(usize, usize, usize, usize, usize)> {
        let (prime, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (r, s, t, f) = self.may_degree()?;
        Ok((prime, r, s, t, f))
    }

    fn prime_profile_page_degree(&mut self) -> Result<(usize, usize, usize, usize, usize, usize)> {
        let (prime, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (n, r, s, t, f) = self.profile_page_degree()?;
        Ok((prime, n, r, s, t, f))
    }

    fn as_operator_word(&mut self) -> Result<Vec<AsOperator>> {
        self.expect(Token::LBracket)?;
        let mut word = Vec::new();
        while self.token() != Token::RBracket {
            let operator = match self.token() {
                Token::Word("F") => AsOperator::F,
                Token::Word("V") => AsOperator::V,
                _ => {
                    return Err(Error::new(
                        self.span(),
                        "expected F or V in an Artin–Schreier word",
                    ))
                }
            };
            word.push(operator);
            self.advance()?;
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBracket)?;
        Ok(word)
    }

    fn as_operator_words(&mut self) -> Result<Vec<Vec<AsOperator>>> {
        self.expect(Token::LBracket)?;
        let mut words = Vec::new();
        while self.token() != Token::RBracket {
            words.push(self.as_operator_word()?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBracket)?;
        Ok(words)
    }

    fn bp_vpage_spec(&mut self) -> Result<BpVPageSpec> {
        let (p, n, r) = self.profile_bidegree()?;
        self.expect(Token::Comma)?;
        let (a, h, t, f) = self.prime_profile_bidegree()?;
        Ok((p, n, r, a, h, t, f))
    }

    fn ty(&mut self, nest: usize) -> Result<ParsedType> {
        check_nesting(nest, self.span())?;
        let mut head = self.type_atom(nest)?;
        while self.eat(Token::Star)? {
            let right = self.type_atom(nest)?;
            let span = head.span;
            let depth = head.depth.max(right.depth) + 1;
            head = ty(
                span,
                TypeForm::Product(Box::new(head), Box::new(right)),
                depth,
            )?;
        }
        if self.eat(Token::Arrow)? {
            let right = self.ty(nest + 1)?;
            let span = head.span;
            let depth = head.depth.max(right.depth) + 1;
            head = ty(
                span,
                TypeForm::Arrow(Box::new(head), Box::new(right)),
                depth,
            )?;
        }
        Ok(head)
    }

    fn type_atom(&mut self, nest: usize) -> Result<ParsedType> {
        let span = self.span();
        if self.eat(Token::LParen)? {
            let ty = self.ty(nest + 1)?;
            self.expect(Token::RParen)?;
            return Ok(ty);
        }
        if self.eat(Token::Word("Star"))? {
            self.expect(Token::LParen)?;
            let base = self.ty(nest + 1)?;
            self.expect(Token::Comma)?;
            let point = self.label()?;
            self.expect(Token::RParen)?;
            let depth = base.depth + 1;
            return ty(span, TypeForm::Star(Box::new(base), point), depth);
        }
        let name = self.label()?;
        if matches!(name.text.as_str(), "PCobarAn" | "PExtAn") && self.eat(Token::LParen)? {
            let (p, n, s, t) = self.prime_profile_bidegree()?;
            self.expect(Token::RParen)?;
            let kind = if name.text == "PCobarAn" {
                TypeForm::PCobarAn(p, n, s, t)
            } else {
                TypeForm::PExtAn(p, n, s, t)
            };
            return ty(span, kind, 0);
        }
        if name.text == "RhoBock" && self.eat(Token::LParen)? {
            let grade = self.rho_grade()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RhoBock(grade), 0);
        }
        if matches!(
            name.text.as_str(),
            "MCobar" | "MExt" | "RMCobar" | "RMExt" | "RMTriple" | "RhoLoc" | "RhoMin"
        ) && self.eat(Token::LParen)?
        {
            let grade = self.motivic_grade()?;
            self.expect(Token::RParen)?;
            let form = match name.text.as_str() {
                "MCobar" => TypeForm::MCobar(grade),
                "RMCobar" => TypeForm::RMCobar(grade),
                "RMExt" => TypeForm::RMExt(grade),
                "RMTriple" => TypeForm::RMTriple(grade),
                "RhoLoc" => TypeForm::RhoLoc(grade),
                "RhoMin" => TypeForm::RhoMin(grade),
                _ => TypeForm::MExt(grade),
            };
            return ty(span, form, 0);
        }
        if matches!(name.text.as_str(), "PCobar" | "PExt") && self.eat(Token::LParen)? {
            let (p, s, t) = self.prime_bidegree()?;
            self.expect(Token::RParen)?;
            let kind = if name.text == "PCobar" {
                TypeForm::PCobar(p, s, t)
            } else {
                TypeForm::PExt(p, s, t)
            };
            return ty(span, kind, 0);
        }
        if matches!(name.text.as_str(), "CobarAn" | "ExtAn") && self.eat(Token::LParen)? {
            let (n, s, t) = self.profile_bidegree()?;
            self.expect(Token::RParen)?;
            let kind = if name.text == "CobarAn" {
                TypeForm::CobarAn(n, s, t)
            } else {
                TypeForm::ExtAn(n, s, t)
            };
            return ty(span, kind, 0);
        }
        if name.text == "XExt" && self.eat(Token::LParen)? {
            let (base, s, t) = self.space_bidegree(nest)?;
            self.expect(Token::RParen)?;
            let depth = base.depth + 1;
            return ty(span, TypeForm::XExt(Box::new(base), s, t), depth);
        }
        if matches!(name.text.as_str(), "MapCobar" | "MapExt") && self.eat(Token::LParen)? {
            let (source, target, s, t) = self.map_bidegree(nest)?;
            self.expect(Token::RParen)?;
            let depth = source.depth.max(target.depth) + 1;
            let source = Box::new(source);
            let target = Box::new(target);
            let kind = if name.text == "MapCobar" {
                TypeForm::MapCobar(source, target, s, t)
            } else {
                TypeForm::MapExt(source, target, s, t)
            };
            return ty(span, kind, depth);
        }
        if name.text == "JetPage" && self.eat(Token::LParen)? {
            let spec = self.jet_page(nest)?;
            self.expect(Token::RParen)?;
            let depth = spec.chain.source.depth.max(spec.chain.target.depth) + 1;
            return ty(span, TypeForm::JetPage(spec), depth);
        }
        if name.text == "JetHom" && self.eat(Token::LParen)? {
            let (source, target) = self.map_endpoints(nest)?;
            let path = self.label()?;
            self.expect(Token::Comma)?;
            let (parity, _) = self.number()?;
            self.expect(Token::RParen)?;
            let depth = source.depth.max(target.depth) + 1;
            return ty(
                span,
                TypeForm::JetHom(JetHomSpec {
                    source: Box::new(source),
                    target: Box::new(target),
                    path,
                    parity,
                }),
                depth,
            );
        }
        if matches!(name.text.as_str(), "QKoszul" | "QExt") && self.eat(Token::LParen)? {
            let spec = self.q_spec(nest)?;
            self.expect(Token::RParen)?;
            let depth = spec.source.depth.max(spec.target.depth) + 1;
            let kind = if name.text == "QKoszul" {
                TypeForm::QKoszul(spec)
            } else {
                TypeForm::QExt(spec)
            };
            return ty(span, kind, depth);
        }
        if matches!(name.text.as_str(), "VKoszul" | "VExt") && self.eat(Token::LParen)? {
            let spec = self.v_spec(nest)?;
            self.expect(Token::RParen)?;
            let depth = spec.source.depth.max(spec.target.depth) + 1;
            let kind = if name.text == "VKoszul" {
                TypeForm::VKoszul(spec)
            } else {
                TypeForm::VExt(spec)
            };
            return ty(span, kind, depth);
        }
        if name.text == "BPVPage" && self.eat(Token::LParen)? {
            let spec = self.bp_vpage_spec()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::BPVPage(spec), 0);
        }
        if matches!(name.text.as_str(), "BPBock" | "BPNov") && self.eat(Token::LParen)? {
            let mode = if name.text == "BPBock" {
                BpPageMode::Bockstein
            } else {
                BpPageMode::Novikov
            };
            let spec = self.prime_profile_page_degree()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::BpPage(mode, spec), 0);
        }
        if name.text == "BPHeight" && self.eat(Token::LParen)? {
            let spec = self.bp_height_spec()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::BPHeight(spec), 0);
        }
        if name.text == "BPRelExt" && self.eat(Token::LParen)? {
            let shape = self.bp_relative_shape()?;
            self.expect(Token::Comma)?;
            let (t, _) = self.number()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::BPRelExt(shape, t), 0);
        }
        if name.text == "BPRel" && self.eat(Token::LParen)? {
            let shape = self.bp_relative_shape()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::BPRel(shape), 0);
        }
        if matches!(name.text.as_str(), "BPExt" | "BPTriple") && self.eat(Token::LParen)? {
            let (spec, t) = self.bp_ext_grade()?;
            self.expect(Token::RParen)?;
            let kind = if name.text == "BPExt" {
                TypeForm::BPExt(spec, t)
            } else {
                TypeForm::BPTriple(spec, t)
            };
            return ty(span, kind, 0);
        }
        if matches!(name.text.as_str(), "BPIExt" | "BPThickExt") && self.eat(Token::LParen)? {
            let grade = self.bp_mod_ext_grade(name.text == "BPThickExt")?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::BPModExt(grade), 0);
        }
        if name.text == "BPThick" && self.eat(Token::LParen)? {
            let (prime, height, arrows) = self.prime_bidegree()?;
            self.expect(Token::Comma)?;
            let (h, _) = self.number()?;
            self.expect(Token::Comma)?;
            let (power, _) = self.number()?;
            self.expect(Token::RParen)?;
            return ty(
                span,
                TypeForm::BPThick(BpSpec(prime, height, arrows), h, power),
                0,
            );
        }
        if matches!(name.text.as_str(), "BP" | "BPI") && self.eat(Token::LParen)? {
            let (prime, height, arrows) = self.prime_bidegree()?;
            let depth = if name.text == "BPI" {
                self.expect(Token::Comma)?;
                Some(self.number()?.0)
            } else {
                None
            };
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::BP(BpSpec(prime, height, arrows), depth), 0);
        }
        let cyclic = match name.text.as_str() {
            "RVHN" | "RVHNChain" => Some(crate::rv_cyclic_homology::Mode::Negative),
            "RVHP" | "RVHPChain" => Some(crate::rv_cyclic_homology::Mode::Periodic),
            "RVHC" | "RVHCChain" => Some(crate::rv_cyclic_homology::Mode::Ordinary),
            _ => None,
        };
        if let Some(mode) = cyclic {
            if self.eat(Token::LParen)? {
                let spec = self.rv_hh_spec()?;
                self.expect(Token::RParen)?;
                return ty(
                    span,
                    TypeForm::RVUntruncatedCyclic(spec, mode, !name.text.ends_with("Chain")),
                    0,
                );
            }
        }
        if matches!(name.text.as_str(), "RVNCChain" | "RVNC") && self.eat(Token::LParen)? {
            let spec = self.rv_nc_spec()?;
            self.expect(Token::RParen)?;
            return ty(
                span,
                TypeForm::RVNegativeCyclic(spec, name.text == "RVNC"),
                0,
            );
        }
        if matches!(name.text.as_str(), "RVHHChain" | "RVHH") && self.eat(Token::LParen)? {
            let spec = self.rv_hh_spec()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVHochschild(spec, name.text == "RVHH"), 0);
        }
        if matches!(name.text.as_str(), "MixedMin" | "MixedMinU") && self.eat(Token::LParen)? {
            let total = name.text == "MixedMinU";
            let name = self.label()?;
            self.expect(Token::Comma)?;
            let b = self.signed_degree()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::MixedMinimal(name, b, total), 0);
        }
        if matches!(name.text.as_str(), "MixedChain" | "MixedU") && self.eat(Token::LParen)? {
            let complex = self.label()?;
            self.expect(Token::Comma)?;
            let degree = self.signed_degree()?;
            self.expect(Token::RParen)?;
            return ty(
                span,
                TypeForm::MixedChain(complex, degree, name.text == "MixedU"),
                0,
            );
        }
        if name.text == "MixedE" && self.eat(Token::LParen)? {
            let spec = self.mixed_spec()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::MixedPage(spec), 0);
        }
        if name.text == "RVBarCobar" && self.eat(Token::LParen)? {
            let grade = self.rv_bar_spec()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVBarCobar(grade), 0);
        }
        if name.text == "RVBarE" && self.eat(Token::LParen)? {
            let grade = self.rv_bar_page_spec()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVBarE(grade), 0);
        }
        if matches!(name.text.as_str(), "RVBarDual" | "RVPathExt") && self.eat(Token::LParen)? {
            let spec = self.rv_bar_spec()?;
            self.expect(Token::RParen)?;
            let kind = if name.text == "RVBarDual" {
                TypeForm::RVBarDual(spec)
            } else {
                TypeForm::RVPathExt(spec)
            };
            return ty(span, kind, 0);
        }
        if name.text == "RVBar" && self.eat(Token::LParen)? {
            let grade = self.rv_bar_spec()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVBar(grade), 0);
        }
        if matches!(
            name.text.as_str(),
            "RVKoszul" | "RVTor" | "RVMKoszul" | "RVMTor"
        ) && self.eat(Token::LParen)?
        {
            let grade = self.rv_grade(name.text.starts_with("RVM"))?;
            self.expect(Token::RParen)?;
            return ty(
                span,
                TypeForm::RVComplex(grade, name.text.ends_with("Tor")),
                0,
            );
        }
        if name.text == "ShLocal" && self.eat(Token::LParen)? {
            let (name, point, mode) = self.sheaf_local_site()?;
            self.expect(Token::Comma)?;
            let degree = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::SheafLocal(name, point, mode, degree), 0);
        }
        if name.text == "ShE" && self.eat(Token::LParen)? {
            let name = self.label()?;
            self.expect(Token::Comma)?;
            let grade = self.sheaf_page_grade()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::SheafPage(name, grade), 0);
        }
        if matches!(name.text.as_str(), "ShTot" | "ShHyper") && self.eat(Token::LParen)? {
            let form = if name.text == "ShTot" {
                SheafForm::Cochain
            } else {
                SheafForm::Class
            };
            let name = self.label()?;
            self.expect(Token::Comma)?;
            let degree = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::SheafComplex(form, name, degree), 0);
        }
        if matches!(name.text.as_str(), "ShHomC" | "ShExt") && self.eat(Token::LParen)? {
            let form = if name.text == "ShHomC" {
                SheafForm::Cochain
            } else {
                SheafForm::Class
            };
            let (source, target, degree) = self.sheaf_hom_grade()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::SheafHom(form, source, target, degree), 0);
        }
        if matches!(name.text.as_str(), "ShC" | "ShH") && self.eat(Token::LParen)? {
            let form = if name.text == "ShC" {
                SheafForm::Cochain
            } else {
                SheafForm::Class
            };
            let sheaf = self.label()?;
            self.expect(Token::Comma)?;
            let degree = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::Sheaf(form, sheaf, degree), 0);
        }
        if name.text == "ASDeRham" && self.eat(Token::LParen)? {
            let cover = self.label()?;
            self.expect(Token::Comma)?;
            let modulus = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::ASDeRham(cover, modulus), 0);
        }
        if name.text == "EOType" && self.eat(Token::LParen)? {
            let genus = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::EOType(genus), 0);
        }
        if name.text == "BinaryScalar" && self.eat(Token::LParen)? {
            let modulus = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::BinaryScalar(modulus), 0);
        }
        if name.text == "ASOddTate" && self.eat(Token::LParen)? {
            let cover = self.label()?;
            self.expect(Token::Comma)?;
            let degree = self.signed_number("Tate degree")?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::ASOddTate(cover, degree), 0);
        }
        if name.text == "ASOddDeRham" && self.eat(Token::LParen)? {
            let cover = self.label()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::ASOddDeRham(cover), 0);
        }
        if name.text == "PrimeScalar" && self.eat(Token::LParen)? {
            let p = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::PrimeScalar(p), 0);
        }
        if name.text == "ASOddDiff" && self.eat(Token::LParen)? {
            let cover = self.label()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::ASOddDiff(cover), 0);
        }
        if matches!(name.text.as_str(), "ASDiff" | "WittDiff") && self.eat(Token::LParen)? {
            let family = if name.text == "WittDiff" {
                DiffFamily::Witt
            } else {
                DiffFamily::ArtinSchreier
            };
            let cover = self.label()?;
            self.expect(Token::Comma)?;
            let modulus = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::CurveDiff(family, cover, modulus), 0);
        }
        if name.text == "RVPoly" && self.eat(Token::LParen)? {
            let selected = self.exponents()?;
            self.expect(Token::Comma)?;
            let (degree, _) = self.number()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVPoly(selected, degree), 0);
        }
        if name.text == "VBock" && self.eat(Token::LParen)? {
            let spec = self.vbock_spec(nest)?;
            self.expect(Token::RParen)?;
            let depth = spec.localized.source.depth.max(spec.localized.target.depth) + 1;
            return ty(span, TypeForm::VBock(spec), depth);
        }
        if name.text == "RVWzw" && self.eat(Token::LParen)? {
            let (rank, level) = self.fusion_grade()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVWzw(rank, level), 0);
        }
        if name.text == "RVCycloQ" && self.eat(Token::LParen)? {
            let (order, _) = self.number()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVCycloQ(order), 0);
        }
        if name.text == "RVCyclo" && self.eat(Token::LParen)? {
            let (order, _) = self.number()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVCyclo(order), 0);
        }
        if name.text == "RVFusion" && self.eat(Token::LParen)? {
            let (rank, _) = self.number()?;
            self.expect(Token::Comma)?;
            let (level, _) = self.number()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVFusion(rank, level), 0);
        }
        if name.text == "RVOre" {
            return ty(span, TypeForm::RVOre, 0);
        }
        if name.text == "RVQGrP" && self.eat(Token::LParen)? {
            let (prime, rank, ambient) = self.quantum_mod_grade()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVQGrP(prime, rank, ambient), 0);
        }
        if name.text == "RVQGr" && self.eat(Token::LParen)? {
            let (rank, ambient) = self.quantum_grade()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVQGr(rank, ambient), 0);
        }
        if name.text == "RVHive" && self.eat(Token::LParen)? {
            let weights = self.hive_spec()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVHive(weights), 0);
        }
        if name.text == "RVSatake" && self.eat(Token::LParen)? {
            let (rank, _) = self.number()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVSatake(rank), 0);
        }
        if matches!(name.text.as_str(), "RVIwHecke" | "RVSphHecke") && self.eat(Token::LParen)? {
            let (rank, _) = self.number()?;
            self.expect(Token::RParen)?;
            let kind = if name.text == "RVIwHecke" {
                TypeForm::RVIwHecke(rank)
            } else {
                TypeForm::RVSphHecke(rank)
            };
            return ty(span, kind, 0);
        }
        if matches!(name.text.as_str(), "RVHom" | "RVLaurent") && self.eat(Token::LParen)? {
            let spec = self.rank_variety_spec(nest, true)?;
            self.expect(Token::RParen)?;
            let depth = spec.source.depth.max(spec.target.depth) + 1;
            let kind = if name.text == "RVHom" {
                TypeForm::RVHom(spec)
            } else {
                TypeForm::RVLaurent(spec)
            };
            return ty(span, kind, depth);
        }
        if name.text == "MapMay" && self.eat(Token::LParen)? {
            let spec = self.map_may_spec(nest)?;
            self.expect(Token::RParen)?;
            let depth = spec.source.depth.max(spec.target.depth) + 1;
            return ty(span, TypeForm::MapMay(spec), depth);
        }
        if matches!(name.text.as_str(), "MayAn" | "CEAn") && self.eat(Token::LParen)? {
            let (n, r, s, t, f) = self.profile_page_degree()?;
            self.expect(Token::RParen)?;
            let kind = if name.text == "CEAn" {
                TypeForm::CeAn(n, r, s, t, f)
            } else {
                TypeForm::MayAn(n, r, s, t, f)
            };
            return ty(span, kind, 0);
        }
        if name.text == "May" && self.eat(Token::LParen)? {
            let (r, s, t, f) = self.may_degree()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::May(r, s, t, f), 0);
        }
        if name.text == "PMay" && self.eat(Token::LParen)? {
            let (p, r, s, t, f) = self.pmay_degree()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::PMay(p, r, s, t, f), 0);
        }
        if matches!(name.text.as_str(), "PMayAn" | "PCEAn") && self.eat(Token::LParen)? {
            let (p, n, r, s, t, f) = self.prime_profile_page_degree()?;
            self.expect(Token::RParen)?;
            let kind = if name.text == "PCEAn" {
                TypeForm::PCeAn(p, n, r, s, t, f)
            } else {
                TypeForm::PMayAn(p, n, r, s, t, f)
            };
            return ty(span, kind, 0);
        }
        if matches!(
            name.text.as_str(),
            "A" | "Bar" | "H2" | "H" | "Susp" | "Stable" | "Cobar" | "Ext"
        ) && self.eat(Token::LParen)?
        {
            let (kind, depth) = if matches!(name.text.as_str(), "A" | "Bar" | "H2") {
                let (degree, _) = self.number()?;
                let kind = match name.text.as_str() {
                    "Bar" => TypeForm::Bar(degree),
                    "H2" => TypeForm::CyclicH2(degree),
                    _ => TypeForm::Operations(degree),
                };
                (kind, 0)
            } else if matches!(name.text.as_str(), "Cobar" | "Ext") {
                let (s, t) = self.bidegree()?;
                let kind = match name.text.as_str() {
                    "Cobar" => TypeForm::Cobar(s, t),
                    _ => TypeForm::Ext(s, t),
                };
                (kind, 0)
            } else {
                let base = self.ty(nest + 1)?;
                let depth = base.depth + 1;
                let kind = if name.text == "Susp" {
                    TypeForm::Suspension(Box::new(base))
                } else {
                    self.expect(Token::Comma)?;
                    let (first, _) = self.number()?;
                    if name.text == "H" {
                        TypeForm::Cohomology(Box::new(base), first)
                    } else {
                        self.expect(Token::Comma)?;
                        let (degree, _) = self.number()?;
                        TypeForm::Stable(Box::new(base), first, degree)
                    }
                };
                (kind, depth)
            };
            self.expect(Token::RParen)?;
            return ty(span, kind, depth);
        }
        ty(span, TypeForm::Name(name), 0)
    }

    fn expr(&mut self, nest: usize) -> Result<ParsedExpr> {
        check_nesting(nest, self.span())?;
        let mut head = self.expr_atom(nest)?;
        while self.eat(Token::LParen)? {
            if let ExprForm::Name(name) = &head.kind {
                if matches!(name.text.as_str(), "massey" | "powerdef" | "bracket")
                    || TransferOp::from_name(&name.text).is_some()
                    || RvTBracketOp::from_name(&name.text).is_some()
                    || MapBracketOp::from_name(&name.text).is_some()
                    || RvTTransferOp::from_name(&name.text).is_some()
                    || RvTSystemOp::from_name(&name.text).is_some()
                    || MapSystemOp::from_name(&name.text).is_some()
                    || OddSystemOp::from_name(&name.text).is_some()
                    || MapTransferOp::from_name(&name.text).is_some()
                {
                    head = self.variadic(head, nest)?;
                    continue;
                }
            }
            let arg = self.expr(nest + 1)?;
            if self.eat(Token::Comma)? {
                let primitive = match &head.kind {
                    ExprForm::Name(name) => BinaryOp::from_name(&name.text),
                    _ => None,
                };
                let Some(primitive) = primitive else {
                    return Err(Error::new(
                        head.span,
                        "a function call takes one argument; use nested calls or pass a pair",
                    ));
                };
                let right = self.expr(nest + 1)?;
                self.expect(Token::RParen)?;
                let depth = arg.depth.max(right.depth) + 1;
                head = expr(
                    head.span,
                    ExprForm::Binary(primitive, Box::new(arg), Box::new(right)),
                    depth,
                )?;
                continue;
            }
            self.expect(Token::RParen)?;
            let span = head.span;
            let depth = head.depth.max(arg.depth) + 1;
            head = expr(span, ExprForm::Call(Box::new(head), Box::new(arg)), depth)?;
        }
        Ok(head)
    }

    fn arguments(&mut self, nest: usize) -> Result<Vec<ParsedExpr>> {
        let mut arguments = Vec::new();
        while self.token() != Token::RParen {
            arguments.push(self.expr(nest + 1)?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RParen)?;
        Ok(arguments)
    }

    fn bp_triple_arguments(
        &mut self,
        nest: usize,
        span: Span,
        name: &str,
    ) -> Result<Box<[ParsedExpr; 3]>> {
        self.exact_arguments(nest, span, name, "three BPExt classes")
    }

    fn exact_arguments<const N: usize>(
        &mut self,
        nest: usize,
        span: Span,
        name: &str,
        expected: &str,
    ) -> Result<Box<[ParsedExpr; N]>> {
        let arguments = self
            .arguments(nest)?
            .try_into()
            .map_err(|_| Error::new(span, format!("{name} requires {expected}")))?;
        Ok(Box::new(arguments))
    }

    fn variadic(&mut self, head: ParsedExpr, nest: usize) -> Result<ParsedExpr> {
        let mut arguments = self.arguments(nest)?;
        let ExprForm::Name(name) = &head.kind else {
            unreachable!();
        };
        let massey = name.text == "massey";
        let rv_t_system = RvTSystemOp::from_name(&name.text);
        let map_system = MapSystemOp::from_name(&name.text);
        let odd_system = OddSystemOp::from_name(&name.text);
        let has_witnesses = massey
            || map_system.is_some()
            || rv_t_system.is_some()
            || odd_system.is_some_and(OddSystemOp::explicit);
        let span = head.span;
        let mut depth = arguments
            .iter()
            .map(|argument| argument.depth + 1)
            .max()
            .unwrap_or(0);
        if arguments.len() == 1 && (self.token() != Token::LBrace || !has_witnesses) {
            return expr(
                span,
                ExprForm::Call(Box::new(head), Box::new(arguments.pop().unwrap())),
                depth,
            );
        }
        if let Some(operation) = TransferOp::from_name(&name.text) {
            return expr(span, ExprForm::Transfer(operation, arguments), depth);
        }
        if let Some(operation) = RvTTransferOp::from_name(&name.text) {
            return expr(span, ExprForm::RvTTransfer(operation, arguments), depth);
        }
        if let Some(operation) = RvTBracketOp::from_name(&name.text) {
            return expr(span, ExprForm::RvTBracket(operation, arguments), depth);
        }
        if let Some(operation) = MapBracketOp::from_name(&name.text) {
            return expr(span, ExprForm::MapBracket(operation, arguments), depth);
        }
        if let Some(operation) = MapTransferOp::from_name(&name.text) {
            return expr(span, ExprForm::MapTransfer(operation, arguments), depth);
        }
        if let Some(operation) = rv_t_system {
            let witnesses = self.witnesses(nest, &mut depth)?;
            return expr(
                span,
                ExprForm::RvTSystem(operation, arguments, witnesses),
                depth,
            );
        }
        if let Some(operation) = map_system {
            let witnesses = self.witnesses(nest, &mut depth)?;
            return expr(
                span,
                ExprForm::MapSystem(operation, arguments, witnesses),
                depth,
            );
        }
        if let Some(operation) = odd_system {
            let witnesses = if operation.explicit() {
                self.witnesses(nest, &mut depth)?
            } else {
                Vec::new()
            };
            return expr(
                span,
                ExprForm::OddSystem(operation, arguments, witnesses),
                depth,
            );
        }
        if name.text == "bracket" {
            return expr(span, ExprForm::Bracket(arguments), depth);
        }
        if !massey {
            return expr(span, ExprForm::PowerDef(arguments), depth);
        }

        let witnesses = self.witnesses(nest, &mut depth)?;
        expr(span, ExprForm::Massey(arguments, witnesses), depth)
    }

    fn witnesses(&mut self, nest: usize, depth: &mut usize) -> Result<Vec<Witness>> {
        self.expect(Token::LBrace)?;
        let mut witnesses = Vec::new();
        while self.token() != Token::RBrace {
            let (from, span) = self.number()?;
            self.expect(Token::Less)?;
            let (to, _) = self.number()?;
            self.expect(Token::FatArrow)?;
            let value = self.expr(nest + 1)?;
            *depth = (*depth).max(value.depth + 1);
            witnesses.push(Witness {
                from,
                to,
                span,
                value,
            });
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(witnesses)
    }

    fn expr_atom(&mut self, nest: usize) -> Result<ParsedExpr> {
        let span = self.span();
        match self.token() {
            Token::Slash => {
                self.advance()?;
                let parameter = self.label()?;
                self.expect(Token::Colon)?;
                let domain = self.ty(nest + 1)?;
                self.expect(Token::FatArrow)?;
                let body = self.expr(nest + 1)?;
                let depth = body.depth + 1;
                expr(
                    span,
                    ExprForm::Lambda(parameter, domain, Box::new(body)),
                    depth,
                )
            }
            Token::LParen => {
                self.advance()?;
                let left = self.expr(nest + 1)?;
                if !self.eat(Token::Comma)? {
                    self.expect(Token::RParen)?;
                    return Ok(left);
                }
                let right = self.expr(nest + 1)?;
                self.expect(Token::RParen)?;
                let depth = left.depth.max(right.depth) + 1;
                expr(span, ExprForm::Pair(Box::new(left), Box::new(right)), depth)
            }
            Token::Word("fst") | Token::Word("snd") => {
                let first = self.token() == Token::Word("fst");
                self.advance()?;
                let value = self.paren_expr(nest)?;
                let depth = value.depth + 1;
                let kind = if first {
                    ExprForm::Fst(Box::new(value))
                } else {
                    ExprForm::Snd(Box::new(value))
                };
                expr(span, kind, depth)
            }
            Token::Word("res") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let point = self.label()?;
                self.expect(Token::RBracket)?;
                let base = self.paren_expr(nest)?;
                let depth = base.depth + 1;
                expr(span, ExprForm::Restrict(point, Box::new(base)), depth)
            }
            Token::Word("tab") | Token::Word("sew") => {
                let is_table = self.token() == Token::Word("tab");
                let keyword = self.token().description();
                self.advance()?;
                let ty = self.ty(nest + 1)?;
                let TypeForm::Arrow(domain, codomain) = ty.kind else {
                    return Err(Error::new(
                        ty.span,
                        format!("{keyword} requires a type of the form X -> Y"),
                    ));
                };
                self.expect(Token::LBrace)?;
                let mut entries = Vec::new();
                while self.token() != Token::RBrace {
                    let label = self.label()?;
                    self.expect(Token::FatArrow)?;
                    entries.push((label, self.expr(nest + 1)?));
                    if !self.eat(Token::Comma)? {
                        break;
                    }
                }
                self.expect(Token::RBrace)?;
                let depth = entries
                    .iter()
                    .map(|(_, value)| value.depth + 1)
                    .max()
                    .unwrap_or(0);
                let kind = if is_table {
                    ExprForm::Table(*domain, *codomain, entries)
                } else {
                    ExprForm::Sew(*domain, *codomain, entries)
                };
                expr(span, kind, depth)
            }
            Token::Word("Star") => {
                let ty = self.type_atom(nest)?;
                self.expect(Token::Dot)?;
                expr(span, ExprForm::Point(ty, self.label()?), 0)
            }
            Token::Word(_) => {
                let name = self.label()?;
                if matches!(name.text.as_str(), "rvHNgen" | "rvHPgen" | "rvHCgen")
                    && self.eat(Token::LBracket)?
                {
                    use crate::rv_cyclic_homology::Mode;
                    let mode = match name.text.as_str() {
                        "rvHNgen" => Mode::Negative,
                        "rvHPgen" => Mode::Periodic,
                        _ => Mode::Ordinary,
                    };
                    let spec = self.rv_hh_spec()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RvCyclicGenerator(spec, mode, index, at), 0);
                }
                if matches!(name.text.as_str(), "rvHNTerm" | "rvHPTerm" | "rvHCTerm")
                    && self.eat(Token::LBracket)?
                {
                    use crate::rv_cyclic_homology::Mode;
                    let mode = match name.text.as_str() {
                        "rvHNTerm" => Mode::Negative,
                        "rvHPTerm" => Mode::Periodic,
                        _ => Mode::Ordinary,
                    };
                    let power = self.signed_degree()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RvCyclicTerm(mode, power, Box::new(input)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "rvCyclicU" | "rvCyclicPart")
                    && self.eat(Token::LBracket)?
                {
                    let power = self.signed_degree()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    let expression = if name.text == "rvCyclicPart" {
                        ExprForm::RvCyclicPart(power, Box::new(input))
                    } else {
                        ExprForm::RvCyclicShift(
                            crate::rv_cyclic_homology::Shift::Periodic(power),
                            Box::new(input),
                        )
                    };
                    return expr(span, expression, depth);
                }
                if name.text == "rvCyclicS" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RvCyclicShift(
                            crate::rv_cyclic_homology::Shift::Suspension,
                            Box::new(input),
                        ),
                        depth,
                    );
                }
                if name.text == "rvNCgen" && self.eat(Token::LBracket)? {
                    let spec = self.rv_nc_spec()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RvNCGenerator(spec, index, at), 0);
                }
                if name.text == "rvNCTerm" && self.eat(Token::LBracket)? {
                    let length = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let power = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RvNCTerm(length, power, Box::new(input)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "rvNCU" | "rvNCPart" | "rvNCTruncate")
                    && self.eat(Token::LBracket)?
                {
                    let index = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    let op = match name.text.as_str() {
                        "rvNCPart" => ExprForm::RvNCPart(index, Box::new(input)),
                        "rvNCU" => ExprForm::RvNCTransfer(
                            crate::rv_negative::Transfer::Shift(index),
                            Box::new(input),
                        ),
                        _ => ExprForm::RvNCTransfer(
                            crate::rv_negative::Transfer::Truncate(index),
                            Box::new(input),
                        ),
                    };
                    return expr(span, op, depth);
                }
                if name.text == "rvHHgen" && self.eat(Token::LBracket)? {
                    let spec = self.rv_hh_spec()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RvHHGenerator(spec, index, at), 0);
                }
                if name.text == "rvpathextgen" && self.eat(Token::LBracket)? {
                    let spec = self.rv_bar_spec()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RvPathGenerator(spec, index, at), 0);
                }
                if name.text == "rvHHToMixed" && self.eat(Token::LBracket)? {
                    let name = self.label()?;
                    let f = if self.eat(Token::Comma)? {
                        Some(self.number()?.0)
                    } else {
                        None
                    };
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RvHHToMixed(name, f, Box::new(input)), depth);
                }
                if name.text == "mixunit" && self.eat(Token::LBracket)? {
                    let complex = self.label()?;
                    self.expect(Token::Comma)?;
                    let source = self.label()?;
                    self.expect(Token::Comma)?;
                    let target = self.label()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::MixedUnit(complex, source, target), 0);
                }
                if name.text == "mixId" && self.eat(Token::LBracket)? {
                    let complex = self.label()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::MixedIdentity(complex), 0);
                }
                if matches!(name.text.as_str(), "mixComp" | "mixTensor")
                    && self.eat(Token::LBracket)?
                {
                    let complex = self.label()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LParen)?;
                    let left = self.expr(nest + 1)?;
                    self.expect(Token::Comma)?;
                    let right = self.expr(nest + 1)?;
                    self.expect(Token::RParen)?;
                    let depth = 1 + left.depth.max(right.depth);
                    let form = if name.text == "mixComp" {
                        ExprForm::MixedComposition(complex, Box::new(left), Box::new(right))
                    } else {
                        ExprForm::MixedTensorProduct(complex, Box::new(left), Box::new(right))
                    };
                    return expr(span, form, depth);
                }
                if matches!(name.text.as_str(), "mixConeI" | "mixConeLift" | "mixConeP")
                    && self.eat(Token::LBracket)?
                {
                    use crate::mixed_cone::Operation;
                    let operation = match name.text.as_str() {
                        "mixConeI" => Operation::Include,
                        "mixConeLift" => Operation::Lift,
                        _ => Operation::Project,
                    };
                    let cone = self.label()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::MixedConeMap(cone, operation, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "mixConeShear" && self.eat(Token::LBracket)? {
                    let cone = self.label()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LParen)?;
                    let homotopy = self.expr(nest + 1)?;
                    self.expect(Token::Comma)?;
                    let input = self.expr(nest + 1)?;
                    self.expect(Token::RParen)?;
                    let depth = 1 + homotopy.depth.max(input.depth);
                    return expr(
                        span,
                        ExprForm::MixedConeShear(cone, Box::new(homotopy), Box::new(input)),
                        depth,
                    );
                }
                if name.text == "mixSwap" && self.eat(Token::LBracket)? {
                    let target = self.label()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::MixedTensorSwap(target, Box::new(input)),
                        depth,
                    );
                }
                if matches!(
                    name.text.as_str(),
                    "mixTensorCycle" | "mixTensorNull" | "mixTensorTor" | "mixTensorTorNull"
                ) && self.eat(Token::LBracket)?
                {
                    use crate::mixed_tensor::Certificate;
                    let kind = match name.text.as_str() {
                        "mixTensorCycle" => Certificate::Product,
                        "mixTensorNull" => Certificate::ProductNull,
                        "mixTensorTor" => Certificate::Secondary,
                        _ => Certificate::SecondaryNull,
                    };
                    let complex = self.label()?;
                    self.expect(Token::Comma)?;
                    let i = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let j = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::MixedTensorCertificate(complex, i, j, kind),
                        0,
                    );
                }
                if name.text == "mixchain" && self.eat(Token::LBracket)? {
                    let name = self.label()?;
                    self.expect(Token::Comma)?;
                    let b = self.signed_degree()?;
                    self.expect(Token::RBracket)?;
                    let terms = self.mixed_terms()?;
                    return expr(span, ExprForm::MixedLiteral(name, b, terms), 0);
                }
                if matches!(name.text.as_str(), "mixBarRep" | "mixBarNull")
                    && self.eat(Token::LBracket)?
                {
                    let complex = self.label()?;
                    self.expect(Token::Comma)?;
                    let i = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::MixedBarTerm(complex, i, name.text == "mixBarNull"),
                        0,
                    );
                }
                if name.text == "mixMinGen" && self.eat(Token::LBracket)? {
                    let name = self.label()?;
                    self.expect(Token::Comma)?;
                    let b = self.signed_degree()?;
                    self.expect(Token::Comma)?;
                    let (i, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::MixedMinGenerator(name, b, i, at), 0);
                }
                if matches!(name.text.as_str(), "mixMinP" | "mixMinH")
                    && self.eat(Token::LBracket)?
                {
                    use crate::mixed_minimal::Operation;
                    let op = if name.text == "mixMinP" {
                        Operation::Projection
                    } else {
                        Operation::Homotopy
                    };
                    let model = self.label()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::MixedMinMap(model, op, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "mixDelta" && self.eat(Token::LBracket)? {
                    let order = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::MixedDelta(order, Box::new(input)), depth);
                }
                if name.text == "mixgen" && self.eat(Token::LBracket)? {
                    let spec = self.mixed_spec()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::MixedGenerator(spec, index, at), 0);
                }
                if name.text == "mixU" && self.eat(Token::LBracket)? {
                    let exponent = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::MixedEmbed(exponent, Box::new(input)), depth);
                }
                if name.text == "mixPageClass" && self.eat(Token::LBracket)? {
                    let r = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let f = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::MixedClass(r, f, Box::new(input)), depth);
                }
                if name.text == "rvbaregen" && self.eat(Token::LBracket)? {
                    let spec = self.rv_bar_page_spec()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RvBarEGenerator(spec, index, at), 0);
                }
                if name.text == "rvBarClass" && self.eat(Token::LBracket)? {
                    let (r, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let length = self.signed_degree()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RvBarClass(r, length, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvHHWord" && self.eat(Token::LParen)? {
                    let inputs = self.arguments(nest)?;
                    let depth = inputs
                        .iter()
                        .map(|input| input.depth + 1)
                        .max()
                        .unwrap_or(0);
                    return expr(span, ExprForm::RvHHWord(inputs), depth);
                }
                if name.text == "rvHHLength" && self.eat(Token::LBracket)? {
                    let (length, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RvHHLength(length, Box::new(input)), depth);
                }
                if name.text == "rvBarWord" && self.eat(Token::LBracket)? {
                    let path = self.label()?;
                    self.expect(Token::Comma)?;
                    let cuts = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LParen)?;
                    let inputs = self.arguments(nest)?;
                    let depth = inputs
                        .iter()
                        .map(|value| value.depth + 1)
                        .max()
                        .unwrap_or(0);
                    return expr(span, ExprForm::RvBarWord(path, cuts, inputs), depth);
                }
                if name.text == "rvBarLength" && self.eat(Token::LBracket)? {
                    let (length, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RvBarLength(length, Box::new(input)), depth);
                }
                if name.text == "rvhface" && self.eat(Token::LBracket)? {
                    let ideal = self.label()?;
                    self.expect(Token::Comma)?;
                    let vertices = self.exponents()?;
                    self.expect(Token::Comma)?;
                    let face = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RvHFace(ideal, vertices, face), 0);
                }
                if name.text == "rvMComponent" && self.eat(Token::LBracket)? {
                    let multidegree = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RvMComponent(multidegree, Box::new(input)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "rvkterm" | "rvmkterm")
                    && self.eat(Token::LBracket)?
                {
                    let ideal = self.label()?;
                    self.expect(Token::Comma)?;
                    let wedge = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    let kind = if name.text == "rvmkterm" {
                        ExprForm::RvMTerm(ideal, wedge, Box::new(input))
                    } else {
                        ExprForm::RvKoszulTerm(ideal, wedge, Box::new(input))
                    };
                    return expr(span, kind, depth);
                }
                if matches!(name.text.as_str(), "rvtorgen" | "rvmtorgen")
                    && self.eat(Token::LBracket)?
                {
                    let grade = self.rv_grade(name.text == "rvmtorgen")?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RvTorGenerator(grade, index, at), 0);
                }
                if name.text == "rvmon" && self.eat(Token::LBracket)? {
                    let selected = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LBrace)?;
                    let powers = self.unsigned_list(Token::RBrace)?;
                    return expr(span, ExprForm::RvMonomial(selected, powers), 0);
                }
                if name.text == "bpVgen" && self.eat(Token::LBracket)? {
                    let spec = self.bp_vpage_spec()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::BpVPageGenerator(spec, index, at), 0);
                }
                if name.text == "bpVClass" && self.eat(Token::LBracket)? {
                    let spec = self.prime_bidegree()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::BpVPageClass(spec, Box::new(input)), depth);
                }
                if let Some(operation) = BpVPageOp::from_name(&name.text) {
                    if self.token() == Token::LParen {
                        let input = self.paren_expr(nest)?;
                        let depth = input.depth + 1;
                        return expr(
                            span,
                            ExprForm::BpVPageOperation(name, operation, Box::new(input)),
                            depth,
                        );
                    }
                }
                if matches!(name.text.as_str(), "bpBockgen" | "bpNovgen")
                    && self.eat(Token::LBracket)?
                {
                    let mode = if name.text == "bpBockgen" {
                        BpPageMode::Bockstein
                    } else {
                        BpPageMode::Novikov
                    };
                    let spec = self.prime_profile_page_degree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::BpPageGenerator(mode, spec, index, at), 0);
                }
                if matches!(name.text.as_str(), "bpBockClass" | "bpNovClass")
                    && self.eat(Token::LBracket)?
                {
                    let mode = if name.text == "bpBockClass" {
                        BpPageMode::Bockstein
                    } else {
                        BpPageMode::Novikov
                    };
                    let specification = self.prime_bidegree()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::BpPageClass(mode, specification, Box::new(input)),
                        depth,
                    );
                }
                if let Some((mode, operation)) = BpPageOp::from_name(&name.text) {
                    if self.token() == Token::LParen {
                        let input = self.paren_expr(nest)?;
                        let depth = input.depth + 1;
                        return expr(
                            span,
                            ExprForm::BpPageOperation(name, mode, operation, Box::new(input)),
                            depth,
                        );
                    }
                }
                if matches!(name.text.as_str(), "bpIExtgen" | "bpThickExtgen")
                    && self.eat(Token::LBracket)?
                {
                    let grade = self.bp_mod_ext_grade(name.text == "bpThickExtgen")?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::BpModExtGen(grade, index, at), 0);
                }
                if matches!(name.text.as_str(), "bpIExtClass" | "bpThickExtClass")
                    && self.eat(Token::LBracket)?
                {
                    let (t, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::BpModExtClass(name.text == "bpThickExtClass", t, Box::new(input)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "bpIExtRep" | "bpThickExtRep")
                    && self.token() == Token::LParen
                {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::BpModExtRep(name.text == "bpThickExtRep", Box::new(input)),
                        depth,
                    );
                }
                if name.text == "bpRaise" && self.eat(Token::LBracket)? {
                    let (height, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::BpRaise(height, Box::new(input)), depth);
                }
                if name.text == "bpTriple" && self.eat(Token::LParen)? {
                    let inputs = self.bp_triple_arguments(nest, span, "bpTriple")?;
                    let depth = inputs.iter().map(|input| input.depth + 1).max().unwrap();
                    return expr(span, ExprForm::BpTriple(inputs), depth);
                }
                if name.text == "bpTripleBase" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::BpTripleBase(Box::new(input)), depth);
                }
                if name.text == "bpHigherSystem" && self.eat(Token::LBracket)? {
                    let (arity, at) = self.number()?;
                    let count = crate::bp_higher::interval_count(arity).map_err(|e| e.at(at))?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LParen)?;
                    let inputs = self.arguments(nest)?;
                    if inputs.len() != count - 1 {
                        return Err(Error::new(span, format!("bpHigherSystem requires {arity} BPExt leaves and {} proper witnesses", count-1-arity)));
                    }
                    let depth = inputs.iter().map(|input| input.depth + 1).max().unwrap();
                    return expr(span, ExprForm::BpHigher(Some(arity), inputs), depth);
                }
                if name.text == "bpChosenHigher" && self.eat(Token::LParen)? {
                    let inputs = self.arguments(nest)?;
                    crate::bp_higher::interval_count(inputs.len()).map_err(|e| e.at(span))?;
                    let depth = inputs.iter().map(|input| input.depth + 1).max().unwrap();
                    return expr(span, ExprForm::BpHigher(None, inputs), depth);
                }
                if name.text == "bpMassey" && self.eat(Token::LParen)? {
                    let arguments = self.arguments(nest)?;
                    if arguments.len() != 5 {
                        return Err(Error::new(
                            span,
                            "bpMassey requires five cochains A, B, C, U, V",
                        ));
                    }
                    let depth = arguments.iter().map(|arg| arg.depth + 1).max().unwrap();
                    return expr(span, ExprForm::BpMassey(arguments), depth);
                }
                if name.text == "bpHeightgen" && self.eat(Token::LBracket)? {
                    let spec = self.bp_height_spec()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::BpHeightGenerator(spec, index, at), 0);
                }
                if name.text == "bpHeightClass" && self.eat(Token::LBracket)? {
                    let (r, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let h = self.signed_wide("a signed BP height", "BP height filtration")?;
                    self.expect(Token::Comma)?;
                    let (t, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    let op = BpHeightOp::Class(r, h, t);
                    return expr(
                        span,
                        ExprForm::BpHeightOperation(op, Box::new(input)),
                        depth,
                    );
                }
                if self.token() == Token::LParen {
                    let operation = match name.text.as_str() {
                        "bpHeightRep" => Some(BpHeightOp::Representative),
                        "bpHeightD" => Some(BpHeightOp::Differential),
                        "bpHeightNext" => Some(BpHeightOp::Next),
                        "bpHeightPrev" => Some(BpHeightOp::Previous),
                        _ => None,
                    };
                    if let Some(op) = operation {
                        let input = self.paren_expr(nest)?;
                        let depth = input.depth + 1;
                        return expr(
                            span,
                            ExprForm::BpHeightOperation(op, Box::new(input)),
                            depth,
                        );
                    }
                }
                if name.text == "bprelgen" && self.eat(Token::LBracket)? {
                    let shape = self.bp_relative_shape()?;
                    self.expect(Token::Comma)?;
                    let (t, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::BpRelGenerator(shape, t, index, at), 0);
                }
                if name.text == "bpRelClass" && self.eat(Token::LBracket)? {
                    let (t, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::BpRelOperation(BpRelOp::Class(t), Box::new(input)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "bpRelRep" | "bpHeightDelta")
                    && self.token() == Token::LParen
                {
                    let operation = if name.text == "bpRelRep" {
                        BpRelOp::Representative
                    } else {
                        BpRelOp::HeightDelta
                    };
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::BpRelOperation(operation, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "bpRel" && self.eat(Token::LBracket)? {
                    let (height, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::BpRelOperation(BpRelOp::Quotient(height), Box::new(input)),
                        depth,
                    );
                }
                if name.text == "bpRelLift" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::BpRelOperation(BpRelOp::Lift, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "bpextgen" && self.eat(Token::LBracket)? {
                    let (spec, t) = self.bp_ext_grade()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::BpExtGenerator(spec, t, index, at), 0);
                }
                if matches!(name.text.as_str(), "bpExtClass" | "bpScale")
                    && self.eat(Token::LBracket)?
                {
                    let operation = if name.text == "bpExtClass" {
                        BpExtOp::Class(self.number()?.0)
                    } else {
                        BpExtOp::Scale(self.bp_coefficient()?)
                    };
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::BpExtOperation(operation, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "bpExtRep" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::BpExtOperation(BpExtOp::Representative, Box::new(input)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "bpint" | "bpv" | "bpt")
                    && self.eat(Token::LBracket)?
                {
                    let (prime, height, arrows) = self.prime_bidegree()?;
                    let shape = BpSpec(prime, height, arrows);
                    self.expect(Token::Comma)?;
                    let kind = if name.text == "bpint" {
                        ExprForm::BpInteger(shape, self.bp_coefficient()?)
                    } else {
                        let slot = if name.text == "bpt" {
                            let (slot, _) = self.number()?;
                            self.expect(Token::Comma)?;
                            Some(slot)
                        } else {
                            None
                        };
                        let (index, _) = self.number()?;
                        ExprForm::BpVariable(shape, slot, index)
                    };
                    self.expect(Token::RBracket)?;
                    return expr(span, kind, 0);
                }
                if matches!(
                    name.text.as_str(),
                    "bpPow"
                        | "bpCoface"
                        | "bpNull"
                        | "bpReduce"
                        | "bpBockP"
                        | "bpDegeneracy"
                        | "bpThick"
                ) && self.eat(Token::LBracket)?
                {
                    let (index, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let operation = match name.text.as_str() {
                        "bpPow" => BpOperation::Power(index),
                        "bpCoface" => BpOperation::Coface(index),
                        "bpReduce" => BpOperation::Reduce(index),
                        "bpBockP" => BpOperation::BockP(index),
                        "bpDegeneracy" => BpOperation::Degeneracy(index),
                        "bpNull" => BpOperation::Null(index),
                        "bpThick" => BpOperation::Thick(index),
                        _ => unreachable!(),
                    };
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::BpOperation(operation, Box::new(input)),
                        depth,
                    );
                }
                if matches!(
                    name.text.as_str(),
                    "bpNeg"
                        | "bpD"
                        | "bpLift"
                        | "bpDelta"
                        | "bpNormalize"
                        | "bpThickLift"
                        | "bpThickDelta"
                ) && self.token() == Token::LParen
                {
                    let operation = match name.text.as_str() {
                        "bpNeg" => BpOperation::Negate,
                        "bpD" => BpOperation::Differential,
                        "bpLift" => BpOperation::Lift,
                        "bpDelta" => BpOperation::Delta,
                        "bpNormalize" => BpOperation::Normalize,
                        "bpThickLift" => BpOperation::ThickLift,
                        "bpThickDelta" => BpOperation::ThickDelta,
                        _ => unreachable!(),
                    };
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::BpOperation(operation, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvMinor" && self.eat(Token::LBracket)? {
                    let spec = self.rv_matrix(nest)?;
                    let rows = self.rv_coordinates()?;
                    self.expect(Token::Comma)?;
                    let columns = self.rv_coordinates()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RvMinor(spec, rows, columns), 0);
                }
                if name.text == "rvHasse" && self.eat(Token::LBracket)? {
                    let orders = self.unsigned_list(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RvHasse(orders, Box::new(input)), depth);
                }
                if name.text == "rvVanish" && self.eat(Token::LBracket)? {
                    let (modulus, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let point = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RvVanish(modulus, point, Box::new(input)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "rvNormal" | "rvInIdeal" | "rvWitness")
                    && self.eat(Token::LBracket)?
                {
                    let ideal = self.label()?;
                    let op = match name.text.as_str() {
                        "rvNormal" => RvIdealOp::Normal,
                        "rvInIdeal" => RvIdealOp::Contains,
                        "rvWitness" => {
                            self.expect(Token::Comma)?;
                            RvIdealOp::Witness(self.number()?.0)
                        }
                        _ => unreachable!(),
                    };
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RvIdeal(op, ideal, Box::new(input)), depth);
                }
                if name.text == "eoType" && self.eat(Token::LBrace)? {
                    let sequence = self.unsigned_list(Token::RBrace)?;
                    return expr(span, ExprForm::EOType(sequence), 0);
                }
                if matches!(name.text.as_str(), "rvASEOType" | "rvASOddEOType")
                    && self.eat(Token::LBracket)?
                {
                    let cover = self.label()?;
                    self.expect(Token::RBracket)?;
                    let form = if name.text == "rvASEOType" {
                        ExprForm::EOCover(cover)
                    } else {
                        ExprForm::EOOddCover(cover)
                    };
                    return expr(span, form, 0);
                }
                if name.text == "rvASDeRhamGen" && self.eat(Token::LBracket)? {
                    let cover = self.label()?;
                    self.expect(Token::Comma)?;
                    let modulus = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::AsDeRhamGen(cover, modulus, index, at), 0);
                }
                if name.text == "rvFieldElt" && self.eat(Token::LBracket)? {
                    let (modulus, scalar) = self.bidegree()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::BinaryScalar(modulus, scalar), 0);
                }
                if name.text == "rvASDeRhamScale" && self.eat(Token::LBracket)? {
                    let scalar = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    let operation = crate::rv_as_derham::Operation::Scale(scalar);
                    return expr(
                        span,
                        ExprForm::AsDeRhamOperation(operation, Box::new(input)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "rvASF" | "rvASV" | "rvASDeRhamIncl")
                    && self.token() == Token::LParen
                {
                    use crate::rv_as_derham::Operation;
                    let operation = match name.text.as_str() {
                        "rvASF" => Operation::Frobenius,
                        "rvASV" => Operation::Verschiebung,
                        "rvASDeRhamIncl" => Operation::Include,
                        _ => unreachable!(),
                    };
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::AsDeRhamOperation(operation, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvASOddTategen" && self.eat(Token::LBracket)? {
                    let cover = self.label()?;
                    self.expect(Token::Comma)?;
                    let degree = self.signed_number("Tate degree")?;
                    self.expect(Token::Comma)?;
                    let index = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::AsOddTateGen(cover, degree, index), 0);
                }
                if matches!(name.text.as_str(), "rvASOddTateD" | "rvASOddTateClass")
                    && self.eat(Token::LBracket)?
                {
                    let degree = self.signed_number("Tate degree")?;
                    self.expect(Token::RBracket)?;
                    let operation = if name.text == "rvASOddTateD" {
                        AsOddTateOp::Differential(degree)
                    } else {
                        AsOddTateOp::Class(degree)
                    };
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::AsOddTateOperation(operation, Box::new(input)),
                        depth,
                    );
                }
                if matches!(
                    name.text.as_str(),
                    "rvASOddTateRep" | "rvASOddTateF" | "rvASOddTateV"
                ) && self.token() == Token::LParen
                {
                    let operation = match name.text.as_str() {
                        "rvASOddTateRep" => AsOddTateOp::Representative,
                        "rvASOddTateF" => AsOddTateOp::Frobenius,
                        "rvASOddTateV" => AsOddTateOp::Verschiebung,
                        _ => unreachable!(),
                    };
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::AsOddTateOperation(operation, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvASOddDeRhamGen" && self.eat(Token::LBracket)? {
                    let cover = self.label()?;
                    self.expect(Token::Comma)?;
                    let index = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::AsOddDeRhamGen(cover, index), 0);
                }
                if name.text == "rvPrimeElt" && self.eat(Token::LBracket)? {
                    let prime = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let coefficient =
                        self.signed_wide("a prime field coefficient", "prime field coefficient")?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::PrimeScalar(prime, coefficient), 0);
                }
                if name.text == "rvASOddDeRhamTranslate" && self.eat(Token::LBracket)? {
                    let scalar = self.signed_wide(
                        "a signed translation",
                        "odd Artin–Schreier de Rham translation",
                    )?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    let operation = AsOddDeRhamOp::Translate(scalar);
                    return expr(
                        span,
                        ExprForm::AsOddDeRhamOperation(operation, Box::new(input)),
                        depth,
                    );
                }
                if matches!(
                    name.text.as_str(),
                    "rvASOddF" | "rvASOddV" | "rvASOddDeRhamIncl"
                ) && self.token() == Token::LParen
                {
                    let operation = match name.text.as_str() {
                        "rvASOddF" => AsOddDeRhamOp::Frobenius,
                        "rvASOddV" => AsOddDeRhamOp::Verschiebung,
                        "rvASOddDeRhamIncl" => AsOddDeRhamOp::Include,
                        _ => unreachable!(),
                    };
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::AsOddDeRhamOperation(operation, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvASOddOmega" && self.eat(Token::LBracket)? {
                    let cover = self.label()?;
                    self.expect(Token::Comma)?;
                    let (b, r) = self.bidegree()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::AsOddOmega(cover, b, r), 0);
                }
                if matches!(
                    name.text.as_str(),
                    "rvASOddCartierIter" | "rvASOddTranslate"
                ) && self.eat(Token::LBracket)?
                {
                    let operation = if name.text == "rvASOddTranslate" {
                        AsOddDiffOp::Translate(
                            self.signed_wide("a signed translation", "Artin–Schreier translation")?,
                        )
                    } else {
                        AsOddDiffOp::Iterate(self.number()?.0)
                    };
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::AsOddDiffOperation(operation, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvASOddCartier" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::AsOddDiffOperation(AsOddDiffOp::Cartier, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvASOddCoeff" && self.eat(Token::LBracket)? {
                    let cover = self.label()?;
                    self.expect(Token::Comma)?;
                    let (character, coefficient) = self.bidegree()?;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::AsOddCoefficient(cover, character, coefficient),
                        0,
                    );
                }
                if matches!(name.text.as_str(), "rvASOmega" | "rvWittOmega")
                    && self.eat(Token::LBracket)?
                {
                    let family = if name.text == "rvWittOmega" {
                        DiffFamily::Witt
                    } else {
                        DiffFamily::ArtinSchreier
                    };
                    let cover = self.label()?;
                    self.expect(Token::Comma)?;
                    let modulus = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let index = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RvOmega(family, cover, modulus, index), 0);
                }
                if matches!(name.text.as_str(), "rvASDiffScale" | "rvWittDiffScale")
                    && self.eat(Token::LBracket)?
                {
                    let family = if name.text == "rvWittDiffScale" {
                        DiffFamily::Witt
                    } else {
                        DiffFamily::ArtinSchreier
                    };
                    let scalar = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RvDiffScale(family, scalar, Box::new(input)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "rvASCartier" | "rvWittCartier")
                    && self.token() == Token::LParen
                {
                    let family = if name.text == "rvWittCartier" {
                        DiffFamily::Witt
                    } else {
                        DiffFamily::ArtinSchreier
                    };
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RvCartier(family, Box::new(input)), depth);
                }
                if name.text == "rvFrobWronskian" && self.eat(Token::LBracket)? {
                    let curve = self.label()?;
                    self.expect(Token::Comma)?;
                    let (exponent, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RvFrobeniusWronskian(curve, exponent), 0);
                }
                if name.text == "rvOscMinor" && self.eat(Token::LBracket)? {
                    let curve = self.label()?;
                    self.expect(Token::Comma)?;
                    let (k, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let columns = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RvOscMinor(curve, k, columns), 0);
                }
                if name.text == "rvGenerator" && self.eat(Token::LBracket)? {
                    let ideal = self.label()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RvGenerator(ideal, index, at), 0);
                }
                if name.text == "rvInMinors" && self.eat(Token::LBracket)? {
                    let spec = self.rv_matrix(nest)?;
                    let (size, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RvInMinors(spec, size, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "milnor" && self.eat(Token::LBrace)? {
                    let exponents = self.exponents()?;
                    self.expect(Token::RBrace)?;
                    return expr(span, ExprForm::Milnor(exponents), 0);
                }
                if name.text == "cword" && self.eat(Token::LBrace)? {
                    let mut factors = Vec::new();
                    while self.token() != Token::RBrace {
                        let at = self.span();
                        factors.push((self.exponents()?, at));
                        if !self.eat(Token::Comma)? {
                            break;
                        }
                    }
                    self.expect(Token::RBrace)?;
                    return expr(span, ExprForm::CobarWord(factors), 0);
                }
                if name.text == "tauShift" && self.eat(Token::LBracket)? {
                    let (power, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::TauShift(power, Box::new(input)), depth);
                }
                if matches!(name.text.as_str(), "rhoMin" | "rhoMinCoef")
                    && self.eat(Token::LBracket)?
                {
                    let (f, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    let op = if name.text == "rhoMin" {
                        RhoMinOp::Inject(f)
                    } else {
                        RhoMinOp::Coefficient(f)
                    };
                    return expr(span, ExprForm::RhoMinOperation(op, Box::new(input)), depth);
                }
                let rho_minimal = match name.text.as_str() {
                    "rhoMinD" => Some(RhoMinOp::Differential),
                    "rhoMinI" => Some(RhoMinOp::Lift),
                    "rhoMinP" => Some(RhoMinOp::Project),
                    "rhoMinH" => Some(RhoMinOp::Homotopy),
                    _ => None,
                };
                if let Some(op) = rho_minimal.filter(|_| self.token() == Token::LParen) {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RhoMinOperation(op, Box::new(input)), depth);
                }
                if name.text == "rhoShift" && self.eat(Token::LBracket)? {
                    let power = self.signed_wide("a signed rho exponent", "rho exponent")?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RhoShift(power, Box::new(input)), depth);
                }
                if name.text == "rholocword" && self.eat(Token::LBracket)? {
                    let t = self.signed_wide(
                        "a signed topological degree",
                        "rho-localized topological degree",
                    )?;
                    self.expect(Token::Comma)?;
                    let (vertex, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LBrace)?;
                    let mut factors = Vec::new();
                    while self.token() != Token::RBrace {
                        let at = self.span();
                        self.expect(Token::LParen)?;
                        let (u, _) = self.number()?;
                        self.expect(Token::Comma)?;
                        let xi = self.exponents()?;
                        self.expect(Token::RParen)?;
                        factors.push((u, xi, at));
                        if !self.eat(Token::Comma)? {
                            break;
                        }
                    }
                    self.expect(Token::RBrace)?;
                    return expr(span, ExprForm::RhoLocWord(t, vertex, factors), 0);
                }
                if name.text == "rhoNull" && self.eat(Token::LBracket)? {
                    let (power, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RhoNull(power, Box::new(input)), depth);
                }
                if matches!(name.text.as_str(), "rmMassey" | "rmTriple")
                    && self.eat(Token::LParen)?
                {
                    let inputs = self.exact_arguments(nest, span, &name.text, "three inputs")?;
                    let depth = inputs.iter().map(|x| x.depth).max().unwrap_or(0) + 1;
                    return expr(
                        span,
                        ExprForm::RealMassey(name.text == "rmTriple", inputs),
                        depth,
                    );
                }
                if name.text == "rmMasseyNull" && self.eat(Token::LParen)? {
                    let inputs = self.exact_arguments(nest, span, &name.text, "two inputs")?;
                    let depth = inputs.iter().map(|x| x.depth).max().unwrap_or(0) + 1;
                    return expr(span, ExprForm::RealMasseyNull(inputs), depth);
                }
                if name.text == "rmTripleBase" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RealTripleClass(None, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rmTripleDir" && self.eat(Token::LBracket)? {
                    let (index, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RealTripleClass(Some(index), Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rhoSection" && self.eat(Token::LBracket)? {
                    let t = self.signed_wide(
                        "a signed topological degree",
                        "rho-localized topological degree",
                    )?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RhoLocOperation(RhoLocOp::Section(t), Box::new(input)),
                        depth,
                    );
                }
                if self.token() == Token::LParen {
                    let operation = match name.text.as_str() {
                        "rhoLocalize" => Some(RhoLocOp::Localize),
                        "rhoProject" => Some(RhoLocOp::Project),
                        "rhoHomotopy" => Some(RhoLocOp::Homotopy),
                        _ => None,
                    };
                    if let Some(operation) = operation {
                        let input = self.paren_expr(nest)?;
                        let depth = input.depth + 1;
                        return expr(
                            span,
                            ExprForm::RhoLocOperation(operation, Box::new(input)),
                            depth,
                        );
                    }
                }
                if name.text == "rhobockgen" && self.eat(Token::LBracket)? {
                    let grade = self.rho_grade()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RhoGenerator(grade, index, at), 0);
                }
                if name.text == "rhoClass" && self.eat(Token::LBracket)? {
                    let (r, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let (f, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RhoOperation(RhoOp::Class(r, f), Box::new(input)),
                        depth,
                    );
                }
                if self.token() == Token::LParen {
                    let operation = match name.text.as_str() {
                        "rhoRep" => Some(RhoOp::Representative),
                        "rhoD" => Some(RhoOp::Differential),
                        "rhoNext" => Some(RhoOp::Next),
                        "rhoPrev" => Some(RhoOp::Previous),
                        _ => None,
                    };
                    if let Some(operation) = operation {
                        let input = self.paren_expr(nest)?;
                        let depth = input.depth + 1;
                        return expr(
                            span,
                            ExprForm::RhoOperation(operation, Box::new(input)),
                            depth,
                        );
                    }
                }
                if name.text == "realComplex" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RealComplex(Box::new(input)), depth);
                }
                if name.text == "realize" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::Realize(Box::new(input)), depth);
                }
                if matches!(name.text.as_str(), "pword" | "mcword" | "rmcword")
                    && self.eat(Token::LBracket)?
                {
                    let (coefficient, _) = self.number()?;
                    let second = if name.text == "rmcword" {
                        self.expect(Token::Comma)?;
                        self.number()?.0
                    } else {
                        0
                    };
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LBrace)?;
                    let mut factors = Vec::new();
                    while self.token() != Token::RBrace {
                        let at = self.span();
                        self.expect(Token::LParen)?;
                        let exponents = self.exponents()?;
                        self.expect(Token::Comma)?;
                        let taus = self.exponents()?;
                        self.expect(Token::RParen)?;
                        factors.push((exponents, taus, at));
                        if !self.eat(Token::Comma)? {
                            break;
                        }
                    }
                    self.expect(Token::RBrace)?;
                    let form = match name.text.as_str() {
                        "mcword" => ExprForm::MCobarWord(coefficient, factors),
                        "rmcword" => ExprForm::RMCobarWord(coefficient, second, factors),
                        _ => ExprForm::PCobarWord(coefficient, factors),
                    };
                    return expr(span, form, 0);
                }
                if name.text == "rvScale" && self.eat(Token::LBracket)? {
                    let (scalar, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::RVScale(scalar, Box::new(value)), depth);
                }
                if name.text == "shCellTerm" && self.eat(Token::LBracket)? {
                    let context = self.label()?;
                    self.expect(Token::Comma)?;
                    let q = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let cell = self.label()?;
                    self.expect(Token::Comma)?;
                    let index = self.number()?.0;
                    let coefficient = if self.eat(Token::Comma)? {
                        self.signed_number("signed sheaf coefficient")?
                    } else {
                        1
                    };
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::SheafCellTerm(context, q, cell, index, coefficient),
                        0,
                    );
                }
                if name.text == "shHomTerm" && self.eat(Token::LBracket)? {
                    let source = self.label()?;
                    self.expect(Token::Comma)?;
                    let target = self.label()?;
                    self.expect(Token::Comma)?;
                    let chain = self.sheaf_chain()?;
                    self.expect(Token::Comma)?;
                    let (column, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let (row, _) = self.number()?;
                    let coefficient = if self.eat(Token::Comma)? {
                        self.signed_number("signed sheaf coefficient")?
                    } else {
                        1
                    };
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::SheafHomTerm(source, target, chain, column, row, coefficient),
                        0,
                    );
                }
                if name.text == "shExtCompose" && self.token() == Token::LParen {
                    let arguments = self.paren_pair(nest)?;
                    let depth = 1 + arguments[0].depth.max(arguments[1].depth);
                    return expr(span, ExprForm::SheafCompose(Box::new(arguments)), depth);
                }
                if name.text == "shExtGen" && self.eat(Token::LBracket)? {
                    let (source, target, degree) = self.sheaf_hom_grade()?;
                    self.expect(Token::Comma)?;
                    let (index, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::SheafExtGenerator(source, target, degree, index),
                        0,
                    );
                }
                if matches!(name.text.as_str(), "shMapC" | "shHomId")
                    && self.eat(Token::LBracket)?
                {
                    let operand = self.label()?;
                    self.expect(Token::RBracket)?;
                    let form = if name.text == "shMapC" {
                        ExprForm::SheafHomMap(operand)
                    } else {
                        ExprForm::SheafHomIdentity(operand)
                    };
                    return expr(span, form, 0);
                }
                if name.text == "shTerm" && self.eat(Token::LBracket)? {
                    let sheaf = self.label()?;
                    self.expect(Token::Comma)?;
                    let chain = self.sheaf_chain()?;
                    self.expect(Token::Comma)?;
                    let index = self.number()?.0;
                    let coefficient = if self.eat(Token::Comma)? {
                        self.signed_number("signed sheaf coefficient")?
                    } else {
                        1
                    };
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::SheafTerm(sheaf, chain, index, coefficient),
                        0,
                    );
                }
                if name.text == "shlocalgen" && self.eat(Token::LBracket)? {
                    let (complex, point, mode) = self.sheaf_local_site()?;
                    self.expect(Token::Comma)?;
                    let degree = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let index = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::SheafLocalGenerator(complex, point, mode, degree, index),
                        0,
                    );
                }
                if name.text == "shLocalClass" && self.eat(Token::LBracket)? {
                    let (complex, point, mode) = self.sheaf_local_site()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::SheafLocalClass(complex, point, mode, Box::new(input)),
                        depth,
                    );
                }
                if self.token() == Token::LParen {
                    let operation = match name.text.as_str() {
                        "shLocalRep" => Some(SupportOp::Rep),
                        "shLocalForget" => Some(SupportOp::Forget),
                        "shLocalRestrict" => Some(SupportOp::Restrict),
                        "shLocalConnect" => Some(SupportOp::Connect),
                        _ => None,
                    };
                    if let Some(operation) = operation {
                        let input = self.paren_expr(nest)?;
                        let depth = input.depth + 1;
                        return expr(
                            span,
                            ExprForm::SheafLocalOperation(operation, Box::new(input)),
                            depth,
                        );
                    }
                }
                if name.text == "shLerayIncl" && self.eat(Token::LBracket)? {
                    let complex = self.label()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::SheafLerayIncl(complex, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "shApply" && self.eat(Token::LBracket)? {
                    let map = self.label()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::SheafApply(map, Box::new(input)), depth);
                }
                if name.text == "shTotEmbed" && self.eat(Token::LBracket)? {
                    let complex = self.label()?;
                    self.expect(Token::Comma)?;
                    let q = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::SheafEmbed(complex, q, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "shTotPart" && self.eat(Token::LBracket)? {
                    let q = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::SheafPart(q, Box::new(input)), depth);
                }
                if name.text == "shegen" && self.eat(Token::LBracket)? {
                    let complex = self.label()?;
                    self.expect(Token::Comma)?;
                    let grade = self.sheaf_page_grade()?;
                    self.expect(Token::Comma)?;
                    let index = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::SheafPageGenerator(complex, grade, index), 0);
                }
                if matches!(name.text.as_str(), "shEClass" | "shELiftTo")
                    && self.eat(Token::LBracket)?
                {
                    let operation = if name.text == "shEClass" {
                        SheafPageLift::Class
                    } else {
                        SheafPageLift::Correct
                    };
                    let r = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let p = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::SheafPageLift(operation, r, p, Box::new(input)),
                        depth,
                    );
                }
                if self.token() == Token::LParen {
                    let operation = match name.text.as_str() {
                        "shED" => Some(SheafPageOp::Differential),
                        "shERep" => Some(SheafPageOp::Representative),
                        "shENext" => Some(SheafPageOp::Next),
                        "shEPrev" => Some(SheafPageOp::Previous),
                        _ => None,
                    };
                    if let Some(operation) = operation {
                        let input = self.paren_expr(nest)?;
                        let depth = input.depth + 1;
                        return expr(
                            span,
                            ExprForm::SheafPageOperation(operation, Box::new(input)),
                            depth,
                        );
                    }
                }
                if name.text == "shhypergen" && self.eat(Token::LBracket)? {
                    let complex = self.label()?;
                    self.expect(Token::Comma)?;
                    let degree = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let index = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::SheafHyperGenerator(complex, degree, index),
                        0,
                    );
                }
                if name.text == "shgen" && self.eat(Token::LBracket)? {
                    let sheaf = self.label()?;
                    self.expect(Token::Comma)?;
                    let degree = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let index = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::SheafGenerator(sheaf, degree, index), 0);
                }
                if let Some(operation) = match name.text.as_str() {
                    "shD" if self.token() == Token::LParen => Some(SheafOp::Differential),
                    "shClass" if self.token() == Token::LParen => Some(SheafOp::Class),
                    "shRep" if self.token() == Token::LParen => Some(SheafOp::Representative),
                    "shScale" if self.eat(Token::LBracket)? => {
                        let coefficient = self.signed_number("signed sheaf coefficient")?;
                        self.expect(Token::RBracket)?;
                        Some(SheafOp::Scale(coefficient))
                    }
                    _ => None,
                } {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::SheafOperation(operation, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "pscale" && self.eat(Token::LBracket)? {
                    let (scalar, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::PScale(scalar, Box::new(value)), depth);
                }
                if name.text == "pextgen" && self.eat(Token::LBracket)? {
                    let (p, s, t) = self.prime_bidegree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::PExtGenerator(p, s, t, index, at), 0);
                }
                if name.text == "pextangen" && self.eat(Token::LBracket)? {
                    let (p, n, s, t) = self.prime_profile_bidegree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::PExtAnGenerator(p, n, s, t, index, at), 0);
                }
                if matches!(name.text.as_str(), "rmbarGen" | "rmbarRep" | "rmbarNull")
                    && self.eat(Token::LBracket)?
                {
                    let (s, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let k = self.signed_wide("a signed degree t-w", "rho bar degree")?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let op = match name.text.as_str() {
                        "rmbarGen" => RealBarOp::Generator,
                        "rmbarRep" => RealBarOp::Representative,
                        _ => RealBarOp::Witness,
                    };
                    return expr(span, ExprForm::RealBar(op, s, k, index, at), 0);
                }
                if name.text == "rmextgen" && self.eat(Token::LBracket)? {
                    let grade = self.motivic_grade()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RMExtGenerator(grade, index, at), 0);
                }
                if name.text == "mextgen" && self.eat(Token::LBracket)? {
                    let grade = self.motivic_grade()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::MExtGenerator(grade, index, at), 0);
                }
                if name.text == "extgen" && self.eat(Token::LBracket)? {
                    let (s, t) = self.bidegree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::ExtGenerator(s, t, index, at), 0);
                }
                if name.text == "extangen" && self.eat(Token::LBracket)? {
                    let (n, s, t) = self.profile_bidegree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::ExtAnGenerator(n, s, t, index, at), 0);
                }
                if name.text == "resA" && self.eat(Token::LBracket)? {
                    let (n, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::RestrictAlgebra(n, Box::new(value)), depth);
                }
                if name.text == "cCup" && self.eat(Token::LBracket)? {
                    let (i, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LParen)?;
                    let left = self.expr(nest + 1)?;
                    self.expect(Token::Comma)?;
                    let right = self.expr(nest + 1)?;
                    self.expect(Token::RParen)?;
                    let depth = left.depth.max(right.depth) + 1;
                    return expr(
                        span,
                        ExprForm::CobarCup(i, Box::new(left), Box::new(right)),
                        depth,
                    );
                }
                if name.text == "extSq" && self.eat(Token::LBracket)? {
                    let (k, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::ExtSquare(k, Box::new(value)), depth);
                }
                if name.text == "cNakamura" && self.eat(Token::LBracket)? {
                    let (k, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::Nakamura(k, Box::new(value)), depth);
                }
                if name.text == "mayLiftTo" && self.eat(Token::LBracket)? {
                    let (r, f) = self.bidegree()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::MayLiftTo(r, f, Box::new(value)), depth);
                }
                if name.text == "resHigher" && self.eat(Token::LBracket)? {
                    let (n, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LParen)?;
                    let inputs = self.arguments(nest)?;
                    let depth = inputs
                        .iter()
                        .map(|input| input.depth + 1)
                        .max()
                        .unwrap_or(0);
                    return expr(span, ExprForm::HigherRestriction(n, inputs), depth);
                }
                if name.text == "xextgen" && self.eat(Token::LBracket)? {
                    let (base, s, t) = self.space_bidegree(nest)?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::XExtGenerator(base, s, t, index, at), 0);
                }
                if name.text == "mapextgen" && self.eat(Token::LBracket)? {
                    let (source, target, s, t) = self.map_bidegree(nest)?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::MapExtGenerator(source, target, s, t, index, at),
                        0,
                    );
                }
                if name.text == "mword" && self.eat(Token::LBracket)? {
                    let source = self.map_generator(nest)?;
                    self.expect(Token::Comma)?;
                    let target = self.map_generator(nest)?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::MapWord(source, target, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "eid" && self.eat(Token::LBracket)? {
                    let base = self.ty(nest + 1)?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::MapExtIdentity(base), 0);
                }
                if name.text == "jetterm" && self.eat(Token::LBracket)? {
                    let source = self.map_generator(nest)?;
                    self.expect(Token::Comma)?;
                    let target = self.map_generator(nest)?;
                    self.expect(Token::Comma)?;
                    let path = self.label()?;
                    self.expect(Token::Comma)?;
                    let (power, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::JetTerm(JetTerm {
                            source,
                            target,
                            path,
                            power,
                        }),
                        0,
                    );
                }
                if name.text == "jetgen" && self.eat(Token::LBracket)? {
                    let spec = self.jet_page(nest)?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::JetGenerator(spec, index, at), 0);
                }
                if name.text == "jetClass" && self.eat(Token::LBracket)? {
                    let (r, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let (f, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::JetOperation(JetOp::Class(r, f), Box::new(input)),
                        depth,
                    );
                }
                if self.token() == Token::LParen {
                    let operation = match name.text.as_str() {
                        "jetD" => Some(JetOp::Differential),
                        "jetRep" => Some(JetOp::Representative),
                        "jetNext" => Some(JetOp::Next),
                        "jetPrev" => Some(JetOp::Prev),
                        _ => None,
                    };
                    if let Some(operation) = operation {
                        let input = self.paren_expr(nest)?;
                        let depth = input.depth + 1;
                        return expr(
                            span,
                            ExprForm::JetOperation(operation, Box::new(input)),
                            depth,
                        );
                    }
                }
                if name.text == "qterm" && self.eat(Token::LBracket)? {
                    let source = self.map_generator(nest)?;
                    self.expect(Token::Comma)?;
                    let target = self.map_generator(nest)?;
                    self.expect(Token::Comma)?;
                    let selected = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LBrace)?;
                    let mut powers = Vec::new();
                    while self.token() != Token::RBrace {
                        powers.push(self.number()?.0);
                        if !self.eat(Token::Comma)? {
                            break;
                        }
                    }
                    self.expect(Token::RBrace)?;
                    return expr(
                        span,
                        ExprForm::QTerm(QTerm {
                            source,
                            target,
                            selected,
                            powers,
                        }),
                        0,
                    );
                }
                if name.text == "qid" && self.eat(Token::LBracket)? {
                    let base = self.ty(nest + 1)?;
                    self.expect(Token::Comma)?;
                    let selected = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::QIdentity(base, selected), 0);
                }
                if name.text == "rvgen" && self.eat(Token::LBracket)? {
                    let spec = self.rank_variety_spec(nest, true)?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RVGenerator(spec, index, at), 0);
                }
                if name.text == "rvId" && self.eat(Token::LBracket)? {
                    let base = self.ty(nest + 1)?;
                    self.expect(Token::Comma)?;
                    let selected = self.exponents()?;
                    self.expect(Token::Comma)?;
                    let (modulus, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let point = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::RVIdentity(base, selected, modulus, point),
                        0,
                    );
                }
                if matches!(
                    name.text.as_str(),
                    "rvFrob" | "rvExtend" | "rvDescend" | "rvTrace"
                ) && self.eat(Token::LBracket)?
                {
                    use crate::rv_field::Operation;
                    let (parameter, _) = self.number()?;
                    let operation = if name.text == "rvFrob" {
                        Operation::Frobenius(parameter)
                    } else {
                        self.expect(Token::Comma)?;
                        let (root, _) = self.number()?;
                        match name.text.as_str() {
                            "rvExtend" => Operation::Extend(parameter, root),
                            "rvDescend" => Operation::Descend(parameter, root),
                            "rvTrace" => Operation::Trace(parameter, root),
                            _ => unreachable!(),
                        }
                    };
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RVField(operation, Box::new(input)), depth);
                }
                if name.text == "rvInverse" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RVInverse(Box::new(input)), depth);
                }
                if matches!(
                    name.text.as_str(),
                    "rvHilbert90" | "rvTwistedDescend" | "rvTwistedExtend"
                ) && self.eat(Token::LBracket)?
                {
                    let (modulus, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let (root, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    if name.text == "rvHilbert90" {
                        let input = self.paren_expr(nest)?;
                        let depth = input.depth + 1;
                        return expr(
                            span,
                            ExprForm::RVHilbert90(modulus, root, Box::new(input)),
                            depth,
                        );
                    }
                    let operation = if name.text == "rvTwistedDescend" {
                        crate::rv_descent::Twist::Descend
                    } else {
                        crate::rv_descent::Twist::Extend
                    };
                    self.expect(Token::LParen)?;
                    let left = self.expr(nest + 1)?;
                    self.expect(Token::Comma)?;
                    let input = self.expr(nest + 1)?;
                    self.expect(Token::Comma)?;
                    let right = self.expr(nest + 1)?;
                    self.expect(Token::RParen)?;
                    let depth = 1 + left.depth.max(input.depth).max(right.depth);
                    return expr(
                        span,
                        ExprForm::RVTwisted(
                            operation,
                            modulus,
                            root,
                            Box::new([left, input, right]),
                        ),
                        depth,
                    );
                }
                if matches!(
                    name.text.as_str(),
                    "rvFrobNorm" | "rvShintani" | "rvShintaniFrame"
                ) && self.eat(Token::LBracket)?
                {
                    use crate::rv_shintani::Operation;
                    let operation = match name.text.as_str() {
                        "rvFrobNorm" => Operation::Norm,
                        "rvShintani" => Operation::Companion,
                        "rvShintaniFrame" => Operation::Frame,
                        _ => unreachable!(),
                    };
                    let (modulus, root) = self.embedding_parameters()?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RVShintani(operation, modulus, root, Box::new(input)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "rvInIwADLV" | "rvInIwADLVClosure")
                    && self.eat(Token::LBracket)?
                {
                    let boundary = if name.text == "rvInIwADLV" {
                        crate::rv_iwahori::Boundary::Exact
                    } else {
                        crate::rv_iwahori::Boundary::Closure
                    };
                    let (modulus, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let (root, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let window = self.affine_window()?;
                    self.expect(Token::RBracket)?;
                    let arguments = self.paren_pair(nest)?;
                    let depth = 1 + arguments[0].depth.max(arguments[1].depth);
                    return expr(
                        span,
                        ExprForm::RVIwAdlv(boundary, modulus, root, window, Box::new(arguments)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "rvInADLV" | "rvInADLVClosure")
                    && self.eat(Token::LBracket)?
                {
                    let boundary = if name.text == "rvInADLV" {
                        crate::rv_hodge::Boundary::Exact
                    } else {
                        crate::rv_hodge::Boundary::Closure
                    };
                    let (modulus, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let (root, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let mu = self.coweight()?;
                    self.expect(Token::RBracket)?;
                    let arguments = self.paren_pair(nest)?;
                    let depth = 1 + arguments[0].depth.max(arguments[1].depth);
                    return expr(
                        span,
                        ExprForm::RVAdlv(boundary, modulus, root, mu, Box::new(arguments)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "rvHeckeT" | "rvHeckeKL")
                    && self.eat(Token::LBracket)?
                {
                    let canonical = name.text == "rvHeckeKL";
                    let window = self.affine_window()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::HeckeBasis(canonical, window), 0);
                }
                if matches!(name.text.as_str(), "rvHeckeId" | "rvHeckeQ")
                    && self.eat(Token::LBracket)?
                {
                    let degree = usize::from(name.text == "rvHeckeQ");
                    let (rank, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::HeckeScalar(rank, degree), 0);
                }
                if name.text == "neg" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::HeckeNegate(Box::new(input)), depth);
                }
                if name.text == "rvHeckeMul" && self.token() == Token::LParen {
                    let arguments = self.paren_pair(nest)?;
                    let depth = 1 + arguments[0].depth.max(arguments[1].depth);
                    return expr(span, ExprForm::HeckeProduct(Box::new(arguments)), depth);
                }
                if matches!(name.text.as_str(), "rvHive" | "rvHiveWitness")
                    && self.eat(Token::LBracket)?
                {
                    let weights = self.hive_spec()?;
                    self.expect(Token::RBracket)?;
                    let interior = if name.text == "rvHive" {
                        self.expect(Token::LBrace)?;
                        Some(self.signed_list(
                            Token::RBrace,
                            "a signed hive vertex",
                            "hive vertex",
                        )?)
                    } else {
                        None
                    };
                    return expr(span, ExprForm::Hive(weights, interior), 0);
                }
                if name.text == "rvFusionBasis" && self.eat(Token::LBracket)? {
                    let (level, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let weight = self.signed_tuple("a signed fusion weight", "fusion weight")?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::FusionBasis(level, weight), 0);
                }
                if matches!(name.text.as_str(), "rvFusionOne" | "rvFusionHandle")
                    && self.eat(Token::LBracket)?
                {
                    let (rank, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let (level, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let form = if name.text == "rvFusionHandle" {
                        ExprForm::FusionHandle(rank, level)
                    } else {
                        ExprForm::FusionOne(rank, level)
                    };
                    return expr(span, form, 0);
                }
                if matches!(
                    name.text.as_str(),
                    "rvWzwBasis" | "rvWzwSEntry" | "rvWzwTEntry"
                ) && self.eat(Token::LBracket)?
                {
                    let (rank, level) = self.fusion_grade()?;
                    self.expect(Token::Comma)?;
                    let weight = self.signed_tuple("a signed WZW weight", "WZW weight")?;
                    let form = match name.text.as_str() {
                        "rvWzwSEntry" => {
                            self.expect(Token::Comma)?;
                            let right = self.signed_tuple("a signed WZW weight", "WZW weight")?;
                            ExprForm::WzwSEntry(rank, level, weight, right)
                        }
                        "rvWzwTEntry" => {
                            let power = if self.eat(Token::Comma)? {
                                self.signed_wide("a signed twist power", "twist power")?
                            } else {
                                1
                            };
                            ExprForm::WzwTEntry(rank, level, weight, power)
                        }
                        _ => ExprForm::WzwBasis(rank, level, weight),
                    };
                    self.expect(Token::RBracket)?;
                    return expr(span, form, 0);
                }
                if name.text == "rvLensTV" && self.eat(Token::LBracket)? {
                    let (rank, level) = self.fusion_grade()?;
                    self.expect(Token::Comma)?;
                    let (p, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let q = self.signed_wide("a signed lens parameter", "lens parameter")?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::LensTv(rank, level, p, q), 0);
                }
                if name.text == "rvWzwTV" && self.eat(Token::LBracket)? {
                    let (rank, level) = self.fusion_grade()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LBrace)?;
                    let mut steps = Vec::new();
                    while self.token() != Token::RBrace {
                        use crate::rv_lens::{Step, MAX_WORD};
                        if steps.len() == MAX_WORD {
                            return Err(Error::new(
                                self.span(),
                                "WZW gluing word exceeds the limit of 256 letters",
                            ));
                        }
                        let step = match self.token() {
                            Token::Word("S") => {
                                self.advance()?;
                                Step::S
                            }
                            Token::Word("T") => {
                                self.advance()?;
                                self.expect(Token::LBracket)?;
                                let power =
                                    self.signed_wide("a signed twist power", "twist power")?;
                                self.expect(Token::RBracket)?;
                                Step::T(power)
                            }
                            _ => return Err(self.expected("S or T[power] in the WZW gluing word")),
                        };
                        steps.push(step);
                        if !self.eat(Token::Comma)? {
                            break;
                        }
                    }
                    self.expect(Token::RBrace)?;
                    return expr(span, ExprForm::WzwTv(rank, level, steps), 0);
                }
                if name.text == "rvCycloRat" && self.eat(Token::LBracket)? {
                    let Token::Number(decimal) = self.token() else {
                        return Err(self.expected("a positive cyclotomic rational denominator"));
                    };
                    let denominator = decimal.parse::<u128>().map_err(|_| {
                        Error::new(
                            self.span(),
                            "cyclotomic rational denominator exceeds the u128 range",
                        )
                    })?;
                    self.advance()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::CycloRatio(denominator, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvWzwCoord" && self.eat(Token::LBracket)? {
                    let weight = self.signed_tuple("a signed WZW weight", "WZW weight")?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::WzwCoordinate(weight, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvWzwScale" && self.token() == Token::LParen {
                    let arguments = self.paren_pair(nest)?;
                    let depth = 1 + arguments[0].depth.max(arguments[1].depth);
                    return expr(span, ExprForm::WzwScale(Box::new(arguments)), depth);
                }
                if let Some(operation) = match name.text.as_str() {
                    "rvWzwS" if self.token() == Token::LParen => Some(WzwOp::S(false)),
                    "rvWzwSAdjoint" if self.token() == Token::LParen => Some(WzwOp::S(true)),
                    "rvWzwCharge" if self.token() == Token::LParen => Some(WzwOp::Charge),
                    "rvWzwT" if self.eat(Token::LBracket)? => {
                        let power = self.signed_wide("a signed twist power", "twist power")?;
                        self.expect(Token::RBracket)?;
                        Some(WzwOp::T(power))
                    }
                    _ => None,
                } {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::WzwUnary(operation, Box::new(input)), depth);
                }
                if name.text == "rvCycloEmbed" && self.eat(Token::LBracket)? {
                    let (order, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::CycloEmbed(order, Box::new(input)), depth);
                }
                if name.text == "rvCyclo" && self.eat(Token::LBracket)? {
                    let (order, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let power =
                        self.signed_wide("a signed cyclotomic exponent", "cyclotomic exponent")?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::CycloRoot(order, power), 0);
                }
                if name.text == "rvCycloMul" && self.token() == Token::LParen {
                    let arguments = self.paren_pair(nest)?;
                    let depth = 1 + arguments[0].depth.max(arguments[1].depth);
                    return expr(span, ExprForm::CycloProduct(Box::new(arguments)), depth);
                }
                if name.text == "rvCycloScale" && self.eat(Token::LBracket)? {
                    let factor = self
                        .signed_wide("a signed cyclotomic coefficient", "cyclotomic coefficient")?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::CycloScale(factor, Box::new(input)), depth);
                }
                if name.text == "rvCycloConj" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::CycloConjugate(Box::new(input)), depth);
                }
                if name.text == "rvQDim" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::Verlinde(None, Box::new(input)), depth);
                }
                if name.text == "rvVerlinde" && self.eat(Token::LBracket)? {
                    let point = self.signed_tuple("a signed torus weight", "torus weight")?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::Verlinde(Some(point), Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvFusionMul" && self.token() == Token::LParen {
                    let arguments = self.paren_pair(nest)?;
                    let depth = 1 + arguments[0].depth.max(arguments[1].depth);
                    return expr(span, ExprForm::FusionProduct(Box::new(arguments)), depth);
                }
                if name.text == "rvFusionScale" && self.eat(Token::LBracket)? {
                    let factor =
                        self.signed_wide("a signed fusion coefficient", "fusion coefficient")?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::FusionScale(factor, Box::new(input)), depth);
                }
                if name.text == "rvFusionDual" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::FusionDual(Box::new(input)), depth);
                }
                if matches!(name.text.as_str(), "rvQGrBasis" | "rvQGrOne")
                    && self.eat(Token::LBracket)?
                {
                    let (rank, ambient) = self.quantum_grade()?;
                    let shape = if name.text == "rvQGrBasis" {
                        self.expect(Token::Comma)?;
                        self.exponents()?
                    } else {
                        Vec::new()
                    };
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::QuantumBasis(rank, ambient, shape), 0);
                }
                if name.text == "rvOre" && self.eat(Token::LBracket)? {
                    let (q, theta) = self.quantum_grade()?;
                    let coefficient = if self.eat(Token::Comma)? {
                        self.signed_wide("a signed Ore coefficient", "Ore coefficient")?
                    } else {
                        1
                    };
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::OreTerm(q, theta, coefficient), 0);
                }
                if matches!(name.text.as_str(), "rvOreMul" | "rvOreApply")
                    && self.token() == Token::LParen
                {
                    let arguments = self.paren_pair(nest)?;
                    let depth = 1 + arguments[0].depth.max(arguments[1].depth);
                    let form = if name.text == "rvOreMul" {
                        ExprForm::OreProduct(Box::new(arguments))
                    } else {
                        ExprForm::OreApply(Box::new(arguments))
                    };
                    return expr(span, form, depth);
                }
                if name.text == "rvQGrMod" && self.eat(Token::LBracket)? {
                    let (prime, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::QuantumReduce(prime, Box::new(input)), depth);
                }
                if matches!(name.text.as_str(), "rvQGrPBasis" | "rvQGrCurvatureCol")
                    && self.eat(Token::LBracket)?
                {
                    let (prime, rank, ambient) = self.quantum_mod_grade()?;
                    self.expect(Token::Comma)?;
                    let shape = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    let column = name.text == "rvQGrCurvatureCol";
                    return expr(
                        span,
                        ExprForm::QuantumPSection(prime, rank, ambient, shape, column),
                        0,
                    );
                }
                if matches!(
                    name.text.as_str(),
                    "rvQGrCurvature" | "rvQGrCurvatureSatake"
                ) && self.token() == Token::LParen
                {
                    let satake = name.text == "rvQGrCurvatureSatake";
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::QuantumCurvature(satake, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvQGrNabla" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::QuantumNabla(Box::new(input)), depth);
                }
                if name.text == "rvQGrAnnihilator" && self.eat(Token::LBracket)? {
                    let (rank, ambient) = self.quantum_grade()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::QuantumAnnihilator(rank, ambient), 0);
                }
                if name.text == "rvQGrMul" && self.token() == Token::LParen {
                    let arguments = self.paren_pair(nest)?;
                    let depth = 1 + arguments[0].depth.max(arguments[1].depth);
                    return expr(span, ExprForm::QuantumProduct(Box::new(arguments)), depth);
                }
                if matches!(name.text.as_str(), "rvQGrScale" | "rvQGrShift")
                    && self.eat(Token::LBracket)?
                {
                    let operation =
                        if name.text == "rvQGrScale" {
                            QuantumOp::Scale(self.signed_wide(
                                "a signed Schubert coefficient",
                                "Schubert coefficient",
                            )?)
                        } else {
                            QuantumOp::Shift(self.number()?.0)
                        };
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::QuantumUnary(operation, Box::new(input)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "rvQGrTranspose" | "rvQGrClassical")
                    && self.token() == Token::LParen
                {
                    let operation = if name.text == "rvQGrTranspose" {
                        QuantumOp::Transpose
                    } else {
                        QuantumOp::Classical
                    };
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::QuantumUnary(operation, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvSphBasis" && self.eat(Token::LBracket)? {
                    let coweight = self.affine_window()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::SphericalBasis(coweight), 0);
                }
                if matches!(name.text.as_str(), "rvSphId" | "rvSphQ")
                    && self.eat(Token::LBracket)?
                {
                    let degree = usize::from(name.text == "rvSphQ");
                    let (rank, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::SphericalScalar(rank, degree), 0);
                }
                if name.text == "rvSphMul" && self.token() == Token::LParen {
                    let arguments = self.paren_pair(nest)?;
                    let depth = 1 + arguments[0].depth.max(arguments[1].depth);
                    return expr(span, ExprForm::SphericalProduct(Box::new(arguments)), depth);
                }
                if name.text == "rvSphIwahori" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::SphericalEmbedding(Box::new(input)), depth);
                }
                if name.text == "rvSatakeInverse" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::SatakeInverse(Box::new(input)), depth);
                }
                if matches!(name.text.as_str(), "rvSatSchur" | "rvSphIC")
                    && self.eat(Token::LBracket)?
                {
                    let intersection = name.text == "rvSphIC";
                    let coweight =
                        self.signed_tuple("a signed coweight entry", "coweight entry")?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::SatakeCharacter(intersection, coweight), 0);
                }
                if name.text == "rvSatake" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::SatakeTransform(Box::new(input)), depth);
                }
                if name.text == "rvSatV" && self.eat(Token::LBracket)? {
                    let (rank, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let power =
                        self.signed_wide("a signed Laurent v exponent", "Laurent v exponent")?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::SatakeScalar(rank, power), 0);
                }
                if name.text == "rvSatakeMul" && self.token() == Token::LParen {
                    let arguments = self.paren_pair(nest)?;
                    let depth = 1 + arguments[0].depth.max(arguments[1].depth);
                    return expr(span, ExprForm::SatakeProduct(Box::new(arguments)), depth);
                }
                if name.text == "rvU" && self.eat(Token::LBracket)? {
                    let power = self.signed_number("Laurent exponent")?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RVLaurentMonomial(power, Box::new(input)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "rvUNorm" | "rvNewtonCoeff")
                    && self.eat(Token::LBracket)?
                {
                    use crate::rv_newton::Operation;
                    let (modulus, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let (root, _) = self.number()?;
                    let operation = if name.text == "rvNewtonCoeff" {
                        self.expect(Token::Comma)?;
                        Operation::Coefficient(self.number()?.0)
                    } else {
                        Operation::Norm
                    };
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RVNewton(operation, modulus, root, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvSpec" && self.eat(Token::LBracket)? {
                    let (modulus, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let point = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RVSpecialize(modulus, point, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "qextgen" && self.eat(Token::LBracket)? {
                    let spec = self.q_spec(nest)?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::QExtGenerator(spec, index, at), 0);
                }
                if name.text == "vextgen" && self.eat(Token::LBracket)? {
                    let spec = self.v_spec(nest)?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::VExtGenerator(spec, index, at), 0);
                }
                if name.text == "vbockgen" && self.eat(Token::LBracket)? {
                    let spec = self.vbock_spec(nest)?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::VBockGenerator(spec, index, at), 0);
                }
                if name.text == "vbClass" && self.eat(Token::LBracket)? {
                    let (r, f) = self.bidegree()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::VBockClass(r, f, Box::new(value)), depth);
                }
                if name.text == "vShift" && self.eat(Token::LBracket)? {
                    let k = self.signed_degree()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::VShift(k, Box::new(input)), depth);
                }
                if name.text == "vRatio" && self.eat(Token::LBracket)? {
                    let (k, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::VRatio(k, Box::new(input)), depth);
                }
                if matches!(name.text.as_str(), "vbarRep" | "vbarNull")
                    && self.eat(Token::LBracket)?
                {
                    let kind = if name.text == "vbarRep" {
                        BarCertificate::Representative
                    } else {
                        BarCertificate::Nullhomotopy
                    };
                    let spec = self.bar_family(nest)?;
                    self.expect(Token::Comma)?;
                    let degree = self.signed_degree()?;
                    self.expect(Token::Comma)?;
                    let (exponent, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let (index, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::VBarCertificate(kind, spec, degree, exponent, index),
                        0,
                    );
                }
                if name.text == "qRes" && self.token() == Token::LBracket {
                    let selected = self.exponents()?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::QRestrict(selected, Box::new(input)), depth);
                }
                if name.text == "mapmaygen" && self.eat(Token::LBracket)? {
                    let spec = self.map_may_spec(nest)?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::MapMayGenerator(spec, index, at), 0);
                }
                if matches!(name.text.as_str(), "mayangen" | "ceangen")
                    && self.eat(Token::LBracket)?
                {
                    let (n, r, s, t, f) = self.profile_page_degree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let kind = if name.text == "ceangen" {
                        ExprForm::CeGenerator(n, r, s, t, f, index, at)
                    } else {
                        ExprForm::MayAnGenerator(n, r, s, t, f, index, at)
                    };
                    return expr(span, kind, 0);
                }
                if name.text == "maygen" && self.eat(Token::LBracket)? {
                    let (r, s, t, f) = self.may_degree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::MayGenerator(r, s, t, f, index, at), 0);
                }
                if name.text == "pmaygen" && self.eat(Token::LBracket)? {
                    let (p, r, s, t, f) = self.pmay_degree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::PMayGenerator(p, r, s, t, f, index, at), 0);
                }
                if matches!(name.text.as_str(), "pmayangen" | "pceangen")
                    && self.eat(Token::LBracket)?
                {
                    let (p, n, r, s, t, f) = self.prime_profile_page_degree()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let kind = if name.text == "pceangen" {
                        ExprForm::PCeAnGenerator(p, n, r, s, t, f, index, at)
                    } else {
                        ExprForm::PMayAnGenerator(p, n, r, s, t, f, index, at)
                    };
                    return expr(span, kind, 0);
                }
                if matches!(name.text.as_str(), "mayClass" | "ceClass")
                    && self.eat(Token::LBracket)?
                {
                    let (r, f) = self.bidegree()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    let kind = if name.text == "ceClass" {
                        ExprForm::CeClass(r, f, Box::new(value))
                    } else {
                        ExprForm::MayClass(r, f, Box::new(value))
                    };
                    return expr(span, kind, depth);
                }
                if matches!(name.text.as_str(), "bar" | "binom" | "carry")
                    && self.eat(Token::LBracket)?
                {
                    let (modulus, at) = self.number()?;
                    let form = match name.text.as_str() {
                        "binom" => {
                            self.expect(Token::Comma)?;
                            let (index, index_at) = self.number()?;
                            self.expect(Token::RBracket)?;
                            CyclicForm::Binomial(index, index_at)
                        }
                        "bar" => {
                            self.expect(Token::RBracket)?;
                            self.expect(Token::LBrace)?;
                            let mut support = Vec::new();
                            while self.token() != Token::RBrace {
                                support.push(self.number()?);
                                if !self.eat(Token::Comma)? {
                                    break;
                                }
                            }
                            self.expect(Token::RBrace)?;
                            CyclicForm::Support(support)
                        }
                        _ => {
                            self.expect(Token::RBracket)?;
                            CyclicForm::Carry
                        }
                    };
                    return expr(span, ExprForm::Cyclic(modulus, at, form), 0);
                }
                if name.text == "sq" && self.eat(Token::LBracket)? {
                    let (degree, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::Square(degree), 0);
                }
                if name.text == "coact" && self.eat(Token::LBracket)? {
                    let (n, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::Coact(n, index, at, Box::new(value)), depth);
                }
                if name.text == "shift" && self.eat(Token::LBracket)? {
                    let (shift, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    let value = self.paren_expr(nest)?;
                    let depth = value.depth + 1;
                    return expr(span, ExprForm::Shift(shift, Box::new(value)), depth);
                }
                let suspension = match name.text.as_str() {
                    "north" => Some(SuspensionForm::North),
                    "south" => Some(SuspensionForm::South),
                    "equator" => Some(SuspensionForm::Equator),
                    _ => None,
                };
                if let Some(form) = suspension {
                    if self.eat(Token::LBracket)? {
                        let base = self.ty(nest + 1)?;
                        self.expect(Token::RBracket)?;
                        return expr(span, ExprForm::Suspension(base, form), 0);
                    }
                }
                if matches!(name.text.as_str(), "zero" | "gen" | "cocycle")
                    && self.eat(Token::LBracket)?
                {
                    let base = self.ty(nest + 1)?;
                    if name.text == "zero" && self.eat(Token::RBracket)? {
                        return expr(span, ExprForm::Zero(base), 0);
                    }
                    self.expect(Token::Comma)?;
                    let (degree, _) = self.number()?;
                    let form = match name.text.as_str() {
                        "gen" => {
                            self.expect(Token::Comma)?;
                            let (index, at) = self.number()?;
                            self.expect(Token::RBracket)?;
                            ClassForm::Generator(index, at)
                        }
                        "cocycle" => {
                            self.expect(Token::RBracket)?;
                            ClassForm::Cocycle(self.chains()?)
                        }
                        _ => {
                            self.expect(Token::RBracket)?;
                            ClassForm::Zero
                        }
                    };
                    return expr(span, ExprForm::Class(base, degree, form), 0);
                }
                if self.eat(Token::Dot)? {
                    let ty = ty(span, TypeForm::Name(name), 0)?;
                    return expr(span, ExprForm::Point(ty, self.label()?), 0);
                }
                expr(span, ExprForm::Name(name), 0)
            }
            _ => Err(self.expected("an expression")),
        }
    }

    fn rv_bar_page_spec(&mut self) -> Result<RvBarPageSpec> {
        let path = self.label()?;
        self.expect(Token::Comma)?;
        let (start, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (end, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (r, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (b, _) = self.number()?;
        self.expect(Token::Comma)?;
        let length = self.signed_degree()?;
        Ok(RvBarPageSpec {
            chain: RvBarSpec {
                path,
                start,
                end,
                b,
            },
            r,
            length,
        })
    }

    fn rv_bar_spec(&mut self) -> Result<RvBarSpec> {
        let path = self.label()?;
        self.expect(Token::Comma)?;
        let (start, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (end, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (b, _) = self.number()?;
        Ok(RvBarSpec {
            path,
            start,
            end,
            b,
        })
    }

    fn rv_grade(&mut self, multigraded: bool) -> Result<RvGradeSpec> {
        let ideal = self.label()?;
        self.expect(Token::Comma)?;
        let (s, _) = self.number()?;
        self.expect(Token::Comma)?;
        let degree = if multigraded {
            RvDegree::Multi(self.exponents()?)
        } else {
            RvDegree::Total(self.number()?.0)
        };
        Ok(RvGradeSpec { ideal, s, degree })
    }

    fn affine_window(&mut self) -> Result<Vec<i128>> {
        self.signed_tuple("a signed affine window entry", "affine window entry")
    }

    fn quantum_mod_grade(&mut self) -> Result<(usize, usize, usize)> {
        let (prime, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (rank, ambient) = self.quantum_grade()?;
        Ok((prime, rank, ambient))
    }

    fn quantum_grade(&mut self) -> Result<(usize, usize)> {
        let (rank, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (ambient, _) = self.number()?;
        Ok((rank, ambient))
    }

    fn fusion_grade(&mut self) -> Result<(usize, usize)> {
        let (rank, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (level, _) = self.number()?;
        Ok((rank, level))
    }

    fn hive_spec(&mut self) -> Result<HiveSpec> {
        let lambda = self.signed_tuple("a signed hive weight", "hive weight")?;
        self.expect(Token::Comma)?;
        let mu = self.signed_tuple("a signed hive weight", "hive weight")?;
        self.expect(Token::Comma)?;
        let nu = self.signed_tuple("a signed hive weight", "hive weight")?;
        Ok([lambda, mu, nu])
    }

    fn signed_tuple(&mut self, expected: &str, label: &str) -> Result<Vec<i128>> {
        self.expect(Token::LBracket)?;
        self.signed_list(Token::RBracket, expected, label)
    }

    fn signed_list(&mut self, end: Token<'src>, expected: &str, label: &str) -> Result<Vec<i128>> {
        let mut values = Vec::new();
        while self.token() != end {
            values.push(self.signed_wide(expected, label)?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(end)?;
        Ok(values)
    }

    fn coweight(&mut self) -> Result<Vec<i64>> {
        self.expect(Token::LBracket)?;
        let mut mu = Vec::new();
        while self.token() != Token::RBracket {
            mu.push(self.signed_number("ADLV coweight entry")?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBracket)?;
        Ok(mu)
    }

    fn paren_pair(&mut self, nest: usize) -> Result<[ParsedExpr; 2]> {
        self.expect(Token::LParen)?;
        let first = self.expr(nest + 1)?;
        self.expect(Token::Comma)?;
        let second = self.expr(nest + 1)?;
        self.expect(Token::RParen)?;
        Ok([first, second])
    }

    fn embedding_parameters(&mut self) -> Result<(usize, usize)> {
        let (modulus, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (root, _) = self.number()?;
        self.expect(Token::RBracket)?;
        Ok((modulus, root))
    }

    fn paren_expr(&mut self, nest: usize) -> Result<ParsedExpr> {
        self.expect(Token::LParen)?;
        let value = self.expr(nest + 1)?;
        self.expect(Token::RParen)?;
        Ok(value)
    }
}

fn check_nesting(depth: usize, span: Span) -> Result<()> {
    if depth > MAX_NESTING {
        return Err(Error::new(span, "syntax nesting exceeds 64 levels"));
    }
    Ok(())
}

fn ty(span: Span, kind: TypeForm, depth: usize) -> Result<ParsedType> {
    check_nesting(depth, span)?;
    Ok(ParsedType { span, kind, depth })
}

fn expr(span: Span, kind: ExprForm, depth: usize) -> Result<ParsedExpr> {
    check_nesting(depth, span)?;
    Ok(ParsedExpr { span, kind, depth })
}
