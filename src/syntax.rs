use crate::bp::Operation as BpOperation;
use crate::bp_page::{Mode as BpPageMode, Operation as BpPageOp};
use crate::bp_vpage::Operation as BpVPageOp;
use crate::comodule::Generator;
use crate::dependent_sets::Operation as SetOperation;
use crate::eo::Report as EoReport;
use crate::error::{Error, Result, Span};
use crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::eigenpacket;
use crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::hecke;
use crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::period_polynomial as period;
use crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::quadratic_hecke;
use crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::residual_hecke;
use crate::map_bracket::Operation as MapBracketOp;
use crate::map_system::Operation as MapSystemOp;
use crate::map_transfer::Operation as MapTransferOp;
use crate::model::DiffFamily;
use crate::motivic::Grade as MotivicGrade;
use crate::odd_cobar::FactorLiteral;
use crate::odd_massey::Operation as OddSystemOp;
use crate::real_extensions::Operation as RealExtOp;
use crate::real_rho::{self, Operation as RealRhoOp};
use crate::rho_bock;
use crate::rho_localized;
use crate::rho_minimal::Operation as RhoMinOp;
use crate::rv_arithmetic::Series as RankSeries;
use crate::rv_as_characters::Query as AsCharacterQuery;
use crate::rv_as_cover::Query as AsQuery;
use crate::rv_as_derham::{Operator as AsOperator, WordQuery as AsWordQuery};
use crate::rv_as_odd::Query as AsOddQuery;
use crate::rv_as_odd_cartier::Operation as AsOddDiffOp;
use crate::rv_as_odd_derham::Operation as AsOddDeRhamOp;
use crate::rv_as_odd_tate::Operation as AsOddTateOp;
use crate::rv_as_two::Query as AsTwoQuery;
use crate::rv_as_two_cartier::{Operation as AsTwoDiffOp, Query as AsTwoCartierQuery};
use crate::rv_crystal::{Entries as CrystalEntries, Query as CrystalQuery};
use crate::rv_gauss::Operation as GaussOp;
use crate::rv_phodge::{QuadraticSpec, Query as PHodgeQuery, Specification as PHodgeSpec};
use crate::rv_qp_quad::Unary as QuadUnary;
use crate::rv_rees::Operation as ReesOp;
use crate::rv_t_bracket::Operation as RvTBracketOp;
use crate::rv_t_system::Operation as RvTSystemOp;
use crate::rv_t_transfer::Operation as RvTTransferOp;
use crate::rv_witt_cover::Query as WittQuery;
use crate::rv_witt_tower::Query as TowerQuery;
use crate::set_kan::Operation as KanOp;
use crate::set_logic::Operation as SetLogic;
use crate::sheaf::{self, Form as SheafForm, Operation as SheafOp};
use crate::sheaf_higher::Operation as SheafHigherOp;
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
    SetNat(Label, Label),
    Star(Box<ParsedType>, Label),
    Cohomology(Box<ParsedType>, usize),
    MixedHomology(Label, i64),
    MixedChain(Label, i64, bool),
    MixedMinimal(Label, i64, bool),
    MixedPage(MixedSpec),
    Sheaf(SheafForm, Label, usize),
    SheafDeformation(Label, usize),
    Convolution(crate::sheaf_disk::convolution::Form, Label),
    Parabolic(crate::sheaf_disk::parabolic::Form, Label),
    DiskHom(Label, Label),
    SheafHom(SheafForm, Label, Label, usize),
    SheafComplex(SheafForm, Label, usize),
    SheafPage(Label, SheafPageGrade),
    SheafLocal(Label, Label, SupportMode, usize),
    Operations(usize),
    Cobar(usize, usize),
    MCobar(MotivicGrade),
    RMCobar(MotivicGrade),
    RealRho(bool, real_rho::Grade),
    RhoBock(rho_bock::Grade),
    RhoLoc(MotivicGrade),
    RhoMin(MotivicGrade),
    MExt(MotivicGrade),
    RMExt(MotivicGrade),
    RMBar(MotivicGrade),
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
    RVWittSpectrum,
    RVCTState(Label, usize, usize),
    RVCTResidue(usize, usize),
    RVFusion(usize, usize),
    RVCyclo(usize),
    Jacobian(Label),
    Kummer(Label),
    KummerK3(Label),
    Enriques(Label),
    K3PencilBase(Label),
    K3PencilSection(Label),
    K3PencilKummer2(Label),
    KummerDisc(Label),
    KummerWeil(Label),
    KummerDiv(Label),
    K3Lattice(Label),
    Nakajima(Label, i128, i128),
    HilbertCohomology(Label, i128, i128),
    K3Clifford(Label, u8),
    OG6Clifford(Label, u8),
    KSExterior(Label, i128),
    KSEnd(Label),
    K3Spinor(Label, u8),
    Mukai(Label),
    MukaiQuotient(Label),
    OG10(Label),
    AbMukai(Label),
    OG6(Label),
    HilbertSquare(Label, i128),
    SiegelWeight,
    ModularHecke(u128),
    HeckeEigen(Label),
    PeriodPolynomial(period::Grade),
    ECMap(Label, Label),
    ECCurve(Label),
    ECPoint(Label),
    ECKernel(Label, usize),
    ECIso(Label),
    HWMatrix(usize, usize, usize),
    HWWMatrix(Label, usize),
    HWScalar(usize, usize),
    RVUnram(Label),
    RVCycloQ(usize),
    RVQpQuad(Label),
    RVResSpace(usize, usize),
    RVValMat(usize, usize, usize),
    HyperArithmetic(usize, usize),
    HyperCover(usize, usize),
    StableTree(usize),
    StableMarks(usize, usize),
    Lee(Label, crate::lee::Grade),
    LeeDual(Label, i64),
    Khovanov(Label, crate::khovanov::Grade, crate::khovanov::Form),
    KhovanovMod(Label, usize, Option<usize>, crate::khovanov::Grade, i64),
    CMOrder(Label, usize),
    ECEndOrder(Label),
    CMClass(Label),
    RvTl(usize, usize, usize),
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
    ASTwoDiff(Label),
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

pub(crate) struct WzMirror {
    pub(crate) period: Label,
    pub(crate) columns: Option<Vec<Vec<i128>>>,
}

#[derive(Clone)]
pub(crate) enum TlOp {
    Generator(usize),
    Braid(Vec<i128>),
    Project(usize, usize),
    Negate,
}

#[derive(Clone, Copy)]
pub(crate) enum WzwOp {
    S(bool),
    T(i128),
    Charge,
    Galois(i128),
    CoeffGalois(i128),
    Negate,
}

impl WzwOp {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::S(false) => "rvWzwS",
            Self::S(true) => "rvWzwSAdjoint",
            Self::T(_) => "rvWzwT",
            Self::Charge => "rvWzwCharge",
            Self::Galois(_) => "rvWzwGalois",
            Self::CoeffGalois(_) => "rvWzwCoeffGalois",
            Self::Negate => "neg",
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum OreOp {
    Factorial(usize),
    EulerFactor,
    EulerReduce,
}

impl OreOp {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Factorial(_) => "rvOreFactorial",
            Self::EulerFactor => "rvOreEulerFactor",
            Self::EulerReduce => "rvOreEulerReduce",
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
    MixInduce,
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
    UnramMultiply,
    BpCup,
    ResSum,
    ResIntersection,
    QpQuadMultiply,
    CMMultiply,
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
            "mixInduce" => Some(Self::MixInduce),
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
            "rvUnramMul" => Some(Self::UnramMultiply),
            "bpCup" => Some(Self::BpCup),
            "rvResSum" => Some(Self::ResSum),
            "rvResIntersection" => Some(Self::ResIntersection),
            "rvQpQuadMul" => Some(Self::QpQuadMultiply),
            "rvCMMul" => Some(Self::CMMultiply),
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
    SetTopology(
        crate::set_topology::Operation,
        Label,
        Option<Box<ParsedExpr>>,
    ),
    SetLogic(SetLogic, Vec<Label>, Vec<ParsedExpr>),
    SetJUnit(Label, bool),
    SetJExtend(Label, Box<ParsedExpr>),
    SetProjection(Label),
    SetSection(Label, usize),
    SetGluing(crate::set_gluing::Operation, [Label; 2], Vec<ParsedExpr>),
    SetIdentity(Label),
    SetCompose(Box<[ParsedExpr; 2]>),
    SetTotal(Label, Label, Box<ParsedExpr>),
    SetDependent(SetOperation, Label, Option<Label>, Option<Box<ParsedExpr>>),
    SetKan(KanOp, Label, Option<Label>, Option<Box<ParsedExpr>>),
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
    MixedHGenerator(Label, i64, usize, Span),
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
    SheafKuranishi(Box<crate::sheaf_kuranishi::Query>),
    Lee(crate::lee::Query),
    LeeDt(Box<crate::lee::dt::Query>),
    Khovanov(Box<crate::khovanov::Query>),
    SheafDefStart(Label),
    SheafDefOperation(crate::sheaf_deform::Operation, Box<ParsedExpr>),
    SheafDefAppend(Box<[ParsedExpr; 2]>),
    SheafDefGauge(usize, bool, Box<[ParsedExpr; 2]>),
    Convolution(crate::sheaf_disk::convolution::Expression),
    Parabolic(crate::sheaf_disk::parabolic::Expression),
    Disk(crate::sheaf_disk::Expression),
    SheafConeApex(Label, usize),
    SheafConeTail(Label, bool, Box<ParsedExpr>),
    SheafConeTransfer(Label, bool, Box<ParsedExpr>),
    SheafBCApply(Label, Box<ParsedExpr>),
    SheafCompose(Box<[ParsedExpr; 2]>),
    SheafHigher(SheafHigherOp, Vec<ParsedExpr>, Vec<Witness>),
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
    RVCTInitial(Label, usize, usize),
    RVCTQuery(Label, usize, usize, u128),
    RVCTLiteral(usize, usize, u128),
    RVCTSection(usize, Box<ParsedExpr>),
    RVCTRead(Box<ParsedExpr>),
    OreTerm(usize, usize, i128),
    OreProduct(Box<[ParsedExpr; 2]>),
    OreUnary(OreOp, Box<ParsedExpr>),
    WittSpectrum(Box<crate::rv_witt_spectrum::Expression>),
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
    CycloGalois(i128, Box<ParsedExpr>),
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
    CycloRatioProduct(Box<[ParsedExpr; 2]>),
    CycloRatioGalois(i128, Box<ParsedExpr>),
    TlPath(usize, Vec<usize>),
    TlOperation(TlOp, Box<ParsedExpr>),
    TlCoordinate(Vec<usize>, Box<ParsedExpr>),
    TlScale(Box<[ParsedExpr; 2]>),
    TlTrace(usize, usize, usize, Vec<i128>),
    TlJones(usize, usize, Vec<i128>),
    TlColoredJones(usize, usize, usize, Vec<i128>),
    TlFramed(usize, usize, Vec<i128>, Vec<usize>, Vec<i128>),
    TlSurgery(crate::rv_tl::Surgery, usize, usize, Vec<i128>, Vec<i128>),
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
    RealRhoGenerator(real_rho::Grade, usize, Span),
    RealBar(RealBarOp, usize, i128, usize, Span),
    RealBarTerm(MotivicGrade, usize, Span),
    RealExtension(RealExtOp, Box<ParsedExpr>),
    RhoProductNull(usize, Box<[ParsedExpr; 2]>),
    RealRhoOperation(RealRhoOp, Box<ParsedExpr>),
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
    CrystalEO(Label),
    CrystalCongruent(Label, Label, usize),
    EOOddCover(Label),
    QpQuadSpecialize(Label, Box<ParsedExpr>),
    CMPi(Label, usize),
    CMTerm(Label, usize, i128, i128),
    CMUnary(crate::legendre_cm::Unary, Box<ParsedExpr>),
    UnramElement(Label, Vec<i128>),
    UnramUnary(crate::rv_unram::Unary, Box<ParsedExpr>),
    ECMap(Box<crate::ec_morphism::Expression>),
    Elliptic(Box<crate::elliptic::Expression>),
    EllipticCM(Box<crate::elliptic_cm_recognition::Expression>),
    EllipticLocal(Box<crate::elliptic_local::Expression>),
    HasseWitt(Box<crate::hyperelliptic::Expression>),
    Jacobian(Box<crate::hyperelliptic_jacobian::Expression>),
    JacobianRichelotKernel(Label, Label, usize, usize),
    JacobianLevelCoordinate(Label, usize, Box<ParsedExpr>),
    Kummer(Box<crate::hyperelliptic_kummer_surface::Expression>),
    Enriques(Box<crate::hyperelliptic_kummer_surface::enriques::points::Expression>),
    K3Pencil(Box<crate::hyperelliptic_kummer_surface::pencil::Expression>),
    K3PencilSection(Box<crate::hyperelliptic_kummer_surface::pencil::sections::Expression>),
    K3Mordell(Box<crate::hyperelliptic_kummer_surface::mordell::Expression>),
    KummerDisc(Box<crate::hyperelliptic_kummer_surface::nikulin::Expression>),
    KummerWeil(Box<crate::hyperelliptic_kummer_surface::weil::Expression>),
    KummerDiv(Box<crate::hyperelliptic_kummer_surface::polarized::Expression>),
    K3Lattice(Box<crate::hyperelliptic_kummer_surface::unimodular::Expression>),
    Nakajima(Box<crate::hyperelliptic_kummer_surface::nakajima::Expression>),
    HilbertCohomology(Box<crate::hyperelliptic_kummer_surface::nakajima::cohomology::Expression>),
    K3Clifford(Box<crate::hyperelliptic_kummer_surface::clifford::Expression>),
    OG6Clifford(
        Box<crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::Expression>,
    ),
    KSExterior(
        Box<crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::exterior::Expression>,
    ),
    KSEnd(
        Box<crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::projector::Expression>,
    ),
    K3Spinor(Box<crate::hyperelliptic_kummer_surface::spinor::Expression>),
    K3Period(Box<crate::hyperelliptic_kummer_surface::period::Expression>),
    CMClass(Box<crate::hyperelliptic_kummer_surface::cm_ring::Expression>),
    CMNeighbor(crate::hyperelliptic_kummer_surface::cm_neighbors::Expression),
    EllipticModular(Box<crate::elliptic_modular::Expression>),
    CMClassPolynomial(Box<crate::cm_class_poly::Expression>),
    CMArtin(Box<crate::cm_artin::Expression>),
    ECEndOrder(Box<crate::elliptic_end_order::Expression>),
    ECEndCanonical(Box<crate::elliptic_end_order::canonical::Expression>),
    ECHorizontal(Box<crate::elliptic_end_order::horizontal::Address>),
    ECCMCycle(Box<crate::elliptic_end_order::cycle::Expression>),
    Mukai(Box<crate::hyperelliptic_kummer_surface::mukai::Expression>),
    MukaiModuli(crate::hyperelliptic_kummer_surface::mukai::moduli::Expression),
    OG10(Box<crate::hyperelliptic_kummer_surface::mukai::moduli::og10::Expression>),
    AbMukai(Box<crate::hyperelliptic_kummer_surface::mukai::abelian::Expression>),
    OG6(Box<crate::hyperelliptic_kummer_surface::mukai::abelian::og6::Expression>),
    OG6Cusp(Box<crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::Expression>),
    OG6Nilpotent(
        Box<crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::nilpotent::Expression>,
    ),
    K3AutoDiv(Box<crate::hyperelliptic_kummer_surface::dynamics::Expression>),
    K3AutoLattice(Box<crate::hyperelliptic_kummer_surface::dynamics::Expression>),
    HilbertSquare(Box<crate::hyperelliptic_kummer_surface::hilbert_square::Expression>),
    K3AutoPoint(Box<crate::hyperelliptic_kummer_surface::dynamics::PointExpression>),
    SiegelWeight(crate::siegel::Weight),
    ModularHecke(hecke::Expression),
    HeckeEigen(eigenpacket::Expression),
    HeckeResidual(residual_hecke::Expression),
    PeriodPolynomial(period::Expression),
    SiegelReflect(usize, Box<ParsedExpr>),
    HWWClosed(Box<crate::hw_closed::Expression>),
    LegendreClosed(crate::legendre_closed::Operation, Box<ParsedExpr>),
    QpQuadElement(Label, i128, i128, i128),
    ResSpace(usize, usize, Vec<Vec<i128>>),
    ResSpecialize(usize, usize, Vec<i64>, Vec<Vec<(i128, i128)>>),
    ResSpecializeQuad(Label, usize, Vec<i64>, Vec<Vec<ParsedExpr>>),
    ResFil(Label, i64, Vec<i64>),
    ResOrthogonal(Box<ParsedExpr>),
    ValMat(usize, usize, Vec<Vec<i128>>),
    ValMatQuad(Label, usize, Vec<Vec<ParsedExpr>>),
    ValFil(Label, i64),
    ValDual(Box<ParsedExpr>),
    HyperArithmetic(
        usize,
        crate::hyperelliptic_graph::arith::Rational,
        Vec<crate::hyperelliptic_graph::arith::Rational>,
    ),
    HyperGeometric(Box<ParsedExpr>),
    HyperCover(Box<ParsedExpr>),
    StableLine(usize, Vec<[i128; 2]>),
    StableLineQuad(Label, Vec<[ParsedExpr; 2]>),
    StableProjection(Box<ParsedExpr>),
    QpQuadUnary(QuadUnary, Box<ParsedExpr>),
    AsDeRhamOperation(crate::rv_as_derham::Operation, Box<ParsedExpr>),
    AsOddCoefficient(Label, usize, usize),
    AsTwoCoefficient(Label, usize, usize),
    AsTwoOmega(Label, i128, usize),
    AsTwoOrdinaryGen(Label, usize),
    AsTwoDiffOperation(AsTwoDiffOp, Box<ParsedExpr>),
    AsOddOmega(Label, usize, usize),
    AsOddDiffOperation(AsOddDiffOp, Box<ParsedExpr>),
    AsOddDeRhamGen(Label, usize),
    AsOddTateGen(Label, i64, usize),
    AsOddTateOperation(AsOddTateOp, Box<ParsedExpr>),
    AsOddDeRhamOperation(AsOddDeRhamOp, Box<ParsedExpr>),
    PrimeScalar(usize, i128),
    KummerHilbertCoefficient(
        Box<crate::hyperelliptic_kummer_surface::hilbert::characteristic::Query>,
    ),
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
    SetFamily(
        Label,
        Label,
        Vec<(Label, usize)>,
        Vec<crate::set_family::Transition>,
    ),
    SetElements(Label, Label),
    SetNatural(Label, Label, Label, Vec<crate::set_family::Component>),
    SetSections(Label),
    SetRecollement(Label, Box<[ParsedExpr; 2]>),
    SetRestrict(Label, Label, Label, crate::set_gluing::Side),
    SetGlueLimit(Label, Label, Label),
    SetGlue(Label, Label, ParsedExpr),
    SetPullback(Label, Label, ParsedExpr),
    SetKan(Label, Label, ParsedExpr, bool),
    SetTopology(crate::set_topology::Command),
    SetOmega(Label, Label),
    SetSubobject(Label, ParsedExpr),
    SetSieve(Label, Label, usize, bool),
    SetJ(Label, Label),
    SetJCheck(Label),
    SetDependent(Label, Label, Label, bool),
    SetMonodromy(Label, Label, Label, Label),
    SetMonoQuery(Label, crate::set_monodromy::Query),
    SetLinearize(Label, Label, usize),
    SetOrbitalMap(
        Label,
        Label,
        Label,
        Label,
        crate::set_monodromy::Coefficients,
    ),
    SheafMapRank(Label, Label),
    MixedComplex(
        Label,
        usize,
        Vec<(Label, i64)>,
        crate::mixed::Rules,
        crate::mixed::Rules,
    ),
    RvHHMixed(Label, Label, Vec<usize>),
    MixedHDimension(Label, i64),
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
    Convolution(crate::sheaf_disk::convolution::Command),
    Parabolic(crate::sheaf_disk::parabolic::Command),
    Disk(crate::sheaf_disk::Command),
    SheafConeGlue(SheafInputForm, [Label; 4], usize, Vec<ParsedExpr>),
    SheafConeHomotopy([Label; 3], Vec<ParsedExpr>),
    SheafBCSquare(Label, [Label; 4]),
    SheafBCBuild(SheafInputForm, Label, Label, Label),
    SheafBCMap(Label, Label, usize, bool),
    SheafBCLinear(Label, Label, usize, crate::sheaf_base_change::LinearQuery),
    SheafBCUniversal(Label, usize),
    SheafBCComma(Label, usize, Label, Label, Option<(i64, usize)>),
    SheafMap(Label, Label, Label, Vec<(Label, Vec<i64>)>),
    SheafKuranishi(Label, Label, usize, Option<Vec<ParsedExpr>>),
    SheafKuranishiEquation(Label, usize),
    SheafKuranishiDimensions(Label),
    SheafExtDimension(Label, Label, usize),
    SheafComplex(Label, Vec<Label>),
    SheafResolution(Label, Label, usize, crate::sheaf_cellular::Presentation),
    SheafCellular(Label, Label, Label),
    SheafHyperDimension(Label, usize),
    SheafPageDimension(Label, SheafPageGrade),
    SheafExtra(SheafInputForm, Label, Label, Label),
    SheafLocalDimension(Label, Label, SupportMode, usize),
    SheafDimension(Label, usize),
    SheafSerre(SheafInputForm, Label, Label),
    SheafSerreShift(Label),
    SheafStalkCohomology(Label, Label, usize),
    SheafLeray(SheafInputForm, Label, Label, Label),
    SheafTruncate(Label, Label, usize),
    SheafDeligne(SheafInputForm, Label, Label, Label, usize),
    SheafImage(SheafImageOp, Label, Label, Label),
    Lee(crate::lee::Command),
    LeeDt(Box<crate::lee::dt::Command>),
    Khovanov(Box<crate::khovanov::Command>),
    ECMapReport(crate::ec_morphism::Report, Box<ParsedExpr>),
    Elliptic(Box<crate::elliptic::Command>),
    EllipticCM(Box<crate::elliptic_cm_recognition::Command>),
    EllipticLocal(Box<crate::elliptic_local::Command>),
    HasseWitt(Box<crate::hyperelliptic::Command>),
    Jacobian(Box<crate::hyperelliptic_jacobian::Command>),
    KummerK3(Box<crate::hyperelliptic_kummer_surface::k3::Command>),
    K3Pencil(Box<crate::hyperelliptic_kummer_surface::pencil::Command>),
    K3FiberEquation(Box<crate::hyperelliptic_kummer_surface::pencil::fiber::Command>),
    K3Mordell(Box<crate::hyperelliptic_kummer_surface::mordell::Command>),
    KummerEnriques(Box<crate::hyperelliptic_kummer_surface::enriques::Command>),
    DeclareEnriques(Box<crate::hyperelliptic_kummer_surface::enriques::points::Declare>),
    KummerHilbert(Box<crate::hyperelliptic_kummer_surface::hilbert::Command>),
    HilbertLlv(Box<crate::hyperelliptic_kummer_surface::llv_character::Command>),
    HilbertStratum(Box<crate::hyperelliptic_kummer_surface::stratum::Command>),
    DeclareSpecht(Box<crate::hyperelliptic_kummer_surface::specht::Declaration>),
    CrystallineNewton(crate::hyperelliptic_kummer_surface::newton::Command),
    GeneralizedKummer(Box<crate::generalized_kummer::Command>),
    Hodge(crate::generalized_kummer::hodge::Command),
    KummerDiv(Box<crate::hyperelliptic_kummer_surface::polarized::Command>),
    JacobianLevel(Box<crate::jacobian_level::Command>),
    KummerNodeLattice(Box<crate::hyperelliptic_kummer_surface::nikulin::Command>),
    KummerWeilCoefficient(crate::hyperelliptic_kummer_surface::weil::Coefficient),
    K3Lattice(Box<crate::hyperelliptic_kummer_surface::unimodular::Command>),
    NakajimaPair(Box<[ParsedExpr; 2]>),
    HilbertCohomology(Box<crate::hyperelliptic_kummer_surface::nakajima::cohomology::Report>),
    K3CliffordCoefficient(crate::hyperelliptic_kummer_surface::clifford::Coefficient),
    OG6Clifford(Box<crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::Command>),
    KSExterior(
        Box<
            crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::exterior::Command,
        >,
    ),
    KSEnd(
        Box<crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::projector::Command>,
    ),
    K3SpinorCoefficient(crate::hyperelliptic_kummer_surface::spinor::Coefficient),
    NakajimaPrimary(Box<crate::hyperelliptic_kummer_surface::nakajima::primary::Command>),
    K3Period(Box<crate::hyperelliptic_kummer_surface::period::Statement>),
    HilbertTautological(Box<crate::hyperelliptic_kummer_surface::nakajima::tautological::Command>),
    CMClass(Box<crate::hyperelliptic_kummer_surface::cm_ring::Statement>),
    CMNeighbors(crate::hyperelliptic_kummer_surface::cm_neighbors::Statement),
    EllipticModular(crate::elliptic_modular::Statement),
    CMClassPolynomial(Box<crate::cm_class_poly::Statement>),
    CMArtin(Box<crate::cm_artin::Statement>),
    ECEndOrder(Box<crate::elliptic_end_order::Statement>),
    ECHorizontal(Box<crate::elliptic_end_order::horizontal::Statement>),
    ECCMCycle(Box<crate::elliptic_end_order::cycle::Statement>),
    Mukai(Box<crate::hyperelliptic_kummer_surface::mukai::Command>),
    MukaiModuli(Box<crate::hyperelliptic_kummer_surface::mukai::moduli::Command>),
    OG10(Box<crate::hyperelliptic_kummer_surface::mukai::moduli::og10::Command>),
    AbMukaiPair(Box<[ParsedExpr; 2]>),
    OG6(Box<crate::hyperelliptic_kummer_surface::mukai::abelian::og6::Command>),
    OG6Cusp(Box<crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::Command>),
    OG6Nilpotent(
        Box<crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::nilpotent::Command>,
    ),
    K3Automorphism(Box<crate::hyperelliptic_kummer_surface::dynamics::Command>),
    HilbertSquareIntegral(Box<ParsedExpr>),
    K3AutoPeriod(Box<crate::hyperelliptic_kummer_surface::dynamics::Expression>),
    HilbertSquareNakajimaDen(Box<ParsedExpr>),
    SiegelPair(Box<[ParsedExpr; 2]>),
    SiegelRoot(Box<ParsedExpr>),
    SiegelDelta5At(Box<ParsedExpr>),
    SiegelMaass(Box<crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::maass::Command>),
    SiegelChamber(
        Box<crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::farey::Command>,
    ),
    ModularHecke(hecke::Report, Box<ParsedExpr>),
    HeckeEigen(eigenpacket::Statement),
    HeckeResidual(residual_hecke::Statement),
    HeckeQuadratic(quadratic_hecke::Query),
    PeriodPolynomial(period::Report),
    Hurwitz(Box<crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::hurwitz::Query>),
    Richelot(Box<crate::richelot::Command>),
    Bielliptic(Box<crate::bielliptic::Command>),
    EllipticGlue(Box<crate::elliptic_gluing::Command>),
    HWWClosed(Box<crate::hw_closed::Command>),
    Legendre(crate::legendre::Query),
    UnramRing(Label, usize, usize, Vec<i128>),
    UnramValuation(ParsedExpr),
    LegendreClosedReport(crate::legendre_closed::Query, ParsedExpr),
    LegendreCM(Label, ParsedExpr),
    CMReport(crate::legendre_cm::Query, Label, usize),
    QModule(Label, Vec<usize>, crate::q_module::Presentation),
    Let(Label, ParsedExpr),
    Show(ParsedExpr),
    RVInvariantFactors(usize, usize, ParsedExpr),
    RVNewton(usize, usize, ParsedExpr),
    RVHodge(ParsedExpr),
    TlLinking(usize, Vec<i128>, Vec<i128>),
    RVIw(bool, ParsedExpr),
    RVIwKL(Vec<i128>, Vec<i128>),
    HeckeCoefficient(bool, Option<i128>, Vec<i128>, ParsedExpr),
    SatakeQuery(Option<i128>, Vec<i128>, ParsedExpr),
    HiveReport(HiveSpec, bool),
    HiveVertex(usize, usize, ParsedExpr),
    QuantumCoefficient(usize, Vec<usize>, ParsedExpr),
    GromovWitten(usize, usize, usize, [Vec<usize>; 3]),
    OreCoefficient(usize, usize, ParsedExpr),
    WittSpectrumQuery(crate::rv_witt_spectrum::Query, Box<ParsedExpr>),
    WzCheck(Option<WzMirror>, crate::rv_wz::Certificate, ParsedExpr),
    QuantumPeriod(usize, usize, usize),
    QuantumRegularized(usize, usize, usize),
    QuantumRegularizedRhs(usize, usize),
    QuantumCurvature(crate::rv_pcurv::Report, usize, usize, usize),
    FusionCoefficient(Vec<i128>, ParsedExpr),
    CycloCoefficient(usize, ParsedExpr),
    ConformalDimension(usize, usize, usize, Vec<Vec<i128>>),
    ConformalSeries(crate::rv_conformal::Report, usize, usize, Vec<Vec<i128>>),
    ConformalModulo(usize, usize, u128, u128, Vec<Vec<i128>>),
    ConformalExact(usize, usize, u128, Vec<Vec<i128>>),
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
    RealRhoDimension(real_rho::Grade),
    RhoOrder(Box<ParsedExpr>),
    RhoDepth(Box<ParsedExpr>),
    RhoHidden(Box<[ParsedExpr; 2]>),
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
    RvPHodge(Label, PHodgeSpec),
    RvPHodgeQuad(Label, Label, QuadraticSpec<ParsedExpr>),
    RvQpQuad(Label, usize, i128, i128),
    RvQpQuadReport(Option<usize>, ParsedExpr),
    ResReport(crate::rv_res_space::Query, ParsedExpr),
    ResLeq(ParsedExpr, ParsedExpr),
    ValMatReport(crate::tropical_matroid::Query, ParsedExpr),
    HyperArithReport(crate::hyperelliptic_graph::arith::Query, ParsedExpr),
    HyperReport(crate::hyperelliptic_graph::Query, ParsedExpr),
    StableReport(crate::stable_marks::Query, ParsedExpr),
    RvPHodgeQuery(Label, PHodgeQuery),
    RvPeriod(Label, usize, Vec<crate::rv_period::Literal>),
    RVCTExact(Label, u128),
    RvAsOdd(Label, usize, Vec<(usize, i128)>),
    RvAsOddQuery(AsOddQuery, Label),
    RvAsTwo(Label, usize, Vec<(i128, i128)>),
    RvAsTwoQuery(AsTwoQuery, Label),
    EOReport(EoReport, ParsedExpr),
    EOWitness(ParsedExpr, ParsedExpr),
    RvCrystal(Label, usize, usize, CrystalEntries),
    RvCrystalQuery(Label, CrystalQuery),
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

    fn wide_unsigned(&mut self, label: &str) -> Result<u128> {
        let span = self.span();
        let Token::Number(text) = self.token() else {
            return Err(self.expected("a nonnegative decimal integer"));
        };
        let value = text
            .parse()
            .map_err(|_| Error::new(span, format!("{label} exceeds u128")))?;
        self.advance()?;
        Ok(value)
    }

    fn cartier_level(&mut self) -> Result<(usize, usize)> {
        let prime = self.number()?.0;
        self.expect(Token::Comma)?;
        Ok((prime, self.number()?.0))
    }

    fn rho_grade(&mut self) -> Result<rho_bock::Grade> {
        let (r, _) = self.number()?;
        self.expect(Token::Comma)?;
        let chain = self.motivic_grade()?;
        self.expect(Token::Comma)?;
        let (f, _) = self.number()?;
        Ok(rho_bock::Grade { chain, r, f })
    }

    fn real_rho_grade(&mut self) -> Result<real_rho::Grade> {
        let (exponent, at) = self.number()?;
        self.expect(Token::Comma)?;
        let chain = self.motivic_grade()?;
        real_rho::Grade::new(exponent, chain).map_err(|error| error.at(at))
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

    fn conformal_size(&mut self, label: &str) -> Result<u128> {
        let Token::Number(text) = self.token() else {
            return Err(self.expected("a nonnegative decimal integer"));
        };
        let value = text.parse().map_err(|_| {
            Error::new(
                self.span(),
                format!("conformal {label} exceeds the u128 range"),
            )
        })?;
        self.advance()?;
        Ok(value)
    }

    fn conformal_punctures(&mut self) -> Result<Vec<Vec<i128>>> {
        self.expect(Token::LBracket)?;
        let mut punctures = Vec::new();
        while self.token() != Token::RBracket {
            punctures.push(self.signed_tuple("a signed fusion weight", "fusion weight")?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBracket)?;
        Ok(punctures)
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

    fn phodge_quad_spec(&mut self) -> Result<QuadraticSpec<ParsedExpr>> {
        self.expect(Token::LBrace)?;
        self.expect(Token::LBrace)?;
        let mut eigenvalues = Vec::new();
        while self.token() != Token::RBrace {
            eigenvalues.push(self.expr(0)?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        self.expect(Token::Comma)?;
        self.expect(Token::LBrace)?;
        let mut monodromy = Vec::new();
        while self.token() != Token::RBrace {
            self.expect(Token::LParen)?;
            let (i, j) = self.bidegree()?;
            self.expect(Token::Comma)?;
            monodromy.push((i, j, self.expr(0)?));
            self.expect(Token::RParen)?;
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        self.expect(Token::Comma)?;
        self.expect(Token::LBrace)?;
        let mut filtration = Vec::new();
        while self.token() != Token::RBrace {
            self.expect(Token::LParen)?;
            let weight = self.signed_number("Hodge weight")?;
            self.expect(Token::Comma)?;
            self.expect(Token::LBracket)?;
            let mut column = Vec::new();
            while self.token() != Token::RBracket {
                column.push(self.expr(0)?);
                if !self.eat(Token::Comma)? {
                    break;
                }
            }
            self.expect(Token::RBracket)?;
            self.expect(Token::RParen)?;
            filtration.push((weight, column));
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        self.expect(Token::RBrace)?;
        Ok(QuadraticSpec {
            eigenvalues,
            monodromy,
            filtration,
        })
    }

    fn phodge_rational(&mut self) -> Result<(i128, i128)> {
        let numerator = self.signed_wide("an integer numerator", "filtered module numerator")?;
        self.expect(Token::Comma)?;
        let denominator =
            self.signed_wide("a positive denominator", "filtered module denominator")?;
        Ok((numerator, denominator))
    }

    fn phodge_spec(&mut self, prime: usize) -> Result<PHodgeSpec> {
        self.expect(Token::LBrace)?;
        self.expect(Token::LBrace)?;
        let mut eigenvalues = Vec::new();
        while self.token() != Token::RBrace {
            self.expect(Token::LParen)?;
            eigenvalues.push(self.phodge_rational()?);
            self.expect(Token::RParen)?;
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        self.expect(Token::Comma)?;
        self.expect(Token::LBrace)?;
        let mut monodromy = Vec::new();
        while self.token() != Token::RBrace {
            self.expect(Token::LParen)?;
            let (row, column) = self.bidegree()?;
            self.expect(Token::Comma)?;
            let (a, b) = self.phodge_rational()?;
            self.expect(Token::RParen)?;
            monodromy.push((row, column, a, b));
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        self.expect(Token::Comma)?;
        self.expect(Token::LBrace)?;
        let mut filtration = Vec::new();
        while self.token() != Token::RBrace {
            self.expect(Token::LParen)?;
            let weight = self.signed_number("Hodge weight")?;
            self.expect(Token::Comma)?;
            self.expect(Token::LBracket)?;
            let mut column = Vec::new();
            while self.token() != Token::RBracket {
                column
                    .push(self.signed_wide("an integer column entry", "filtration column entry")?);
                if !self.eat(Token::Comma)? {
                    break;
                }
            }
            self.expect(Token::RBracket)?;
            self.expect(Token::RParen)?;
            filtration.push((weight, column));
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        self.expect(Token::RBrace)?;
        Ok(PHodgeSpec {
            prime,
            eigenvalues,
            monodromy,
            filtration,
        })
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

    fn set_named_rows<T>(&mut self, mut row: impl FnMut(&mut Self) -> Result<T>) -> Result<Vec<T>> {
        self.expect(Token::LBrace)?;
        let mut rows = Vec::new();
        while self.token() != Token::RBrace {
            rows.push(row(self)?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(rows)
    }

    fn legendre_query(&mut self, name: &str) -> Result<crate::legendre::Query> {
        use crate::legendre::Query;
        self.expect(Token::LBracket)?;
        let prime = self.number()?.0;
        self.expect(Token::Comma)?;
        let query = match name {
            "legendreTrace" => Query::Trace {
                prime,
                lambda: self.number()?.0,
            },
            "legendreHensel" | "legendreDwork" => {
                let lambda = self.number()?.0;
                self.expect(Token::Comma)?;
                let precision = self.number()?.0;
                if name == "legendreDwork" {
                    Query::Dwork {
                        prime,
                        lambda,
                        precision,
                    }
                } else {
                    Query::Hensel {
                        prime,
                        lambda,
                        precision,
                    }
                }
            }
            "dworkPeriod" => {
                let z = self.signed_wide("a signed Dwork argument", "Dwork argument")?;
                self.expect(Token::Comma)?;
                let order = self.number()?.0;
                self.expect(Token::Comma)?;
                let precision = self.number()?.0;
                Query::Period {
                    prime,
                    z,
                    order,
                    precision,
                }
            }
            "dworkCongruence" => {
                let precision = self.number()?.0;
                let z = if self.eat(Token::Comma)? {
                    Some(self.signed_wide("a signed Dwork argument", "Dwork argument")?)
                } else {
                    None
                };
                Query::Congruence {
                    prime,
                    precision,
                    z,
                }
            }
            _ => unreachable!(),
        };
        self.expect(Token::RBracket)?;
        Ok(query)
    }

    fn set_topology_command(&mut self, command: &str) -> Result<crate::set_topology::Command> {
        use crate::set_topology::Command;
        self.advance()?;
        let name = match command {
            "settopology" | "settopologyMap" | "settopologymeet" | "settopologyjoin"
            | "setreflect" => Some(self.label()?),
            _ => None,
        };
        if command == "settopologyMap" {
            return Ok(Command::FromMap(
                name.expect("topology name"),
                Box::new(self.paren_expr(0)?),
            ));
        }
        self.expect(Token::LBracket)?;
        let first = self.label()?;
        let parsed = match command {
            "settopology" => {
                self.expect(Token::RBracket)?;
                let points = self.set_named_rows(Self::label)?;
                return Ok(Command::Subset(name.expect("topology name"), first, points));
            }
            "settopologysupport" => Command::Support(first),
            _ => {
                self.expect(Token::Comma)?;
                let second = self.label()?;
                match command {
                    "settopologymeet" | "settopologyjoin" => Command::Lattice(
                        name.expect("topology name"),
                        first,
                        second,
                        command == "settopologymeet",
                    ),
                    "settopologyjoinsteps" => Command::JoinSteps(first, second),
                    "setreflect" => Command::Reflect(name.expect("reflection name"), first, second),
                    "setreflectcheck" => Command::Check(first, second),
                    _ => unreachable!(),
                }
            }
        };
        self.expect(Token::RBracket)?;
        Ok(parsed)
    }

    fn set_topology_expression(
        &mut self,
        operation: crate::set_topology::Operation,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        self.advance()?;
        self.expect(Token::LBracket)?;
        let name = self.label()?;
        self.expect(Token::RBracket)?;
        let (input, depth) = if matches!(operation, crate::set_topology::Operation::Extend) {
            let input = self.paren_expr(nest)?;
            let depth = input.depth + 1;
            (Some(Box::new(input)), depth)
        } else {
            (None, 0)
        };
        expr(span, ExprForm::SetTopology(operation, name, input), depth)
    }

    fn sheaf_kuranishi_statement(&mut self, name: &str) -> Result<StatementForm> {
        self.advance()?;
        if name == "shkuranishi" {
            let name = self.label()?;
            self.expect(Token::LBracket)?;
            let sheaf = self.label()?;
            self.expect(Token::Comma)?;
            let order = self.number()?.0;
            self.expect(Token::RBracket)?;
            let directions = if self.token() == Token::LBrace {
                Some(self.set_named_rows(|parser| parser.expr(0))?)
            } else {
                None
            };
            return Ok(StatementForm::SheafKuranishi(
                name, sheaf, order, directions,
            ));
        }
        self.expect(Token::LBracket)?;
        let chart = self.label()?;
        let statement = if name == "shkuranishiEq" {
            self.expect(Token::Comma)?;
            let index = self.number()?.0;
            StatementForm::SheafKuranishiEquation(chart, index)
        } else {
            StatementForm::SheafKuranishiDimensions(chart)
        };
        self.expect(Token::RBracket)?;
        Ok(statement)
    }

    fn sheaf_kuranishi_expression(&mut self, name: &str, span: Span) -> Result<ParsedExpr> {
        use crate::sheaf_kuranishi::{Coefficient, Query};
        self.expect(Token::LBracket)?;
        let chart = self.label()?;
        self.expect(Token::Comma)?;
        let query = if matches!(name, "shKurArc" | "shKurArcKappa") {
            let amplitudes = self.sheaf_matrix()?;
            self.expect(Token::Comma)?;
            let weights = self.exponents()?;
            self.expect(Token::Comma)?;
            let order = self.number()?.0;
            Query::Arc(chart, amplitudes, weights, order, name == "shKurArcKappa")
        } else {
            let coefficient = match name {
                "shKurA" => Coefficient::Lift,
                "shKurR" => Coefficient::Residual,
                _ => Coefficient::Obstruction,
            };
            Query::Coefficient(chart, self.exponents()?, coefficient)
        };
        self.expect(Token::RBracket)?;
        expr(span, ExprForm::SheafKuranishi(Box::new(query)), 0)
    }

    fn kh_grade(&mut self) -> Result<crate::khovanov::Grade> {
        self.expect(Token::Comma)?;
        let i = self.signed_number("Khovanov homological degree")?;
        self.expect(Token::Comma)?;
        let j = self.signed_number("Khovanov quantum degree")?;
        Ok(crate::khovanov::Grade { i, j })
    }

    fn kh_k(&mut self) -> Result<i64> {
        self.expect(Token::Comma)?;
        self.signed_number("Khovanov annular degree")
    }

    fn kh_page_indices(&mut self) -> Result<(usize, usize, crate::khovanov::Grade, i64)> {
        self.expect(Token::Comma)?;
        let p = self.number()?.0;
        self.expect(Token::Comma)?;
        let r = self.number()?.0;
        let grade = self.kh_grade()?;
        let k = self.kh_k()?;
        Ok((p, r, grade, k))
    }

    fn kh_statement(&mut self, command: &str) -> Result<StatementForm> {
        use crate::khovanov::Command;
        self.advance()?;
        let command = if command == "khbraid" {
            let name = self.label()?;
            self.expect(Token::LBracket)?;
            let strands = self.number()?.0;
            self.expect(Token::Comma)?;
            let word = self.signed_tuple("a signed braid index", "Khovanov braid index")?;
            Command::Braid(name, strands, word)
        } else {
            self.expect(Token::LBracket)?;
            let name = self.label()?;
            match command {
                "khGroup" => Command::Group(name, self.kh_grade()?),
                "khAnnularGroup" => Command::AnnularGroup(name, self.kh_grade()?, self.kh_k()?),
                "khAnnularEuler" => Command::AnnularEuler(name),
                "khAnnularCheck" => Command::AnnularCheck(name),
                "khActionCheck" => Command::ActionCheck(name),
                "khPageDim" => {
                    let (p, r, grade, k) = self.kh_page_indices()?;
                    Command::PageDimension(name, p, r, grade, k)
                }
                "khIrreps" => {
                    let grade = self.kh_grade()?;
                    Command::Irreps(name, grade.i, grade.j)
                }
                "khKappa" => Command::Kappa(name),
                _ => Command::Euler(name),
            }
        };
        self.expect(Token::RBracket)?;
        Ok(StatementForm::Khovanov(Box::new(command)))
    }

    fn lee_dt_parameter(&mut self) -> Result<crate::lee::dt::Parameter> {
        let span = self.span();
        let p = self.number()?.0;
        self.expect(Token::Comma)?;
        let q = self.number()?.0;
        crate::lee::dt::Parameter::new(p, q).map_err(|error| error.at(span))
    }

    fn lee_dt_command(&mut self, command: &str) -> Result<StatementForm> {
        use crate::lee::dt::Command;
        self.advance()?;
        self.expect(Token::LBracket)?;
        let query = match command {
            "leeDtWeight" => {
                let t = self.lee_dt_parameter()?;
                self.expect(Token::RBracket)?;
                Command::Weight(t, Box::new(self.paren_expr(0)?))
            }
            _ => {
                let braid = self.label()?;
                let result = if command == "leeDt" {
                    self.expect(Token::Comma)?;
                    Command::Value(braid, self.lee_dt_parameter()?)
                } else {
                    Command::Profile(braid)
                };
                self.expect(Token::RBracket)?;
                result
            }
        };
        Ok(StatementForm::LeeDt(Box::new(query)))
    }

    fn lee_dt_expression(&mut self, command: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::lee::dt::{Query, Witness};
        let (query, depth) = match command {
            "leeDtDualD" => {
                self.expect(Token::LParen)?;
                let input = self.expr(nest)?;
                self.expect(Token::RParen)?;
                let depth = 1 + input.depth;
                (Query::DualD(Box::new(input)), depth)
            }
            "leeDtPair" => {
                self.expect(Token::LParen)?;
                let left = self.expr(nest)?;
                self.expect(Token::Comma)?;
                let right = self.expr(nest)?;
                self.expect(Token::RParen)?;
                let depth = 1 + left.depth.max(right.depth);
                (Query::Pair(Box::new([left, right])), depth)
            }
            _ => {
                self.expect(Token::LBracket)?;
                let braid = self.label()?;
                let operation = if command == "leeOrientation" {
                    Query::Orientation(braid)
                } else {
                    self.expect(Token::Comma)?;
                    let t = self.lee_dt_parameter()?;
                    let mode = match command {
                        "leeDtRep" => Witness::Representative,
                        "leeDtNull" => Witness::Primitive,
                        _ => Witness::Dual,
                    };
                    Query::Witness(braid, t, mode)
                };
                self.expect(Token::RBracket)?;
                (operation, 0)
            }
        };
        expr(span, ExprForm::LeeDt(Box::new(query)), depth)
    }

    fn lee_grade(&mut self) -> Result<crate::lee::Grade> {
        let i = self.signed_number("Lee homological degree")?;
        self.expect(Token::Comma)?;
        let f = self.signed_number("Lee filtration")?;
        Ok(crate::lee::Grade { i, f })
    }

    fn lee_command(&mut self, command: &str) -> Result<StatementForm> {
        use crate::lee::Command;
        let command = match command {
            "leeFiltration" => Command::Filtration(Box::new(self.paren_expr(0)?)),
            "leeIndependent" => {
                self.expect(Token::LParen)?;
                let mut inputs = Vec::new();
                while self.token() != Token::RParen {
                    inputs.push(self.expr(1)?);
                    if !self.eat(Token::Comma)? {
                        break;
                    }
                }
                self.expect(Token::RParen)?;
                Command::Independent(inputs)
            }
            _ => {
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                let operation = match command {
                    "leeS" => Command::Rasmussen(name),
                    "leeLevels" => {
                        self.expect(Token::Comma)?;
                        Command::Levels(name, self.signed_number("Lee homological degree")?)
                    }
                    _ => {
                        self.expect(Token::Comma)?;
                        Command::Rank(name, self.lee_grade()?)
                    }
                };
                self.expect(Token::RBracket)?;
                operation
            }
        };
        Ok(StatementForm::Lee(command))
    }

    fn lee_expression(&mut self, command: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::lee::{Operation, Query};
        let (query, depth) = if command == "leeRep" {
            self.expect(Token::LBracket)?;
            let name = self.label()?;
            self.expect(Token::Comma)?;
            let grade = self.lee_grade()?;
            self.expect(Token::Comma)?;
            let index = self.number()?.0;
            self.expect(Token::RBracket)?;
            (Query::Representative(name, grade, index), 0)
        } else {
            let operation = match command {
                "leeLift" => Operation::Lift,
                "leePhi" => Operation::Perturbation,
                "leeD" => Operation::Differential,
                "leeNull" => Operation::Null,
                "leeFilter" => {
                    self.expect(Token::LBracket)?;
                    let f = self.signed_number("Lee filtration")?;
                    self.expect(Token::RBracket)?;
                    Operation::Filter(f)
                }
                _ => {
                    self.expect(Token::LBracket)?;
                    let numerator =
                        self.signed_wide("a signed integer numerator", "Lee scalar numerator")?;
                    self.expect(Token::Comma)?;
                    let denominator = self
                        .signed_wide("a positive integer denominator", "Lee scalar denominator")?;
                    self.expect(Token::RBracket)?;
                    Operation::Scale(crate::lee_linear::Rational::new(numerator, denominator)?)
                }
            };
            let input = self.paren_expr(nest)?;
            let depth = input.depth + 1;
            (Query::Operation(operation, Box::new(input)), depth)
        };
        expr(span, ExprForm::Lee(query), depth)
    }

    fn kh_expression(&mut self, command: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::khovanov::{
            page::{Map, Operation as Page},
            Action, Generator, Operation, Query,
        };
        let operation = match command {
            "khD" => Some(Operation::Differential),
            "khMod" | "khPageClass" => {
                self.expect(Token::LBracket)?;
                let index = self.number()?.0;
                let k = self.kh_k()?;
                self.expect(Token::RBracket)?;
                Some(Operation::Page(if command == "khMod" {
                    Page::Reduce(index, k)
                } else {
                    Page::Class(index, k)
                }))
            }
            "khPageRep" => Some(Operation::Page(Page::Representative)),
            "khPageD" => Some(Operation::Page(Page::Map(Map::Differential))),
            "khPageNext" => Some(Operation::Page(Page::Map(Map::Next))),
            "khPagePrev" => Some(Operation::Page(Page::Map(Map::Previous))),
            "khPageBoundary" => Some(Operation::Page(Page::Map(Map::Boundary))),
            "khE" => Some(Operation::Action(Action::E)),
            "khF" => Some(Operation::Action(Action::F)),
            "khH" => Some(Operation::Action(Action::H)),
            "khVMinus" => Some(Operation::Action(Action::Minus)),
            "khVZero" => Some(Operation::Action(Action::Zero)),
            "khWPlus" => Some(Operation::Action(Action::PlusTwice)),
            "khLee0" => Some(Operation::Action(Action::LeeZero)),
            "khLeePlus" => Some(Operation::Action(Action::LeePlus)),
            "khD0" => Some(Operation::AnnularDifferential(0)),
            "khDMinus2" => Some(Operation::AnnularDifferential(-2)),
            "khForgetK" => Some(Operation::ForgetK),
            "khClass" => Some(Operation::Project),
            "khRep" => Some(Operation::Representative),
            "khScale" => {
                self.expect(Token::LBracket)?;
                let scale = self.signed_wide("an integer scalar", "Khovanov scalar")?;
                self.expect(Token::RBracket)?;
                Some(Operation::Scale(scale))
            }
            "khAnnular" | "khToF2" => {
                self.expect(Token::LBracket)?;
                let k = self.signed_number("Khovanov annular degree")?;
                self.expect(Token::RBracket)?;
                Some(if command == "khAnnular" {
                    Operation::AnnularPart(k)
                } else {
                    Operation::ToF2(k)
                })
            }
            _ => None,
        };
        let (query, depth) = if let Some(operation) = operation {
            self.expect(Token::LParen)?;
            let input = self.expr(nest)?;
            self.expect(Token::RParen)?;
            let depth = 1 + input.depth;
            (Query::Operation(operation, Box::new(input)), depth)
        } else if command == "khPair" {
            self.expect(Token::LParen)?;
            let left = self.expr(nest)?;
            self.expect(Token::Comma)?;
            let right = self.expr(nest)?;
            self.expect(Token::RParen)?;
            let depth = 1 + left.depth.max(right.depth);
            (Query::Pair(Box::new([left, right])), depth)
        } else {
            self.expect(Token::LBracket)?;
            let name = self.label()?;
            let query = match command {
                "khTerm" => {
                    self.expect(Token::Comma)?;
                    let ones = self.exponents()?;
                    self.expect(Token::Comma)?;
                    let xs = self.exponents()?;
                    let coefficient = if self.eat(Token::Comma)? {
                        self.signed_wide("an integer coefficient", "Khovanov coefficient")?
                    } else {
                        1
                    };
                    Query::Term(name, ones, xs, coefficient)
                }
                "khEulerAt" => {
                    self.expect(Token::Comma)?;
                    Query::EulerAt(name, self.number()?.0)
                }
                "khPsi" => Query::Psi(name),
                "khKappaNull" => Query::KappaNull(name),
                "khKappaDual" => Query::KappaDual(name),
                "khPageGen" => {
                    let (p, r, grade, k) = self.kh_page_indices()?;
                    self.expect(Token::Comma)?;
                    Query::PageGenerator(name, p, r, grade, k, self.number()?.0)
                }
                "khHighest" => {
                    let grade = self.kh_grade()?;
                    let weight = self.kh_k()?;
                    self.expect(Token::Comma)?;
                    Query::Highest(name, grade, weight, self.number()?.0)
                }
                _ => {
                    let grade = self.kh_grade()?;
                    let annular = matches!(command, "khABasis" | "khAGen" | "khANull")
                        .then(|| self.kh_k())
                        .transpose()?;
                    self.expect(Token::Comma)?;
                    let index = self.number()?.0;
                    let kind = match command {
                        "khBasis" | "khABasis" => Generator::Basis,
                        "khGen" | "khAGen" => Generator::Class,
                        _ => Generator::Null,
                    };
                    Query::Generator(name, grade, annular, index, kind)
                }
            };
            self.expect(Token::RBracket)?;
            (query, 0)
        };
        expr(span, ExprForm::Khovanov(Box::new(query)), depth)
    }

    fn hw_shape(&mut self) -> Result<(usize, usize, usize)> {
        let prime = self.number()?.0;
        self.expect(Token::Comma)?;
        let genus = self.number()?.0;
        self.expect(Token::Comma)?;
        let precision = self.number()?.0;
        Ok((prime, genus, precision))
    }

    fn hw_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic::{Command, Report};
        self.advance()?;
        if name == "rvHWVal" {
            return Ok(StatementForm::HasseWitt(Box::new(Command::Valuation(
                Box::new(self.paren_expr(0)?),
            ))));
        }
        if name == "rvhyperelliptic" {
            let curve = self.label()?;
            self.expect(Token::LBracket)?;
            if matches!(self.token(), Token::Word(_)) {
                let ring = self.label()?;
                self.expect(Token::RBracket)?;
                self.expect(Token::LBrace)?;
                let mut coefficients = Vec::new();
                while !self.eat(Token::RBrace)? {
                    coefficients.push(self.expr(0)?);
                    if !self.eat(Token::Comma)? {
                        self.expect(Token::RBrace)?;
                        break;
                    }
                }
                return Ok(StatementForm::HWWClosed(Box::new(
                    crate::hw_closed::Command::Define(curve, ring, coefficients),
                )));
            }
            let p = self.number()?.0;
            self.expect(Token::RBracket)?;
            self.expect(Token::LBrace)?;
            let coefficients = self.signed_list(
                Token::RBrace,
                "an integer lift coefficient",
                "hyperelliptic lift coefficient",
            )?;
            return Ok(StatementForm::HasseWitt(Box::new(Command::Define(
                curve,
                p,
                coefficients,
            ))));
        }
        self.expect(Token::LBracket)?;
        let curve = self.label()?;
        let report = match name {
            "rvHWStatus" => Report::Status,
            "rvHWWeil" => Report::Weil,
            _ => {
                self.expect(Token::Comma)?;
                let index = self.number()?.0;
                match name {
                    "rvHWPoints" => Report::Points(index),
                    "rvHWConverges" => Report::Convergence(index),
                    _ => Report::Certificate(index),
                }
            }
        };
        self.expect(Token::RBracket)?;
        Ok(StatementForm::HasseWitt(Box::new(Command::Report(
            curve, report,
        ))))
    }

    fn hw_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic::{Expression, Unary};
        let mut depth = 0;
        let expression = match name {
            "rvHWCompose" => {
                let pair = self.paren_pair(nest)?;
                depth = 1 + pair[0].depth.max(pair[1].depth);
                Expression::Compose(Box::new(pair))
            }
            "rvHWInv" | "rvHWEntry" | "rvHWChar" | "rvHWReduce" => {
                let operation = if name == "rvHWInv" {
                    Unary::Invert
                } else {
                    self.expect(Token::LBracket)?;
                    let i = self.number()?.0;
                    let operation = match name {
                        "rvHWEntry" => {
                            self.expect(Token::Comma)?;
                            Unary::Entry(i, self.number()?.0)
                        }
                        "rvHWChar" => Unary::Coefficient(i),
                        _ => Unary::Reduce(i),
                    };
                    self.expect(Token::RBracket)?;
                    operation
                };
                let input = self.paren_expr(nest)?;
                depth = input.depth + 1;
                Expression::Unary(operation, Box::new(input))
            }
            "rvHWMat" => {
                self.expect(Token::LBracket)?;
                let (p, g, s) = self.hw_shape()?;
                self.expect(Token::RBracket)?;
                self.expect(Token::LBrace)?;
                let entries = self.signed_list(
                    Token::RBrace,
                    "an integer matrix entry",
                    "Hasse–Witt matrix entry",
                )?;
                Expression::Matrix(p, g, s, entries)
            }
            "rvHWScalar" => {
                self.expect(Token::LBracket)?;
                let p = self.number()?.0;
                self.expect(Token::Comma)?;
                let s = self.number()?.0;
                self.expect(Token::Comma)?;
                let c = self.signed_wide("an integer scalar", "Hasse–Witt scalar")?;
                self.expect(Token::RBracket)?;
                Expression::Scalar(p, s, c)
            }
            _ => {
                self.expect(Token::LBracket)?;
                let curve = self.label()?;
                self.expect(Token::Comma)?;
                let index = self.number()?.0;
                let query = if name == "rvHWUnit" {
                    Expression::Unit(curve, index)
                } else {
                    self.expect(Token::Comma)?;
                    let precision = self.number()?.0;
                    Expression::Raw(curve, index, precision)
                };
                self.expect(Token::RBracket)?;
                query
            }
        };
        expr(span, ExprForm::HasseWitt(Box::new(expression)), depth)
    }

    fn jacobian_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_jacobian::Command;
        self.advance()?;
        let command = match name {
            "rvjacbase" => {
                let target = self.label()?;
                self.expect(Token::LBracket)?;
                let base = self.label()?;
                self.expect(Token::Comma)?;
                let field = self.label()?;
                self.expect(Token::RBracket)?;
                Command::Base(target, base, field)
            }
            "rvjacchart" => {
                let target = self.label()?;
                self.expect(Token::LBracket)?;
                let source = self.label()?;
                self.expect(Token::RBracket)?;
                Command::Chart(target, source, Box::new(self.paren_expr(0)?))
            }
            "rvJacCompose" => Command::Compose(Box::new(self.paren_pair(0)?)),
            "rvJacOrder" => Command::Order(Box::new(self.paren_expr(0)?)),
            "rvJacKummer" => Command::Kummer(Box::new(self.paren_expr(0)?)),
            "rvJacHalfFiber" => {
                self.expect(Token::LBracket)?;
                let curve = self.label()?;
                self.expect(Token::RBracket)?;
                Command::HalfFiber(curve, Box::new(self.paren_pair(0)?))
            }
            _ => {
                self.expect(Token::LBracket)?;
                let curve = self.label()?;
                self.expect(Token::RBracket)?;
                Command::Count(curve, name == "rvJacEnumerate")
            }
        };
        Ok(StatementForm::Jacobian(Box::new(command)))
    }

    fn jacobian_level_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::jacobian_level::Command;
        self.advance()?;
        let command = if name == "rvjactorsionbasis" {
            let basis = self.label()?;
            self.expect(Token::LBracket)?;
            let curve = self.label()?;
            self.expect(Token::Comma)?;
            let level = self.number()?.0;
            self.expect(Token::RBracket)?;
            let root = self.paren_expr(0)?;
            let generators = self.set_named_rows(|parser| parser.expr(0))?;
            Command::Define(basis, curve, level, Box::new(root), generators)
        } else {
            self.expect(Token::LBracket)?;
            let basis = self.label()?;
            self.expect(Token::Comma)?;
            let endpoint = self.label()?;
            if name == "rvJacRichelotMatrix" {
                self.expect(Token::RBracket)?;
                Command::RichelotMatrix(basis, endpoint)
            } else {
                self.expect(Token::Comma)?;
                let power = self.number()?.0;
                self.expect(Token::RBracket)?;
                Command::Frobenius(basis, endpoint, power)
            }
        };
        Ok(StatementForm::JacobianLevel(Box::new(command)))
    }

    fn abelian_mukai_expression(
        &mut self,
        operation: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::mukai::abelian::Expression;
        let (expression, depth) = if operation == "rvAbMukaiKummerH2" {
            let input = self.paren_expr(nest)?;
            let depth = 1 + input.depth;
            (Expression::Kummer(Box::new(input)), depth)
        } else {
            self.expect(Token::LBracket)?;
            let marking = self.label()?;
            self.expect(Token::Comma)?;
            if operation == "rvAbMukaiBasis" {
                let index = self.number()?.0;
                self.expect(Token::RBracket)?;
                (Expression::Basis(marking, index), 0)
            } else {
                let r = self.signed_wide("a signed rank", "abelian Mukai rank")?;
                self.expect(Token::Comma)?;
                let s = self.signed_wide("a signed integer", "abelian Mukai last coordinate")?;
                self.expect(Token::RBracket)?;
                self.expect(Token::LBrace)?;
                let wedge = self.signed_list(
                    Token::RBrace,
                    "a signed wedge coordinate",
                    "abelian Mukai wedge coordinate",
                )?;
                let wedge = wedge
                    .try_into()
                    .map_err(|_| Error::new(span, "rvAbMukai requires six wedge coordinates"))?;
                (Expression::Literal(marking, r, s, wedge), 0)
            }
        };
        expr(span, ExprForm::AbMukai(Box::new(expression)), depth)
    }

    fn ks_end_statement(&mut self, operation: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::projector::Command;
        self.advance()?;
        let command = match operation {
            "rvKSVolumeSquare" => {
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::RBracket)?;
                Command::VolumeSquare(name)
            }
            "rvKSEndCoeff" => {
                self.expect(Token::LBracket)?;
                let i = self.number()?.0;
                self.expect(Token::Comma)?;
                let j = self.number()?.0;
                self.expect(Token::RBracket)?;
                Command::Coefficient(i, j, Box::new(self.paren_expr(0)?))
            }
            "rvKSEndTrace" => Command::Trace(Box::new(self.paren_expr(0)?)),
            _ => Command::TracePair(Box::new(self.paren_pair(0)?)),
        };
        Ok(StatementForm::KSEnd(Box::new(command)))
    }

    fn ks_end_expression(
        &mut self,
        operation: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::projector::Expression;
        let (expression, depth) = match operation {
            "rvKSEndUnit" => {
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::Comma)?;
                let i = self.number()?.0;
                self.expect(Token::Comma)?;
                let j = self.number()?.0;
                self.expect(Token::RBracket)?;
                (Expression::Unit(name, i, j), 0)
            }
            "rvKSVolume" | "rvKSEmbed" => {
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::RBracket)?;
                if operation == "rvKSVolume" {
                    (Expression::Volume(name), 0)
                } else {
                    let input = self.paren_expr(nest)?;
                    let depth = 1 + input.depth;
                    (Expression::Embed(name, Box::new(input)), depth)
                }
            }
            "rvKSEndCompose" => {
                let operands = self.paren_pair(nest)?;
                let depth = 1 + operands[0].depth.max(operands[1].depth);
                (Expression::Compose(Box::new(operands)), depth)
            }
            "rvKSEndScale" => {
                self.expect(Token::LBracket)?;
                let numerator =
                    self.signed_wide("a signed numerator", "Kuga–Satake endomorphism numerator")?;
                self.expect(Token::Comma)?;
                let denominator = self.wide_unsigned("Kuga–Satake endomorphism denominator")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (
                    Expression::Scale(numerator, denominator, Box::new(input)),
                    depth,
                )
            }
            _ => {
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                let expression = if operation == "rvKSEndLeft" {
                    Expression::Left(Box::new(input))
                } else {
                    Expression::Project(Box::new(input))
                };
                (expression, depth)
            }
        };
        expr(span, ExprForm::KSEnd(Box::new(expression)), depth)
    }

    fn ks_exterior_mask(&mut self) -> Result<u64> {
        let span = self.span();
        self.signed_wide(
            "an unsigned 64-bit subset mask",
            "Kuga–Satake exterior mask",
        )?
        .try_into()
        .map_err(|_| {
            Error::new(
                span,
                "Kuga–Satake exterior mask must be an unsigned 64-bit integer",
            )
        })
    }

    fn ks_exterior_statement(&mut self, operation: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::exterior::Command;
        self.advance()?;
        self.expect(Token::LBracket)?;
        let command = if matches!(operation, "rvKSExteriorData" | "rvKSExteriorSmith") {
            let name = self.label()?;
            self.expect(Token::Comma)?;
            let degree = self.number()?.0;
            let command = if operation == "rvKSExteriorSmith" {
                self.expect(Token::Comma)?;
                Command::Smith(name, degree, self.number()?.0)
            } else {
                Command::Data(name, degree)
            };
            self.expect(Token::RBracket)?;
            command
        } else {
            let mask = self.ks_exterior_mask()?;
            self.expect(Token::RBracket)?;
            Command::Coefficient(mask, Box::new(self.paren_expr(0)?))
        };
        Ok(StatementForm::KSExterior(Box::new(command)))
    }

    fn ks_exterior_expression(
        &mut self,
        operation: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::exterior::{
            Expression, Operation,
        };
        let (expression, depth) = match operation {
            "rvKSExteriorCup" => {
                let operands = self.paren_pair(nest)?;
                let depth = 1 + operands[0].depth.max(operands[1].depth);
                (Expression::Cup(Box::new(operands)), depth)
            }
            "rvKSExteriorEmbed" | "rvKSExteriorToSpin" => {
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                let expression = if operation == "rvKSExteriorEmbed" {
                    Expression::Embed(Box::new(input))
                } else {
                    Expression::ToSpin(Box::new(input))
                };
                (expression, depth)
            }
            "rvKSExteriorScale" => {
                self.expect(Token::LBracket)?;
                let k = self.signed_wide("a signed integer", "Kuga–Satake exterior scalar")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::Scale(k, Box::new(input)), depth)
            }
            "rvKSExteriorBasis" => {
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::Comma)?;
                let mask = self.ks_exterior_mask()?;
                self.expect(Token::RBracket)?;
                (Expression::Basis(name, mask), 0)
            }
            _ => {
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::Comma)?;
                let action = match operation {
                    "rvKSExteriorPower" => Operation::Power(self.number()?.0),
                    "rvKSExteriorDivided" => Operation::Divided(self.number()?.0),
                    _ => Operation::Exp(
                        self.signed_wide("a signed integer", "Kuga–Satake exterior exponent")?,
                    ),
                };
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::Apply(name, action, Box::new(input)), depth)
            }
        };
        expr(span, ExprForm::KSExterior(Box::new(expression)), depth)
    }

    fn og6_spin_statement(&mut self, operation: &str, span: Span) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::Command;
        self.advance()?;
        self.expect(Token::LBracket)?;
        let command = if operation == "rvOG6SpinData" {
            let name = self.label()?;
            self.expect(Token::RBracket)?;
            Command::Data(name)
        } else {
            let mask = self.number()?.0.try_into().map_err(|_| {
                Error::new(span, "OG6 Clifford blade mask must be between 0 and 127")
            })?;
            self.expect(Token::RBracket)?;
            Command::Coefficient(mask, Box::new(self.paren_expr(0)?))
        };
        Ok(StatementForm::OG6Clifford(Box::new(command)))
    }

    fn og6_clifford_expression(
        &mut self,
        operation: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::Expression;
        let (expression, depth) = if operation == "rvOG6CliffordMul" {
            let operands = self.paren_pair(nest)?;
            let depth = 1 + operands[0].depth.max(operands[1].depth);
            (Expression::Multiply(Box::new(operands)), depth)
        } else if operation == "rvOG6CliffordScale" {
            self.expect(Token::LBracket)?;
            let k = self.signed_wide("a signed integer", "OG6 Clifford scalar")?;
            self.expect(Token::RBracket)?;
            let input = self.paren_expr(nest)?;
            let depth = 1 + input.depth;
            (Expression::Scale(k, Box::new(input)), depth)
        } else {
            self.expect(Token::LBracket)?;
            let name = self.label()?;
            match operation {
                "rvOG6CliffordBasis" | "rvOG6PolarizedBasis" => {
                    self.expect(Token::Comma)?;
                    let index = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let expression = if operation == "rvOG6CliffordBasis" {
                        let mask = index.try_into().map_err(|_| {
                            Error::new(span, "OG6 Clifford blade mask must be between 0 and 127")
                        })?;
                        Expression::Basis(name, mask)
                    } else {
                        Expression::PolarizedBasis(name, index)
                    };
                    (expression, 0)
                }
                "rvOG6SpinLog" => {
                    self.expect(Token::RBracket)?;
                    (Expression::Log(name), 0)
                }
                "rvOG6SpinExp" => {
                    self.expect(Token::Comma)?;
                    let m = self.signed_wide("a signed integer", "OG6 spin exponent")?;
                    self.expect(Token::RBracket)?;
                    (Expression::Exp(name, m), 0)
                }
                _ => {
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = 1 + input.depth;
                    (Expression::Vector(name, Box::new(input)), depth)
                }
            }
        };
        expr(span, ExprForm::OG6Clifford(Box::new(expression)), depth)
    }

    fn og6_nilpotent_statement(&mut self, operation: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::nilpotent::Command;
        self.advance()?;
        let command = if matches!(operation, "rvog6nilpotentii" | "rvog6nilpotentiii") {
            let name = self.label()?;
            self.expect(Token::LBracket)?;
            let cusp = self.label()?;
            self.expect(Token::Comma)?;
            if operation == "rvog6nilpotentii" {
                let k = self.signed_wide("a signed integer", "OG6 nilpotent coefficient")?;
                self.expect(Token::RBracket)?;
                Command::TypeII(name, cusp, k)
            } else {
                let index = self.number()?.0;
                self.expect(Token::RBracket)?;
                Command::TypeIII(name, cusp, index, Box::new(self.paren_expr(0)?))
            }
        } else {
            self.expect(Token::LBracket)?;
            let name = self.label()?;
            if operation == "rvOG6NilpotentData" {
                self.expect(Token::RBracket)?;
                Command::Data(name)
            } else {
                self.expect(Token::Comma)?;
                let weight = self.signed_wide("a signed integer", "OG6 filtration weight")?;
                self.expect(Token::RBracket)?;
                Command::InWeight(name, weight, Box::new(self.paren_expr(0)?))
            }
        };
        Ok(StatementForm::OG6Nilpotent(Box::new(command)))
    }

    fn og6_nilpotent_expression(
        &mut self,
        operation: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::nilpotent::{
            Expression, Operation,
        };
        self.expect(Token::LBracket)?;
        let name = self.label()?;
        self.expect(Token::Comma)?;
        let operation = if operation == "rvOG6NilpotentPower" {
            Operation::Power(self.number()?.0)
        } else {
            Operation::Exp(self.signed_wide("a signed integer", "OG6 nilpotent exponent")?)
        };
        self.expect(Token::RBracket)?;
        let input = self.paren_expr(nest)?;
        let depth = 1 + input.depth;
        expr(
            span,
            ExprForm::OG6Nilpotent(Box::new(Expression(name, operation, Box::new(input)))),
            depth,
        )
    }

    fn og6_cusp_statement(&mut self, operation: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::Command;
        self.advance()?;
        let command = if operation == "rvog6cusp" {
            let name = self.label()?;
            self.expect(Token::LParen)?;
            let h = self.expr(1)?;
            self.expect(Token::Comma)?;
            let e = self.expr(1)?;
            self.expect(Token::Comma)?;
            let f = self.expr(1)?;
            self.expect(Token::RParen)?;
            Command::Define(name, Box::new([h, e, f]))
        } else {
            self.expect(Token::LBracket)?;
            let name = self.label()?;
            self.expect(Token::RBracket)?;
            Command::Gram(name)
        };
        Ok(StatementForm::OG6Cusp(Box::new(command)))
    }

    fn og6_cusp_expression(
        &mut self,
        operation: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::Expression;
        self.expect(Token::LBracket)?;
        let name = self.label()?;
        self.expect(Token::Comma)?;
        let (expression, depth) = match operation {
            "rvOG6CuspBasis" => {
                let index = self.number()?.0;
                self.expect(Token::RBracket)?;
                (Expression::Basis(name, index), 0)
            }
            "rvOG6CuspTransvection" => {
                let index = self.number()?.0;
                if index > 1 {
                    return Err(Error::new(
                        span,
                        "rvOG6CuspTransvection axis must be 0 or 1",
                    ));
                }
                self.expect(Token::RBracket)?;
                let inputs = self.paren_pair(nest)?;
                let depth = 1 + inputs[0].depth.max(inputs[1].depth);
                (
                    Expression::Transvection(name, index, Box::new(inputs)),
                    depth,
                )
            }
            _ => {
                let k = self.signed_wide("a signed integer", "OG6 cusp central coefficient")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::Central(name, k, Box::new(input)), depth)
            }
        };
        expr(span, ExprForm::OG6Cusp(Box::new(expression)), depth)
    }

    fn og6_statement(&mut self, operation: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::Command;
        self.advance()?;
        let command = match operation {
            "rvog6" => Command::Define(self.label()?, Box::new(self.paren_expr(0)?)),
            "rvOG6Lattice" => {
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::RBracket)?;
                Command::Lattice(name)
            }
            "rvOG6Pair" => Command::Pair(Box::new(self.paren_pair(0)?)),
            "rvOG6SameOrbit" => Command::SameOrbit(Box::new(self.paren_pair(0)?)),
            "rvOG6Orbit" => Command::Orbit(Box::new(self.paren_expr(0)?)),
            "rvOG6Wall" => Command::Wall(Box::new(self.paren_expr(0)?)),
            "rvOG6Fujiki" => {
                self.expect(Token::LBracket)?;
                let exponent = self.number()?.0;
                self.expect(Token::RBracket)?;
                Command::Fujiki(exponent, Box::new(self.paren_pair(0)?))
            }
            "rvOG6Div" => Command::Divisibility(Box::new(self.paren_expr(0)?)),
            _ => Command::RiemannRoch(Box::new(self.paren_expr(0)?)),
        };
        Ok(StatementForm::OG6(Box::new(command)))
    }

    fn og6_expression(&mut self, operation: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::Expression;
        if operation == "rvOG6Reflect" {
            let inputs = self.paren_pair(nest)?;
            let depth = 1 + inputs[0].depth.max(inputs[1].depth);
            return expr(
                span,
                ExprForm::OG6(Box::new(Expression::Reflect(Box::new(inputs)))),
                depth,
            );
        }
        self.expect(Token::LBracket)?;
        let name = self.label()?;
        let (expression, depth) = match operation {
            "rvOG6Basis" => {
                self.expect(Token::Comma)?;
                let index = self.number()?.0;
                self.expect(Token::RBracket)?;
                (Expression::Basis(name, index), 0)
            }
            "rvOG6Exc" | "rvOG6ExcHalf" => {
                self.expect(Token::RBracket)?;
                (
                    Expression::Exceptional(name, operation == "rvOG6ExcHalf"),
                    0,
                )
            }
            _ => {
                self.expect(Token::Comma)?;
                let k = self.signed_wide("a signed integer", "OG6 half-exceptional coefficient")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::Build(name, k, Box::new(input)), depth)
            }
        };
        expr(span, ExprForm::OG6(Box::new(expression)), depth)
    }

    fn og10_statement(&mut self, operation: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::mukai::moduli::og10::Command;
        self.advance()?;
        let command = match operation {
            "rvOG10Lattice" => {
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::RBracket)?;
                Command::Lattice(name)
            }
            "rvOG10Pair" => Command::Pair(Box::new(self.paren_pair(0)?)),
            "rvOG10Wall" => Command::Wall(Box::new(self.paren_expr(0)?)),
            "rvOG10Walls" => Command::Walls(Box::new(self.paren_pair(0)?)),
            "rvOG10Fujiki" => {
                self.expect(Token::LBracket)?;
                let exponent = self.number()?.0;
                self.expect(Token::RBracket)?;
                Command::Fujiki(exponent, Box::new(self.paren_pair(0)?))
            }
            "rvOG10Div" => Command::Divisibility(Box::new(self.paren_expr(0)?)),
            _ => Command::RiemannRoch(Box::new(self.paren_expr(0)?)),
        };
        Ok(StatementForm::OG10(Box::new(command)))
    }

    fn og10_expression(&mut self, operation: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::mukai::moduli::og10::Expression;
        if matches!(operation, "rvOG10Reflect" | "rvOG10WallAt") {
            let index = if operation == "rvOG10WallAt" {
                self.expect(Token::LBracket)?;
                let index = self.number()?.0;
                self.expect(Token::RBracket)?;
                Some(index)
            } else {
                None
            };
            let inputs = self.paren_pair(nest)?;
            let depth = 1 + inputs[0].depth.max(inputs[1].depth);
            let inputs = Box::new(inputs);
            let expression = match index {
                Some(index) => Expression::WallAt(index, inputs),
                None => Expression::Reflect(inputs),
            };
            return expr(span, ExprForm::OG10(Box::new(expression)), depth);
        }
        self.expect(Token::LBracket)?;
        let name = self.label()?;
        let (expression, depth) = match operation {
            "rvOG10Basis" => {
                self.expect(Token::Comma)?;
                let index = self.number()?.0;
                self.expect(Token::RBracket)?;
                (Expression::Basis(name, index), 0)
            }
            "rvOG10Exc" => {
                self.expect(Token::RBracket)?;
                (Expression::Exceptional(name), 0)
            }
            _ => {
                self.expect(Token::Comma)?;
                let generator = self.signed_wide("a signed integer", "OG10 gluing coefficient")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::Build(name, generator, Box::new(input)), depth)
            }
        };
        expr(span, ExprForm::OG10(Box::new(expression)), depth)
    }

    fn mukai_moduli_statement(&mut self, operation: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::mukai::moduli::Command;
        self.advance()?;
        if operation == "rvmukaimoduli" {
            let name = self.label()?;
            let input = self.paren_expr(0)?;
            return Ok(StatementForm::MukaiModuli(Box::new(Command::Define(
                name,
                Box::new(input),
            ))));
        }
        if operation == "rvMukaiQuotientPair" {
            return Ok(StatementForm::MukaiModuli(Box::new(Command::QuotientPair(
                Box::new(self.paren_pair(0)?),
            ))));
        }
        self.expect(Token::LBracket)?;
        let name = self.label()?;
        let command = if operation == "rvMukaiModuliFujiki" {
            self.expect(Token::Comma)?;
            let exponent = self.number()?.0;
            self.expect(Token::RBracket)?;
            Command::Fujiki(name, exponent, Box::new(self.paren_pair(0)?))
        } else {
            self.expect(Token::RBracket)?;
            match operation {
                "rvMukaiModuliLattice" => Command::Lattice(name),
                "rvMukaiModuliDiv" => Command::Divisibility(name, Box::new(self.paren_expr(0)?)),
                _ => Command::RiemannRoch(name, Box::new(self.paren_expr(0)?)),
            }
        };
        Ok(StatementForm::MukaiModuli(Box::new(command)))
    }

    fn mukai_moduli_expression(
        &mut self,
        operation: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::mukai::moduli::Expression;
        if operation == "rvMukaiQuotientLift" {
            let input = self.paren_expr(nest)?;
            let depth = 1 + input.depth;
            return expr(
                span,
                ExprForm::MukaiModuli(Expression::Lift(Box::new(input))),
                depth,
            );
        }
        self.expect(Token::LBracket)?;
        let name = self.label()?;
        let index = if operation == "rvMukaiModuliBasis" {
            self.expect(Token::Comma)?;
            Some(self.number()?.0)
        } else {
            None
        };
        self.expect(Token::RBracket)?;
        let (expression, depth) = match operation {
            "rvMukaiModuliBasis" => (Expression::Basis(name, index.unwrap()), 0),
            "rvMukaiModuliDual" => (Expression::Dual(name), 0),
            "rvMukaiQuotient" => {
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::Quotient(name, Box::new(input)), depth)
            }
            _ => (Expression::Glue(name), 0),
        };
        expr(span, ExprForm::MukaiModuli(expression), depth)
    }

    fn mukai_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::mukai::Command;
        self.advance()?;
        let command = if name == "rvMukaiPair" {
            Command::Pair(Box::new(self.paren_pair(0)?))
        } else {
            self.expect(Token::LBracket)?;
            let n = self.number()?.0;
            self.expect(Token::Comma)?;
            let j = self.number()?.0;
            self.expect(Token::RBracket)?;
            Command::Fujiki(n, j, Box::new(self.paren_pair(0)?))
        };
        Ok(StatementForm::Mukai(Box::new(command)))
    }

    fn mukai_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::mukai::Expression;
        let (expression, depth) = match name {
            "rvMukai" => {
                self.expect(Token::LBracket)?;
                let r = self.signed_wide("a signed Mukai rank", "Mukai rank")?;
                self.expect(Token::Comma)?;
                let s = self.signed_wide("a signed Mukai degree", "Mukai degree")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Vector(r, s, Box::new(input)), depth)
            }
            "rvMukaiIdeal" => {
                self.expect(Token::LBracket)?;
                let basis = self.label()?;
                self.expect(Token::Comma)?;
                let length = self.wide_unsigned("ideal sheaf length")?;
                self.expect(Token::RBracket)?;
                (Expression::Ideal(basis, length), 0)
            }
            "rvMukaiTensor" | "rvMukaiTwist" => {
                let inputs = self.paren_pair(nest)?;
                let depth = 1 + inputs[0].depth.max(inputs[1].depth);
                let expression = if name == "rvMukaiTensor" {
                    Expression::Tensor(Box::new(inputs))
                } else {
                    Expression::Twist(Box::new(inputs))
                };
                (expression, depth)
            }
            _ => {
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                let input = Box::new(input);
                let expression = match name {
                    "rvK3FromKummerDiv" => Expression::FromDivisor(input),
                    "rvK3ToKummerDiv" => Expression::ToDivisor(input),
                    "rvMukaiC" => Expression::Center(input),
                    _ => Expression::Line(input),
                };
                (expression, depth)
            }
        };
        expr(span, ExprForm::Mukai(Box::new(expression)), depth)
    }

    fn k3_automorphism_statement(&mut self, operation: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::dynamics::{Command, Expression, Generator};
        self.advance()?;
        let command = if operation == "rvk3autom" {
            let name = self.label()?;
            self.expect(Token::LBracket)?;
            let basis = self.label()?;
            self.expect(Token::RBracket)?;
            let generators = self.set_named_rows(|parser| {
                let name = parser.label()?;
                let generator = match name.text.as_str() {
                    "translate" => Generator::Translate,
                    "switch" => Generator::Switch,
                    "project" => Generator::Project,
                    _ => {
                        return Err(Error::new(
                            name.span,
                            "expected translate, switch, or project",
                        ))
                    }
                };
                parser.expect(Token::LParen)?;
                let (mask, span) = parser.number()?;
                if mask >= 16 {
                    return Err(Error::new(
                        span,
                        "K3 automorphism mask must be between 0 and 15",
                    ));
                }
                parser.expect(Token::RParen)?;
                Ok(generator(mask as u8))
            })?;
            Command::Define(name, basis, generators)
        } else {
            self.expect(Token::LBracket)?;
            let name = self.label()?;
            let positive_iteration = |parser: &mut Self| {
                let span = parser.span();
                let power = parser.wide_unsigned("K3 Frobenius iteration count")?;
                if power == 0 {
                    return Err(Error::new(
                        span,
                        "K3 Frobenius iteration count must be positive",
                    ));
                }
                Ok(power)
            };
            let command = match operation {
                "rvK3StratumFrobTrace"
                | "rvK3StratumEulerFactor"
                | "rvK3SpechtFrobTrace"
                | "rvK3SpechtEulerFactor" => {
                    self.expect(Token::Comma)?;
                    let base = self.label()?;
                    self.expect(Token::Comma)?;
                    let shape = self.label()?;
                    self.expect(Token::Comma)?;
                    let r = positive_iteration(self)?;
                    self.expect(Token::Comma)?;
                    let i = self.wide_unsigned("cohomological degree")?;
                    let coefficient = if matches!(
                        operation,
                        "rvK3StratumEulerFactor" | "rvK3SpechtEulerFactor"
                    ) {
                        self.expect(Token::Comma)?;
                        Some(self.wide_unsigned("Euler-factor coefficient index")?)
                    } else {
                        None
                    };
                    if matches!(operation, "rvK3SpechtFrobTrace" | "rvK3SpechtEulerFactor") {
                        Command::SpechtTrace(name, base, shape, r, i, coefficient)
                    } else {
                        Command::StratumTrace(name, base, shape, r, i, coefficient)
                    }
                }
                "rvK3FrobHilbOrbits" => {
                    self.expect(Token::Comma)?;
                    let base = self.label()?;
                    self.expect(Token::Comma)?;
                    let length = self.wide_unsigned("K3 Hilbert-scheme length")?;
                    self.expect(Token::Comma)?;
                    let period = self.wide_unsigned("Hilbert-scheme primitive period")?;
                    Command::HilbertOrbits(name, base, length, period)
                }
                "rvK3FrobFixed" | "rvK3FrobZeta" | "rvK3FrobHilbFixed" | "rvK3FrobHilbTrace" => {
                    self.expect(Token::Comma)?;
                    let base = self.label()?;
                    if operation == "rvK3FrobZeta" {
                        Command::FrobZeta(name, base)
                    } else {
                        self.expect(Token::Comma)?;
                        let length = if operation == "rvK3FrobFixed" {
                            None
                        } else {
                            let n = self.wide_unsigned("K3 Hilbert-scheme length")?;
                            self.expect(Token::Comma)?;
                            Some(n)
                        };
                        let power = positive_iteration(self)?;
                        if let Some(n) = length {
                            let grade = if operation == "rvK3FrobHilbTrace" {
                                self.expect(Token::Comma)?;
                                Some(self.wide_unsigned("cohomological degree")?)
                            } else {
                                None
                            };
                            Command::FrobHilbert(name, base, n, power, grade)
                        } else {
                            Command::FrobFixed(name, base, power)
                        }
                    }
                }
                "rvK3AutoPeriod" => {
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(0)?;
                    return Ok(StatementForm::K3AutoPeriod(Box::new(Expression(
                        name,
                        Box::new(input),
                    ))));
                }
                "rvK3AutoCycles" => Command::Cycles(name),
                "rvK3AutoNikulin" => Command::Nikulin(name),
                _ => Command::Entropy(name),
            };
            self.expect(Token::RBracket)?;
            command
        };
        Ok(StatementForm::K3Automorphism(Box::new(command)))
    }

    fn ec_end_order_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::elliptic_end_order::Statement;
        self.advance()?;
        let command = if matches!(name, "rvECEndTrace" | "rvECEndNorm") {
            Statement::Invariant(name == "rvECEndNorm", Box::new(self.paren_expr(0)?))
        } else {
            let certificate = if name == "rvecendorder" {
                Some(self.label()?)
            } else {
                None
            };
            self.expect(Token::LBracket)?;
            let first = self.label()?;
            let command = if name == "rvECEndData" {
                Statement::Data(first)
            } else {
                self.expect(Token::Comma)?;
                let second = self.label()?;
                if let Some(certificate) = certificate {
                    Statement::Define(certificate, first, second)
                } else {
                    Statement::Check(first, second)
                }
            };
            self.expect(Token::RBracket)?;
            command
        };
        Ok(StatementForm::ECEndOrder(Box::new(command)))
    }

    fn ec_horizontal_address(
        &mut self,
        certificate: Label,
    ) -> Result<crate::elliptic_end_order::horizontal::Address> {
        self.expect(Token::Comma)?;
        let ell = self.number()?.0;
        self.expect(Token::Comma)?;
        let (orientation, span) = self.number()?;
        if orientation > 1 {
            return Err(Error::new(span, "CM horizontal orientation must be 0 or 1"));
        }
        self.expect(Token::RBracket)?;
        Ok(crate::elliptic_end_order::horizontal::Address(
            certificate,
            ell,
            orientation == 1,
        ))
    }

    fn ec_cm_cycle_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::elliptic_end_order::cycle::Statement;
        self.advance()?;
        let command = if name == "rveccmcycle" {
            let cycle = self.label()?;
            self.expect(Token::LBracket)?;
            let certificate = self.label()?;
            Statement::Define(cycle, self.ec_horizontal_address(certificate)?)
        } else {
            self.expect(Token::LBracket)?;
            let cycle = self.label()?;
            let command = match name {
                "rvECCMCycleData" => Statement::Data(cycle),
                "rvECCMCycleStep" => {
                    self.expect(Token::Comma)?;
                    Statement::Step(cycle, self.number()?.0)
                }
                "rvECCMCycleVerify" => Statement::Verify(cycle),
                _ => unreachable!(),
            };
            self.expect(Token::RBracket)?;
            command
        };
        Ok(StatementForm::ECCMCycle(Box::new(command)))
    }

    fn ec_horizontal_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::elliptic_end_order::horizontal::Statement;
        self.advance()?;
        let command = if name == "rveccmhorizontal" {
            let map = self.label()?;
            self.expect(Token::LBracket)?;
            let certificate = self.label()?;
            self.expect(Token::Comma)?;
            let target = self.label()?;
            Statement::Define(map, target, self.ec_horizontal_address(certificate)?)
        } else {
            self.expect(Token::LBracket)?;
            let certificate = self.label()?;
            Statement::Data(self.ec_horizontal_address(certificate)?)
        };
        Ok(StatementForm::ECHorizontal(Box::new(command)))
    }

    fn ec_end_canonical_expression(
        &mut self,
        name: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::elliptic_end_order::canonical::Expression;
        self.expect(Token::LBracket)?;
        let certificate = self.label()?;
        self.expect(Token::Comma)?;
        let ring = self.label()?;
        self.expect(Token::RBracket)?;
        let (expression, depth) = match name {
            "rvECEndCanonicalJ" => (Expression::JInvariant(certificate, ring), 0),
            "rvECEndCanonicalPhi" => (Expression::ModularResidual(certificate, ring), 0),
            "rvECEndUnitRoot" | "rvECEndNonunitRoot" => (
                Expression::Root(certificate, ring, name == "rvECEndNonunitRoot"),
                0,
            ),
            _ => {
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (
                    Expression::Embed(
                        certificate,
                        ring,
                        name == "rvECEndEmbedConj",
                        Box::new(input),
                    ),
                    depth,
                )
            }
        };
        expr(span, ExprForm::ECEndCanonical(Box::new(expression)), depth)
    }

    fn ec_end_order_expression(
        &mut self,
        name: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::elliptic_end_order::Expression;
        let (expression, depth) = match name {
            "rvECEndElem" | "rvECEndFrobenius" | "rvECEndJ" => {
                self.expect(Token::LBracket)?;
                let label = self.label()?;
                let expression = match name {
                    "rvECEndElem" => {
                        self.expect(Token::Comma)?;
                        let a =
                            self.signed_wide("a signed coordinate", "CM endomorphism coordinate")?;
                        self.expect(Token::Comma)?;
                        let b =
                            self.signed_wide("a signed coordinate", "CM endomorphism coordinate")?;
                        Expression::Literal(label, a, b)
                    }
                    "rvECEndFrobenius" => Expression::Frobenius(label),
                    _ => Expression::JInvariant(label),
                };
                self.expect(Token::RBracket)?;
                (expression, 0)
            }
            "rvECEndMul" => {
                let pair = self.paren_pair(nest)?;
                let depth = 1 + pair[0].depth.max(pair[1].depth);
                (Expression::Multiply(Box::new(pair)), depth)
            }
            _ => {
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::Conjugate(Box::new(input)), depth)
            }
        };
        expr(span, ExprForm::ECEndOrder(Box::new(expression)), depth)
    }

    fn cm_artin_address(&mut self) -> Result<crate::cm_artin::Address> {
        let order = self.label()?;
        self.expect(Token::Comma)?;
        let p = self.number()?.0;
        self.expect(Token::Comma)?;
        let (orientation, span) = self.number()?;
        if orientation > 1 {
            return Err(Error::new(span, "CM Artin orientation must be 0 or 1"));
        }
        Ok(crate::cm_artin::Address(order, p, orientation == 1))
    }

    fn cm_artin_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::cm_artin::Statement;
        self.advance()?;
        self.expect(Token::LBracket)?;
        let statement = if matches!(name, "rvCMArtinPrime" | "rvCMArtinOrbit") {
            let address = self.cm_artin_address()?;
            self.expect(Token::RBracket)?;
            if name == "rvCMArtinPrime" {
                Statement::Prime(address)
            } else {
                Statement::Orbit(address, Box::new(self.paren_expr(0)?))
            }
        } else {
            let first = self.label()?;
            self.expect(Token::Comma)?;
            let field = self.label()?;
            self.expect(Token::RBracket)?;
            if name == "rvCMArtinFactors" {
                Statement::Factors(first, field)
            } else {
                Statement::Lift(first, field)
            }
        };
        Ok(StatementForm::CMArtin(Box::new(statement)))
    }

    fn cm_artin_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::cm_artin::Expression;
        self.expect(Token::LBracket)?;
        let address = self.cm_artin_address()?;
        let (expression, depth) = if name == "rvCMArtinClass" {
            self.expect(Token::RBracket)?;
            (Expression::Class(address), 0)
        } else {
            self.expect(Token::Comma)?;
            let exponent = self.signed_wide("a signed exponent", "CM Artin exponent")?;
            self.expect(Token::RBracket)?;
            let input = self.paren_expr(nest)?;
            let depth = input.depth + 1;
            (
                Expression::Translate(address, exponent, Box::new(input)),
                depth,
            )
        };
        expr(span, ExprForm::CMArtin(Box::new(expression)), depth)
    }

    fn cm_class_poly_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::cm_class_poly::Statement;
        self.advance()?;
        let command = if matches!(name, "rvcmclasspoly" | "rvcmclasspolylift") {
            let target = self.label()?;
            self.expect(Token::LBracket)?;
            if name == "rvcmclasspoly" {
                Statement::Seed(
                    target,
                    self.signed_wide(
                        "a negative discriminant",
                        "CM class polynomial discriminant",
                    )?,
                )
            } else {
                let source = self.label()?;
                self.expect(Token::Comma)?;
                Statement::Lift(target, source, self.number()?.0)
            }
        } else {
            self.expect(Token::LBracket)?;
            let polynomial = self.label()?;
            if name == "rvCMClassPolyData" {
                Statement::Data(polynomial)
            } else {
                self.expect(Token::Comma)?;
                Statement::Coefficient(polynomial, self.number()?.0)
            }
        };
        self.expect(Token::RBracket)?;
        Ok(StatementForm::CMClassPolynomial(Box::new(command)))
    }

    fn cm_class_poly_expression(
        &mut self,
        name: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::cm_class_poly::Expression;
        self.expect(Token::LBracket)?;
        let polynomial = self.label()?;
        let (expression, depth) = if name == "rvCMClassPolyEval" {
            self.expect(Token::RBracket)?;
            let input = self.paren_expr(nest)?;
            let depth = input.depth + 1;
            (Expression::Evaluate(polynomial, Box::new(input)), depth)
        } else {
            self.expect(Token::Comma)?;
            let ring = self.label()?;
            self.expect(Token::Comma)?;
            let index = self.number()?.0;
            self.expect(Token::RBracket)?;
            (Expression::ModularCoefficient(polynomial, ring, index), 0)
        };
        expr(
            span,
            ExprForm::CMClassPolynomial(Box::new(expression)),
            depth,
        )
    }

    fn elliptic_modular_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::elliptic_modular::Statement;
        self.advance()?;
        self.expect(Token::LBracket)?;
        let level = self.number()?.0;
        let statement = if name == "rvModularPolyDegree" {
            Statement::Degree(level)
        } else {
            self.expect(Token::Comma)?;
            let i = self.number()?.0;
            self.expect(Token::Comma)?;
            let j = self.number()?.0;
            Statement::Coefficient(level, i, j)
        };
        self.expect(Token::RBracket)?;
        Ok(StatementForm::EllipticModular(statement))
    }

    fn elliptic_modular_expression(
        &mut self,
        name: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::elliptic_modular::Expression;
        self.expect(Token::LBracket)?;
        let level = self.number()?.0;
        let (expression, depth) = if matches!(name, "rvModularPolyEval" | "rvModularPolyWittEval") {
            self.expect(Token::RBracket)?;
            let inputs = self.paren_pair(nest)?;
            let depth = inputs[0].depth.max(inputs[1].depth) + 1;
            let inputs = Box::new(inputs);
            let expression = if name == "rvModularPolyWittEval" {
                Expression::WittEvaluate(level, inputs)
            } else {
                Expression::Evaluate(level, inputs)
            };
            (expression, depth)
        } else {
            self.expect(Token::Comma)?;
            let j = self.number()?.0;
            self.expect(Token::RBracket)?;
            let input = self.paren_expr(nest)?;
            let depth = input.depth + 1;
            (
                Expression::SpecializeCoefficient(level, j, Box::new(input)),
                depth,
            )
        };
        expr(span, ExprForm::EllipticModular(Box::new(expression)), depth)
    }

    fn cm_neighbor_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::cm_neighbors::Statement;
        self.advance()?;
        self.expect(Token::LBracket)?;
        let ell = self.number()?.0;
        let address = if name == "rvCMNeighborData" {
            self.expect(Token::Comma)?;
            Some(self.number()?.0)
        } else {
            None
        };
        self.expect(Token::RBracket)?;
        let input = Box::new(self.paren_expr(0)?);
        let statement = match address {
            Some(address) => Statement::Address(ell, address, input),
            None if name == "rvCMNeighborCounts" => Statement::Counts(ell, input),
            None => Statement::Neighbors(ell, input),
        };
        Ok(StatementForm::CMNeighbors(statement))
    }

    fn cm_neighbor_expression(&mut self, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::cm_neighbors::Expression;
        self.expect(Token::LBracket)?;
        let target = self.label()?;
        self.expect(Token::Comma)?;
        let ell = self.number()?.0;
        self.expect(Token::Comma)?;
        let address = self.number()?.0;
        self.expect(Token::RBracket)?;
        let input = self.paren_expr(nest)?;
        let depth = input.depth + 1;
        expr(
            span,
            ExprForm::CMNeighbor(Expression(target, ell, address, Box::new(input))),
            depth,
        )
    }

    fn cm_class_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::cm_ring::Statement;
        self.advance()?;
        let statement = match name {
            "rvcmorder" => {
                let order = self.label()?;
                self.expect(Token::LBracket)?;
                let discriminant =
                    self.signed_wide("a negative discriminant", "CM order discriminant")?;
                self.expect(Token::RBracket)?;
                Statement::Define(order, discriminant)
            }
            "rvCMClassGroup" => {
                self.expect(Token::LBracket)?;
                let order = self.label()?;
                self.expect(Token::RBracket)?;
                Statement::Group(order)
            }
            "rvCMClassOrder" => Statement::Order(Box::new(self.paren_expr(0)?)),
            _ => Statement::SameGenus(Box::new(self.paren_pair(0)?)),
        };
        Ok(StatementForm::CMClass(Box::new(statement)))
    }

    fn cm_class_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::cm_ring::Expression;
        let (expression, depth) = match name {
            "rvCMClass" | "rvCMPrincipal" | "rvK3CMClass" => {
                self.expect(Token::LBracket)?;
                let order = self.label()?;
                let expression = match name {
                    "rvCMClass" => {
                        let mut coefficients = [0; 3];
                        for coefficient in &mut coefficients {
                            self.expect(Token::Comma)?;
                            *coefficient = self
                                .signed_wide("a signed form coefficient", "CM class coefficient")?;
                        }
                        Expression::Literal(order, coefficients)
                    }
                    "rvK3CMClass" => {
                        self.expect(Token::Comma)?;
                        Expression::Period(order, self.label()?)
                    }
                    _ => Expression::Principal(order),
                };
                self.expect(Token::RBracket)?;
                (expression, 0)
            }
            "rvCMPow" => {
                self.expect(Token::LBracket)?;
                let power = self.signed_wide("a signed exponent", "CM class exponent")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Power(power, Box::new(input)), depth)
            }
            "rvCMInverse" => {
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Inverse(Box::new(input)), depth)
            }
            _ => {
                let inputs = self.paren_pair(nest)?;
                let depth = 1 + inputs[0].depth.max(inputs[1].depth);
                (Expression::Compose(Box::new(inputs)), depth)
            }
        };
        expr(span, ExprForm::CMClass(Box::new(expression)), depth)
    }

    fn k3_period_statement(&mut self, command: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::cm::Command;
        use crate::hyperelliptic_kummer_surface::period::Statement;
        self.advance()?;
        let statement = if command == "rvk3period" {
            let period = self.label()?;
            Statement::Define(period, Box::new(self.paren_pair(0)?))
        } else {
            self.expect(Token::LBracket)?;
            let period = self.label()?;
            if command == "rvK3CMEquivalent" {
                self.expect(Token::Comma)?;
                let other = self.label()?;
                self.expect(Token::RBracket)?;
                Statement::CM(Command::Equivalent(period, other))
            } else {
                self.expect(Token::RBracket)?;
                match command {
                    "rvK3PeriodReport" => Statement::Report(period),
                    "rvK3CMForm" => Statement::CM(Command::Form(period)),
                    "rvK3CMPeriods" => Statement::CM(Command::Periods(period)),
                    "rvK3CMClassNumber" => Statement::CM(Command::ClassNumber(period)),
                    _ => Statement::Is11(period, Box::new(self.paren_expr(0)?)),
                }
            }
        };
        Ok(StatementForm::K3Period(Box::new(statement)))
    }

    fn k3_period_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::period::Expression;
        self.expect(Token::LBracket)?;
        let period = self.label()?;
        if name == "rvK3CMBasis" {
            self.expect(Token::Comma)?;
            let (index, at) = self.number()?;
            if index > 1 {
                return Err(Error::new(at, "K3 CM basis index must be 0 or 1"));
            }
            self.expect(Token::RBracket)?;
            return expr(
                span,
                ExprForm::K3Period(Box::new(Expression::CMBasis(period, index))),
                0,
            );
        }
        self.expect(Token::RBracket)?;
        let (expression, depth) = if name == "rvK3PeriodBlade" {
            (Expression::Blade(period), 0)
        } else {
            let input = self.paren_expr(nest)?;
            let depth = 1 + input.depth;
            (Expression::CM(period, Box::new(input)), depth)
        };
        expr(span, ExprForm::K3Period(Box::new(expression)), depth)
    }

    fn spinor_coefficient_statement(&mut self) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::spinor::Coefficient;
        self.advance()?;
        self.expect(Token::LBracket)?;
        let index = self.spinor_index()?;
        self.expect(Token::RBracket)?;
        Ok(StatementForm::K3SpinorCoefficient(Coefficient(
            index,
            Box::new(self.paren_expr(0)?),
        )))
    }

    fn spinor_index(&mut self) -> Result<u32> {
        use crate::hyperelliptic_kummer_surface::spinor::DIMENSION;
        let (index, at) = self.number()?;
        if index >= DIMENSION {
            return Err(Error::new(at, "K3 spinor index must be between 0 and 1023"));
        }
        Ok(index as u32)
    }

    fn spinor_chirality(&mut self) -> Result<u8> {
        let (chirality, at) = self.number()?;
        if chirality > 1 {
            return Err(Error::new(at, "K3 spinor chirality must be 0 or 1"));
        }
        Ok(chirality as u8)
    }

    fn spinor_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::spinor::Expression;
        let (expression, depth) = match name {
            "rvK3SpinorBasis" => {
                self.expect(Token::LBracket)?;
                let marking = self.label()?;
                self.expect(Token::Comma)?;
                let chirality = self.spinor_chirality()?;
                self.expect(Token::Comma)?;
                let index = self.spinor_index()?;
                self.expect(Token::RBracket)?;
                (Expression::Basis(marking, chirality, index), 0)
            }
            "rvK3SpinorScale" => {
                self.expect(Token::LBracket)?;
                let numerator =
                    self.signed_wide("a signed numerator", "K3 spinor scalar numerator")?;
                self.expect(Token::Comma)?;
                let denominator = self.wide_unsigned("K3 spinor scalar denominator")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (
                    Expression::Scale(numerator, denominator, Box::new(input)),
                    depth,
                )
            }
            _ => {
                let inputs = self.paren_pair(nest)?;
                let depth = 1 + inputs[0].depth.max(inputs[1].depth);
                (Expression::Act(Box::new(inputs)), depth)
            }
        };
        expr(span, ExprForm::K3Spinor(Box::new(expression)), depth)
    }

    fn clifford_coefficient_statement(&mut self) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::clifford::{Coefficient, VOLUME_MASK};
        self.advance()?;
        self.expect(Token::LBracket)?;
        let (mask, at) = self.number()?;
        if mask > VOLUME_MASK as usize {
            return Err(Error::new(at, "Clifford monomial mask must fit 22 bits"));
        }
        self.expect(Token::RBracket)?;
        Ok(StatementForm::K3CliffordCoefficient(Coefficient(
            mask as u32,
            Box::new(self.paren_expr(0)?),
        )))
    }

    fn clifford_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::clifford::Expression;
        let (expression, depth) = match name {
            "rvCliffordOne" | "rvCliffordVolume" => {
                self.expect(Token::LBracket)?;
                let basis = self.label()?;
                self.expect(Token::RBracket)?;
                let expression = if name == "rvCliffordOne" {
                    Expression::One(basis)
                } else {
                    Expression::Volume(basis)
                };
                (expression, 0)
            }
            "rvCliffordScale" => {
                self.expect(Token::LBracket)?;
                let numerator =
                    self.signed_wide("a signed numerator", "Clifford scalar numerator")?;
                self.expect(Token::Comma)?;
                let denominator = self.wide_unsigned("Clifford scalar denominator")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (
                    Expression::Scale(numerator, denominator, Box::new(input)),
                    depth,
                )
            }
            "rvCliffordMul" | "rvCliffordAct" => {
                let inputs = self.paren_pair(nest)?;
                let depth = 1 + inputs[0].depth.max(inputs[1].depth);
                let expression = if name == "rvCliffordMul" {
                    Expression::Multiply(Box::new(inputs))
                } else {
                    Expression::Act(Box::new(inputs))
                };
                (expression, depth)
            }
            _ => {
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                let expression = if name == "rvCliffordVector" {
                    Expression::Vector(Box::new(input))
                } else {
                    Expression::Reverse(Box::new(input))
                };
                (expression, depth)
            }
        };
        expr(span, ExprForm::K3Clifford(Box::new(expression)), depth)
    }

    fn k3_lattice_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::unimodular::Command;
        self.advance()?;
        let command = if name == "rvK3LatticePair" {
            Command::Pair(Box::new(self.paren_pair(0)?))
        } else {
            self.expect(Token::LBracket)?;
            let basis = self.label()?;
            self.expect(Token::RBracket)?;
            Command::Gram(basis)
        };
        Ok(StatementForm::K3Lattice(Box::new(command)))
    }

    fn hilbert_square_expression(
        &mut self,
        name: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::hilbert_square::{llv::Operation, Expression};
        let (expression, depth) = match name {
            "rvHilbSquareCup" => {
                let inputs = self.paren_pair(nest)?;
                let depth = 1 + inputs[0].depth.max(inputs[1].depth);
                (Expression::Cup(Box::new(inputs)), depth)
            }
            "rvHilbSquareLower" | "rvHilbSquareLambda" | "rvHilbSquarePrimitive" => {
                let operation = match name {
                    "rvHilbSquareLower" => Operation::Lower,
                    "rvHilbSquareLambda" => Operation::Lambda,
                    _ => Operation::Primitive,
                };
                let inputs = self.paren_pair(nest)?;
                let depth = 1 + inputs[0].depth.max(inputs[1].depth);
                (Expression::Llv(operation, Box::new(inputs)), depth)
            }
            "rvHilbSquareRotate" => {
                let [a, b, class] = self.paren_triple(nest)?;
                let depth = 1 + a.depth.max(b.depth).max(class.depth);
                (Expression::Rotate(Box::new([a, b, class])), depth)
            }
            "rvNakajimaHilbSquare" => {
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::FromNakajima(Box::new(input)), depth)
            }
            "rvHilbSquareNakajima" => {
                self.expect(Token::LBracket)?;
                let scale = self.signed_wide("a signed integer", "Fock conversion scale")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::ToNakajima(scale, Box::new(input)), depth)
            }
            "rvHilbSquareBasis" => {
                self.expect(Token::LBracket)?;
                let basis = self.label()?;
                self.expect(Token::Comma)?;
                let degree =
                    self.signed_wide("a signed cohomological degree", "Hilbert-square degree")?;
                self.expect(Token::Comma)?;
                let index = self.number()?.0;
                self.expect(Token::RBracket)?;
                (Expression::Basis(basis, degree, index), 0)
            }
            "rvHilbSquareC2" | "rvHilbSquareC4" => {
                self.expect(Token::LBracket)?;
                let basis = self.label()?;
                self.expect(Token::RBracket)?;
                (Expression::Chern(basis, name == "rvHilbSquareC4"), 0)
            }
            _ => {
                self.expect(Token::LBracket)?;
                let numerator =
                    self.signed_wide("a signed integer", "Hilbert-square coefficient")?;
                let denominator = if name == "rvHilbSquareScale" {
                    self.expect(Token::Comma)?;
                    Some(self.signed_wide("a signed integer", "Hilbert-square denominator")?)
                } else {
                    None
                };
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                let expression = if let Some(denominator) = denominator {
                    Expression::Scale(numerator, denominator, Box::new(input))
                } else {
                    Expression::Divisor(numerator, Box::new(input))
                };
                (expression, depth)
            }
        };
        expr(span, ExprForm::HilbertSquare(Box::new(expression)), depth)
    }

    fn siegel_weight_coordinates(&mut self) -> Result<crate::siegel::Weight> {
        let mut coordinates = [0; 3];
        for (index, coordinate) in coordinates.iter_mut().enumerate() {
            if index != 0 {
                self.expect(Token::Comma)?;
            }
            *coordinate = self.signed_wide("a signed doubled coordinate", "Siegel weight")?;
        }
        self.expect(Token::RBracket)?;
        Ok(crate::siegel::Weight(coordinates))
    }

    fn siegel_maass_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::maass::{
            Command, Family, Monomial,
        };
        self.advance()?;
        let command = if name == "rvsiegelmaass" {
            let name = self.label()?;
            self.expect(Token::LBracket)?;
            let weight = self.wide_unsigned("Maaß weight")?;
            self.expect(Token::RBracket)?;
            self.expect(Token::LBrace)?;
            let mut terms = Vec::new();
            while self.token() != Token::RBrace {
                let family = self.label()?;
                let family = match family.text.as_str() {
                    "minus" => Family::MinusTwo,
                    "zero" => Family::Zero,
                    _ => {
                        return Err(Error::new(
                            family.span,
                            "expected a minus or zero Jacobi monomial",
                        ))
                    }
                };
                self.expect(Token::LParen)?;
                let mut exponents = [0; 3];
                for exponent in &mut exponents {
                    *exponent = self.wide_unsigned("Jacobi modular-factor exponent")?;
                    self.expect(Token::Comma)?;
                }
                let scale =
                    self.signed_wide("a signed Jacobi coefficient", "Jacobi coefficient")?;
                self.expect(Token::RParen)?;
                terms.push((Monomial(family, exponents), scale));
                if !self.eat(Token::Comma)? {
                    break;
                }
            }
            self.expect(Token::RBrace)?;
            Command::Declare(name, weight, terms)
        } else {
            self.expect(Token::LBracket)?;
            let query = if name == "rvIgusaChi12" {
                Command::Chi12(self.elliptic_genus_index("P exponent")?)
            } else {
                let form = self.label()?;
                self.expect(Token::Comma)?;
                if name == "rvSiegelMaassJacobi" {
                    let n = self.wide_unsigned("q exponent")?;
                    self.expect(Token::Comma)?;
                    let r = self.signed_wide("a signed Laurent exponent", "y exponent")?;
                    Command::Jacobi(form, n, r)
                } else {
                    Command::Lift(form, self.elliptic_genus_index("P exponent")?)
                }
            };
            self.expect(Token::RBracket)?;
            query
        };
        Ok(StatementForm::SiegelMaass(Box::new(command)))
    }

    fn period_grade(&mut self) -> Result<period::Grade> {
        let weight = self.wide_unsigned("period-polynomial weight")?;
        self.expect(Token::Comma)?;
        let name = self.label()?;
        let sign = match name.text.as_str() {
            "plus" => period::Sign::Plus,
            "minus" => period::Sign::Minus,
            _ => {
                return Err(Error::new(
                    name.span,
                    "period-polynomial sign must be plus or minus",
                ))
            }
        };
        period::Grade::new(weight, sign).map_err(|error| error.at(name.span))
    }

    fn period_statement(&mut self, name: &str) -> Result<StatementForm> {
        self.advance()?;
        self.expect(Token::LBracket)?;
        let report = if name == "rvPeriodCoeff" {
            let index = self.wide_unsigned("period-polynomial coefficient index")?;
            self.expect(Token::RBracket)?;
            period::Report::Coefficient(index, Box::new(self.paren_expr(0)?))
        } else {
            let grade = self.period_grade()?;
            self.expect(Token::RBracket)?;
            if name == "rvPeriodDim" {
                period::Report::Dimension(grade)
            } else {
                period::Report::Matrix(grade, Box::new(self.paren_expr(0)?))
            }
        };
        Ok(StatementForm::PeriodPolynomial(report))
    }

    fn period_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use period::Expression;
        let (operation, depth) = if name == "rvPeriodHecke" {
            let [operator, polynomial] = self.paren_pair(nest)?;
            let depth = operator.depth.max(polynomial.depth) + 1;
            (Expression::Hecke(Box::new([operator, polynomial])), depth)
        } else {
            self.expect(Token::LBracket)?;
            if name == "rvPeriodScale" {
                let numerator =
                    self.signed_wide("an integer numerator", "period-polynomial numerator")?;
                self.expect(Token::Comma)?;
                let denominator =
                    self.signed_wide("a positive denominator", "period-polynomial denominator")?;
                self.expect(Token::RBracket)?;
                let scalar =
                    period::ratio(numerator, denominator).map_err(|error| error.at(span))?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Scale(scalar, Box::new(input)), depth)
            } else {
                let grade = self.period_grade()?;
                if name == "rvPeriodBasis" {
                    self.expect(Token::Comma)?;
                    let (index, _) = self.number()?;
                    self.expect(Token::RBracket)?;
                    (Expression::Basis(grade, index), 0)
                } else {
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LParen)?;
                    let mut coefficients = Vec::new();
                    while self.token() != Token::RParen {
                        let fraction = self.eat(Token::LParen)?;
                        let numerator = self
                            .signed_wide("an integer numerator", "period-polynomial numerator")?;
                        let denominator = if fraction {
                            self.expect(Token::Comma)?;
                            let denominator = self.signed_wide(
                                "a positive denominator",
                                "period-polynomial denominator",
                            )?;
                            self.expect(Token::RParen)?;
                            denominator
                        } else {
                            1
                        };
                        coefficients.push(
                            period::ratio(numerator, denominator)
                                .map_err(|error| error.at(span))?,
                        );
                        if !self.eat(Token::Comma)? {
                            break;
                        }
                    }
                    self.expect(Token::RParen)?;
                    (Expression::Polynomial(grade, coefficients), 0)
                }
            }
        };
        expr(span, ExprForm::PeriodPolynomial(operation), depth)
    }

    fn hecke_quadratic_statement(&mut self, name: &str) -> Result<StatementForm> {
        use quadratic_hecke::Query;
        self.advance()?;
        let query = if name == "rvHeckeQuadraticIntegral" {
            Query::Integral(Box::new(self.paren_expr(0)?))
        } else {
            self.expect(Token::LBracket)?;
            let name_of_packet = self.label()?;
            let query = match name {
                "rvHeckeQuadraticOrder" => Query::Data(name_of_packet),
                "rvHeckeQuadraticCoeff" => {
                    self.expect(Token::Comma)?;
                    let index = self.wide_unsigned("quadratic Hecke Fourier index")?;
                    Query::Fourier(name_of_packet, index)
                }
                _ => {
                    self.expect(Token::Comma)?;
                    let span = self.span();
                    let prime = self.wide_unsigned("quadratic Hecke local prime")?;
                    let prime = u64::try_from(prime)
                        .map_err(|_| Error::new(span, "quadratic Hecke local prime exceeds u64"))?;
                    Query::Local(name_of_packet, prime)
                }
            };
            self.expect(Token::RBracket)?;
            query
        };
        Ok(StatementForm::HeckeQuadratic(query))
    }

    fn hecke_residual_statement(&mut self, name: &str) -> Result<StatementForm> {
        use residual_hecke::Statement;
        self.advance()?;
        let statement = if name == "rvheckeresidual" {
            let name = self.label()?;
            self.expect(Token::LBracket)?;
            let source = self.label()?;
            self.expect(Token::Comma)?;
            let (prime, _) = self.number()?;
            self.expect(Token::Comma)?;
            let root = self.signed_wide("a residual generator root", "residual generator root")?;
            self.expect(Token::RBracket)?;
            Statement::Declare(name, source, prime, root)
        } else {
            self.expect(Token::LBracket)?;
            let name = self.label()?;
            self.expect(Token::RBracket)?;
            Statement::Sturm(name)
        };
        Ok(StatementForm::HeckeResidual(statement))
    }

    fn hecke_residual_expression(
        &mut self,
        name: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use residual_hecke::Expression;
        self.expect(Token::LBracket)?;
        let (operation, depth) = match name {
            "rvHeckeResidualReduce" => {
                let branch = self.label()?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Reduce(branch, Box::new(input)), depth)
            }
            "rvHeckeResidualCoeff" => {
                let branch = self.label()?;
                self.expect(Token::Comma)?;
                let index = self.wide_unsigned("residual Fourier index")?;
                self.expect(Token::RBracket)?;
                (Expression::Fourier(branch, index), 0)
            }
            _ => {
                let weight = self.wide_unsigned("Eisenstein weight")?;
                self.expect(Token::Comma)?;
                let (prime, _) = self.number()?;
                self.expect(Token::Comma)?;
                let index = self.wide_unsigned("Eisenstein Fourier index")?;
                self.expect(Token::RBracket)?;
                (Expression::Eisenstein(weight, prime, index), 0)
            }
        };
        expr(span, ExprForm::HeckeResidual(operation), depth)
    }

    fn hecke_eigen_statement(&mut self, name: &str) -> Result<StatementForm> {
        use eigenpacket::{Report, Statement};
        self.advance()?;
        let statement = match name {
            "rvheckepacket" => Statement::Declare(self.label()?, Box::new(self.paren_expr(0)?)),
            "rvHeckeEigenData" => {
                self.expect(Token::LBracket)?;
                let packet = self.label()?;
                self.expect(Token::RBracket)?;
                Statement::Data(packet)
            }
            _ => {
                let report = match name {
                    "rvHeckeEigenTrace" => Report::Trace,
                    "rvHeckeEigenNorm" => Report::Norm,
                    _ => {
                        self.expect(Token::LBracket)?;
                        let j = self.wide_unsigned("Hecke eigenpacket coordinate")?;
                        self.expect(Token::RBracket)?;
                        Report::Coordinate(j)
                    }
                };
                Statement::Report(report, Box::new(self.paren_expr(0)?))
            }
        };
        Ok(StatementForm::HeckeEigen(statement))
    }

    fn hecke_eigen_expression(
        &mut self,
        name: &Label,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use eigenpacket::Expression;
        if name.text == "rvHeckeEigenMul" {
            let [a, b] = self.paren_pair(nest)?;
            let depth = a.depth.max(b.depth) + 1;
            return expr(
                span,
                ExprForm::HeckeEigen(Expression::Multiply(Box::new([a, b]))),
                depth,
            );
        }
        self.expect(Token::LBracket)?;
        let packet = self.label()?;
        let (operation, depth) = match name.text.as_str() {
            "rvHeckeEigenImage" => {
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Image(packet, Box::new(input)), depth)
            }
            "rvHeckeEigenGenerator" => {
                self.expect(Token::RBracket)?;
                (Expression::Generator(packet), 0)
            }
            "rvHeckeEigenCoeff" => {
                self.expect(Token::Comma)?;
                let n = self.wide_unsigned("Hecke Fourier index")?;
                self.expect(Token::RBracket)?;
                (Expression::Coefficient(packet, n), 0)
            }
            _ => {
                self.expect(Token::Comma)?;
                let numerator = self.signed_wide(
                    "a signed Hecke eigenpacket numerator",
                    "Hecke eigenpacket scalar numerator",
                )?;
                self.expect(Token::Comma)?;
                let denominator = self.wide_unsigned("Hecke eigenpacket scalar denominator")?;
                self.expect(Token::RBracket)?;
                (Expression::Scalar(packet, numerator, denominator), 0)
            }
        };
        expr(span, ExprForm::HeckeEigen(operation), depth)
    }

    fn hecke_operator_statement(&mut self, name: &str) -> Result<StatementForm> {
        self.advance()?;
        let query = match name {
            "rvHeckeActionTrace" => hecke::Report::Trace,
            "rvHeckeChar" => hecke::Report::Characteristic,
            _ => {
                self.expect(Token::LBracket)?;
                let index = self.wide_unsigned("Hecke characteristic coefficient index")?;
                self.expect(Token::RBracket)?;
                hecke::Report::CharacteristicCoefficient(index)
            }
        };
        Ok(StatementForm::ModularHecke(
            query,
            Box::new(self.paren_expr(0)?),
        ))
    }

    fn hurwitz_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::hurwitz::Query;
        self.advance()?;
        self.expect(Token::LBracket)?;
        let query = if name == "rvHeckeTrace" {
            let weight = self.wide_unsigned("cusp form weight")?;
            self.expect(Token::Comma)?;
            Query::Trace(weight, self.wide_unsigned("Hecke index")?)
        } else {
            let n = self.signed_wide("a signed Hurwitz index", "Hurwitz index")?;
            if name == "rvHurwitzClass" {
                Query::Class(n)
            } else {
                self.expect(Token::Comma)?;
                let r = self.signed_wide("a signed Laurent exponent", "y exponent")?;
                Query::Jacobi(n, r)
            }
        };
        self.expect(Token::RBracket)?;
        Ok(StatementForm::Hurwitz(Box::new(query)))
    }

    fn igusa_inverse_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::hilbert::{
            elliptic_genus::{
                inverse::{Chamber, Query},
                Query as EllipticQuery,
            },
            Command,
        };
        self.advance()?;
        self.expect(Token::LBracket)?;
        let m = self.signed_wide("a signed P exponent", "P exponent")?;
        self.expect(Token::Comma)?;
        let n = self.signed_wide("a signed q exponent", "q exponent")?;
        let query = if name == "rvIgusaResidue" {
            Query::Residue(m, n)
        } else {
            self.expect(Token::Comma)?;
            let r = self.signed_wide("a signed Laurent exponent", "y exponent")?;
            match name {
                "rvIgusaWallJump" => Query::WallJump(m, n, r),
                "rvIgusaFinite" => Query::Finite(m, n, r),
                "rvIgusaMockWeight2" => Query::MockWeightTwo(m, n, r),
                _ => {
                    self.expect(Token::Comma)?;
                    let side = self.label()?;
                    let chamber = match side.text.as_str() {
                        "inside" => Chamber::Inside,
                        "outside" => Chamber::Outside,
                        _ => {
                            return Err(Error::new(
                                side.span,
                                "Igusa reciprocal chamber must be inside or outside",
                            ))
                        }
                    };
                    if name == "rvIgusaPolar" {
                        Query::Polar(m, n, r, chamber)
                    } else {
                        Query::Coefficient(m, n, r, chamber)
                    }
                }
            }
        };
        self.expect(Token::RBracket)?;
        Ok(StatementForm::KummerHilbert(Box::new(
            Command::EllipticGenus(EllipticQuery::IgusaInverse(query)),
        )))
    }

    fn siegel_chamber_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::farey::Command;
        self.advance()?;
        let command = if name == "rvsiegelchamber" {
            let name = self.label()?;
            self.expect(Token::LBracket)?;
            let mut coordinates = [0; 3];
            for (i, coordinate) in coordinates.iter_mut().enumerate() {
                if i != 0 {
                    self.expect(Token::Comma)?;
                }
                *coordinate =
                    self.signed_wide("a signed period coordinate", "Siegel period coordinate")?;
            }
            Command::Declare(name, coordinates)
        } else {
            self.expect(Token::LBracket)?;
            let chamber = self.label()?;
            if name == "rvSiegelChamberData" {
                Command::Data(chamber)
            } else {
                let target = if name == "rvIgusaChamberPath" {
                    self.expect(Token::Comma)?;
                    Some(self.label()?)
                } else {
                    None
                };
                let mut indices = [0; 3];
                for index in &mut indices {
                    self.expect(Token::Comma)?;
                    *index = self.signed_wide(
                        "a signed Fourier exponent",
                        "Igusa chamber Fourier exponent",
                    )?;
                }
                match target {
                    Some(target) => Command::Path(chamber, target, indices),
                    None => Command::Coefficient(chamber, indices),
                }
            }
        };
        self.expect(Token::RBracket)?;
        Ok(StatementForm::SiegelChamber(Box::new(command)))
    }

    fn k3_lattice_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::unimodular::{Expression, RANK};
        let (expression, depth) = match name {
            "rvK3Lattice" => {
                self.expect(Token::LBracket)?;
                let basis = self.label()?;
                self.expect(Token::RBracket)?;
                self.expect(Token::LBrace)?;
                let w = self.signed_tuple("a signed wedge coordinate", "K3 lattice coordinate")?;
                self.expect(Token::Comma)?;
                let b =
                    self.signed_tuple("a signed exceptional coordinate", "K3 lattice coordinate")?;
                self.expect(Token::RBrace)?;
                let w = w
                    .try_into()
                    .map_err(|_| Error::new(span, "rvK3Lattice requires six wedge coordinates"))?;
                let b = b.try_into().map_err(|_| {
                    Error::new(span, "rvK3Lattice requires sixteen exceptional coordinates")
                })?;
                (Expression::Literal(basis, Box::new((w, b))), 0)
            }
            "rvK3LatticeBasis" => {
                self.expect(Token::LBracket)?;
                let basis = self.label()?;
                self.expect(Token::Comma)?;
                let (index, at) = self.number()?;
                if index >= RANK {
                    return Err(Error::new(
                        at,
                        "rvK3LatticeBasis index must be between 0 and 21",
                    ));
                }
                self.expect(Token::RBracket)?;
                (Expression::Basis(basis, index), 0)
            }
            "rvK3Reflect" => {
                let inputs = self.paren_pair(nest)?;
                let depth = 1 + inputs[0].depth.max(inputs[1].depth);
                (Expression::Reflect(Box::new(inputs)), depth)
            }
            _ => {
                self.expect(Token::LBracket)?;
                let (mask, at) = self.number()?;
                if mask >= 16 {
                    return Err(Error::new(
                        at,
                        "rvK3Transvection mask must be between 0 and 15",
                    ));
                }
                self.expect(Token::Comma)?;
                let amount =
                    self.signed_wide("a signed transvection amount", "K3 transvection amount")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (
                    Expression::Transvection(mask as u8, amount, Box::new(input)),
                    depth,
                )
            }
        };
        expr(span, ExprForm::K3Lattice(Box::new(expression)), depth)
    }

    fn nakajima_grade_args(&mut self) -> Result<(Label, i128, i128)> {
        let basis = self.label()?;
        self.expect(Token::Comma)?;
        let length = self.signed_wide("a signed Nakajima length grade", "Nakajima length grade")?;
        self.expect(Token::Comma)?;
        let degree = self.signed_wide(
            "a signed Nakajima cohomological degree",
            "Nakajima cohomological degree",
        )?;
        Ok((basis, length, degree))
    }

    fn nakajima_primary_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::nakajima::primary::Command;
        self.advance()?;
        let command = if name == "rvNakajimaPrimaryDim" {
            self.expect(Token::LBracket)?;
            let (basis, length, degree) = self.nakajima_grade_args()?;
            self.expect(Token::RBracket)?;
            Command::Dimension(basis, length, degree)
        } else {
            Command::IsPrimary(Box::new(self.paren_expr(0)?))
        };
        Ok(StatementForm::NakajimaPrimary(Box::new(command)))
    }

    fn nakajima_mode(&mut self) -> Result<crate::hyperelliptic_kummer_surface::nakajima::Mode> {
        let span = self.span();
        let length = self.wide_unsigned("Nakajima mode length")?;
        self.expect(Token::Comma)?;
        let color = self.number()?.0;
        crate::hyperelliptic_kummer_surface::nakajima::Mode::new(length, color)
            .map_err(|error| error.at(span))
    }

    fn hilbert_cohomology_expression(
        &mut self,
        name: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::nakajima::cohomology::Expression;
        let (expression, depth) = if name == "rvHilbCohRotate" {
            let [a, b, class] = self.paren_triple(nest)?;
            let depth = 1 + a.depth.max(b.depth).max(class.depth);
            (Expression::Rotate(Box::new([a, b, class])), depth)
        } else if matches!(
            name,
            "rvHilbCohCup" | "rvHilbCohLefschetzInverse" | "rvHilbCohLower" | "rvHilbCohLLVLower"
        ) {
            let inputs = self.paren_pair(nest)?;
            let depth = 1 + inputs[0].depth.max(inputs[1].depth);
            use crate::hyperelliptic_kummer_surface::nakajima::cohomology::lefschetz::Operation;
            let expression = match name {
                "rvHilbCohLefschetzInverse" => {
                    Expression::Lefschetz(Operation::Inverse, Box::new(inputs))
                }
                "rvHilbCohLower" => Expression::Lefschetz(Operation::Lower, Box::new(inputs)),
                "rvHilbCohLLVLower" => Expression::LlvLower(Box::new(inputs)),
                _ => Expression::Cup(Box::new(inputs)),
            };
            (expression, depth)
        } else {
            self.expect(Token::LBracket)?;
            match name {
                "rvHilbCohUnit" | "rvHilbCohPoint" => {
                    let basis = self.label()?;
                    self.expect(Token::Comma)?;
                    let n = self.wide_unsigned("Hilbert length")?;
                    self.expect(Token::RBracket)?;
                    let expression = if name == "rvHilbCohUnit" {
                        Expression::Unit(basis, n)
                    } else {
                        Expression::Point(basis, n)
                    };
                    (expression, 0)
                }
                "rvHilbCohChern" | "rvHilbCohCh" | "rvHilbCohSchur" => {
                    let n = self.wide_unsigned("Hilbert length")?;
                    self.expect(Token::Comma)?;
                    let j = self.wide_unsigned("characteristic class index")?;
                    self.expect(Token::RBracket)?;
                    let shape = if name == "rvHilbCohSchur" {
                        Some(
                            self.set_named_rows(|parser| parser.wide_unsigned("Schur row length"))?,
                        )
                    } else {
                        None
                    };
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    (
                        if let Some(shape) = shape {
                            Expression::Schur(n, j, shape, Box::new(input))
                        } else {
                            Expression::Chern(n, j, name == "rvHilbCohCh", Box::new(input))
                        },
                        depth,
                    )
                }
                _ => {
                    let k = self.signed_wide("a signed integer", "Hilbert cohomology numerator")?;
                    let denominator = if name == "rvHilbCohFock" {
                        None
                    } else {
                        self.expect(Token::Comma)?;
                        Some(self.wide_unsigned("Hilbert cohomology denominator")?)
                    };
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    let input = Box::new(input);
                    let expression = match denominator {
                        None => Expression::ToFock(k, input),
                        Some(d) if name == "rvHilbCohQ" => Expression::FromFock(k, d, input),
                        Some(d) => Expression::Scale(k, d, input),
                    };
                    (expression, depth)
                }
            }
        };
        expr(
            span,
            ExprForm::HilbertCohomology(Box::new(expression)),
            depth,
        )
    }

    fn hilbert_tautological_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::nakajima::tautological::Command;
        self.advance()?;
        self.expect(Token::LBracket)?;
        let length = self.wide_unsigned("Hilbert length")?;
        let codimension = if name == "rvHilbTautChernDen" {
            self.expect(Token::Comma)?;
            Some(self.wide_unsigned("Chern index")?)
        } else {
            None
        };
        self.expect(Token::RBracket)?;
        let input = self.paren_expr(0)?;
        Ok(StatementForm::HilbertTautological(Box::new(
            match codimension {
                Some(j) => Command::ChernDenominator(length, j, input),
                None => Command::TopSegre(length, input),
            },
        )))
    }

    fn nakajima_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::nakajima::{Expression, Factor, Operation};
        if name == "rvHilbTautChern" {
            self.expect(Token::LBracket)?;
            let length = self.wide_unsigned("Hilbert length")?;
            self.expect(Token::Comma)?;
            let codimension = self.wide_unsigned("Chern index")?;
            self.expect(Token::Comma)?;
            let scale = self.signed_wide("a signed integer", "Chern scale")?;
            self.expect(Token::RBracket)?;
            let input = self.paren_expr(nest)?;
            let depth = 1 + input.depth;
            let expression =
                Expression::TautologicalChern(length, codimension, scale, Box::new(input));
            return expr(span, ExprForm::Nakajima(Box::new(expression)), depth);
        }
        if name == "rvNakajimaCup" {
            let inputs = self.paren_pair(nest)?;
            let depth = 1 + inputs[0].depth.max(inputs[1].depth);
            return expr(
                span,
                ExprForm::Nakajima(Box::new(Expression::Cup(Box::new(inputs)))),
                depth,
            );
        }
        if name == "rvNakajimaBoundary" {
            let input = self.paren_expr(nest)?;
            let depth = input.depth + 1;
            let expression = Expression::Apply(Operation::Boundary, Box::new(input));
            return expr(span, ExprForm::Nakajima(Box::new(expression)), depth);
        }
        self.expect(Token::LBracket)?;
        let (expression, depth) = match name {
            "rvNakajimaStratum" => {
                let shape = self.label()?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Stratum(shape, Box::new(input)), depth)
            }
            "rvNakajimaPrimaryGen" => {
                let (basis, length, degree) = self.nakajima_grade_args()?;
                self.expect(Token::Comma)?;
                let index = self.number()?.0;
                self.expect(Token::RBracket)?;
                (
                    Expression::PrimaryGenerator(basis, length, degree, index),
                    0,
                )
            }
            "rvNakajimaVac" | "rvNakajimaMonomial" => {
                let basis = self.label()?;
                self.expect(Token::RBracket)?;
                let expression = if name == "rvNakajimaVac" {
                    Expression::Vacuum(basis)
                } else {
                    let factors = self.set_named_rows(|parser| {
                        parser.expect(Token::LParen)?;
                        let mode = parser.nakajima_mode()?;
                        parser.expect(Token::Comma)?;
                        let exponent = parser.wide_unsigned("Nakajima multiplicity")?;
                        parser.expect(Token::RParen)?;
                        Ok(Factor { mode, exponent })
                    })?;
                    Expression::Monomial(basis, factors)
                };
                (expression, 0)
            }
            _ => {
                let operation = match name {
                    "rvNakajimaCreate" => Operation::Create(self.nakajima_mode()?),
                    "rvNakajimaAnn" => Operation::Annihilate(self.nakajima_mode()?),
                    "rvVirasoro" => Operation::Virasoro(
                        self.signed_wide("a signed Virasoro index", "Virasoro index")?,
                    ),
                    _ => Operation::Scale(
                        self.signed_wide("a signed Nakajima scalar", "Nakajima scalar")?,
                    ),
                };
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Apply(operation, Box::new(input)), depth)
            }
        };
        expr(span, ExprForm::Nakajima(Box::new(expression)), depth)
    }

    fn hilbert_stratum_statement(&mut self, operation: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::{hilbert::SurfaceKind, stratum::Command};
        self.advance()?;
        let command = if operation == "rvhilbstratum" {
            let name = self.label()?;
            let parts = self.set_named_rows(|parser| {
                parser.expect(Token::LParen)?;
                let size = parser.wide_unsigned("Hilbert–Chow part size")?;
                parser.expect(Token::Comma)?;
                let count = parser.wide_unsigned("Hilbert–Chow multiplicity")?;
                parser.expect(Token::RParen)?;
                Ok((size, count))
            })?;
            Command::Declare(name, parts)
        } else {
            self.expect(Token::LBracket)?;
            let name = self.label()?;
            let command = if operation == "rvHilbStratumData" {
                Command::Data(name)
            } else {
                self.expect(Token::Comma)?;
                let order = if operation == "rvGenKummerStratumHodge" {
                    let order = self.wide_unsigned("translation order")?;
                    self.expect(Token::Comma)?;
                    Some(order)
                } else {
                    None
                };
                let p = self.wide_unsigned("first Hodge index")?;
                self.expect(Token::Comma)?;
                let q = self.wide_unsigned("second Hodge index")?;
                if let Some(order) = order {
                    Command::GeneralizedHodge(name, order, p, q)
                } else {
                    let surface = if operation == "rvK3StratumHodge" {
                        SurfaceKind::K3
                    } else {
                        SurfaceKind::Abelian
                    };
                    Command::SurfaceHodge(name, surface, p, q)
                }
            };
            self.expect(Token::RBracket)?;
            command
        };
        Ok(StatementForm::HilbertStratum(Box::new(command)))
    }

    fn specht_stratum_statement(&mut self) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::specht::Declaration;
        self.advance()?;
        let name = self.label()?;
        self.expect(Token::LBracket)?;
        let shape = self.label()?;
        self.expect(Token::RBracket)?;
        let diagrams = self.set_named_rows(|parser| {
            parser.expect(Token::LParen)?;
            let size = parser.wide_unsigned("Hilbert–Chow part size")?;
            parser.expect(Token::Comma)?;
            parser.expect(Token::LBracket)?;
            let mut rows = Vec::new();
            while parser.token() != Token::RBracket {
                rows.push(parser.wide_unsigned("Specht row length")?);
                if !parser.eat(Token::Comma)? {
                    break;
                }
            }
            parser.expect(Token::RBracket)?;
            parser.expect(Token::RParen)?;
            Ok((size, rows))
        })?;
        Ok(StatementForm::DeclareSpecht(Box::new(Declaration(
            name, shape, diagrams,
        ))))
    }

    fn witt_cyclic_power(&mut self, cyclic: bool) -> Result<crate::rv_witt_spectrum::CyclicPower> {
        use crate::rv_witt_spectrum::CyclicPower;
        let n = self.wide_unsigned("cyclic tensor length")?;
        if cyclic {
            self.expect(Token::Comma)?;
            let d = self.wide_unsigned("cyclotomic order")?;
            Ok(CyclicPower::Isotypic(n, d))
        } else {
            Ok(CyclicPower::Lie(n))
        }
    }

    fn witt_spectrum_expression(
        &mut self,
        name: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::rv_witt_spectrum::{Bridge, Expression, HilbertBridge};
        let iteration = |parser: &mut Self| {
            let span = parser.span();
            let power = parser.wide_unsigned("K3 Frobenius iteration count")?;
            if power == 0 {
                return Err(Error::new(
                    span,
                    "K3 Frobenius iteration count must be positive",
                ));
            }
            Ok(power)
        };
        let (form, depth) = match name {
            "rvWittSpectrum" => {
                self.expect(Token::LBrace)?;
                let coefficients = |parser: &mut Self| {
                    parser.expect(Token::LBracket)?;
                    let mut out = Vec::new();
                    while parser.token() != Token::RBracket {
                        out.push(
                            parser.signed_wide("a signed Witt coefficient", "Witt coefficient")?,
                        );
                        if !parser.eat(Token::Comma)? {
                            break;
                        }
                    }
                    parser.expect(Token::RBracket)?;
                    Ok(out)
                };
                let p = coefficients(self)?;
                let q = if self.eat(Token::Comma)? {
                    Some(coefficients(self)?)
                } else {
                    None
                };
                self.expect(Token::RBrace)?;
                (Expression::Literal(p, q), 0)
            }
            "rvWittSpectrumLine" => {
                self.expect(Token::LBracket)?;
                let a = self.signed_wide("a signed eigenvalue", "Witt line eigenvalue")?;
                self.expect(Token::RBracket)?;
                (Expression::Line(a), 0)
            }
            "rvWittSpectrumTensor" => {
                let args = self.paren_pair(nest)?;
                let depth = 1 + args[0].depth.max(args[1].depth);
                (Expression::Tensor(Box::new(args)), depth)
            }
            "rvWittSpectrumAdams"
            | "rvWittSpectrumExterior"
            | "rvWittSpectrumVerschiebung"
            | "rvWittSpectrumNorm" => {
                self.expect(Token::LBracket)?;
                let index = self.wide_unsigned("Witt operation index")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                let input = Box::new(input);
                let form = match name {
                    "rvWittSpectrumAdams" => Expression::Adams(index, input),
                    "rvWittSpectrumExterior" => Expression::Exterior(index, input),
                    "rvWittSpectrumVerschiebung" => Expression::Verschiebung(index, input),
                    _ => Expression::Norm(index, input),
                };
                (form, depth)
            }
            "rvWittSpectrumCyclic" | "rvWittSpectrumLie" => {
                self.expect(Token::LBracket)?;
                let projection = self.witt_cyclic_power(name == "rvWittSpectrumCyclic")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::CyclicPower(projection, Box::new(input)), depth)
            }
            "rvK3FrobHilbSpectrum" => {
                self.expect(Token::LBracket)?;
                let word = self.label()?;
                self.expect(Token::Comma)?;
                let base = self.label()?;
                self.expect(Token::Comma)?;
                let length = self.wide_unsigned("K3 Hilbert-scheme length")?;
                self.expect(Token::Comma)?;
                let r = iteration(self)?;
                self.expect(Token::RBracket)?;
                (
                    Expression::HilbertBridge(HilbertBridge {
                        word,
                        base,
                        length,
                        r,
                    }),
                    0,
                )
            }
            _ => {
                self.expect(Token::LBracket)?;
                let word = self.label()?;
                self.expect(Token::Comma)?;
                let base = self.label()?;
                self.expect(Token::Comma)?;
                let stratum = self.label()?;
                self.expect(Token::Comma)?;
                let r = iteration(self)?;
                self.expect(Token::Comma)?;
                let i = self.wide_unsigned("cohomological degree")?;
                self.expect(Token::RBracket)?;
                let specht = name == "rvK3SpechtSpectrum";
                (
                    Expression::Bridge(Bridge {
                        word,
                        base,
                        stratum,
                        r,
                        i,
                        specht,
                    }),
                    0,
                )
            }
        };
        expr(span, ExprForm::WittSpectrum(Box::new(form)), depth)
    }

    fn kummer_node_lattice_statement(&mut self) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::nikulin::Command;
        self.advance()?;
        self.expect(Token::LBracket)?;
        let basis = self.label()?;
        self.expect(Token::Comma)?;
        let base = self.label()?;
        self.expect(Token::Comma)?;
        let power = self.number()?.0;
        self.expect(Token::RBracket)?;
        Ok(StatementForm::KummerNodeLattice(Box::new(Command(
            basis, base, power,
        ))))
    }

    fn kummer_weil_expression(
        &mut self,
        name: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::weil::{Expression, Operation};
        let (expression, depth) = match name {
            "rvKummerThetaShell" => {
                self.expect(Token::LBracket)?;
                let marking = self.label()?;
                self.expect(Token::Comma)?;
                let norm = self.number()?.0;
                self.expect(Token::RBracket)?;
                (Expression::Shell(marking, norm), 0)
            }
            "rvKummerWeilFrob" => {
                self.expect(Token::LBracket)?;
                let base = self.label()?;
                self.expect(Token::Comma)?;
                let power = self.number()?.0;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Frobenius(base, power, Box::new(input)), depth)
            }
            "rvKummerWeilBasis" => {
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Basis(Box::new(input)), depth)
            }
            _ => {
                let operation = match name {
                    "rvKummerWeilScale" => {
                        self.expect(Token::LBracket)?;
                        let numerator =
                            self.signed_wide("a signed numerator", "Kummer Weil scalar numerator")?;
                        self.expect(Token::Comma)?;
                        let denominator = self.wide_unsigned("Kummer Weil scalar denominator")?;
                        self.expect(Token::RBracket)?;
                        Operation::Scale(numerator, denominator)
                    }
                    "rvKummerWeilS" => Operation::S,
                    _ => Operation::T,
                };
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Apply(operation, Box::new(input)), depth)
            }
        };
        expr(span, ExprForm::KummerWeil(Box::new(expression)), depth)
    }

    fn kummer_disc_expression(
        &mut self,
        name: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::nikulin::{Expression, ORDER};
        let (expression, depth) = match name {
            "rvKummerDisc" => {
                self.expect(Token::LBracket)?;
                let basis = self.label()?;
                self.expect(Token::Comma)?;
                let (mask, at) = self.number()?;
                if mask >= usize::from(ORDER) {
                    return Err(Error::new(at, "rvKummerDisc mask must be between 0 and 63"));
                }
                self.expect(Token::RBracket)?;
                (Expression::Element(basis, mask as u8), 0)
            }
            "rvKummerDiscPair" => {
                let inputs = self.paren_pair(nest)?;
                let depth = 1 + inputs[0].depth.max(inputs[1].depth);
                (Expression::Pair(Box::new(inputs)), depth)
            }
            "rvKummerDiscFrob" => {
                self.expect(Token::LBracket)?;
                let base = self.label()?;
                self.expect(Token::Comma)?;
                let power = self.number()?.0;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Frobenius(base, power, Box::new(input)), depth)
            }
            _ => {
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Quadratic(Box::new(input)), depth)
            }
        };
        expr(span, ExprForm::KummerDisc(Box::new(expression)), depth)
    }

    fn generalized_kummer_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::generalized_kummer::Command;
        self.advance()?;
        self.expect(Token::LBracket)?;
        if matches!(name, "rvGenKummerBetti" | "rvJacHilbBetti") {
            let span = self.span();
            let length = self.wide_unsigned("Jacobian Hilbert-scheme length")?;
            self.expect(Token::Comma)?;
            let i = self.wide_unsigned("cohomological degree")?;
            self.expect(Token::RBracket)?;
            let command = if name == "rvGenKummerBetti" {
                if length == 0 {
                    return Err(Error::new(
                        span,
                        "generalized Kummer length must be positive",
                    ));
                }
                Command::Betti(length, i)
            } else {
                Command::HilbertBetti(length, i)
            };
            return Ok(StatementForm::GeneralizedKummer(Box::new(command)));
        }
        if matches!(
            name,
            "rvGenKummerTranslateChar" | "rvGenKummerTranslateFixed" | "rvGenKummerInversionEuler"
        ) {
            let length = self.wide_unsigned("generalized Kummer length")?;
            self.expect(Token::Comma)?;
            let order = self.wide_unsigned("translation order")?;
            let command = if name == "rvGenKummerTranslateChar" {
                self.expect(Token::Comma)?;
                Command::TranslateChar(length, order, self.wide_unsigned("cohomological degree")?)
            } else if name == "rvGenKummerInversionEuler" {
                Command::InversionEuler(length, order)
            } else {
                Command::TranslateFixed(length, order)
            };
            self.expect(Token::RBracket)?;
            return Ok(StatementForm::GeneralizedKummer(Box::new(command)));
        }
        let first = self.label()?;
        self.expect(Token::Comma)?;
        let positive_degree = |parser: &mut Self| {
            let span = parser.span();
            let degree = parser.wide_unsigned("Jacobian Hilbert-scheme extension degree")?;
            if degree == 0 {
                return Err(Error::new(span, "extension degree must be positive"));
            }
            Ok(degree)
        };
        let command = if name == "rvGenKummerStratumFrob" {
            let base = self.label()?;
            self.expect(Token::Comma)?;
            let shape = self.label()?;
            self.expect(Token::Comma)?;
            let degree = positive_degree(self)?;
            self.expect(Token::Comma)?;
            let i = self.wide_unsigned("cohomological degree")?;
            self.expect(Token::RBracket)?;
            Command::StratumFrob(first, base, shape, degree, i, Box::new(self.paren_expr(0)?))
        } else if matches!(
            name,
            "rvGenKummerLevel"
                | "rvGenKummerTrace"
                | "rvGenKummerTranslateFrob"
                | "rvGenKummerTranslatePoints"
                | "rvGenKummerInversionFrob"
                | "rvGenKummerInversionPoints"
        ) {
            let base = self.label()?;
            self.expect(Token::Comma)?;
            let degree = positive_degree(self)?;
            let i = if matches!(
                name,
                "rvGenKummerTrace" | "rvGenKummerTranslateFrob" | "rvGenKummerInversionFrob"
            ) {
                self.expect(Token::Comma)?;
                Some(self.wide_unsigned("cohomological degree")?)
            } else {
                None
            };
            self.expect(Token::RBracket)?;
            match name {
                "rvGenKummerTrace" => Command::Trace(first, base, degree, i.unwrap()),
                "rvGenKummerTranslateFrob" => Command::TranslateFrob(
                    first,
                    base,
                    degree,
                    i.unwrap(),
                    Box::new(self.paren_expr(0)?),
                ),
                "rvGenKummerTranslatePoints" => {
                    Command::TranslatePoints(first, base, degree, Box::new(self.paren_expr(0)?))
                }
                "rvGenKummerInversionFrob" => Command::InversionFrob(
                    first,
                    base,
                    degree,
                    i.unwrap(),
                    Box::new(self.paren_expr(0)?),
                ),
                "rvGenKummerInversionPoints" => {
                    Command::InversionPoints(first, base, degree, Box::new(self.paren_expr(0)?))
                }
                _ => Command::Level(first, base, degree),
            }
        } else {
            let span = self.span();
            let length = self.wide_unsigned("Jacobian Hilbert-scheme length")?;
            if matches!(name, "rvJacHilbPoints" | "rvJacHilbTrace") {
                self.expect(Token::Comma)?;
                let degree = positive_degree(self)?;
                let command = if name == "rvJacHilbTrace" {
                    self.expect(Token::Comma)?;
                    Command::HilbertTrace(
                        first,
                        length,
                        degree,
                        self.wide_unsigned("cohomological degree")?,
                    )
                } else {
                    Command::Hilbert(first, length, degree)
                };
                self.expect(Token::RBracket)?;
                command
            } else {
                if length == 0 {
                    return Err(Error::new(span, "summation fiber length must be positive"));
                }
                self.expect(Token::RBracket)?;
                if name == "rvGenKummerFiber" {
                    Command::Fiber(first, length, Box::new(self.paren_expr(0)?))
                } else {
                    Command::Zero(first, length)
                }
            }
        };
        Ok(StatementForm::GeneralizedKummer(Box::new(command)))
    }

    fn k3_mordell_statement(&mut self, operation: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::mordell::Command;
        self.advance()?;
        let command = if operation == "rvk3mordell" {
            let name = self.label()?;
            self.expect(Token::LBracket)?;
            let pencil = self.label()?;
            let zero = self.kummer_theta_index("K3 Mordell–Weil zero trope mask", 15)? as u8;
            self.expect(Token::RBracket)?;
            Command::Define(name, pencil, zero)
        } else {
            self.expect(Token::LBracket)?;
            let name = self.label()?;
            if operation == "rvK3MordellBSDFactors" {
                let mask = self.kummer_theta_index("K3 BSD section trope mask", 15)? as u8;
                self.expect(Token::Comma)?;
                let multiple =
                    self.signed_wide("a signed section multiple", "K3 BSD section multiple")?;
                self.expect(Token::RBracket)?;
                return Ok(StatementForm::K3Mordell(Box::new(Command::BsdFactors(
                    name, mask, multiple,
                ))));
            }
            self.expect(Token::RBracket)?;
            match operation {
                "rvK3MordellCoord" => Command::Coord(name, Box::new(self.paren_expr(0)?)),
                "rvK3MordellHeight" => Command::Height(name, Box::new(self.paren_pair(0)?)),
                _ => Command::Lattice(name),
            }
        };
        Ok(StatementForm::K3Mordell(Box::new(command)))
    }

    fn k3_mordell_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::mordell::{Expression, Shift};
        self.expect(Token::LBracket)?;
        let marking = self.label()?;
        self.expect(Token::Comma)?;
        let steps = self.signed_wide(
            "a signed translation count",
            "K3 Mordell–Weil translation count",
        )?;
        self.expect(Token::RBracket)?;
        let input = self.paren_expr(nest)?;
        let depth = 1 + input.depth;
        let shift = if name == "rvK3MordellShiftDiv" {
            Shift::Divisor
        } else {
            Shift::Point
        };
        expr(
            span,
            ExprForm::K3Mordell(Box::new(Expression(marking, shift, steps, Box::new(input)))),
            depth,
        )
    }

    fn k3_fiber_equation_statement(&mut self, operation: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::pencil::fiber::{Command, Query};
        self.advance()?;
        let target = (operation == "rvk3fiberjac")
            .then(|| self.label())
            .transpose()?;
        self.expect(Token::LBracket)?;
        let pencil = self.label()?;
        let query = if let Some(target) = target {
            Query::Register(target)
        } else {
            match operation {
                "rvK3FiberQuartic" => Query::Quartic,
                "rvK3FiberJacobian" => Query::Jacobian,
                _ => {
                    self.expect(Token::Comma)?;
                    let span = self.span();
                    let degree = self.wide_unsigned("K3 fiber extension degree")?;
                    if degree == 0 {
                        return Err(Error::new(
                            span,
                            "K3 fiber extension degree must be positive",
                        ));
                    }
                    Query::Points(degree)
                }
            }
        };
        self.expect(Token::RBracket)?;
        let parameter = self.paren_expr(0)?;
        Ok(StatementForm::K3FiberEquation(Box::new(Command(
            pencil,
            query,
            Box::new(parameter),
        ))))
    }

    fn k3_pencil_statement(&mut self, operation: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::pencil::Command;
        self.advance()?;
        let command = if operation == "rvk3pencil" {
            let name = self.label()?;
            self.expect(Token::LBracket)?;
            let basis = self.label()?;
            let first = self.kummer_theta_index("K3 pencil node mask", 15)? as u8;
            let second = self.kummer_theta_index("K3 pencil node mask", 15)? as u8;
            self.expect(Token::RBracket)?;
            Command::Define(name, basis, first, second)
        } else {
            self.expect(Token::LBracket)?;
            let name = self.label()?;
            self.expect(Token::RBracket)?;
            match operation {
                "rvK3PencilFiber" => Command::Fiber(name, Box::new(self.paren_expr(0)?)),
                "rvK3FiberReduction" => Command::Reduction(name, Box::new(self.paren_expr(0)?)),
                "rvK3FiberEuler" => Command::Euler(name, Box::new(self.paren_expr(0)?)),
                "rvK3PencilLFunction" => Command::LFunction(name),
                "rvK3PencilBSD" => Command::Bsd(name),
                "rvK3PencilBSDFactors" => Command::BsdFactors(name),
                "rvK3PencilTorsion2" => Command::Torsion2(name),
                "rvK3PencilSelmer" => Command::Selmer(name),
                "rvK3PencilWeierstrass" => Command::Weierstrass(name),
                "rvK3PencilDiscriminant" => Command::Discriminant(name),
                _ => Command::Census(name),
            }
        };
        Ok(StatementForm::K3Pencil(Box::new(command)))
    }

    fn k3_pencil_section_expression(
        &mut self,
        name: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::pencil::sections::Expression;
        let (expression, depth) = match name {
            "rvK3PencilSectionPoint" => {
                self.expect(Token::LBracket)?;
                let pencil = self.label()?;
                self.expect(Token::RBracket)?;
                let coefficient = |parser: &mut Self| {
                    parser.signed_tuple(
                        "a signed residue coordinate",
                        "K3 pencil section coefficient",
                    )
                };
                let coordinates = [
                    self.set_named_rows(coefficient)?,
                    self.set_named_rows(coefficient)?,
                    self.set_named_rows(coefficient)?,
                    self.set_named_rows(coefficient)?,
                ];
                (Expression::Literal(pencil, Box::new(coordinates)), 0)
            }
            "rvK3PencilTorsion2At" | "rvK3PencilSelmerBasis" => {
                self.expect(Token::LBracket)?;
                let pencil = self.label()?;
                self.expect(Token::Comma)?;
                let index = self.number()?.0;
                self.expect(Token::RBracket)?;
                let expression = if name == "rvK3PencilSelmerBasis" {
                    Expression::SelmerBasis(pencil, index)
                } else {
                    Expression::Torsion(pencil, index)
                };
                (expression, 0)
            }
            _ => {
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                let expression = if name == "rvK3PencilSelmerContains" {
                    Expression::SelmerContains(Box::new(input))
                } else {
                    Expression::Descent(Box::new(input), name == "rvK3PencilDivisible2")
                };
                (expression, depth)
            }
        };
        expr(span, ExprForm::K3PencilSection(Box::new(expression)), depth)
    }

    fn k3_pencil_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::pencil::{Expression, Operation};
        self.expect(Token::LBracket)?;
        let pencil = self.label()?;
        let mask = if name == "rvK3PencilSection" {
            self.kummer_theta_index("K3 pencil trope mask", 15)? as u8
        } else {
            0
        };
        self.expect(Token::RBracket)?;
        let (operation, depth) = match name {
            "rvK3PencilDiv" => (Operation::Divisor, 0),
            "rvK3PencilParam" => {
                let inputs = self.paren_pair(nest)?;
                let depth = 1 + inputs[0].depth.max(inputs[1].depth);
                (Operation::Param(Box::new(inputs)), depth)
            }
            _ => {
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                let operation = match name {
                    "rvK3PencilSection" => Operation::Section(mask, Box::new(input)),
                    "rvK3FiberJ" => Operation::FiberJ(Box::new(input)),
                    _ => Operation::Map(Box::new(input)),
                };
                (operation, depth)
            }
        };
        expr(
            span,
            ExprForm::K3Pencil(Box::new(Expression(pencil, operation))),
            depth,
        )
    }

    fn declare_enriques_statement(&mut self) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::enriques::points::Declare;
        self.advance()?;
        let name = self.label()?;
        self.expect(Token::LBracket)?;
        let base = self.label()?;
        self.expect(Token::Comma)?;
        let cover = self.label()?;
        self.expect(Token::Comma)?;
        let basis = self.label()?;
        let mask = self.kummer_theta_index("Enriques switch mask", 15)?;
        self.expect(Token::RBracket)?;
        Ok(StatementForm::DeclareEnriques(Box::new(Declare(
            name, base, cover, basis, mask,
        ))))
    }

    fn enriques_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::enriques::points::{Expression, Read};
        let (expression, depth) = if name == "rvEnriquesOrbit" {
            self.expect(Token::LBracket)?;
            let quotient = self.label()?;
            self.expect(Token::RBracket)?;
            let input = self.paren_expr(nest)?;
            let depth = input.depth + 1;
            (Expression::Orbit(quotient, Box::new(input)), depth)
        } else {
            let operation = match name {
                "rvEnriquesSplit" => Read::Split,
                "rvEnriquesRepresentative" => Read::Representative,
                _ => Read::Lift,
            };
            let input = self.paren_expr(nest)?;
            let depth = input.depth + 1;
            (Expression::Read(operation, Box::new(input)), depth)
        };
        expr(span, ExprForm::Enriques(Box::new(expression)), depth)
    }

    fn kummer_enriques_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::enriques::{Command, Query};
        self.advance()?;
        self.expect(Token::LBracket)?;
        let basis = self.label()?;
        self.expect(Token::Comma)?;
        let curve = self.label()?;
        let mask = self.kummer_theta_index("Enriques switch mask", 15)?;
        let query = if name == "rvKummerEnriquesZeta" {
            Query::Zeta
        } else {
            self.expect(Token::Comma)?;
            let span = self.span();
            let degree = self.wide_unsigned("Enriques extension degree")?;
            if degree == 0 {
                return Err(Error::new(
                    span,
                    "Enriques extension degree must be positive",
                ));
            }
            Query::Points(degree)
        };
        self.expect(Token::RBracket)?;
        Ok(StatementForm::KummerEnriques(Box::new(Command(
            basis, curve, mask, query,
        ))))
    }

    fn kummer_k3_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::k3::{Command, Query};
        self.advance()?;
        self.expect(Token::LBracket)?;
        let curve = self.label()?;
        let query = if name == "rvKummerK3Zeta" {
            Query::Zeta
        } else {
            self.expect(Token::Comma)?;
            let span = self.span();
            let degree = self.wide_unsigned("Kummer K3 extension degree")?;
            if degree == 0 {
                return Err(Error::new(
                    span,
                    "Kummer K3 extension degree must be positive",
                ));
            }
            match name {
                "rvKummerK3Points" => Query::Points(degree),
                "rvKummerK3ArtinTate" => Query::ArtinTate(degree),
                _ => Query::Picard(degree),
            }
        };
        self.expect(Token::RBracket)?;
        Ok(StatementForm::KummerK3(Box::new(Command(curve, query))))
    }

    fn crystalline_newton_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::newton::{Command, Kind};
        self.advance()?;
        self.expect(Token::LBracket)?;
        let curve = self.label()?;
        let command = if name == "rvKummerK3Height" {
            Command::Height(curve)
        } else {
            self.expect(Token::Comma)?;
            let kind = match name {
                "rvKummerHilbNewton" => Kind::K3Hilbert,
                "rvJacHilbNewton" => Kind::JacobianHilbert,
                _ => Kind::Generalized,
            };
            let length = self.wide_unsigned("Hilbert-scheme length")?;
            self.expect(Token::Comma)?;
            let i = self.wide_unsigned("cohomological degree")?;
            Command::Polygon(curve, kind, length, i)
        };
        self.expect(Token::RBracket)?;
        Ok(StatementForm::CrystallineNewton(command))
    }

    fn hodge_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::generalized_kummer::hodge::{Command, Operation};
        self.advance()?;
        self.expect(Token::LBracket)?;
        let n = self.wide_unsigned("Hilbert-scheme length")?;
        self.expect(Token::Comma)?;
        let operation = match name {
            "rvK3HilbHodge" => Operation::K3Hilbert,
            "rvJacHilbHodge" => Operation::JacobianHilbert,
            "rvGenKummerHodge" => Operation::Generalized,
            _ => {
                let h = self.wide_unsigned("translation order")?;
                self.expect(Token::Comma)?;
                match name {
                    "rvGenKummerTranslateHodge" => Operation::Translate(h),
                    "rvGenKummerInversionHodge" => Operation::Invert(h),
                    _ => Operation::Spectrum(h),
                }
            }
        };
        let p = self.wide_unsigned("first Hodge index")?;
        self.expect(Token::Comma)?;
        let q = self.wide_unsigned("second Hodge index")?;
        self.expect(Token::RBracket)?;
        Ok(StatementForm::Hodge(Command(operation, n, p, q)))
    }

    fn hilbert_llv_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::llv_character::Command;
        self.advance()?;
        let shape = |parser: &mut Self| {
            parser.set_named_rows(|parser| parser.wide_unsigned("highest weight coordinate"))
        };
        let command = if name == "rvHilbLLVDimension" {
            Command::Dimension(shape(self)?)
        } else {
            self.expect(Token::LBracket)?;
            let n = self.wide_unsigned("Hilbert length")?;
            match name {
                "rvHilbReducedLLVMultiplicity"
                | "rvHilbReducedLLVDecomposition"
                | "rvHilbPrimitiveLLVMultiplicity"
                | "rvHilbPrimitiveLLVDecomposition" => {
                    self.expect(Token::Comma)?;
                    let d = self.signed_wide("a signed degree", "cohomological degree")?;
                    self.expect(Token::RBracket)?;
                    match name {
                        "rvHilbReducedLLVMultiplicity" => {
                            Command::ReducedMultiplicity(n, d, shape(self)?)
                        }
                        "rvHilbReducedLLVDecomposition" => Command::ReducedDecompose(n, d),
                        "rvHilbPrimitiveLLVMultiplicity" => {
                            let weight = self.set_named_rows(|parser| {
                                parser.signed_wide(
                                    "a signed highest weight coordinate",
                                    "highest weight coordinate",
                                )
                            })?;
                            Command::PrimitiveMultiplicity(n, d, weight)
                        }
                        _ => Command::PrimitiveDecompose(n, d),
                    }
                }
                "rvHilbLLVHodge" | "rvHilbLLVHodgeDecomposition" => {
                    self.expect(Token::Comma)?;
                    let p = self.signed_wide("a signed Hodge index", "first Hodge degree")?;
                    self.expect(Token::Comma)?;
                    let q = self.signed_wide("a signed Hodge index", "second Hodge degree")?;
                    self.expect(Token::RBracket)?;
                    if name == "rvHilbLLVHodge" {
                        Command::Hodge(n, p, q, shape(self)?)
                    } else {
                        Command::HodgeDecompose(n, p, q)
                    }
                }
                _ => {
                    self.expect(Token::RBracket)?;
                    match name {
                        "rvHilbLLVMultiplicity" => Command::Multiplicity(n, shape(self)?),
                        "rvHilbLLVHodgeTable" => {
                            let weight = if self.token() == Token::LBrace {
                                Some(shape(self)?)
                            } else {
                                None
                            };
                            Command::HodgeTable(n, weight)
                        }
                        _ => Command::Decompose(n),
                    }
                }
            }
        };
        Ok(StatementForm::HilbertLlv(Box::new(command)))
    }

    fn kummer_hilbert_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::hilbert::{elliptic_genus::Query, Command};
        self.advance()?;
        if name == "rvKummerHilbChar" {
            let query = self.kummer_hilbert_characteristic()?;
            return Ok(StatementForm::KummerHilbert(Box::new(
                Command::Characteristic(query),
            )));
        }
        self.expect(Token::LBracket)?;
        let command = if matches!(
            name,
            "rvK3EllipticGenus" | "rvGenKummerElliptic" | "rvIgusaPhi10"
        ) {
            let first = if name == "rvGenKummerElliptic" {
                "generalized Kummer length"
            } else {
                "P exponent"
            };
            let index = self.elliptic_genus_index(first)?;
            let query = match name {
                "rvK3EllipticGenus" => Query::EllipticGenus(index),
                "rvGenKummerElliptic" => Query::Generalized(index),
                _ => Query::Igusa(index),
            };
            Command::EllipticGenus(query)
        } else if name == "rvK3SymplecticElliptic" {
            let (order, span) = self.number()?;
            if !(1..=8).contains(&order) {
                return Err(Error::new(
                    span,
                    "symplectic K3 order must be between 1 and 8",
                ));
            }
            self.expect(Token::Comma)?;
            Command::EllipticGenus(Query::Symplectic(
                order,
                self.elliptic_genus_index("P exponent")?,
            ))
        } else if name == "rvKummerTranslateElliptic" {
            let basis = self.label()?;
            let mask = self.kummer_theta_index("Kummer translation mask", 15)?;
            self.expect(Token::Comma)?;
            Command::EllipticGenus(Query::Translate(
                basis,
                mask,
                self.elliptic_genus_index("P exponent")?,
            ))
        } else if name == "rvKummerEnriquesElliptic" {
            let basis = self.label()?;
            self.expect(Token::Comma)?;
            let base = self.label()?;
            let mask = self.kummer_theta_index("Enriques switch mask", 15)?;
            self.expect(Token::Comma)?;
            Command::EllipticGenus(Query::Enriques(
                basis,
                base,
                mask,
                self.elliptic_genus_index("P exponent")?,
            ))
        } else if name == "rvKummerHilbBetti" {
            let length = self.wide_unsigned("Kummer Hilbert-scheme length")?;
            self.expect(Token::Comma)?;
            let degree = self.wide_unsigned("cohomological degree")?;
            Command::Betti(length, degree)
        } else {
            let curve = self.label()?;
            self.expect(Token::Comma)?;
            let length = self.wide_unsigned("Kummer Hilbert-scheme length")?;
            self.expect(Token::Comma)?;
            let span = self.span();
            let extension = self.wide_unsigned("Kummer Hilbert-scheme extension degree")?;
            if extension == 0 {
                return Err(Error::new(
                    span,
                    "Kummer Hilbert-scheme extension degree must be positive",
                ));
            }
            if name == "rvKummerHilbPoints" {
                Command::Points(curve, length, extension)
            } else if name == "rvKummerHilbWeilRestrictPoints" {
                self.expect(Token::Comma)?;
                let index = self.wide_unsigned("point-count extension degree")?;
                Command::WeilRestriction(curve, length, extension, index)
            } else {
                self.expect(Token::Comma)?;
                let degree = self.wide_unsigned("cohomological degree")?;
                Command::Trace(curve, length, extension, degree)
            }
        };
        self.expect(Token::RBracket)?;
        Ok(StatementForm::KummerHilbert(Box::new(command)))
    }

    fn elliptic_genus_index(
        &mut self,
        first: &str,
    ) -> Result<crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::Index> {
        use crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::Index;
        let p = self.wide_unsigned(first)?;
        self.expect(Token::Comma)?;
        let q = self.wide_unsigned("q exponent")?;
        self.expect(Token::Comma)?;
        let r = self.signed_wide("a signed Laurent exponent", "y exponent")?;
        Ok(Index(p, q, r))
    }

    fn kummer_divisor_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hyperelliptic_kummer_surface::polarized::{Command, Fujiki};
        self.advance()?;
        let command = if name == "rvKummerFujiki" {
            self.expect(Token::LBracket)?;
            let dimension = self.number()?.0;
            self.expect(Token::Comma)?;
            let right_exponent = self.number()?.0;
            self.expect(Token::Comma)?;
            let left_shift = self.signed_wide("an integer shift", "Kummer Fujiki shift")?;
            self.expect(Token::Comma)?;
            let right_shift = self.signed_wide("an integer shift", "Kummer Fujiki shift")?;
            self.expect(Token::RBracket)?;
            Command::Fujiki(
                Fujiki {
                    dimension,
                    right_exponent,
                    left_shift,
                    right_shift,
                },
                Box::new(self.paren_pair(0)?),
            )
        } else {
            Command::Pair(Box::new(self.paren_pair(0)?))
        };
        Ok(StatementForm::KummerDiv(Box::new(command)))
    }

    fn kummer_divisor_expression(
        &mut self,
        name: &str,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::polarized::{Expression, Generator, COEFFICIENTS};
        self.expect(Token::LBracket)?;
        let (expression, depth) = match name {
            "rvKummerDivScale" => {
                let n = self.signed_wide("an integer multiple", "Kummer divisor scalar")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::Scale(n, Box::new(input)), depth)
            }
            "rvKummerDivFrob" => {
                let base = self.label()?;
                self.expect(Token::Comma)?;
                let power = self.number()?.0;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::Frobenius(base, power, Box::new(input)), depth)
            }
            "rvKummerDiv" => {
                let basis = self.label()?;
                self.expect(Token::RBracket)?;
                self.expect(Token::LBrace)?;
                let mut coefficients = [0; COEFFICIENTS];
                for (i, coefficient) in coefficients.iter_mut().enumerate() {
                    if i != 0 {
                        self.expect(Token::Comma)?;
                    }
                    *coefficient = self.signed_wide(
                        "a signed Kummer divisor coefficient",
                        "Kummer divisor coefficient",
                    )?;
                }
                self.expect(Token::RBrace)?;
                (Expression::Literal(basis, Box::new(coefficients)), 0)
            }
            _ => {
                let basis = self.label()?;
                let generator = if name == "rvKummerHyperplane" {
                    Generator::Hyperplane
                } else {
                    self.expect(Token::Comma)?;
                    let (index, span) = self.number()?;
                    if index >= COEFFICIENTS - 1 {
                        return Err(Error::new(
                            span,
                            "marked Kummer node or trope index must be between 0 and 15",
                        ));
                    }
                    if name == "rvKummerNodeDiv" {
                        Generator::Node(index)
                    } else {
                        Generator::Trope(index)
                    }
                };
                self.expect(Token::RBracket)?;
                (Expression::Generator(basis, generator), 0)
            }
        };
        expr(span, ExprForm::KummerDiv(Box::new(expression)), depth)
    }

    fn kummer_hilbert_characteristic(
        &mut self,
    ) -> Result<crate::hyperelliptic_kummer_surface::hilbert::characteristic::Query> {
        use crate::hyperelliptic_kummer_surface::hilbert::characteristic::Query;
        self.expect(Token::LBracket)?;
        let curve = self.label()?;
        self.expect(Token::Comma)?;
        let length = self.wide_unsigned("Kummer Hilbert-scheme length")?;
        self.expect(Token::Comma)?;
        let span = self.span();
        let extension = self.wide_unsigned("Kummer Hilbert-scheme extension degree")?;
        if extension == 0 {
            return Err(Error::new(
                span,
                "Kummer Hilbert-scheme extension degree must be positive",
            ));
        }
        self.expect(Token::Comma)?;
        let degree = self.wide_unsigned("cohomological degree")?;
        self.expect(Token::Comma)?;
        let prime = self.number()?.0;
        self.expect(Token::Comma)?;
        let order = self.number()?.0;
        self.expect(Token::RBracket)?;
        Ok(Query {
            curve,
            length,
            extension,
            degree,
            prime,
            order,
        })
    }

    fn kummer_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_kummer_surface::{
            BaseOperation, Expression, Operation, ResolvedTheta,
        };
        let (expression, depth) = if name == "rvKummer" {
            self.expect(Token::LBracket)?;
            let curve = self.label()?;
            self.expect(Token::RBracket)?;
            self.expect(Token::LParen)?;
            let k0 = self.expr(nest + 1)?;
            self.expect(Token::Comma)?;
            let k1 = self.expr(nest + 1)?;
            self.expect(Token::Comma)?;
            let k2 = self.expr(nest + 1)?;
            self.expect(Token::Comma)?;
            let k3 = self.expr(nest + 1)?;
            self.expect(Token::RParen)?;
            let inputs = [k0, k1, k2, k3];
            let depth = 1 + inputs.iter().map(|input| input.depth).max().unwrap_or(0);
            (Expression::Coordinates(curve, Box::new(inputs)), depth)
        } else if name == "rvKummerExceptional" {
            self.expect(Token::LParen)?;
            let base = self.expr(nest + 1)?;
            let mut inputs = Vec::with_capacity(4);
            for _ in 0..4 {
                self.expect(Token::Comma)?;
                inputs.push(self.expr(nest + 1)?);
            }
            self.expect(Token::RParen)?;
            let depth = 1 + inputs
                .iter()
                .fold(base.depth, |depth, input| depth.max(input.depth));
            let inputs = inputs.try_into().ok().expect("four tangent coordinates");
            (
                Expression::Exceptional(Box::new(base), Box::new(inputs)),
                depth,
            )
        } else if name == "rvKummerNodeInvolution" {
            let inputs = self.paren_pair(nest)?;
            let depth = 1 + inputs[0].depth.max(inputs[1].depth);
            (Expression::NodeInvolution(Box::new(inputs)), depth)
        } else if matches!(name, "rvKummerRichelot" | "rvKummerRichelotDesc") {
            self.expect(Token::LBracket)?;
            let target = self.label()?;
            self.expect(Token::RBracket)?;
            let input = self.paren_expr(nest)?;
            let depth = 1 + input.depth;
            let expression = if name == "rvKummerRichelot" {
                Expression::Richelot(target, Box::new(input))
            } else {
                Expression::RichelotDescend(target, Box::new(input))
            };
            (expression, depth)
        } else if matches!(
            name,
            "rvKummerExtend"
                | "rvKummerDescend"
                | "rvKummerFrob"
                | "rvKummerK3Extend"
                | "rvKummerK3Descend"
                | "rvKummerK3Frob"
        ) {
            self.expect(Token::LBracket)?;
            let curve = self.label()?;
            let operation = match name {
                "rvKummerExtend" | "rvKummerK3Extend" => BaseOperation::Extend,
                "rvKummerDescend" | "rvKummerK3Descend" => BaseOperation::Descend,
                _ => {
                    self.expect(Token::Comma)?;
                    BaseOperation::Frobenius(self.number()?.0)
                }
            };
            self.expect(Token::RBracket)?;
            let input = self.paren_expr(nest)?;
            let depth = 1 + input.depth;
            let expression = if name.starts_with("rvKummerK3") {
                Expression::ResolvedBase(operation, curve, Box::new(input))
            } else {
                Expression::Base(operation, curve, Box::new(input))
            };
            (expression, depth)
        } else if matches!(
            name,
            "rvKummerK3Translate" | "rvKummerTropeParam" | "rvKummerOnTrope"
        ) {
            self.expect(Token::LBracket)?;
            let basis = self.label()?;
            let mask = self.kummer_theta_index("Kummer theta mask", 15)?;
            self.expect(Token::RBracket)?;
            let (operation, depth) = if name == "rvKummerTropeParam" {
                let pair = self.paren_pair(nest)?;
                let depth = 1 + pair[0].depth.max(pair[1].depth);
                (ResolvedTheta::Parameter(Box::new(pair)), depth)
            } else {
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                let operation = if name == "rvKummerK3Translate" {
                    ResolvedTheta::Translate(Box::new(input))
                } else {
                    ResolvedTheta::OnTrope(Box::new(input))
                };
                (operation, depth)
            };
            (Expression::ResolvedTheta(basis, mask, operation), depth)
        } else if matches!(
            name,
            "rvKummerTranslate"
                | "rvKummerThetaEntry"
                | "rvKummerThetaFactor"
                | "rvKummerTropeCoeff"
        ) {
            use crate::hyperelliptic_kummer_surface::theta::Query;
            self.expect(Token::LBracket)?;
            let basis = self.label()?;
            let mask = self.kummer_theta_index("Kummer theta mask", 15)?;
            if name == "rvKummerTranslate" {
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::Translate(basis, mask, Box::new(input)), depth)
            } else {
                let query = match name {
                    "rvKummerThetaEntry" => Query::Entry(
                        mask,
                        self.kummer_theta_index("Kummer theta row", 3)?,
                        self.kummer_theta_index("Kummer theta column", 3)?,
                    ),
                    "rvKummerThetaFactor" => {
                        Query::Factor(mask, self.kummer_theta_index("Kummer theta mask", 15)?)
                    }
                    _ => Query::Trope(
                        mask,
                        self.kummer_theta_index("Kummer trope coefficient", 3)?,
                    ),
                };
                self.expect(Token::RBracket)?;
                (Expression::ThetaScalar(basis, query), 0)
            }
        } else if name == "rvKummerDiffAdd" {
            self.expect(Token::LParen)?;
            let left = self.expr(nest + 1)?;
            self.expect(Token::Comma)?;
            let right = self.expr(nest + 1)?;
            self.expect(Token::Comma)?;
            let difference = self.expr(nest + 1)?;
            self.expect(Token::RParen)?;
            let depth = 1 + left.depth.max(right.depth).max(difference.depth);
            (
                Expression::DifferentialAdd(Box::new([left, right, difference])),
                depth,
            )
        } else {
            let operation = match name {
                "rvJacKummerPoint" => Operation::Project,
                "rvKummerCoord" | "rvKummerDirection" => {
                    self.expect(Token::LBracket)?;
                    let (coordinate, span) = self.number()?;
                    if coordinate >= 4 {
                        return Err(Error::new(
                            span,
                            format!("{name} index must be between 0 and 3"),
                        ));
                    }
                    self.expect(Token::RBracket)?;
                    if name == "rvKummerCoord" {
                        Operation::Coordinate(coordinate)
                    } else {
                        Operation::Direction(coordinate)
                    }
                }
                "rvKummerScale" => {
                    self.expect(Token::LBracket)?;
                    let n = self.signed_wide("an integer multiple", "Kummer scalar")?;
                    self.expect(Token::RBracket)?;
                    Operation::Scale(n)
                }
                "rvKummerLiftable" => Operation::Liftable,
                "rvKummerSmooth" => Operation::Smooth,
                "rvKummerBlowdown" => Operation::Blowdown,
                "rvKummerIsExceptional" => Operation::IsExceptional,
                "rvKummerSwitch" => Operation::Switch,
                _ => Operation::Lift,
            };
            let input = self.paren_expr(nest)?;
            let depth = 1 + input.depth;
            (Expression::Apply(operation, Box::new(input)), depth)
        };
        expr(span, ExprForm::Kummer(Box::new(expression)), depth)
    }

    fn kummer_theta_index(&mut self, description: &str, maximum: usize) -> Result<usize> {
        self.expect(Token::Comma)?;
        let (index, span) = self.number()?;
        if index > maximum {
            return Err(Error::new(
                span,
                format!("{description} must be between 0 and {maximum}"),
            ));
        }
        Ok(index)
    }

    fn richelot_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::richelot::Command;
        self.advance()?;
        let command = match name {
            "rvRichelotData" | "rvRichelotDescent" => {
                self.expect(Token::LBracket)?;
                let target = self.label()?;
                self.expect(Token::RBracket)?;
                if name == "rvRichelotData" {
                    Command::Data(target)
                } else {
                    Command::Descent(target)
                }
            }
            "rvrichelotdescend" => {
                let target = self.label()?;
                self.expect(Token::LBracket)?;
                let source = self.label()?;
                self.expect(Token::Comma)?;
                let extension_target = self.label()?;
                self.expect(Token::RBracket)?;
                Command::Descend(target, source, extension_target)
            }
            "rvJacRichelotKernels" => {
                self.expect(Token::LBracket)?;
                let (basis, base) = self.richelot_level_names()?;
                self.expect(Token::RBracket)?;
                Command::Kernels(basis, base)
            }
            "rvrichelotfromlevel" => {
                let target = self.label()?;
                self.expect(Token::LBracket)?;
                let (basis, base) = self.richelot_level_names()?;
                self.expect(Token::Comma)?;
                let index = self.number()?.0;
                self.expect(Token::RBracket)?;
                Command::FromLevel(target, basis, base, index)
            }
            _ => {
                let target = self.label()?;
                self.expect(Token::LBracket)?;
                let source = self.label()?;
                self.expect(Token::RBracket)?;
                if name == "rvrichelotdual" {
                    Command::Dual(target, source)
                } else {
                    let coefficient = |parser: &mut Self| {
                        parser.signed_tuple(
                            "a signed residue coordinate",
                            "Richelot coefficient coordinate",
                        )
                    };
                    let factors = [
                        self.set_named_rows(coefficient)?,
                        self.set_named_rows(coefficient)?,
                        self.set_named_rows(coefficient)?,
                    ];
                    Command::Define(target, source, factors)
                }
            }
        };
        Ok(StatementForm::Richelot(Box::new(command)))
    }

    fn richelot_level_names(&mut self) -> Result<(Label, Label)> {
        let basis = self.label()?;
        self.expect(Token::Comma)?;
        let base = self.label()?;
        Ok((basis, base))
    }

    fn elliptic_quotient_statement(&mut self, glue: bool) -> Result<StatementForm> {
        self.advance()?;
        let curve = self.label()?;
        self.expect(Token::LBracket)?;
        let left = self.label()?;
        self.expect(Token::Comma)?;
        let right = self.label()?;
        self.expect(Token::RBracket)?;
        self.expect(Token::LParen)?;
        let a = self.expr(1)?;
        self.expect(Token::Comma)?;
        let b = self.expr(1)?;
        self.expect(Token::Comma)?;
        let c = self.expr(1)?;
        self.expect(Token::Comma)?;
        let d = self.expr(1)?;
        self.expect(Token::RParen)?;
        let matrix = [a, b, c, d];
        let factors = [left, right];
        if glue {
            Ok(StatementForm::EllipticGlue(Box::new(
                crate::elliptic_gluing::Command {
                    name: curve,
                    factors,
                    matrix,
                },
            )))
        } else {
            Ok(StatementForm::Bielliptic(Box::new(
                crate::bielliptic::Command {
                    source: curve,
                    targets: factors,
                    matrix,
                },
            )))
        }
    }

    fn jacobian_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hyperelliptic_jacobian::{Expression, Operation, Pairing};
        if matches!(name, "rvJacWeil2" | "rvJacTate2") {
            let pairing = if name == "rvJacWeil2" {
                Pairing::Weil
            } else {
                Pairing::Tate
            };
            let inputs = self.paren_pair(nest)?;
            let depth = 1 + inputs[0].depth.max(inputs[1].depth);
            let expression = Expression::Pair(pairing, Box::new(inputs));
            return expr(span, ExprForm::Jacobian(Box::new(expression)), depth);
        }
        if name == "rvJacDivisible2" {
            let input = self.paren_expr(nest)?;
            let depth = 1 + input.depth;
            return expr(
                span,
                ExprForm::Jacobian(Box::new(Expression::Apply(
                    Operation::Divisible2,
                    Box::new(input),
                ))),
                depth,
            );
        }
        let coefficient = |parser: &mut Self| {
            let context = if name == "rvJacHalfPoint" {
                "halving witness coefficient coordinate"
            } else {
                "Mumford coefficient coordinate"
            };
            parser.signed_tuple("a signed residue coordinate", context)
        };
        self.expect(Token::LBracket)?;
        let (expression, depth) = match name {
            "rvJacTate" | "rvJacWeil" => {
                let degree_span = self.span();
                let description = format!("{name} order");
                let degree = self.wide_unsigned(&description)?;
                let degree = u64::try_from(degree)
                    .map_err(|_| Error::new(degree_span, format!("{description} exceeds u64")))?;
                self.expect(Token::RBracket)?;
                let inputs = self.paren_pair(nest)?;
                let depth = 1 + inputs[0].depth.max(inputs[1].depth);
                let pairing = if name == "rvJacTate" {
                    Pairing::MillerTate(degree)
                } else {
                    Pairing::MillerWeil(degree)
                };
                (Expression::Pair(pairing, Box::new(inputs)), depth)
            }
            "rvjac" => {
                let curve = self.label()?;
                self.expect(Token::RBracket)?;
                let u = self.set_named_rows(coefficient)?;
                let v = self.set_named_rows(coefficient)?;
                (Expression::Literal(curve, u, v), 0)
            }
            "rvJacPoint" | "rvJacHalfPoint" => {
                let curve = self.label()?;
                self.expect(Token::RBracket)?;
                let xy = self.paren_pair(nest)?;
                let depth = 1 + xy[0].depth.max(xy[1].depth);
                let root = if name == "rvJacHalfPoint" {
                    Some(self.set_named_rows(coefficient)?)
                } else {
                    None
                };
                (Expression::Point(curve, Box::new(xy), root), depth)
            }
            "rvJacHalfAt" => {
                let curve = self.label()?;
                self.expect(Token::Comma)?;
                let Token::Number(index) = self.token() else {
                    return Err(self.expected("a nonnegative decimal integer"));
                };
                let index = index.into();
                self.advance()?;
                self.expect(Token::RBracket)?;
                let xy = self.paren_pair(nest)?;
                let depth = 1 + xy[0].depth.max(xy[1].depth);
                (Expression::HalfAt(curve, index, Box::new(xy)), depth)
            }
            "rvRichelotPoint" | "rvBiellipticPoint" | "rvEllGlueLeft" | "rvEllGlueRight" => {
                let target = self.label()?;
                self.expect(Token::RBracket)?;
                let coordinates = Box::new(self.paren_pair(nest)?);
                let depth = 1 + coordinates[0].depth.max(coordinates[1].depth);
                let expression = match name {
                    "rvRichelotPoint" => Expression::RichelotPoint(target, coordinates),
                    "rvBiellipticPoint" => Expression::BiellipticPoint(target, coordinates),
                    "rvEllGlueLeft" => Expression::EllipticGluePoint(
                        target,
                        crate::elliptic_gluing::Side::Left,
                        coordinates,
                    ),
                    _ => Expression::EllipticGluePoint(
                        target,
                        crate::elliptic_gluing::Side::Right,
                        coordinates,
                    ),
                };
                (expression, depth)
            }
            "rvJacQuadTrace" => {
                let curve = self.label()?;
                self.expect(Token::RBracket)?;
                let radicand = self.paren_expr(nest)?;
                let depth = 1 + radicand.depth;
                let components = [
                    self.set_named_rows(coefficient)?,
                    self.set_named_rows(coefficient)?,
                    self.set_named_rows(coefficient)?,
                    self.set_named_rows(coefficient)?,
                ];
                (
                    Expression::QuadraticTrace(curve, Box::new(radicand), Box::new(components)),
                    depth,
                )
            }
            "rvRichelotPush" | "rvRichelotDescPush" => {
                let target = self.label()?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                let expression = if name == "rvRichelotPush" {
                    Expression::RichelotPush(target, Box::new(input))
                } else {
                    Expression::RichelotDescPush(target, Box::new(input))
                };
                (expression, depth)
            }
            "rvBiellipticPush" | "rvBiellipticPull" => {
                use crate::bielliptic::picard::Direction;
                let target = self.label()?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                let direction = if name == "rvBiellipticPush" {
                    Direction::Push
                } else {
                    Direction::Pull
                };
                (
                    Expression::BiellipticPicard(target, direction, Box::new(input)),
                    depth,
                )
            }
            "rvEllGlueMap" => {
                let target = self.label()?;
                self.expect(Token::RBracket)?;
                let arguments = self.paren_pair(nest)?;
                let depth = 1 + arguments[0].depth.max(arguments[1].depth);
                (
                    Expression::EllipticGlueMap(target, Box::new(arguments)),
                    depth,
                )
            }
            "rvEllGlueDual" => {
                let target = self.label()?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::EllipticGlueDual(target, Box::new(input)), depth)
            }
            "rvJacRechart" | "rvJacUnchart" => {
                use crate::hyperelliptic_chart::Direction;
                let chart = self.label()?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                let direction = if name == "rvJacRechart" {
                    Direction::ToChart
                } else {
                    Direction::ToSource
                };
                (
                    Expression::Transport(chart, direction, Box::new(input)),
                    depth,
                )
            }
            _ => {
                let operation = match name {
                    "rvJacScale" => Operation::Scale(
                        self.signed_wide("an integer multiple", "Jacobian scalar")?,
                    ),
                    "rvJacU" | "rvJacV" => {
                        Operation::Coefficient(name == "rvJacV", self.number()?.0)
                    }
                    "rvJacKummerSupport" => Operation::KummerSupport(self.number()?.0),
                    "rvJacExtend" => Operation::Extend(self.label()?),
                    "rvJacDescend" => Operation::Descend(self.label()?),
                    "rvJacWeilApply" => Operation::Weil(self.label()?),
                    _ => {
                        let base = self.label()?;
                        self.expect(Token::Comma)?;
                        Operation::Frobenius(base, self.number()?.0)
                    }
                };
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::Apply(operation, Box::new(input)), depth)
            }
        };
        expr(span, ExprForm::Jacobian(Box::new(expression)), depth)
    }

    fn hww_statement(&mut self, name: &str) -> Result<StatementForm> {
        use crate::hw_closed::{Command, Query};
        self.advance()?;
        self.expect(Token::LBracket)?;
        let curve = self.label()?;
        let query = match name {
            "rvHWWPoints" => {
                self.expect(Token::Comma)?;
                Query::Points(self.number()?.0)
            }
            "rvHWWWeil" => Query::Weil,
            _ => Query::Verify,
        };
        self.expect(Token::RBracket)?;
        let candidate = if matches!(query, Query::Verify) && self.token() == Token::LParen {
            Some(Box::new(self.paren_expr(0)?))
        } else {
            None
        };
        Ok(StatementForm::HWWClosed(Box::new(Command::Report(
            curve, query, candidate,
        ))))
    }

    fn hww_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::hw_closed::{CurveOp, Expression, Unary};
        let mut depth = 0;
        let expression = match name {
            "rvHWW" | "rvHWWRatio" | "rvHWWUnit" => {
                self.expect(Token::LBracket)?;
                let curve = self.label()?;
                let query = match name {
                    "rvHWW" => {
                        self.expect(Token::Comma)?;
                        CurveOp::Alpha(self.number()?.0)
                    }
                    "rvHWWRatio" => CurveOp::Ratio,
                    _ => CurveOp::Product,
                };
                self.expect(Token::RBracket)?;
                Expression::Curve(curve, query)
            }
            "rvHWWCompose" => {
                let pair = self.paren_pair(nest)?;
                depth = 1 + pair[0].depth.max(pair[1].depth);
                Expression::Compose(Box::new(pair))
            }
            "rvHWWReduce" => {
                self.expect(Token::LBracket)?;
                let ring = self.label()?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                depth = input.depth + 1;
                Expression::Reduce(ring, Box::new(input))
            }
            _ => {
                let operation = match name {
                    "rvHWWInv" => Unary::Inverse,
                    "rvHWWDescend" => Unary::Descend,
                    _ => {
                        self.expect(Token::LBracket)?;
                        let i = self.number()?.0;
                        let operation = match name {
                            "rvHWWEntry" => {
                                self.expect(Token::Comma)?;
                                Unary::Entry(i, self.number()?.0)
                            }
                            "rvHWWChar" => Unary::Characteristic(i),
                            _ => Unary::Frobenius(i),
                        };
                        self.expect(Token::RBracket)?;
                        operation
                    }
                };
                let input = self.paren_expr(nest)?;
                depth = input.depth + 1;
                Expression::Unary(operation, Box::new(input))
            }
        };
        expr(span, ExprForm::HWWClosed(Box::new(expression)), depth)
    }

    fn ec_statement(&mut self, command: &str) -> Result<StatementForm> {
        use crate::elliptic::Command;
        self.advance()?;
        let statement = match command {
            "rvec" => {
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let field = self.label()?;
                self.expect(Token::RBracket)?;
                Command::Curve(name, field, Box::new(self.paren_pair(0)?))
            }
            "rvecisogeny" => {
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let source = self.label()?;
                self.expect(Token::Comma)?;
                let target = self.label()?;
                self.expect(Token::RBracket)?;
                Command::Isogeny(name, source, target, Box::new(self.paren_expr(0)?))
            }
            _ => {
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                let command = match command {
                    "rvECPointCount" => Command::PointCount(name),
                    "rvECCertify" => Command::Certify(name),
                    _ => {
                        self.expect(Token::Comma)?;
                        Command::Kernels(name, self.number()?.0)
                    }
                };
                self.expect(Token::RBracket)?;
                command
            }
        };
        Ok(StatementForm::Elliptic(Box::new(statement)))
    }

    fn ec_cm_field(&mut self) -> Result<Label> {
        self.expect(Token::LBracket)?;
        let field = self.label()?;
        self.expect(Token::RBracket)?;
        Ok(field)
    }

    fn ec_cm_statement(&mut self) -> Result<StatementForm> {
        self.advance()?;
        let family = self.ec_cm_field()?;
        let input = self.paren_expr(0)?;
        Ok(StatementForm::EllipticCM(Box::new(
            crate::elliptic_cm_recognition::Command::Certify(family, Box::new(input)),
        )))
    }

    fn ec_cm_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::elliptic_cm_recognition::Expression;
        let (form, depth) = match name {
            "rvECFrobenius" | "rvECLang" => (
                Expression::Frobenius(self.ec_cm_field()?, name == "rvECLang"),
                0,
            ),
            "rvECSum" => {
                let args = self.paren_pair(nest)?;
                let depth = 1 + args[0].depth.max(args[1].depth);
                (Expression::Sum(Box::new(args)), depth)
            }
            "rvECNegate" => {
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Negate(Box::new(input)), depth)
            }
            "rvECCMRecognize" => {
                let family = self.ec_cm_field()?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Recognize(family, Box::new(input)), depth)
            }
            _ => {
                self.expect(Token::LBracket)?;
                let source = self.label()?;
                self.expect(Token::Comma)?;
                let target = self.label()?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Extend(source, target, Box::new(input)), depth)
            }
        };
        expr(span, ExprForm::EllipticCM(Box::new(form)), depth)
    }

    fn ec_local_test(&mut self) -> Result<(Label, crate::elliptic_local::Test)> {
        self.expect(Token::LBracket)?;
        let curve = self.label()?;
        self.expect(Token::Comma)?;
        let prime = self.number()?.0;
        self.expect(Token::Comma)?;
        let level = self.number()?.0;
        let candidate = if self.eat(Token::Comma)? {
            Some(self.signed_wide("an integer scalar candidate", "Frobenius scalar candidate")?)
        } else {
            None
        };
        self.expect(Token::RBracket)?;
        Ok((
            curve,
            crate::elliptic_local::Test {
                prime,
                level,
                candidate,
            },
        ))
    }

    fn ec_local_statement(&mut self, command: &str) -> Result<StatementForm> {
        use crate::elliptic_local::Command;
        self.advance()?;
        let action = match command {
            "rvEllFrobCert" => {
                let (curve, test) = self.ec_local_test()?;
                Command::Certificate(curve, test)
            }
            "rvellipticcm" => {
                let field = self.label()?;
                self.expect(Token::LBracket)?;
                let curve = self.label()?;
                self.expect(Token::RBracket)?;
                Command::Field(field, curve)
            }
            _ => {
                self.expect(Token::LBracket)?;
                let label = self.label()?;
                let action = if command == "rvEllLocalOrder" {
                    self.expect(Token::Comma)?;
                    Command::Order(label, self.number()?.0)
                } else {
                    Command::Edge(label)
                };
                self.expect(Token::RBracket)?;
                action
            }
        };
        Ok(StatementForm::EllipticLocal(Box::new(action)))
    }

    fn ec_local_expression(&mut self, name: &str, span: Span) -> Result<ParsedExpr> {
        use crate::elliptic_local::Expression;
        let value = if name == "rvEllFrobCoordinates" {
            let (curve, test) = self.ec_local_test()?;
            Expression::Coordinates(curve, test)
        } else {
            self.expect(Token::LBracket)?;
            let field = self.label()?;
            self.expect(Token::Comma)?;
            let curve = self.label()?;
            self.expect(Token::Comma)?;
            let prime = self.number()?.0;
            self.expect(Token::RBracket)?;
            Expression::Generator(field, curve, prime)
        };
        expr(span, ExprForm::EllipticLocal(Box::new(value)), 0)
    }

    fn ec_expression(&mut self, command: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::elliptic::{Expression, MapPolynomial, Unary};
        let mut depth = 0;
        let expression = match command {
            "rvECApply" => {
                let args = self.paren_pair(nest)?;
                depth = 1 + args[0].depth.max(args[1].depth);
                Expression::Apply(Box::new(args))
            }
            "rvECScale" | "rvECX" | "rvECY" | "rvECKCoeff" | "rvECGetKernel" | "rvECXNum"
            | "rvECXDen" | "rvECDNum" | "rvECDDen" => {
                let op = match command {
                    "rvECX" => Unary::PointCoordinate(false),
                    "rvECY" => Unary::PointCoordinate(true),
                    "rvECGetKernel" => Unary::IsogenyKernel,
                    "rvECScale" => {
                        self.expect(Token::LBracket)?;
                        let n =
                            self.signed_wide("a signed point multiple", "elliptic point multiple")?;
                        self.expect(Token::RBracket)?;
                        Unary::Multiple(n)
                    }
                    _ => {
                        self.expect(Token::LBracket)?;
                        let k = self.number()?.0;
                        self.expect(Token::RBracket)?;
                        match command {
                            "rvECKCoeff" => Unary::KernelCoefficient(k),
                            "rvECXNum" => Unary::IsogenyCoefficient(MapPolynomial::XNumerator, k),
                            "rvECXDen" => Unary::IsogenyCoefficient(MapPolynomial::XDenominator, k),
                            "rvECDNum" => {
                                Unary::IsogenyCoefficient(MapPolynomial::DerivativeNumerator, k)
                            }
                            _ => Unary::IsogenyCoefficient(MapPolynomial::DerivativeDenominator, k),
                        }
                    }
                };
                let input = self.paren_expr(nest)?;
                depth = 1 + input.depth;
                Expression::Unary(op, Box::new(input))
            }
            _ => {
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                if command == "rveckernel" {
                    self.expect(Token::Comma)?;
                    let ell = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let coefficients = self.set_named_rows(|parser| {
                        parser.signed_tuple(
                            "a signed residue coordinate",
                            "kernel coefficient coordinate",
                        )
                    })?;
                    Expression::Kernel(name, ell, coefficients)
                } else {
                    self.expect(Token::RBracket)?;
                    match command {
                        "rvECCurve" => Expression::Curve(name),
                        "rvECA" => Expression::CurveCoefficient(name, false),
                        "rvECB" => Expression::CurveCoefficient(name, true),
                        "rvECIso" => Expression::Isogeny(name),
                        "rvECInfinity" => Expression::Infinity(name),
                        _ => {
                            let args = self.paren_pair(nest)?;
                            depth = 1 + args[0].depth.max(args[1].depth);
                            Expression::Point(name, Box::new(args))
                        }
                    }
                }
            }
        };
        expr(span, ExprForm::Elliptic(Box::new(expression)), depth)
    }

    fn ec_map_expression(&mut self, name: &str, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::ec_morphism::{Binary, Coordinate, Expression, Unary};
        let mut depth = 0;
        let expression = match name {
            "rvECMapCompose" | "rvECMapApply" | "rvECPush" => {
                let operation = match name {
                    "rvECMapCompose" => Binary::Compose,
                    "rvECMapApply" => Binary::Apply,
                    _ => Binary::Push,
                };
                let args = self.paren_pair(nest)?;
                depth = 1 + args[0].depth.max(args[1].depth);
                Expression::Binary(operation, Box::new(args))
            }
            "rvECId" | "rvECMultMap" | "rvECChart" => {
                self.expect(Token::LBracket)?;
                let source = self.label()?;
                if name == "rvECId" {
                    self.expect(Token::RBracket)?;
                    Expression::Identity(source)
                } else {
                    self.expect(Token::Comma)?;
                    if name == "rvECMultMap" {
                        let n = self.signed_wide(
                            "a signed integer multiple",
                            "elliptic multiplication index",
                        )?;
                        self.expect(Token::RBracket)?;
                        Expression::Multiple(source, n)
                    } else {
                        let target = self.label()?;
                        self.expect(Token::RBracket)?;
                        let scale = self.paren_expr(nest)?;
                        depth = 1 + scale.depth;
                        Expression::Chart(source, target, Box::new(scale))
                    }
                }
            }
            _ => {
                let operation = match name {
                    "rvECMap" => Unary::FromVelu,
                    "rvECDual" => Unary::Dual,
                    _ => {
                        self.expect(Token::LBracket)?;
                        let index = self.number()?.0;
                        self.expect(Token::RBracket)?;
                        match name {
                            "rvECMapXNum" => Unary::Coefficient(Coordinate::XNumerator, index),
                            "rvECMapXDen" => Unary::Coefficient(Coordinate::XDenominator, index),
                            "rvECMapYNum" => Unary::Coefficient(Coordinate::YNumerator, index),
                            "rvECMapYDen" => Unary::Coefficient(Coordinate::YDenominator, index),
                            _ => Unary::KernelCoefficient(index),
                        }
                    }
                };
                let input = self.paren_expr(nest)?;
                depth = 1 + input.depth;
                Expression::Unary(operation, Box::new(input))
            }
        };
        expr(span, ExprForm::ECMap(Box::new(expression)), depth)
    }

    pub(crate) fn statement(&mut self) -> Result<Option<Statement>> {
        let span = self.span();
        let kind = match self.token() {
            Token::Word(command @ ("shkuranishi" | "shkuranishiEq" | "shkuranishiDims")) => {
                self.sheaf_kuranishi_statement(command)?
            }
            Token::Word(command @ ("leeDt" | "leeDtProfile" | "leeDtWeight")) => {
                self.lee_dt_command(command)?
            }
            Token::Word(
                command @ ("leeRank" | "leeLevels" | "leeS" | "leeFiltration" | "leeIndependent"),
            ) => {
                self.advance()?;
                self.lee_command(command)?
            }
            Token::Word(
                command @ ("khbraid" | "khGroup" | "khEuler" | "khAnnularGroup" | "khAnnularEuler"
                | "khAnnularCheck" | "khActionCheck" | "khIrreps" | "khPageDim"
                | "khKappa"),
            ) => self.kh_statement(command)?,
            Token::Word(
                command @ ("rvhyperelliptic" | "rvHWStatus" | "rvHWPoints" | "rvHWWeil"
                | "rvHWConverges" | "rvHWCertify" | "rvHWVal"),
            ) => self.hw_statement(command)?,
            Token::Word(command @ ("rvHWWWeil" | "rvHWWPoints" | "rvHWWVerify")) => {
                self.hww_statement(command)?
            }
            Token::Word(
                command @ ("rvec" | "rvecisogeny" | "rvECPointCount" | "rvECIsogenies"
                | "rvECCertify"),
            ) => self.ec_statement(command)?,
            Token::Word(
                command @ ("rvEllFrobCert" | "rvEllLocalOrder" | "rvellipticcm" | "rvEllEdge"),
            ) => self.ec_local_statement(command)?,
            Token::Word(name @ ("rvECMapDegree" | "rvECMapCertify")) => {
                use crate::ec_morphism::Report;
                let report = if name == "rvECMapDegree" {
                    Report::Degree
                } else {
                    Report::Certify
                };
                self.advance()?;
                StatementForm::ECMapReport(report, Box::new(self.paren_expr(0)?))
            }
            Token::Word("rvECCMCertify") => self.ec_cm_statement()?,
            Token::Word(
                command @ ("rvjacbase" | "rvJacCount" | "rvJacEnumerate" | "rvJacOrder"
                | "rvJacCompose" | "rvjacchart" | "rvJacKummer" | "rvJacHalfFiber"),
            ) => self.jacobian_statement(command)?,
            Token::Word(command @ ("rvKummerEnriquesPoints" | "rvKummerEnriquesZeta")) => {
                self.kummer_enriques_statement(command)?
            }
            Token::Word("rvenriques") => self.declare_enriques_statement()?,
            Token::Word(
                op @ ("rvk3pencil"
                | "rvK3PencilFiber"
                | "rvK3PencilCensus"
                | "rvK3PencilWeierstrass"
                | "rvK3PencilDiscriminant"
                | "rvK3FiberReduction"
                | "rvK3FiberEuler"
                | "rvK3PencilLFunction"
                | "rvK3PencilBSD"
                | "rvK3PencilBSDFactors"
                | "rvK3PencilTorsion2"
                | "rvK3PencilSelmer"),
            ) => self.k3_pencil_statement(op)?,
            Token::Word(
                op
                @ ("rvk3fiberjac" | "rvK3FiberQuartic" | "rvK3FiberJacobian" | "rvK3FiberPoints"),
            ) => self.k3_fiber_equation_statement(op)?,
            Token::Word(
                op @ ("rvk3mordell"
                | "rvK3MordellCoord"
                | "rvK3MordellHeight"
                | "rvK3MordellLattice"
                | "rvK3MordellBSDFactors"),
            ) => self.k3_mordell_statement(op)?,
            Token::Word(
                command @ ("rvKummerK3Zeta"
                | "rvKummerK3Points"
                | "rvKummerK3Picard"
                | "rvKummerK3ArtinTate"),
            ) => self.kummer_k3_statement(command)?,
            Token::Word(command @ ("rvKummerDivPair" | "rvKummerFujiki")) => {
                self.kummer_divisor_statement(command)?
            }
            Token::Word(
                command @ ("rvKummerHilbPoints"
                | "rvKummerHilbWeilRestrictPoints"
                | "rvKummerHilbTrace"
                | "rvKummerHilbBetti"
                | "rvKummerHilbChar"
                | "rvK3EllipticGenus"
                | "rvGenKummerElliptic"
                | "rvIgusaPhi10"
                | "rvK3SymplecticElliptic"
                | "rvKummerTranslateElliptic"
                | "rvKummerEnriquesElliptic"),
            ) => self.kummer_hilbert_statement(command)?,
            Token::Word(
                command @ ("rvKummerK3Height" | "rvKummerHilbNewton" | "rvJacHilbNewton"
                | "rvGenKummerNewton"),
            ) => self.crystalline_newton_statement(command)?,
            Token::Word(
                command @ ("rvK3HilbHodge"
                | "rvJacHilbHodge"
                | "rvGenKummerHodge"
                | "rvGenKummerTranslateHodge"
                | "rvGenKummerInversionHodge"
                | "rvGenKummerHodgeSpectrum"),
            ) => self.hodge_statement(command)?,
            Token::Word(
                command @ ("rvHilbLLVMultiplicity"
                | "rvHilbLLVDimension"
                | "rvHilbLLVDecomposition"
                | "rvHilbLLVHodge"
                | "rvHilbLLVHodgeDecomposition"
                | "rvHilbLLVHodgeTable"
                | "rvHilbReducedLLVMultiplicity"
                | "rvHilbReducedLLVDecomposition"
                | "rvHilbPrimitiveLLVMultiplicity"
                | "rvHilbPrimitiveLLVDecomposition"),
            ) => self.hilbert_llv_statement(command)?,
            Token::Word("rvKummerNodeLattice") => self.kummer_node_lattice_statement()?,
            Token::Word("rvAbMukaiPair") => {
                self.advance()?;
                StatementForm::AbMukaiPair(Box::new(self.paren_pair(0)?))
            }
            Token::Word(
                operation @ ("rvog6nilpotentii"
                | "rvog6nilpotentiii"
                | "rvOG6NilpotentData"
                | "rvOG6NilpotentInWeight"),
            ) => self.og6_nilpotent_statement(operation)?,
            Token::Word(operation @ ("rvOG6SpinData" | "rvOG6CliffordCoeff")) => {
                self.og6_spin_statement(operation, span)?
            }
            Token::Word(
                operation @ ("rvKSExteriorData" | "rvKSExteriorSmith" | "rvKSExteriorCoeff"),
            ) => self.ks_exterior_statement(operation)?,
            Token::Word(
                operation @ ("rvKSVolumeSquare" | "rvKSEndCoeff" | "rvKSEndTrace"
                | "rvKSEndTracePair"),
            ) => self.ks_end_statement(operation)?,
            Token::Word(operation @ ("rvog6cusp" | "rvOG6CuspGram")) => {
                self.og6_cusp_statement(operation)?
            }
            Token::Word(
                operation @ ("rvog6" | "rvOG6Lattice" | "rvOG6Pair" | "rvOG6Div" | "rvOG6Fujiki"
                | "rvOG6RR" | "rvOG6Orbit" | "rvOG6SameOrbit" | "rvOG6Wall"),
            ) => self.og6_statement(operation)?,
            Token::Word(
                operation @ ("rvOG10Lattice" | "rvOG10Pair" | "rvOG10Div" | "rvOG10Fujiki"
                | "rvOG10RR" | "rvOG10Wall" | "rvOG10Walls"),
            ) => self.og10_statement(operation)?,
            Token::Word("rvspechtstratum") => self.specht_stratum_statement()?,
            Token::Word(
                operation @ ("rvhilbstratum"
                | "rvHilbStratumData"
                | "rvK3StratumHodge"
                | "rvJacStratumHodge"
                | "rvGenKummerStratumHodge"),
            ) => self.hilbert_stratum_statement(operation)?,
            Token::Word(
                operation @ ("rvmukaimoduli"
                | "rvMukaiModuliLattice"
                | "rvMukaiModuliDiv"
                | "rvMukaiModuliFujiki"
                | "rvMukaiModuliRR"
                | "rvMukaiQuotientPair"),
            ) => self.mukai_moduli_statement(operation)?,
            Token::Word(name @ ("rvMukaiPair" | "rvMukaiHilbFujiki")) => {
                self.mukai_statement(name)?
            }
            Token::Word("rvHilbSquareIntegral") => {
                self.advance()?;
                StatementForm::HilbertSquareIntegral(Box::new(self.paren_expr(0)?))
            }
            Token::Word("rvKummerWeilCoeff") => {
                self.advance()?;
                StatementForm::KummerWeilCoefficient(
                    crate::hyperelliptic_kummer_surface::weil::Coefficient(Box::new(
                        self.paren_pair(0)?,
                    )),
                )
            }
            Token::Word("rvHilbSquareNakajimaDen") => {
                self.advance()?;
                StatementForm::HilbertSquareNakajimaDen(Box::new(self.paren_expr(0)?))
            }
            Token::Word(name @ ("rvK3LatticePair" | "rvK3LatticeGram")) => {
                self.k3_lattice_statement(name)?
            }
            Token::Word("rvCliffordCoeff") => self.clifford_coefficient_statement()?,
            Token::Word("rvK3SpinCoeff") => self.spinor_coefficient_statement()?,
            Token::Word(name @ ("rveccmhorizontal" | "rvECCMHorizontalData")) => {
                self.ec_horizontal_statement(name)?
            }
            Token::Word(
                name
                @ ("rveccmcycle" | "rvECCMCycleData" | "rvECCMCycleStep" | "rvECCMCycleVerify"),
            ) => self.ec_cm_cycle_statement(name)?,
            Token::Word(
                name @ ("rvCMArtinPrime" | "rvCMArtinOrbit" | "rvCMArtinFactors" | "rvCMArtinLift"),
            ) => self.cm_artin_statement(name)?,
            Token::Word(
                name @ ("rvecendorder" | "rvECEndData" | "rvECEndCheck" | "rvECEndTrace"
                | "rvECEndNorm"),
            ) => self.ec_end_order_statement(name)?,
            Token::Word(name @ ("rvModularPolyDegree" | "rvModularPolyCoeff")) => {
                self.elliptic_modular_statement(name)?
            }
            Token::Word(
                name @ ("rvcmclasspoly" | "rvcmclasspolylift" | "rvCMClassPolyData"
                | "rvCMClassPolyCoeff"),
            ) => self.cm_class_poly_statement(name)?,
            Token::Word(name @ ("rvCMNeighborData" | "rvCMNeighbors" | "rvCMNeighborCounts")) => {
                self.cm_neighbor_statement(name)?
            }
            Token::Word(
                name @ ("rvcmorder" | "rvCMClassOrder" | "rvCMClassGroup" | "rvCMSameGenus"),
            ) => self.cm_class_statement(name)?,
            Token::Word(name @ ("rvNakajimaPrimaryDim" | "rvNakajimaIsPrimary")) => {
                self.nakajima_primary_statement(name)?
            }
            Token::Word(
                command @ ("rvk3period" | "rvK3PeriodReport" | "rvK3PeriodIs11" | "rvK3CMForm"
                | "rvK3CMPeriods" | "rvK3CMEquivalent" | "rvK3CMClassNumber"),
            ) => self.k3_period_statement(command)?,
            Token::Word(name @ ("rvHilbTautChernDen" | "rvHilbTautSegre")) => {
                self.hilbert_tautological_statement(name)?
            }
            Token::Word("rvHilbCohBeauville") => {
                use crate::hyperelliptic_kummer_surface::nakajima::cohomology::Report;
                self.advance()?;
                StatementForm::HilbertCohomology(Box::new(Report::Beauville(Box::new(
                    self.paren_pair(0)?,
                ))))
            }
            Token::Word("rvHilbExteriorEuler") => {
                use crate::hyperelliptic_kummer_surface::nakajima::cohomology::Report;
                self.advance()?;
                self.expect(Token::LBracket)?;
                let n = self.wide_unsigned("Hilbert length")?;
                self.expect(Token::Comma)?;
                let k = self.wide_unsigned("exterior power")?;
                self.expect(Token::RBracket)?;
                StatementForm::HilbertCohomology(Box::new(Report::ExteriorEuler(
                    n,
                    k,
                    Box::new(self.paren_pair(0)?),
                )))
            }
            Token::Word(name @ ("rvHilbCohDen" | "rvHilbCohIntegral")) => {
                use crate::hyperelliptic_kummer_surface::nakajima::cohomology::Report;
                self.advance()?;
                let input = self.paren_expr(0)?;
                let report = if name == "rvHilbCohIntegral" {
                    Report::Integral(Box::new(input))
                } else {
                    Report::Denominator(Box::new(input))
                };
                StatementForm::HilbertCohomology(Box::new(report))
            }
            Token::Word("rvNakajimaPair") => {
                self.advance()?;
                StatementForm::NakajimaPair(Box::new(self.paren_pair(0)?))
            }
            Token::Word("rvSiegelPair") => {
                self.advance()?;
                StatementForm::SiegelPair(Box::new(self.paren_pair(0)?))
            }
            Token::Word("rvSiegelRoot") => {
                self.advance()?;
                StatementForm::SiegelRoot(Box::new(self.paren_expr(0)?))
            }
            Token::Word("rvSiegelDelta5At") => {
                self.advance()?;
                StatementForm::SiegelDelta5At(Box::new(self.paren_expr(0)?))
            }
            Token::Word(
                name @ ("rvsiegelmaass"
                | "rvSiegelMaassJacobi"
                | "rvSiegelMaassCoeff"
                | "rvIgusaChi12"),
            ) => self.siegel_maass_statement(name)?,
            Token::Word(
                name @ ("rvIgusaInverse" | "rvIgusaWallJump" | "rvIgusaResidue" | "rvIgusaPolar"
                | "rvIgusaFinite" | "rvIgusaMockWeight2"),
            ) => self.igusa_inverse_statement(name)?,
            Token::Word(
                name @ ("rvsiegelchamber"
                | "rvSiegelChamberData"
                | "rvIgusaChamberCoeff"
                | "rvIgusaChamberPath"),
            ) => self.siegel_chamber_statement(name)?,
            Token::Word(
                name @ ("rvheckepacket" | "rvHeckeEigenData" | "rvHeckeEigenTrace"
                | "rvHeckeEigenNorm" | "rvHeckeEigenCoord"),
            ) => self.hecke_eigen_statement(name)?,
            Token::Word(name @ ("rvheckeresidual" | "rvHeckeEisensteinSturm")) => {
                self.hecke_residual_statement(name)?
            }
            Token::Word(
                name @ ("rvHeckeQuadraticOrder"
                | "rvHeckeQuadraticCoeff"
                | "rvHeckeQuadraticLocal"
                | "rvHeckeQuadraticIntegral"),
            ) => self.hecke_quadratic_statement(name)?,
            Token::Word(name @ ("rvPeriodDim" | "rvPeriodCoeff" | "rvPeriodHeckeMatrix")) => {
                self.period_statement(name)?
            }
            Token::Word(name @ ("rvHeckeActionTrace" | "rvHeckeChar" | "rvHeckeCharCoeff")) => {
                self.hecke_operator_statement(name)?
            }
            Token::Word(name @ ("rvHurwitzClass" | "rvHurwitzJacobi" | "rvHeckeTrace")) => {
                self.hurwitz_statement(name)?
            }
            Token::Word("rvSiegelDelta5") => {
                use crate::hyperelliptic_kummer_surface::hilbert::{
                    elliptic_genus::Query, Command,
                };
                self.advance()?;
                self.expect(Token::LBracket)?;
                let weight = self.siegel_weight_coordinates()?;
                StatementForm::KummerHilbert(Box::new(Command::EllipticGenus(Query::Delta5(
                    weight,
                ))))
            }
            Token::Word(
                name @ ("rvk3autom"
                | "rvK3Entropy"
                | "rvK3AutoPeriod"
                | "rvK3AutoCycles"
                | "rvK3FrobFixed"
                | "rvK3FrobZeta"
                | "rvK3FrobHilbFixed"
                | "rvK3FrobHilbOrbits"
                | "rvK3StratumFrobTrace"
                | "rvK3StratumEulerFactor"
                | "rvK3SpechtFrobTrace"
                | "rvK3SpechtEulerFactor"
                | "rvK3FrobHilbTrace"
                | "rvK3AutoNikulin"),
            ) => self.k3_automorphism_statement(name)?,
            Token::Word(
                command @ ("rvJacHilbPoints"
                | "rvGenKummerPoints"
                | "rvGenKummerFiber"
                | "rvGenKummerLevel"
                | "rvGenKummerBetti"
                | "rvGenKummerTrace"
                | "rvGenKummerTranslateChar"
                | "rvGenKummerTranslateFixed"
                | "rvGenKummerTranslateFrob"
                | "rvGenKummerTranslatePoints"
                | "rvGenKummerInversionEuler"
                | "rvGenKummerInversionFrob"
                | "rvGenKummerInversionPoints"
                | "rvGenKummerStratumFrob"
                | "rvJacHilbBetti"
                | "rvJacHilbTrace"),
            ) => self.generalized_kummer_statement(command)?,
            Token::Word(
                command @ ("rvjactorsionbasis" | "rvJacTorsionFrob" | "rvJacRichelotMatrix"),
            ) => self.jacobian_level_statement(command)?,
            Token::Word(
                command @ ("rvrichelot"
                | "rvrichelotdual"
                | "rvrichelotdescend"
                | "rvRichelotData"
                | "rvRichelotDescent"
                | "rvJacRichelotKernels"
                | "rvrichelotfromlevel"),
            ) => self.richelot_statement(command)?,
            Token::Word(name @ ("rvbielliptic" | "rvellglue")) => {
                self.elliptic_quotient_statement(name == "rvellglue")?
            }
            Token::Eof => return Ok(None),
            Token::Word(
                command @ ("settopology"
                | "settopologyMap"
                | "settopologymeet"
                | "settopologyjoin"
                | "settopologysupport"
                | "settopologyjoinsteps"
                | "setreflect"
                | "setreflectcheck"),
            ) => StatementForm::SetTopology(self.set_topology_command(command)?),
            Token::Word("setomega") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let site = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SetOmega(name, site)
            }
            Token::Word("setsub") => {
                self.advance()?;
                let name = self.label()?;
                StatementForm::SetSubobject(name, self.paren_expr(0)?)
            }
            Token::Word("setj") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let family = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SetJ(name, family)
            }
            Token::Word("setjcheck") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let family = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SetJCheck(family)
            }
            Token::Word("setsieve" | "setdense") => {
                let dense = self.token() == Token::Word("setdense");
                self.advance()?;
                self.expect(Token::LBracket)?;
                let omega = self.label()?;
                self.expect(Token::Comma)?;
                let point = self.label()?;
                self.expect(Token::Comma)?;
                let index = self.number()?.0;
                self.expect(Token::RBracket)?;
                StatementForm::SetSieve(omega, point, index, dense)
            }
            Token::Word("setrecollement") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LParen)?;
                let closed = self.expr(1)?;
                self.expect(Token::Comma)?;
                let open = self.expr(1)?;
                self.expect(Token::RParen)?;
                StatementForm::SetRecollement(name, Box::new([closed, open]))
            }
            Token::Word(command @ ("setrestrict" | "setgluelimit" | "setglue")) => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let recollement = self.label()?;
                if command == "setglue" {
                    self.expect(Token::RBracket)?;
                    StatementForm::SetGlue(name, recollement, self.paren_expr(0)?)
                } else {
                    self.expect(Token::Comma)?;
                    let family = self.label()?;
                    if command == "setgluelimit" {
                        self.expect(Token::RBracket)?;
                        StatementForm::SetGlueLimit(name, recollement, family)
                    } else {
                        self.expect(Token::Comma)?;
                        let side = self.label()?;
                        let side = match side.text.as_str() {
                            "closed" => crate::set_gluing::Side::Closed,
                            "open" => crate::set_gluing::Side::Open,
                            _ => return Err(Error::new(side.span, "expected 'closed' or 'open'")),
                        };
                        self.expect(Token::RBracket)?;
                        StatementForm::SetRestrict(name, recollement, family, side)
                    }
                }
            }
            Token::Word(op @ ("setpullback" | "setlan" | "setran")) => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let family = self.label()?;
                self.expect(Token::RBracket)?;
                let map = self.paren_expr(0)?;
                if op == "setpullback" {
                    StatementForm::SetPullback(name, family, map)
                } else {
                    StatementForm::SetKan(name, family, map, op == "setran")
                }
            }
            Token::Word("setsum" | "setproduct") => {
                let product = self.token() == Token::Word("setproduct");
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let total = self.label()?;
                self.expect(Token::Comma)?;
                let body = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SetDependent(name, total, body, product)
            }
            Token::Word("setfamily") => {
                self.advance()?;
                let family = self.label()?;
                self.expect(Token::LBracket)?;
                let site = self.label()?;
                self.expect(Token::RBracket)?;
                let fibers = self.set_named_rows(|parser| {
                    let point = parser.label()?;
                    parser.expect(Token::Colon)?;
                    Ok((point, parser.number()?.0))
                })?;
                let transitions = self.set_named_rows(|parser| {
                    let source = parser.label()?;
                    parser.expect(Token::Less)?;
                    let target = parser.label()?;
                    parser.expect(Token::Colon)?;
                    Ok((source, target, parser.exponents()?))
                })?;
                StatementForm::SetFamily(family, site, fibers, transitions)
            }
            Token::Word("elements") => {
                self.advance()?;
                let space = self.label()?;
                self.expect(Token::LBracket)?;
                let family = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SetElements(space, family)
            }
            Token::Word("setnat") => {
                self.advance()?;
                let map = self.label()?;
                self.expect(Token::LBracket)?;
                let source = self.label()?;
                self.expect(Token::Comma)?;
                let target = self.label()?;
                self.expect(Token::RBracket)?;
                let rows = self.set_named_rows(|parser| {
                    let point = parser.label()?;
                    parser.expect(Token::Colon)?;
                    Ok((point, parser.exponents()?))
                })?;
                StatementForm::SetNatural(map, source, target, rows)
            }
            Token::Word("setsections") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let family = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SetSections(family)
            }
            Token::Word("setmonodromy") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let source = self.label()?;
                self.expect(Token::Comma)?;
                let target = self.label()?;
                self.expect(Token::Comma)?;
                let point = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SetMonodromy(name, source, target, point)
            }
            Token::Word(
                op @ ("setmonosize" | "setmonochars" | "setmonoat" | "setmonostabilizer"
                | "setmonomaps" | "setmonoorbitals"),
            ) => {
                use crate::set_monodromy::Query;
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                let query = match op {
                    "setmonosize" => Query::Size,
                    "setmonochars" => Query::Characters,
                    "setmonoat" => {
                        self.expect(Token::Comma)?;
                        Query::Element(self.number()?.0)
                    }
                    "setmonostabilizer" => {
                        self.expect(Token::Comma)?;
                        Query::Stabilizer(self.number()?.0)
                    }
                    "setmonomaps" => Query::Maps,
                    _ => Query::Orbitals,
                };
                self.expect(Token::RBracket)?;
                StatementForm::SetMonoQuery(name, query)
            }
            Token::Word("shlinearize") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let family = self.label()?;
                self.expect(Token::Comma)?;
                let prime = self.number()?.0;
                self.expect(Token::RBracket)?;
                StatementForm::SetLinearize(name, family, prime)
            }
            Token::Word(op @ ("shorbital" | "shorbitals")) => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let monodromy = self.label()?;
                self.expect(Token::Comma)?;
                let source = self.label()?;
                self.expect(Token::Comma)?;
                let target = self.label()?;
                let weights = if op == "shorbital" {
                    self.expect(Token::Comma)?;
                    let basis = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    vec![(basis, 1)]
                } else {
                    self.expect(Token::RBracket)?;
                    self.set_named_rows(|parser| {
                        let basis = parser.number()?.0;
                        parser.expect(Token::Colon)?;
                        Ok((basis, parser.signed_number("orbital coefficient")?))
                    })?
                };
                StatementForm::SetOrbitalMap(name, monodromy, source, target, weights)
            }
            Token::Word("shmaprank") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::Comma)?;
                let point = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SheafMapRank(name, point)
            }
            Token::Word(
                command @ ("legendreTrace" | "legendreHensel" | "legendreDwork" | "dworkPeriod"
                | "dworkCongruence"),
            ) => {
                self.advance()?;
                StatementForm::Legendre(self.legendre_query(command)?)
            }
            Token::Word("rvunram") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let (p, s) = self.bidegree()?;
                self.expect(Token::RBracket)?;
                self.expect(Token::LBrace)?;
                let h = self.signed_list(
                    Token::RBrace,
                    "a signed polynomial coefficient",
                    "unramified modulus",
                )?;
                StatementForm::UnramRing(name, p, s, h)
            }
            Token::Word("rvUnramVal") => {
                self.advance()?;
                StatementForm::UnramValuation(self.paren_expr(0)?)
            }
            Token::Word("rvLegendreData" | "rvLegendreTrace" | "rvLegendreHensel") => {
                use crate::legendre_closed::Query;
                let query = match self.token() {
                    Token::Word("rvLegendreData") => Query::Data,
                    Token::Word("rvLegendreTrace") => Query::Trace,
                    _ => Query::Hensel,
                };
                self.advance()?;
                StatementForm::LegendreClosedReport(query, self.paren_expr(0)?)
            }
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
            Token::Word("mixHdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::Comma)?;
                let b = self.signed_degree()?;
                self.expect(Token::RBracket)?;
                StatementForm::MixedHDimension(name, b)
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
            Token::Word("shserre" | "shserrecomplex") => {
                let form = if self.token() == Token::Word("shserre") {
                    SheafInputForm::Sheaf
                } else {
                    SheafInputForm::Complex
                };
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let sheaf = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SheafSerre(form, name, sheaf)
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
            Token::Word(
                command @ ("shdisk" | "shdiskshriek" | "shdiskstar" | "shdiskmiddle" | "shdiskdual"
                | "shdiskmap" | "shdiskmat" | "shdiskmapmat" | "shdiskranks"
                | "shdiskhomdim" | "shdiskstalk" | "shdiskcostalk" | "shdiskoperator"
                | "shdiskrealize"),
            ) => StatementForm::Disk(self.disk_command(command)?),
            Token::Word(
                command @ ("shdiskseq" | "shhurwitz" | "shdiskpoint" | "shdiskprod"
                | "shdiskseqsame" | "shdiskmonosame" | "shdiskseqmap"
                | "shdiskseqmapbraid" | "shdiskseqmapat" | "shparrank" | "shparrows"
                | "shparvector"),
            ) => StatementForm::Parabolic(self.parabolic_command(command)?),
            Token::Word(
                command @ ("shconv" | "shconvmap" | "shconvrank" | "shdiskinfinity"
                | "shdiskrigidity" | "shdiskgroup" | "shconvcoords"),
            ) => StatementForm::Convolution(self.convolution_command(command)?),
            Token::Word(op @ ("shconeglue" | "shconegluecomplex")) => {
                let form = if op == "shconeglue" {
                    SheafInputForm::Sheaf
                } else {
                    SheafInputForm::Complex
                };
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let map = self.label()?;
                self.expect(Token::Comma)?;
                let apex = self.label()?;
                self.expect(Token::Comma)?;
                let source = self.label()?;
                self.expect(Token::Comma)?;
                let degree = self.number()?.0;
                self.expect(Token::RBracket)?;
                StatementForm::SheafConeGlue(
                    form,
                    [name, map, apex, source],
                    degree,
                    self.sheaf_cone_columns()?,
                )
            }
            Token::Word("shconehomotopy") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let from = self.label()?;
                self.expect(Token::Comma)?;
                let to = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SheafConeHomotopy([name, from, to], self.sheaf_cone_columns()?)
            }
            Token::Word("shbcsquare") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let f = self.label()?;
                self.expect(Token::Comma)?;
                let u = self.label()?;
                self.expect(Token::Comma)?;
                let v = self.label()?;
                self.expect(Token::Comma)?;
                let fp = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::SheafBCSquare(name, [f, u, v, fp])
            }
            Token::Word(op @ ("shbc" | "shbccomplex")) => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let square = self.label()?;
                self.expect(Token::Comma)?;
                let input = self.label()?;
                self.expect(Token::RBracket)?;
                let form = if op == "shbc" {
                    SheafInputForm::Sheaf
                } else {
                    SheafInputForm::Complex
                };
                StatementForm::SheafBCBuild(form, name, square, input)
            }
            Token::Word(op @ ("shbcmap" | "shbchom")) => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let comparison = self.label()?;
                self.expect(Token::Comma)?;
                let degree = self.number()?.0;
                self.expect(Token::RBracket)?;
                StatementForm::SheafBCMap(name, comparison, degree, op == "shbchom")
            }
            Token::Word(op @ ("shbcrank" | "shbckernel" | "shbccokernel")) => {
                use crate::sheaf_base_change::LinearQuery;
                self.advance()?;
                self.expect(Token::LBracket)?;
                let comparison = self.label()?;
                self.expect(Token::Comma)?;
                let y = self.label()?;
                self.expect(Token::Comma)?;
                let degree = self.number()?.0;
                self.expect(Token::RBracket)?;
                let query = match op {
                    "shbckernel" => LinearQuery::Kernel,
                    "shbccokernel" => LinearQuery::Cokernel,
                    _ => LinearQuery::Rank,
                };
                StatementForm::SheafBCLinear(comparison, y, degree, query)
            }
            Token::Word(op @ ("shbcuniversal" | "shbccomma" | "shbcwitness")) => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let square = self.label()?;
                self.expect(Token::Comma)?;
                let prime = self.number()?.0;
                let statement = if op == "shbcuniversal" {
                    StatementForm::SheafBCUniversal(square, prime)
                } else {
                    self.expect(Token::Comma)?;
                    let y = self.label()?;
                    self.expect(Token::Comma)?;
                    let x = self.label()?;
                    let witness = if op == "shbcwitness" {
                        self.expect(Token::Comma)?;
                        let degree = self.signed_number("augmented comma degree")?;
                        self.expect(Token::Comma)?;
                        Some((degree, self.number()?.0))
                    } else {
                        None
                    };
                    StatementForm::SheafBCComma(square, prime, y, x, witness)
                };
                self.expect(Token::RBracket)?;
                statement
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
            Token::Word("rvTLLinking") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (strands, _) = self.number()?;
                self.expect(Token::Comma)?;
                let word = self.signed_tuple("a signed Artin generator", "Artin generator")?;
                self.expect(Token::Comma)?;
                let framings = self.signed_tuple("a signed framing", "framing")?;
                self.expect(Token::RBracket)?;
                StatementForm::TlLinking(strands, word, framings)
            }
            Token::Word("rvCycloCoeff") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (power, _) = self.number()?;
                self.expect(Token::RBracket)?;
                StatementForm::CycloCoefficient(power, self.paren_expr(0)?)
            }
            Token::Word(
                "rvConformalDim" | "rvConformalSeries" | "rvConformalNum" | "rvConformalDen"
                | "rvConformalMod" | "rvConformalExact",
            ) => {
                use crate::rv_conformal::Report;
                let command = self.token();
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (rank, _) = self.number()?;
                self.expect(Token::Comma)?;
                let (level, _) = self.number()?;
                self.expect(Token::Comma)?;
                let form = match command {
                    Token::Word("rvConformalDim") => {
                        let (genus, _) = self.number()?;
                        self.expect(Token::Comma)?;
                        StatementForm::ConformalDimension(
                            rank,
                            level,
                            genus,
                            self.conformal_punctures()?,
                        )
                    }
                    Token::Word("rvConformalMod") => {
                        let genus = self.conformal_size("genus")?;
                        self.expect(Token::Comma)?;
                        let modulus = self.conformal_size("modulus")?;
                        self.expect(Token::Comma)?;
                        StatementForm::ConformalModulo(
                            rank,
                            level,
                            genus,
                            modulus,
                            self.conformal_punctures()?,
                        )
                    }
                    Token::Word("rvConformalExact") => {
                        let genus = self.conformal_size("genus")?;
                        self.expect(Token::Comma)?;
                        StatementForm::ConformalExact(
                            rank,
                            level,
                            genus,
                            self.conformal_punctures()?,
                        )
                    }
                    _ => {
                        let report = match command {
                            Token::Word("rvConformalNum") => Report::Numerator,
                            Token::Word("rvConformalDen") => Report::Denominator,
                            _ => Report::Series,
                        };
                        StatementForm::ConformalSeries(
                            report,
                            rank,
                            level,
                            self.conformal_punctures()?,
                        )
                    }
                };
                self.expect(Token::RBracket)?;
                form
            }
            Token::Word("rvFusionCoeff") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let weight = self.signed_tuple("a signed fusion weight", "fusion weight")?;
                self.expect(Token::RBracket)?;
                StatementForm::FusionCoefficient(weight, self.paren_expr(0)?)
            }
            Token::Word("rvQGrRegularized" | "rvQGrRegularizedRhs") => {
                let rhs = self.token() == Token::Word("rvQGrRegularizedRhs");
                self.advance()?;
                self.expect(Token::LBracket)?;
                let (rank, ambient) = self.quantum_grade()?;
                let statement = if rhs {
                    StatementForm::QuantumRegularizedRhs(rank, ambient)
                } else {
                    self.expect(Token::Comma)?;
                    let degree = self.number()?.0;
                    StatementForm::QuantumRegularized(rank, ambient, degree)
                };
                self.expect(Token::RBracket)?;
                statement
            }
            Token::Word("rvWZCheck" | "rvCTWZCheck") => {
                let with_period = self.token() == Token::Word("rvCTWZCheck");
                self.advance()?;
                self.expect(Token::LBracket)?;
                let period = if with_period {
                    let name = self.label()?;
                    self.expect(Token::Comma)?;
                    Some(name)
                } else {
                    None
                };
                let (u, v) = self.bidegree()?;
                self.expect(Token::Comma)?;
                let z = self.signed_wide("a nonzero integer", "WZ summand parameter")?;
                self.expect(Token::Comma)?;
                let start = self.number()?.0;
                let mirror = if let Some(period) = period {
                    let columns = if self.eat(Token::Comma)? {
                        self.expect(Token::LBracket)?;
                        let mut columns = Vec::new();
                        while self.token() != Token::RBracket {
                            columns.push(
                                self.signed_tuple("a signed lattice exponent", "lattice exponent")?,
                            );
                            if !self.eat(Token::Comma)? {
                                break;
                            }
                        }
                        self.expect(Token::RBracket)?;
                        Some(columns)
                    } else {
                        None
                    };
                    Some(WzMirror { period, columns })
                } else {
                    None
                };
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(0)?;
                self.expect(Token::LBrace)?;
                let numerator = self.wz_numerator()?;
                self.expect(Token::Comma)?;
                let denominator = self.wz_denominator()?;
                self.expect(Token::RBrace)?;
                StatementForm::WzCheck(
                    mirror,
                    crate::rv_wz::Certificate {
                        u,
                        v,
                        z,
                        start,
                        numerator,
                        denominator,
                    },
                    input,
                )
            }
            Token::Word(
                name @ ("rvWittSpectrumData"
                | "rvWittSpectrumGhost"
                | "rvWittSpectrumCoeff"
                | "rvWittSpectrumNecklace"
                | "rvWittSpectrumNormGhost"
                | "rvWittSpectrumCyclicGhost"
                | "rvWittSpectrumLieGhost"),
            ) => {
                use crate::rv_witt_spectrum::Query;
                self.advance()?;
                let query = if name == "rvWittSpectrumData" {
                    Query::Data
                } else if matches!(name, "rvWittSpectrumCyclicGhost" | "rvWittSpectrumLieGhost") {
                    self.expect(Token::LBracket)?;
                    let projection = self.witt_cyclic_power(name == "rvWittSpectrumCyclicGhost")?;
                    self.expect(Token::Comma)?;
                    let index = self.wide_unsigned("Witt spectrum ghost index")?;
                    self.expect(Token::RBracket)?;
                    Query::CyclicGhost(projection, index)
                } else if name == "rvWittSpectrumNormGhost" {
                    self.expect(Token::LBracket)?;
                    let degree = self.wide_unsigned("tensor-induction degree")?;
                    self.expect(Token::Comma)?;
                    let index = self.wide_unsigned("Witt spectrum ghost index")?;
                    self.expect(Token::RBracket)?;
                    Query::NormGhost(degree, index)
                } else {
                    self.expect(Token::LBracket)?;
                    let index = self.wide_unsigned("Witt spectrum query index")?;
                    self.expect(Token::RBracket)?;
                    match name {
                        "rvWittSpectrumGhost" => Query::Ghost(index),
                        "rvWittSpectrumNecklace" => Query::Necklace(index),
                        _ => Query::Factor(index),
                    }
                };
                StatementForm::WittSpectrumQuery(query, Box::new(self.paren_expr(0)?))
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
            Token::Word("rhodepth") => {
                self.advance()?;
                StatementForm::RhoDepth(Box::new(self.paren_expr(0)?))
            }
            Token::Word("rhohidden") => {
                self.advance()?;
                StatementForm::RhoHidden(self.real_extension_pair(0, span, "rhohidden")?)
            }
            Token::Word("rmrhoextdim") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let grade = self.real_rho_grade()?;
                self.expect(Token::RBracket)?;
                StatementForm::RealRhoDimension(grade)
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
            Token::Word("rvastwo") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let prime = self.number()?.0;
                self.expect(Token::RBracket)?;
                self.expect(Token::LBrace)?;
                let mut terms = Vec::new();
                while self.token() != Token::RBrace {
                    self.expect(Token::LParen)?;
                    let exponent =
                        self.signed_wide("a signed Laurent exponent", "two-pole Laurent exponent")?;
                    self.expect(Token::Comma)?;
                    let coefficient = self.signed_wide(
                        "a signed coefficient",
                        "two-pole Artin–Schreier coefficient",
                    )?;
                    self.expect(Token::RParen)?;
                    terms.push((exponent, coefficient));
                    if !self.eat(Token::Comma)? {
                        break;
                    }
                }
                self.expect(Token::RBrace)?;
                StatementForm::RvAsTwo(name, prime, terms)
            }
            Token::Word(
                "rvASTwoBasis"
                | "rvASTwoCartierRank"
                | "rvASTwoStableRank"
                | "rvASTwoaNumber"
                | "rvASTwoOrdinaryRanks"
                | "rvASTwoHasseWitt"
                | "rvASTwoReduced"
                | "rvASTwoGenus"
                | "rvASTwoPrank"
                | "rvASTwoCharacterL"
                | "rvASTwoCharacterNewton"
                | "rvASTwoLPolynomial"
                | "rvASTwoNewton"
                | "rvASTwoPoints",
            ) => {
                let operation = self.token();
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                let query = match operation {
                    Token::Word("rvASTwoBasis") => AsTwoQuery::Cartier(AsTwoCartierQuery::Basis),
                    Token::Word("rvASTwoCartierRank") => {
                        self.expect(Token::Comma)?;
                        AsTwoQuery::Cartier(AsTwoCartierQuery::Rank(self.number()?.0))
                    }
                    Token::Word("rvASTwoStableRank") => {
                        AsTwoQuery::Cartier(AsTwoCartierQuery::StableRank)
                    }
                    Token::Word("rvASTwoaNumber") => {
                        AsTwoQuery::Cartier(AsTwoCartierQuery::ANumber)
                    }
                    Token::Word("rvASTwoOrdinaryRanks") => {
                        AsTwoQuery::Cartier(AsTwoCartierQuery::OrdinaryRanks)
                    }
                    Token::Word("rvASTwoHasseWitt") => {
                        AsTwoQuery::Cartier(AsTwoCartierQuery::HasseWitt)
                    }
                    Token::Word("rvASTwoReduced") => AsTwoQuery::Reduced,
                    Token::Word("rvASTwoGenus") => AsTwoQuery::Genus,
                    Token::Word("rvASTwoPrank") => AsTwoQuery::PRank,
                    Token::Word("rvASTwoLPolynomial") => {
                        AsTwoQuery::Character(AsCharacterQuery::LPolynomial)
                    }
                    Token::Word("rvASTwoNewton") => AsTwoQuery::Character(AsCharacterQuery::Newton),
                    _ => {
                        self.expect(Token::Comma)?;
                        let index = self.number()?.0;
                        let query = match operation {
                            Token::Word("rvASTwoCharacterL") => AsCharacterQuery::Character(index),
                            Token::Word("rvASTwoCharacterNewton") => {
                                AsCharacterQuery::CharacterNewton(index)
                            }
                            Token::Word("rvASTwoPoints") => AsCharacterQuery::Points(index),
                            _ => unreachable!(),
                        };
                        AsTwoQuery::Character(query)
                    }
                };
                self.expect(Token::RBracket)?;
                StatementForm::RvAsTwoQuery(query, name)
            }
            Token::Word("rvperiod") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let variables = self.number()?.0;
                self.expect(Token::RBracket)?;
                self.expect(Token::LBrace)?;
                let mut terms = Vec::new();
                while self.token() != Token::RBrace {
                    self.expect(Token::LParen)?;
                    let exponents =
                        self.signed_tuple("a signed Laurent exponent", "Laurent exponent")?;
                    self.expect(Token::Comma)?;
                    let value =
                        self.signed_wide("a signed Laurent coefficient", "Laurent coefficient")?;
                    self.expect(Token::RParen)?;
                    terms.push((exponents, value));
                    if !self.eat(Token::Comma)? {
                        break;
                    }
                }
                self.expect(Token::RBrace)?;
                StatementForm::RvPeriod(name, variables, terms)
            }
            Token::Word("rvCTExact") => {
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                self.expect(Token::Comma)?;
                let index = self.wide_unsigned("period index")?;
                self.expect(Token::RBracket)?;
                StatementForm::RVCTExact(name, index)
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
            Token::Word("rvcrystal") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let (prime, genus) = self.bidegree()?;
                self.expect(Token::RBracket)?;
                self.expect(Token::LBrace)?;
                let mut entries = Vec::new();
                while self.token() != Token::RBrace {
                    self.expect(Token::LParen)?;
                    let (row, column) = self.bidegree()?;
                    self.expect(Token::Comma)?;
                    let value =
                        self.signed_wide("a signed matrix entry", "F-crystal matrix entry")?;
                    self.expect(Token::RParen)?;
                    entries.push((row, column, value));
                    if !self.eat(Token::Comma)? {
                        break;
                    }
                }
                self.expect(Token::RBrace)?;
                StatementForm::RvCrystal(name, prime, genus, entries)
            }
            Token::Word("rvCrystalChar" | "rvCrystalNewton" | "rvCrystalMatrices") => {
                let operation = self.token();
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                let query = match operation {
                    Token::Word("rvCrystalChar") => CrystalQuery::Characteristic,
                    Token::Word("rvCrystalNewton") => CrystalQuery::Newton,
                    Token::Word("rvCrystalMatrices") => {
                        self.expect(Token::Comma)?;
                        CrystalQuery::Matrices(self.number()?.0)
                    }
                    _ => unreachable!(),
                };
                self.expect(Token::RBracket)?;
                StatementForm::RvCrystalQuery(name, query)
            }
            Token::Word("rvResDim" | "rvResContains" | "rvResSupport") => {
                use crate::rv_res_space::Query;
                let operation = self.token();
                self.advance()?;
                let query = match operation {
                    Token::Word("rvResDim") => Query::Dimension,
                    Token::Word("rvResSupport") => Query::Support,
                    _ => {
                        self.expect(Token::LBracket)?;
                        let vector =
                            self.signed_tuple("an integer coordinate", "residue coordinate")?;
                        self.expect(Token::RBracket)?;
                        Query::Contains(vector)
                    }
                };
                StatementForm::ResReport(query, self.paren_expr(0)?)
            }
            Token::Word("rvResLeq") => {
                self.advance()?;
                let [first, second] = self.paren_pair(0)?;
                StatementForm::ResLeq(first, second)
            }
            Token::Word("rvqpadic") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let p = self.number()?.0;
                self.expect(Token::Comma)?;
                let d = self.signed_wide("a signed radicand", "quadratic radicand")?;
                self.expect(Token::Comma)?;
                let r = self.signed_wide("a signed residue", "quadratic branch residue")?;
                self.expect(Token::RBracket)?;
                StatementForm::RvQpQuad(name, p, d, r)
            }
            Token::Word(
                op @ ("rvValPlucker" | "rvValSupport" | "rvValInitial" | "rvValExchange"
                | "rvValTie"),
            ) => {
                self.advance()?;
                use crate::tropical_matroid::Query;
                let query = if op == "rvValSupport" {
                    Query::Support
                } else {
                    self.expect(Token::LBracket)?;
                    let query = match op {
                        "rvValPlucker" => Query::Plucker(self.exponents()?),
                        "rvValInitial" => Query::Initial(self.tropical_weights()?),
                        _ => {
                            let left = self.exponents()?;
                            self.expect(Token::Comma)?;
                            let right = self.exponents()?;
                            if op == "rvValExchange" {
                                self.expect(Token::Comma)?;
                                Query::Exchange(left, right, self.number()?.0)
                            } else {
                                Query::Tie(left, right)
                            }
                        }
                    };
                    self.expect(Token::RBracket)?;
                    query
                };
                StatementForm::ValMatReport(query, self.paren_expr(0)?)
            }
            Token::Word(
                op @ ("rvHyperFrobenius"
                | "rvHyperFixed"
                | "rvHyperLocal"
                | "rvHyperLocalFactors"
                | "rvHyperNodalPoints"),
            ) => {
                use crate::hyperelliptic_graph::arith::Query;
                self.advance()?;
                let query = match op {
                    "rvHyperFrobenius" => Query::Frobenius,
                    "rvHyperFixed" => Query::Fixed,
                    "rvHyperLocal" => Query::Local,
                    "rvHyperLocalFactors" => Query::LocalFactors,
                    _ => {
                        self.expect(Token::LBracket)?;
                        let (degree, _) = self.number()?;
                        self.expect(Token::RBracket)?;
                        Query::NodalPoints(degree)
                    }
                };
                StatementForm::HyperArithReport(query, self.paren_expr(0)?)
            }
            Token::Word(
                op @ ("rvHyperGraph" | "rvHyperHarmonic" | "rvHyperGenus" | "rvHyperCycles"
                | "rvHyperMetric" | "rvHyperMonodromy" | "rvHyperComponent"),
            ) => {
                use crate::hyperelliptic_graph::Query;
                self.advance()?;
                let query = match op {
                    "rvHyperGraph" => Query::Graph,
                    "rvHyperHarmonic" => Query::Harmonic,
                    "rvHyperGenus" => Query::Genus,
                    "rvHyperCycles" => Query::Cycles,
                    _ => {
                        let multiplier = if self.eat(Token::LBracket)? {
                            let n = self.number()?.0;
                            self.expect(Token::RBracket)?;
                            n
                        } else {
                            1
                        };
                        match op {
                            "rvHyperMetric" => Query::Metric(multiplier),
                            "rvHyperMonodromy" => Query::Monodromy(multiplier),
                            _ => Query::Component(multiplier),
                        }
                    }
                };
                StatementForm::HyperReport(query, self.paren_expr(0)?)
            }
            Token::Word(
                op @ ("rvStableSplits" | "rvStableDirections" | "rvStableComponents"
                | "rvStableResidue"),
            ) => {
                self.advance()?;
                use crate::stable_marks::Query;
                let query = match op {
                    "rvStableSplits" => Query::Splits,
                    "rvStableComponents" => Query::Components,
                    _ => {
                        self.expect(Token::LBracket)?;
                        let vertex = self.number()?.0;
                        let query = if op == "rvStableDirections" {
                            Query::Directions(vertex)
                        } else {
                            self.expect(Token::Comma)?;
                            Query::Residue(vertex, self.number()?.0)
                        };
                        self.expect(Token::RBracket)?;
                        query
                    }
                };
                StatementForm::StableReport(query, self.paren_expr(0)?)
            }
            Token::Word("rvlegendrecm") => {
                self.advance()?;
                let field = self.label()?;
                StatementForm::LegendreCM(field, self.paren_expr(0)?)
            }
            Token::Word("rvCMOrder" | "rvCMTrace" | "rvCMIndex" | "rvCMNorm") => {
                use crate::legendre_cm::Query;
                let kind = self.token();
                self.advance()?;
                self.expect(Token::LBracket)?;
                let field = self.label()?;
                self.expect(Token::Comma)?;
                let exponent = self.number()?.0;
                let query = match kind {
                    Token::Word("rvCMOrder") => Query::Order,
                    Token::Word("rvCMTrace") => Query::Trace,
                    Token::Word("rvCMNorm") => Query::Norm,
                    _ => {
                        self.expect(Token::Comma)?;
                        Query::Index(self.number()?.0)
                    }
                };
                self.expect(Token::RBracket)?;
                StatementForm::CMReport(query, field, exponent)
            }
            Token::Word("rvQpQuadVal" | "rvQpQuadResidue") => {
                let valuation = self.token() == Token::Word("rvQpQuadVal");
                self.advance()?;
                let precision = if valuation {
                    None
                } else {
                    self.expect(Token::LBracket)?;
                    let precision = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    Some(precision)
                };
                StatementForm::RvQpQuadReport(precision, self.paren_expr(0)?)
            }
            Token::Word("rvphodgequad") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let field = self.label()?;
                self.expect(Token::RBracket)?;
                StatementForm::RvPHodgeQuad(name, field, self.phodge_quad_spec()?)
            }
            Token::Word("rvphodge") => {
                self.advance()?;
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let prime = self.number()?.0;
                self.expect(Token::RBracket)?;
                StatementForm::RvPHodge(name, self.phodge_spec(prime)?)
            }
            Token::Word("rvPHodgeData" | "rvPHodgeWeak" | "rvPHodgeHN") => {
                let operation = self.token();
                self.advance()?;
                self.expect(Token::LBracket)?;
                let name = self.label()?;
                let query = match operation {
                    Token::Word("rvPHodgeData") => {
                        let coordinates = self
                            .eat(Token::Comma)?
                            .then(|| self.exponents())
                            .transpose()?;
                        PHodgeQuery::Data(coordinates)
                    }
                    Token::Word("rvPHodgeWeak") => PHodgeQuery::WeakAdmissibility,
                    _ => PHodgeQuery::HarderNarasimhan,
                };
                self.expect(Token::RBracket)?;
                StatementForm::RvPHodgeQuery(name, query)
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
        if name.text == "SetNat" && self.eat(Token::LParen)? {
            let source = self.label()?;
            self.expect(Token::Comma)?;
            let target = self.label()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::SetNat(source, target), 0);
        }
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
        if matches!(name.text.as_str(), "RMRhoCobar" | "RMRhoExt") && self.eat(Token::LParen)? {
            let grade = self.real_rho_grade()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RealRho(name.text == "RMRhoExt", grade), 0);
        }
        if matches!(
            name.text.as_str(),
            "MCobar" | "MExt" | "RMCobar" | "RMExt" | "RMBar" | "RMTriple" | "RhoLoc" | "RhoMin"
        ) && self.eat(Token::LParen)?
        {
            let grade = self.motivic_grade()?;
            self.expect(Token::RParen)?;
            let form = match name.text.as_str() {
                "MCobar" => TypeForm::MCobar(grade),
                "RMCobar" => TypeForm::RMCobar(grade),
                "RMExt" => TypeForm::RMExt(grade),
                "RMBar" => TypeForm::RMBar(grade),
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
        if matches!(name.text.as_str(), "MixedChain" | "MixedU" | "MixedH")
            && self.eat(Token::LParen)?
        {
            let complex = self.label()?;
            self.expect(Token::Comma)?;
            let degree = self.signed_degree()?;
            self.expect(Token::RParen)?;
            return ty(
                span,
                if name.text == "MixedH" {
                    TypeForm::MixedHomology(complex, degree)
                } else {
                    TypeForm::MixedChain(complex, degree, name.text == "MixedU")
                },
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
        if name.text == "ShDef" && self.eat(Token::LParen)? {
            let sheaf = self.label()?;
            self.expect(Token::Comma)?;
            let order = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::SheafDeformation(sheaf, order), 0);
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
        if let Some(form) = crate::sheaf_disk::convolution::Form::detect(&name.text) {
            if self.eat(Token::LParen)? {
                let sequence = self.label()?;
                self.expect(Token::RParen)?;
                return ty(span, TypeForm::Convolution(form, sequence), 0);
            }
        }
        if let Some(form) = crate::sheaf_disk::parabolic::Form::detect(&name.text) {
            if self.eat(Token::LParen)? {
                let sequence = self.label()?;
                self.expect(Token::RParen)?;
                return ty(span, TypeForm::Parabolic(form, sequence), 0);
            }
        }
        if name.text == "ShDiskHom" && self.eat(Token::LParen)? {
            let from = self.label()?;
            self.expect(Token::Comma)?;
            let to = self.label()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::DiskHom(from, to), 0);
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
        if name.text == "ASTwoDiff" && self.eat(Token::LParen)? {
            let cover = self.label()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::ASTwoDiff(cover), 0);
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
        if name.text == "RVTL" && self.eat(Token::LParen)? {
            let (level, strands, height) = self.tl_grade()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RvTl(level, strands, height), 0);
        }
        if name.text == "RVResSpace" && self.eat(Token::LParen)? {
            let (p, n) = self.bidegree()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVResSpace(p, n), 0);
        }
        if name.text == "RVHyperArithmetic" && self.eat(Token::LParen)? {
            let (p, marks) = self.bidegree()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::HyperArithmetic(p, marks), 0);
        }
        if name.text == "RVHyperCover" && self.eat(Token::LParen)? {
            let (p, marks) = self.bidegree()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::HyperCover(p, marks), 0);
        }
        if name.text == "RVStableTree" && self.eat(Token::LParen)? {
            let count = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::StableTree(count), 0);
        }
        if name.text == "RVStableMarks" && self.eat(Token::LParen)? {
            let p = self.number()?.0;
            self.expect(Token::Comma)?;
            let count = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::StableMarks(p, count), 0);
        }
        if name.text == "RVValMat" && self.eat(Token::LParen)? {
            let p = self.number()?.0;
            self.expect(Token::Comma)?;
            let n = self.number()?.0;
            self.expect(Token::Comma)?;
            let r = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVValMat(p, n, r), 0);
        }
        if name.text == "RVQpQuad" && self.eat(Token::LParen)? {
            let field = self.label()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVQpQuad(field), 0);
        }
        if name.text == "RVCMOrder" && self.eat(Token::LParen)? {
            let field = self.label()?;
            self.expect(Token::Comma)?;
            let exponent = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::CMOrder(field, exponent), 0);
        }
        if name.text == "RVECEnd" && self.eat(Token::LParen)? {
            let certificate = self.label()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::ECEndOrder(certificate), 0);
        }
        if name.text == "RVCMClass" && self.eat(Token::LParen)? {
            let order = self.label()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::CMClass(order), 0);
        }
        if matches!(
            name.text.as_str(),
            "RVECCurve" | "RVECPoint" | "RVECKernel" | "RVECIso"
        ) && self.eat(Token::LParen)?
        {
            let label = self.label()?;
            let form = match name.text.as_str() {
                "RVECCurve" => TypeForm::ECCurve(label),
                "RVECPoint" => TypeForm::ECPoint(label),
                "RVECIso" => TypeForm::ECIso(label),
                _ => {
                    self.expect(Token::Comma)?;
                    TypeForm::ECKernel(label, self.number()?.0)
                }
            };
            self.expect(Token::RParen)?;
            return ty(span, form, 0);
        }
        if name.text == "RVKSEnd" && self.eat(Token::LParen)? {
            let cusp = self.label()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::KSEnd(cusp), 0);
        }
        if name.text == "RVKSExterior" && self.eat(Token::LParen)? {
            let cusp = self.label()?;
            self.expect(Token::Comma)?;
            let degree = self.signed_wide(
                "a signed cohomological degree",
                "Kuga–Satake exterior degree",
            )?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::KSExterior(cusp, degree), 0);
        }
        if name.text == "RVOG6Clifford" && self.eat(Token::LParen)? {
            let cusp = self.label()?;
            self.expect(Token::Comma)?;
            let parity = self.number()?.0;
            self.expect(Token::RParen)?;
            if parity > 1 {
                return Err(Error::new(span, "OG6 Clifford parity must be 0 or 1"));
            }
            return ty(span, TypeForm::OG6Clifford(cusp, parity as u8), 0);
        }
        if name.text == "RVK3Clifford" && self.eat(Token::LParen)? {
            let basis = self.label()?;
            self.expect(Token::Comma)?;
            let (parity, at) = self.number()?;
            if parity > 1 {
                return Err(Error::new(at, "Clifford parity must be 0 or 1"));
            }
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::K3Clifford(basis, parity as u8), 0);
        }
        if name.text == "RVK3Spinor" && self.eat(Token::LParen)? {
            let basis = self.label()?;
            self.expect(Token::Comma)?;
            let chirality = self.spinor_chirality()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::K3Spinor(basis, chirality), 0);
        }
        if name.text == "RVECMap" && self.eat(Token::LParen)? {
            let source = self.label()?;
            self.expect(Token::Comma)?;
            let target = self.label()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::ECMap(source, target), 0);
        }
        if name.text == "RVPeriod" && self.eat(Token::LParen)? {
            let grade = self.period_grade()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::PeriodPolynomial(grade), 0);
        }
        if name.text == "RVHeckeEigen" && self.eat(Token::LParen)? {
            let packet = self.label()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::HeckeEigen(packet), 0);
        }
        if name.text == "RVHecke" && self.eat(Token::LParen)? {
            let weight = self.wide_unsigned("level-one Hecke weight")?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::ModularHecke(weight), 0);
        }
        if name.text == "RVSiegelWeight" {
            return ty(span, TypeForm::SiegelWeight, 0);
        }
        if matches!(
            name.text.as_str(),
            "RVJac"
                | "RVKummer"
                | "RVKummerK3"
                | "RVEnriques"
                | "RVK3PencilBase"
                | "RVK3PencilSection"
                | "RVK3PencilKummer2"
                | "RVKummerDisc"
                | "RVKummerDiv"
                | "RVK3Lattice"
                | "RVMukai"
                | "RVKummerWeil"
                | "RVMukaiQuotient"
                | "RVOG10"
                | "RVAbMukai"
                | "RVOG6"
        ) && self.eat(Token::LParen)?
        {
            let endpoint = self.label()?;
            self.expect(Token::RParen)?;
            let form = match name.text.as_str() {
                "RVJac" => TypeForm::Jacobian(endpoint),
                "RVKummer" => TypeForm::Kummer(endpoint),
                "RVKummerK3" => TypeForm::KummerK3(endpoint),
                "RVEnriques" => TypeForm::Enriques(endpoint),
                "RVK3PencilBase" => TypeForm::K3PencilBase(endpoint),
                "RVK3PencilSection" => TypeForm::K3PencilSection(endpoint),
                "RVK3PencilKummer2" => TypeForm::K3PencilKummer2(endpoint),
                "RVKummerDiv" => TypeForm::KummerDiv(endpoint),
                "RVMukai" => TypeForm::Mukai(endpoint),
                "RVMukaiQuotient" => TypeForm::MukaiQuotient(endpoint),
                "RVOG10" => TypeForm::OG10(endpoint),
                "RVAbMukai" => TypeForm::AbMukai(endpoint),
                "RVOG6" => TypeForm::OG6(endpoint),
                "RVKummerDisc" => TypeForm::KummerDisc(endpoint),
                "RVKummerWeil" => TypeForm::KummerWeil(endpoint),
                _ => TypeForm::K3Lattice(endpoint),
            };
            return ty(span, form, 0);
        }
        if matches!(name.text.as_str(), "RVNakajima" | "RVHilbCoh") && self.eat(Token::LParen)? {
            let (basis, length, degree) = self.nakajima_grade_args()?;
            self.expect(Token::RParen)?;
            let form = if name.text == "RVNakajima" {
                TypeForm::Nakajima(basis, length, degree)
            } else {
                TypeForm::HilbertCohomology(basis, length, degree)
            };
            return ty(span, form, 0);
        }
        if name.text == "RVHilbSquare" && self.eat(Token::LParen)? {
            let basis = self.label()?;
            self.expect(Token::Comma)?;
            let degree =
                self.signed_wide("a signed cohomological degree", "Hilbert-square degree")?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::HilbertSquare(basis, degree), 0);
        }
        if name.text == "RVHWMatrix" && self.eat(Token::LParen)? {
            let (p, g, s) = self.hw_shape()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::HWMatrix(p, g, s), 0);
        }
        if name.text == "RVHWWMatrix" && self.eat(Token::LParen)? {
            let ring = self.label()?;
            self.expect(Token::Comma)?;
            let genus = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::HWWMatrix(ring, genus), 0);
        }
        if name.text == "RVHWScalar" && self.eat(Token::LParen)? {
            let p = self.number()?.0;
            self.expect(Token::Comma)?;
            let s = self.number()?.0;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::HWScalar(p, s), 0);
        }
        if name.text == "RVUnram" && self.eat(Token::LParen)? {
            let ring = self.label()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVUnram(ring), 0);
        }
        if name.text == "LeeDual" && self.eat(Token::LParen)? {
            let braid = self.label()?;
            self.expect(Token::Comma)?;
            let i = self.signed_number("Lee dual homological degree")?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::LeeDual(braid, i), 0);
        }
        if name.text == "LeeC" && self.eat(Token::LParen)? {
            let braid = self.label()?;
            self.expect(Token::Comma)?;
            let grade = self.lee_grade()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::Lee(braid, grade), 0);
        }
        if matches!(name.text.as_str(), "KhPC" | "KhPage") && self.eat(Token::LParen)? {
            let braid = self.label()?;
            self.expect(Token::Comma)?;
            let p = self.number()?.0;
            let r = if name.text == "KhPage" {
                self.expect(Token::Comma)?;
                Some(self.number()?.0)
            } else {
                None
            };
            let grade = self.kh_grade()?;
            let k = self.kh_k()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::KhovanovMod(braid, p, r, grade, k), 0);
        }
        if matches!(
            name.text.as_str(),
            "KhC" | "KhH" | "KhAC" | "KhAH" | "KhFC" | "KhFD"
        ) && self.eat(Token::LParen)?
        {
            use crate::khovanov::Form;
            let braid = self.label()?;
            let grade = self.kh_grade()?;
            let form = match name.text.as_str() {
                "KhC" => Form::Chain,
                "KhH" => Form::Class,
                "KhAC" => Form::AnnularChain(self.kh_k()?),
                "KhAH" => Form::AnnularClass(self.kh_k()?),
                "KhFC" => Form::Filtered(self.kh_k()?),
                _ => Form::Dual(self.kh_k()?),
            };
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::Khovanov(braid, grade, form), 0);
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
        if name.text == "RVCTState" && self.eat(Token::LParen)? {
            let period = self.label()?;
            self.expect(Token::Comma)?;
            let (prime, exponent) = self.cartier_level()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVCTState(period, prime, exponent), 0);
        }
        if name.text == "RVCTResidue" && self.eat(Token::LParen)? {
            let (prime, exponent) = self.cartier_level()?;
            self.expect(Token::RParen)?;
            return ty(span, TypeForm::RVCTResidue(prime, exponent), 0);
        }
        if name.text == "RVOre" {
            return ty(span, TypeForm::RVOre, 0);
        }
        if name.text == "RVWittSpectrum" {
            return ty(span, TypeForm::RVWittSpectrum, 0);
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
            Token::Word(name) if crate::set_topology::Operation::detect(name).is_some() => self
                .set_topology_expression(
                    crate::set_topology::Operation::detect(name).unwrap(),
                    span,
                    nest,
                ),
            Token::Word("setJUnit" | "setJInverse" | "setJExtend") => {
                let operation = self.label()?;
                self.expect(Token::LBracket)?;
                let family = self.label()?;
                self.expect(Token::RBracket)?;
                if operation.text == "setJExtend" {
                    let argument = self.paren_expr(nest)?;
                    let depth = argument.depth + 1;
                    expr(
                        span,
                        ExprForm::SetJExtend(family, Box::new(argument)),
                        depth,
                    )
                } else {
                    expr(
                        span,
                        ExprForm::SetJUnit(family, operation.text == "setJInverse"),
                        0,
                    )
                }
            }
            Token::Word(name) if SetLogic::detect(name).is_some() => {
                let operation = SetLogic::detect(name).unwrap();
                self.advance()?;
                let mut labels = Vec::new();
                let label_count = operation.label_count();
                if label_count != 0 {
                    self.expect(Token::LBracket)?;
                    for _ in 0..label_count {
                        if !labels.is_empty() {
                            self.expect(Token::Comma)?;
                        }
                        labels.push(self.label()?);
                    }
                    self.expect(Token::RBracket)?;
                }
                let mut inputs = Vec::new();
                let arity = operation.arity();
                let mut depth = 0;
                if arity != 0 {
                    self.expect(Token::LParen)?;
                    for _ in 0..arity {
                        if !inputs.is_empty() {
                            self.expect(Token::Comma)?;
                        }
                        let input = self.expr(nest + 1)?;
                        depth = depth.max(input.depth + 1);
                        inputs.push(input);
                    }
                    self.expect(Token::RParen)?;
                }
                expr(span, ExprForm::SetLogic(operation, labels, inputs), depth)
            }
            Token::Word(name) if crate::set_gluing::Operation::detect(name).is_some() => {
                let operation = crate::set_gluing::Operation::detect(name).unwrap();
                self.advance()?;
                self.expect(Token::LBracket)?;
                let first = self.label()?;
                self.expect(Token::Comma)?;
                let second = self.label()?;
                self.expect(Token::RBracket)?;
                let mut inputs = Vec::new();
                let mut depth = 0;
                if operation.arity() != 0 {
                    self.expect(Token::LParen)?;
                    for i in 0..operation.arity() {
                        if i != 0 {
                            self.expect(Token::Comma)?;
                        }
                        let input = self.expr(nest + 1)?;
                        depth = depth.max(input.depth + 1);
                        inputs.push(input);
                    }
                    self.expect(Token::RParen)?;
                }
                expr(
                    span,
                    ExprForm::SetGluing(operation, [first, second], inputs),
                    depth,
                )
            }
            Token::Word(name) if KanOp::detect(name).is_some() => {
                let operation = KanOp::detect(name).unwrap();
                self.advance()?;
                self.expect(Token::LBracket)?;
                let first = self.label()?;
                let second = if operation.second_label() {
                    self.expect(Token::Comma)?;
                    Some(self.label()?)
                } else {
                    None
                };
                self.expect(Token::RBracket)?;
                let input = if operation.takes_input() {
                    Some(Box::new(self.paren_expr(nest)?))
                } else {
                    None
                };
                let depth = input.as_ref().map_or(0, |expr| expr.depth + 1);
                expr(
                    span,
                    ExprForm::SetKan(operation, first, second, input),
                    depth,
                )
            }
            Token::Word(name) if SetOperation::detect(name).is_some() => {
                let operation = SetOperation::detect(name).unwrap();
                self.advance()?;
                self.expect(Token::LBracket)?;
                let first = self.label()?;
                let second = match operation {
                    SetOperation::PullNatural | SetOperation::Evaluation => {
                        self.expect(Token::Comma)?;
                        Some(self.label()?)
                    }
                    _ => None,
                };
                self.expect(Token::RBracket)?;
                let input = if matches!(operation, SetOperation::Evaluation) {
                    None
                } else {
                    Some(Box::new(self.paren_expr(nest)?))
                };
                let depth = input.as_ref().map_or(0, |expr| expr.depth + 1);
                expr(
                    span,
                    ExprForm::SetDependent(operation, first, second, input),
                    depth,
                )
            }
            Token::Word("setProjection" | "setsect" | "setId" | "setTotal") => {
                let operation = self.label()?;
                self.expect(Token::LBracket)?;
                let first = self.label()?;
                let (form, depth) = match operation.text.as_str() {
                    "setProjection" => {
                        self.expect(Token::RBracket)?;
                        (ExprForm::SetProjection(first), 0)
                    }
                    "setId" => {
                        self.expect(Token::RBracket)?;
                        (ExprForm::SetIdentity(first), 0)
                    }
                    "setsect" => {
                        self.expect(Token::Comma)?;
                        let index = self.number()?.0;
                        self.expect(Token::RBracket)?;
                        (ExprForm::SetSection(first, index), 0)
                    }
                    "setTotal" => {
                        self.expect(Token::Comma)?;
                        let second = self.label()?;
                        self.expect(Token::RBracket)?;
                        let map = self.paren_expr(nest)?;
                        let depth = map.depth + 1;
                        (ExprForm::SetTotal(first, second, Box::new(map)), depth)
                    }
                    _ => unreachable!(),
                };
                expr(span, form, depth)
            }
            Token::Word("setComp") => {
                self.advance()?;
                self.expect(Token::LParen)?;
                let outer = self.expr(nest + 1)?;
                self.expect(Token::Comma)?;
                let inner = self.expr(nest + 1)?;
                self.expect(Token::RParen)?;
                let depth = outer.depth.max(inner.depth) + 1;
                expr(span, ExprForm::SetCompose(Box::new([outer, inner])), depth)
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
                if matches!(
                    name.text.as_str(),
                    "mixConeI" | "mixConeLift" | "mixConeP" | "mixConeConnect"
                ) && self.eat(Token::LBracket)?
                {
                    use crate::mixed_cone::Operation;
                    let operation = match name.text.as_str() {
                        "mixConeI" => Operation::Include,
                        "mixConeLift" => Operation::Lift,
                        "mixConeP" => Operation::Project,
                        _ => Operation::Connect,
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
                if name.text == "mixHgen" && self.eat(Token::LBracket)? {
                    return self.mixed_h_generator(span);
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
                if name.text == "rvCrystalEOType" && self.eat(Token::LBracket)? {
                    let crystal = self.label()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::CrystalEO(crystal), 0);
                }
                if name.text == "rvCrystalCongruent" && self.eat(Token::LBracket)? {
                    let left = self.label()?;
                    self.expect(Token::Comma)?;
                    let right = self.label()?;
                    self.expect(Token::Comma)?;
                    let level = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::CrystalCongruent(left, right, level), 0);
                }
                if name.text == "rvResSpace" && self.eat(Token::LBracket)? {
                    let (p, n) = self.bidegree()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LBrace)?;
                    let mut rows = Vec::new();
                    while self.token() != Token::RBrace {
                        rows.push(self.signed_tuple("an integer entry", "residue matrix entry")?);
                        if !self.eat(Token::Comma)? {
                            break;
                        }
                    }
                    self.expect(Token::RBrace)?;
                    return expr(span, ExprForm::ResSpace(p, n, rows), 0);
                }
                if name.text == "rvResSpecialize" && self.eat(Token::LBracket)? {
                    let (p, n) = self.bidegree()?;
                    self.expect(Token::Comma)?;
                    let weights = self.tropical_weights()?;
                    self.expect(Token::RBracket)?;
                    let rows = self.residue_rational_rows()?;
                    return expr(span, ExprForm::ResSpecialize(p, n, weights, rows), 0);
                }
                if name.text == "rvResSpecializeQuad" && self.eat(Token::LBracket)? {
                    let field = self.label()?;
                    self.expect(Token::Comma)?;
                    let n = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let weights = self.tropical_weights()?;
                    self.expect(Token::RBracket)?;
                    let (rows, depth) = self.tropical_rows(nest)?;
                    return expr(
                        span,
                        ExprForm::ResSpecializeQuad(field, n, weights, rows),
                        depth,
                    );
                }
                if name.text == "rvResFil" && self.eat(Token::LBracket)? {
                    let module = self.label()?;
                    self.expect(Token::Comma)?;
                    let cutoff = self.signed_degree()?;
                    self.expect(Token::Comma)?;
                    let weights = self.tropical_weights()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::ResFil(module, cutoff, weights), 0);
                }
                if name.text == "rvResOrthogonal" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::ResOrthogonal(Box::new(input)), depth);
                }
                if name.text == "rvStableLine" && self.eat(Token::LBracket)? {
                    let p = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LBrace)?;
                    let mut marks = Vec::new();
                    while self.token() != Token::RBrace {
                        let row = self
                            .signed_tuple("an integer marking coordinate", "marking coordinate")?;
                        marks.push(row.try_into().map_err(|_| {
                            Error::new(
                                span,
                                "a projective mark requires two homogeneous coordinates",
                            )
                        })?);
                        if !self.eat(Token::Comma)? {
                            break;
                        }
                    }
                    self.expect(Token::RBrace)?;
                    return expr(span, ExprForm::StableLine(p, marks), 0);
                }
                if name.text == "rvStableLineQuad" && self.eat(Token::LBracket)? {
                    let field = self.label()?;
                    self.expect(Token::RBracket)?;
                    let (rows, depth) = self.tropical_rows(nest)?;
                    let marks = rows
                        .into_iter()
                        .map(|row| {
                            row.try_into().map_err(|_| {
                                Error::new(
                                    span,
                                    "a projective mark requires two homogeneous coordinates",
                                )
                            })
                        })
                        .collect::<Result<Vec<_>>>()?;
                    return expr(span, ExprForm::StableLineQuad(field, marks), depth);
                }
                if name.text == "rvHyperArithmetic" && self.eat(Token::LBracket)? {
                    let (p, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let leading = self.rational_entry(
                        "hyperelliptic leading numerator",
                        "hyperelliptic leading denominator",
                    )?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LBrace)?;
                    let mut roots = Vec::new();
                    while self.token() != Token::RBrace {
                        roots.push(self.rational_entry(
                            "hyperelliptic root numerator",
                            "hyperelliptic root denominator",
                        )?);
                        if !self.eat(Token::Comma)? {
                            break;
                        }
                    }
                    self.expect(Token::RBrace)?;
                    return expr(span, ExprForm::HyperArithmetic(p, leading, roots), 0);
                }
                if name.text == "rvHyperGeometric" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::HyperGeometric(Box::new(input)), depth);
                }
                if name.text == "rvHyperCover" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::HyperCover(Box::new(input)), depth);
                }
                if name.text == "rvStableTree" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::StableProjection(Box::new(input)), depth);
                }
                if name.text == "rvValMat" && self.eat(Token::LBracket)? {
                    let p = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let n = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LBrace)?;
                    let mut rows = Vec::new();
                    while self.token() != Token::RBrace {
                        rows.push(self.signed_tuple("an integer row entry", "tropical row entry")?);
                        if !self.eat(Token::Comma)? {
                            break;
                        }
                    }
                    self.expect(Token::RBrace)?;
                    return expr(span, ExprForm::ValMat(p, n, rows), 0);
                }
                if name.text == "rvValMatQuad" && self.eat(Token::LBracket)? {
                    let field = self.label()?;
                    self.expect(Token::Comma)?;
                    let n = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let (rows, depth) = self.tropical_rows(nest)?;
                    return expr(span, ExprForm::ValMatQuad(field, n, rows), depth);
                }
                if name.text == "rvValFil" && self.eat(Token::LBracket)? {
                    let module = self.label()?;
                    self.expect(Token::Comma)?;
                    let cutoff = self.signed_degree()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::ValFil(module, cutoff), 0);
                }
                if name.text == "rvValDual" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::ValDual(Box::new(input)), depth);
                }
                if name.text == "rvQpQuad" && self.eat(Token::LBracket)? {
                    let field = self.label()?;
                    self.expect(Token::Comma)?;
                    let a = self.signed_wide("a signed coefficient", "quadratic coefficient")?;
                    self.expect(Token::Comma)?;
                    let b = self.signed_wide("a signed coefficient", "quadratic coefficient")?;
                    self.expect(Token::Comma)?;
                    let c = self.signed_wide("a positive denominator", "quadratic denominator")?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::QpQuadElement(field, a, b, c), 0);
                }
                if name.text == "rvQpQuadSpecialize" && self.eat(Token::LBracket)? {
                    let ring = self.label()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::QpQuadSpecialize(ring, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvCMPi" && self.eat(Token::LBracket)? {
                    let field = self.label()?;
                    let exponent = if self.eat(Token::Comma)? {
                        self.number()?.0
                    } else {
                        1
                    };
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::CMPi(field, exponent), 0);
                }
                if name.text == "rvCMTerm" && self.eat(Token::LBracket)? {
                    let field = self.label()?;
                    self.expect(Token::Comma)?;
                    let exponent = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let constant =
                        self.signed_wide("a signed coefficient", "CM order coefficient")?;
                    self.expect(Token::Comma)?;
                    let generator =
                        self.signed_wide("a signed coefficient", "CM order coefficient")?;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::CMTerm(field, exponent, constant, generator),
                        0,
                    );
                }
                if name.text == "rvCMInclude" && self.eat(Token::LBracket)? {
                    let exponent = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::CMUnary(
                            crate::legendre_cm::Unary::Include(exponent),
                            Box::new(input),
                        ),
                        depth,
                    );
                }
                if name.text == "rvCMLift" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::CMUnary(crate::legendre_cm::Unary::Lift, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvQpQuadPow" && self.eat(Token::LBracket)? {
                    let power = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::QpQuadUnary(QuadUnary::Power(power), Box::new(input)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "rvQpQuadConj" | "rvQpQuadInv")
                    && self.token() == Token::LParen
                {
                    let operation = if name.text == "rvQpQuadConj" {
                        QuadUnary::Conjugate
                    } else {
                        QuadUnary::Inverse
                    };
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::QpQuadUnary(operation, Box::new(input)),
                        depth,
                    );
                }
                if (matches!(
                    name.text.as_str(),
                    "rvHilbCohUnit"
                        | "rvHilbCohPoint"
                        | "rvHilbCohQ"
                        | "rvHilbCohScale"
                        | "rvHilbCohFock"
                        | "rvHilbCohChern"
                        | "rvHilbCohCh"
                        | "rvHilbCohSchur"
                ) && self.token() == Token::LBracket)
                    || (matches!(
                        name.text.as_str(),
                        "rvHilbCohCup"
                            | "rvHilbCohLefschetzInverse"
                            | "rvHilbCohLower"
                            | "rvHilbCohLLVLower"
                            | "rvHilbCohRotate"
                    ) && self.token() == Token::LParen)
                {
                    return self.hilbert_cohomology_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvHW"
                        | "rvHWUnit"
                        | "rvHWMat"
                        | "rvHWScalar"
                        | "rvHWEntry"
                        | "rvHWChar"
                        | "rvHWReduce"
                ) && self.token() == Token::LBracket)
                    || (matches!(name.text.as_str(), "rvHWCompose" | "rvHWInv")
                        && self.token() == Token::LParen)
                {
                    return self.hw_expression(&name.text, span, nest);
                }
                if matches!(
                    name.text.as_str(),
                    "rvHWW"
                        | "rvHWWRatio"
                        | "rvHWWUnit"
                        | "rvHWWCompose"
                        | "rvHWWFrob"
                        | "rvHWWInv"
                        | "rvHWWEntry"
                        | "rvHWWChar"
                        | "rvHWWReduce"
                        | "rvHWWDescend"
                ) {
                    return self.hww_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rveckernel"
                        | "rvECCurve"
                        | "rvECA"
                        | "rvECB"
                        | "rvECIso"
                        | "rvECPoint"
                        | "rvECInfinity"
                        | "rvECScale"
                        | "rvECKCoeff"
                        | "rvECXNum"
                        | "rvECXDen"
                        | "rvECDNum"
                        | "rvECDDen"
                ) && self.token() == Token::LBracket)
                    || (matches!(
                        name.text.as_str(),
                        "rvECApply" | "rvECX" | "rvECY" | "rvECGetKernel"
                    ) && self.token() == Token::LParen)
                {
                    return self.ec_expression(&name.text, span, nest);
                }
                if (matches!(name.text.as_str(), "rvECSum" | "rvECNegate")
                    && self.token() == Token::LParen)
                    || (matches!(
                        name.text.as_str(),
                        "rvECFrobenius" | "rvECLang" | "rvECCMRecognize" | "rvECExtend"
                    ) && self.token() == Token::LBracket)
                {
                    return self.ec_cm_expression(&name.text, span, nest);
                }
                if matches!(name.text.as_str(), "rvEllFrobCoordinates" | "rvEllLocalGen")
                    && self.token() == Token::LBracket
                {
                    return self.ec_local_expression(&name.text, span);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvECId"
                        | "rvECMultMap"
                        | "rvECChart"
                        | "rvECMapXNum"
                        | "rvECMapXDen"
                        | "rvECMapYNum"
                        | "rvECMapYDen"
                        | "rvECMapKCoeff"
                ) && self.token() == Token::LBracket)
                    || (matches!(
                        name.text.as_str(),
                        "rvECMapCompose" | "rvECMapApply" | "rvECPush" | "rvECMap" | "rvECDual"
                    ) && self.token() == Token::LParen)
                {
                    return self.ec_map_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvjac"
                        | "rvJacPoint"
                        | "rvJacHalfPoint"
                        | "rvJacHalfAt"
                        | "rvRichelotPoint"
                        | "rvRichelotPush"
                        | "rvRichelotDescPush"
                        | "rvBiellipticPoint"
                        | "rvBiellipticPush"
                        | "rvBiellipticPull"
                        | "rvEllGlueLeft"
                        | "rvEllGlueRight"
                        | "rvEllGlueMap"
                        | "rvEllGlueDual"
                        | "rvJacScale"
                        | "rvJacU"
                        | "rvJacV"
                        | "rvJacKummerSupport"
                        | "rvJacExtend"
                        | "rvJacDescend"
                        | "rvJacRechart"
                        | "rvJacUnchart"
                        | "rvJacQuadTrace"
                        | "rvJacFrob"
                        | "rvJacWeilApply"
                        | "rvJacTate"
                        | "rvJacWeil"
                ) && self.token() == Token::LBracket)
                    || (matches!(
                        name.text.as_str(),
                        "rvJacWeil2" | "rvJacTate2" | "rvJacDivisible2"
                    ) && self.token() == Token::LParen)
                {
                    return self.jacobian_expression(&name.text, name.span, nest);
                }
                if matches!(
                    name.text.as_str(),
                    "rvK3PencilDiv"
                        | "rvK3PencilParam"
                        | "rvK3PencilMap"
                        | "rvK3PencilSection"
                        | "rvK3FiberJ"
                ) && self.token() == Token::LBracket
                {
                    return self.k3_pencil_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvK3PencilSectionPoint" | "rvK3PencilTorsion2At" | "rvK3PencilSelmerBasis"
                ) && self.token() == Token::LBracket)
                    || (matches!(
                        name.text.as_str(),
                        "rvK3PencilDescent2" | "rvK3PencilDivisible2" | "rvK3PencilSelmerContains"
                    ) && self.token() == Token::LParen)
                {
                    return self.k3_pencil_section_expression(&name.text, span, nest);
                }
                if matches!(
                    name.text.as_str(),
                    "rvK3MordellShiftDiv" | "rvK3MordellShift"
                ) && self.token() == Token::LBracket
                {
                    return self.k3_mordell_expression(&name.text, span, nest);
                }
                if (name.text == "rvEnriquesOrbit" && self.token() == Token::LBracket)
                    || (matches!(
                        name.text.as_str(),
                        "rvEnriquesSplit" | "rvEnriquesRepresentative" | "rvEnriquesLift"
                    ) && self.token() == Token::LParen)
                {
                    return self.enriques_expression(&name.text, span, nest);
                }
                if name.text == "rvJacTorsionCoord" && self.eat(Token::LBracket)? {
                    let basis = self.label()?;
                    self.expect(Token::Comma)?;
                    let index = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::JacobianLevelCoordinate(basis, index, Box::new(input)),
                        depth,
                    );
                }
                if (matches!(
                    name.text.as_str(),
                    "rvKummer"
                        | "rvKummerCoord"
                        | "rvKummerScale"
                        | "rvKummerRichelot"
                        | "rvKummerRichelotDesc"
                        | "rvKummerExtend"
                        | "rvKummerDescend"
                        | "rvKummerFrob"
                        | "rvKummerK3Extend"
                        | "rvKummerK3Descend"
                        | "rvKummerK3Frob"
                        | "rvKummerK3Translate"
                        | "rvKummerTropeParam"
                        | "rvKummerOnTrope"
                        | "rvKummerDirection"
                        | "rvKummerTranslate"
                        | "rvKummerThetaEntry"
                        | "rvKummerThetaFactor"
                        | "rvKummerTropeCoeff"
                ) && self.token() == Token::LBracket)
                    || (matches!(
                        name.text.as_str(),
                        "rvJacKummerPoint"
                            | "rvKummerLiftable"
                            | "rvKummerLift"
                            | "rvKummerDiffAdd"
                            | "rvKummerSmooth"
                            | "rvKummerExceptional"
                            | "rvKummerBlowdown"
                            | "rvKummerIsExceptional"
                            | "rvKummerSwitch"
                            | "rvKummerNodeInvolution"
                    ) && self.token() == Token::LParen)
                {
                    return self.kummer_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvKummerWeilScale" | "rvKummerWeilFrob" | "rvKummerThetaShell"
                ) && self.token() == Token::LBracket)
                    || (matches!(
                        name.text.as_str(),
                        "rvKummerWeilBasis" | "rvKummerWeilS" | "rvKummerWeilT"
                    ) && self.token() == Token::LParen)
                {
                    return self.kummer_weil_expression(&name.text, span, nest);
                }
                if (matches!(name.text.as_str(), "rvKummerDisc" | "rvKummerDiscFrob")
                    && self.token() == Token::LBracket)
                    || (matches!(
                        name.text.as_str(),
                        "rvKummerDiscQuadratic" | "rvKummerDiscPair"
                    ) && self.token() == Token::LParen)
                {
                    return self.kummer_disc_expression(&name.text, span, nest);
                }
                if matches!(
                    name.text.as_str(),
                    "rvKummerDiv"
                        | "rvKummerHyperplane"
                        | "rvKummerNodeDiv"
                        | "rvKummerTropeDiv"
                        | "rvKummerDivScale"
                        | "rvKummerDivFrob"
                ) && self.token() == Token::LBracket
                {
                    return self.kummer_divisor_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvCliffordOne" | "rvCliffordVolume" | "rvCliffordScale"
                ) && self.token() == Token::LBracket)
                    || (matches!(
                        name.text.as_str(),
                        "rvCliffordVector"
                            | "rvCliffordReverse"
                            | "rvCliffordMul"
                            | "rvCliffordAct"
                    ) && self.token() == Token::LParen)
                {
                    return self.clifford_expression(&name.text, span, nest);
                }
                if (matches!(name.text.as_str(), "rvK3SpinorBasis" | "rvK3SpinorScale")
                    && self.token() == Token::LBracket)
                    || (name.text == "rvK3SpinAct" && self.token() == Token::LParen)
                {
                    return self.spinor_expression(&name.text, span, nest);
                }
                if matches!(
                    name.text.as_str(),
                    "rvK3PeriodBlade" | "rvK3PeriodCM" | "rvK3CMBasis"
                ) && self.token() == Token::LBracket
                {
                    return self.k3_period_expression(&name.text, span, nest);
                }
                if name.text == "rvCMNeighbor" && self.token() == Token::LBracket {
                    return self.cm_neighbor_expression(span, nest);
                }
                if matches!(
                    name.text.as_str(),
                    "rvModularPolyEval" | "rvModularPolySpecCoeff" | "rvModularPolyWittEval"
                ) && self.token() == Token::LBracket
                {
                    return self.elliptic_modular_expression(&name.text, span, nest);
                }
                if matches!(
                    name.text.as_str(),
                    "rvCMClassPolyEval" | "rvCMClassPolyModCoeff"
                ) && self.token() == Token::LBracket
                {
                    return self.cm_class_poly_expression(&name.text, span, nest);
                }
                if matches!(name.text.as_str(), "rvCMArtinClass" | "rvCMArtinTranslate")
                    && self.token() == Token::LBracket
                {
                    return self.cm_artin_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvECEndElem" | "rvECEndFrobenius" | "rvECEndJ"
                ) && self.token() == Token::LBracket)
                    || (matches!(name.text.as_str(), "rvECEndMul" | "rvECEndConj")
                        && self.token() == Token::LParen)
                {
                    return self.ec_end_order_expression(&name.text, span, nest);
                }
                if matches!(
                    name.text.as_str(),
                    "rvECEndCanonicalJ"
                        | "rvECEndCanonicalPhi"
                        | "rvECEndUnitRoot"
                        | "rvECEndNonunitRoot"
                        | "rvECEndEmbed"
                        | "rvECEndEmbedConj"
                ) && self.token() == Token::LBracket
                {
                    return self.ec_end_canonical_expression(&name.text, span, nest);
                }
                if name.text == "rvECCMKernel" && self.eat(Token::LBracket)? {
                    let certificate = self.label()?;
                    let address = self.ec_horizontal_address(certificate)?;
                    return expr(span, ExprForm::ECHorizontal(Box::new(address)), 0);
                }
                if matches!(
                    name.text.as_str(),
                    "rvECCMCycleApply" | "rvECCMCycleElement"
                ) && self.eat(Token::LBracket)?
                {
                    use crate::elliptic_end_order::cycle::Expression;
                    let cycle = self.label()?;
                    self.expect(Token::RBracket)?;
                    let (expression, depth) = if name.text == "rvECCMCycleApply" {
                        let input = self.paren_expr(nest)?;
                        let depth = input.depth + 1;
                        (Expression::Apply(cycle, Box::new(input)), depth)
                    } else {
                        (Expression::Element(cycle), 0)
                    };
                    return expr(span, ExprForm::ECCMCycle(Box::new(expression)), depth);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvCMClass" | "rvCMPrincipal" | "rvK3CMClass" | "rvCMPow"
                ) && self.token() == Token::LBracket)
                    || (matches!(name.text.as_str(), "rvCMCompose" | "rvCMInverse")
                        && self.token() == Token::LParen)
                {
                    return self.cm_class_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvK3Lattice" | "rvK3LatticeBasis" | "rvK3Transvection"
                ) && self.token() == Token::LBracket)
                    || (name.text == "rvK3Reflect" && self.token() == Token::LParen)
                {
                    return self.k3_lattice_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvNakajimaVac"
                        | "rvNakajimaMonomial"
                        | "rvNakajimaPrimaryGen"
                        | "rvNakajimaStratum"
                        | "rvNakajimaCreate"
                        | "rvNakajimaAnn"
                        | "rvNakajimaScale"
                        | "rvVirasoro"
                        | "rvHilbTautChern"
                ) && self.token() == Token::LBracket)
                    || (matches!(name.text.as_str(), "rvNakajimaBoundary" | "rvNakajimaCup")
                        && self.token() == Token::LParen)
                {
                    return self.nakajima_expression(&name.text, span, nest);
                }
                if (matches!(name.text.as_str(), "rvAbMukai" | "rvAbMukaiBasis")
                    && self.token() == Token::LBracket)
                    || (name.text == "rvAbMukaiKummerH2" && self.token() == Token::LParen)
                {
                    return self.abelian_mukai_expression(&name.text, span, nest);
                }
                if matches!(
                    name.text.as_str(),
                    "rvOG6NilpotentPower" | "rvOG6NilpotentExp"
                ) && self.token() == Token::LBracket
                {
                    return self.og6_nilpotent_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvOG6CliffordBasis"
                        | "rvOG6PolarizedBasis"
                        | "rvOG6CliffordVector"
                        | "rvOG6CliffordScale"
                        | "rvOG6SpinLog"
                        | "rvOG6SpinExp"
                ) && self.token() == Token::LBracket)
                    || (name.text == "rvOG6CliffordMul" && self.token() == Token::LParen)
                {
                    return self.og6_clifford_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvKSExteriorBasis"
                        | "rvKSExteriorScale"
                        | "rvKSExteriorPower"
                        | "rvKSExteriorDivided"
                        | "rvKSExteriorExp"
                ) && self.token() == Token::LBracket)
                    || (matches!(
                        name.text.as_str(),
                        "rvKSExteriorCup" | "rvKSExteriorEmbed" | "rvKSExteriorToSpin"
                    ) && self.token() == Token::LParen)
                {
                    return self.ks_exterior_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvKSEndUnit" | "rvKSVolume" | "rvKSEmbed" | "rvKSEndScale"
                ) && self.token() == Token::LBracket)
                    || (matches!(
                        name.text.as_str(),
                        "rvKSEndLeft" | "rvKSProject" | "rvKSEndCompose"
                    ) && self.token() == Token::LParen)
                {
                    return self.ks_end_expression(&name.text, span, nest);
                }
                if matches!(
                    name.text.as_str(),
                    "rvOG6CuspBasis" | "rvOG6CuspTransvection" | "rvOG6CuspCentral"
                ) && self.token() == Token::LBracket
                {
                    return self.og6_cusp_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvOG6" | "rvOG6Basis" | "rvOG6Exc" | "rvOG6ExcHalf"
                ) && self.token() == Token::LBracket)
                    || (name.text == "rvOG6Reflect" && self.token() == Token::LParen)
                {
                    return self.og6_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvOG10" | "rvOG10Basis" | "rvOG10Exc" | "rvOG10WallAt"
                ) && self.token() == Token::LBracket)
                    || (name.text == "rvOG10Reflect" && self.token() == Token::LParen)
                {
                    return self.og10_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvMukaiModuliBasis"
                        | "rvMukaiModuliGlue"
                        | "rvMukaiModuliDual"
                        | "rvMukaiQuotient"
                ) && self.token() == Token::LBracket)
                    || (name.text == "rvMukaiQuotientLift" && self.token() == Token::LParen)
                {
                    return self.mukai_moduli_expression(&name.text, span, nest);
                }
                if (matches!(name.text.as_str(), "rvMukai" | "rvMukaiIdeal")
                    && self.token() == Token::LBracket)
                    || (matches!(
                        name.text.as_str(),
                        "rvK3FromKummerDiv"
                            | "rvK3ToKummerDiv"
                            | "rvMukaiC"
                            | "rvMukaiLine"
                            | "rvMukaiTensor"
                            | "rvMukaiTwist"
                    ) && self.token() == Token::LParen)
                {
                    return self.mukai_expression(&name.text, span, nest);
                }
                if matches!(name.text.as_str(), "rvK3AutoDiv" | "rvK3AutoLattice")
                    && self.eat(Token::LBracket)?
                {
                    let word = self.label()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = 1 + input.depth;
                    let expression = crate::hyperelliptic_kummer_surface::dynamics::Expression(
                        word,
                        Box::new(input),
                    );
                    let form = if name.text == "rvK3AutoLattice" {
                        ExprForm::K3AutoLattice(Box::new(expression))
                    } else {
                        ExprForm::K3AutoDiv(Box::new(expression))
                    };
                    return expr(span, form, depth);
                }
                if matches!(name.text.as_str(), "rvK3AutoPoint" | "rvK3AutoIter")
                    && self.eat(Token::LBracket)?
                {
                    let word = self.label()?;
                    let power = if name.text == "rvK3AutoIter" {
                        self.expect(Token::Comma)?;
                        self.signed_wide(
                            "a signed iteration count",
                            "K3 automorphism iteration count",
                        )?
                    } else {
                        1
                    };
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = 1 + input.depth;
                    let expression = crate::hyperelliptic_kummer_surface::dynamics::PointExpression(
                        word,
                        power,
                        Box::new(input),
                    );
                    return expr(span, ExprForm::K3AutoPoint(Box::new(expression)), depth);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvHilbSquareBasis"
                        | "rvHilbSquareDiv"
                        | "rvHilbSquareScale"
                        | "rvHilbSquareC2"
                        | "rvHilbSquareC4"
                        | "rvHilbSquareNakajima"
                ) && self.token() == Token::LBracket)
                    || (matches!(
                        name.text.as_str(),
                        "rvHilbSquareCup"
                            | "rvHilbSquareLower"
                            | "rvHilbSquareLambda"
                            | "rvHilbSquarePrimitive"
                            | "rvHilbSquareRotate"
                            | "rvNakajimaHilbSquare"
                    ) && self.token() == Token::LParen)
                {
                    return self.hilbert_square_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvPeriodBasis" | "rvPeriodPoly" | "rvPeriodScale"
                ) && self.token() == Token::LBracket)
                    || (name.text == "rvPeriodHecke" && self.token() == Token::LParen)
                {
                    return self.period_expression(&name.text, span, nest);
                }
                if matches!(
                    name.text.as_str(),
                    "rvHeckeResidualReduce" | "rvHeckeResidualCoeff" | "rvHeckeEisensteinCoeff"
                ) && self.token() == Token::LBracket
                {
                    return self.hecke_residual_expression(&name.text, span, nest);
                }
                if (matches!(
                    name.text.as_str(),
                    "rvHeckeEigenImage"
                        | "rvHeckeEigenGenerator"
                        | "rvHeckeEigenCoeff"
                        | "rvHeckeEigenScalar"
                ) && self.token() == Token::LBracket)
                    || (name.text == "rvHeckeEigenMul" && self.token() == Token::LParen)
                {
                    return self.hecke_eigen_expression(&name, span, nest);
                }
                if matches!(
                    name.text.as_str(),
                    "rvHeckeT" | "rvHeckeUnit" | "rvHeckeScale"
                ) && self.eat(Token::LBracket)?
                {
                    if name.text == "rvHeckeT" && self.token() == Token::LBracket {
                        let window = self.affine_window()?;
                        self.expect(Token::RBracket)?;
                        return expr(span, ExprForm::HeckeBasis(false, window), 0);
                    }
                    let (operation, depth) = if name.text == "rvHeckeScale" {
                        let scalar =
                            self.signed_wide("a signed Hecke coefficient", "Hecke coefficient")?;
                        self.expect(Token::RBracket)?;
                        let input = self.paren_expr(nest)?;
                        let depth = input.depth + 1;
                        (hecke::Expression::Scale(scalar, Box::new(input)), depth)
                    } else {
                        let weight = self.wide_unsigned("level-one Hecke weight")?;
                        let index = if name.text == "rvHeckeT" {
                            self.expect(Token::Comma)?;
                            self.wide_unsigned("Hecke index")?
                        } else {
                            1
                        };
                        self.expect(Token::RBracket)?;
                        (hecke::Expression::Basis(weight, index), 0)
                    };
                    return expr(span, ExprForm::ModularHecke(operation), depth);
                }
                if name.text == "rvHeckeSameAction" && self.token() == Token::LParen {
                    let [left, right] = self.paren_pair(nest)?;
                    let depth = left.depth.max(right.depth) + 1;
                    let operation = hecke::Expression::SameAction(Box::new([left, right]));
                    return expr(span, ExprForm::ModularHecke(operation), depth);
                }
                if name.text == "rvSiegelWeight" && self.eat(Token::LBracket)? {
                    let weight = self.siegel_weight_coordinates()?;
                    return expr(span, ExprForm::SiegelWeight(weight), 0);
                }
                if name.text == "rvSiegelReflect" && self.eat(Token::LBracket)? {
                    let (index, at) = self.number()?;
                    if index >= 3 {
                        return Err(Error::new(
                            at,
                            "Siegel reflection index must be between 0 and 2",
                        ));
                    }
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::SiegelReflect(index, Box::new(input)), depth);
                }
                if name.text == "rvJacRichelotKernel" && self.eat(Token::LBracket)? {
                    let (basis, base) = self.richelot_level_names()?;
                    self.expect(Token::Comma)?;
                    let plane = self.number()?.0;
                    self.expect(Token::Comma)?;
                    let generator = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::JacobianRichelotKernel(basis, base, plane, generator),
                        0,
                    );
                }
                if name.text == "rvUnramElem" && self.eat(Token::LBracket)? {
                    let ring = self.label()?;
                    self.expect(Token::RBracket)?;
                    self.expect(Token::LBrace)?;
                    let coefficients = self.signed_list(
                        Token::RBrace,
                        "a signed coefficient",
                        "unramified element",
                    )?;
                    return expr(span, ExprForm::UnramElement(ring, coefficients), 0);
                }
                if matches!(name.text.as_str(), "rvUnramFrob" | "rvUnramPow")
                    && self.eat(Token::LBracket)?
                {
                    use crate::rv_unram::Unary;
                    let k = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    let operation = if name.text == "rvUnramFrob" {
                        Unary::Frobenius(k)
                    } else {
                        Unary::Power(k as u64)
                    };
                    return expr(
                        span,
                        ExprForm::UnramUnary(operation, Box::new(input)),
                        depth,
                    );
                }
                if matches!(
                    name.text.as_str(),
                    "rvUnramInv" | "rvUnramTeich" | "rvUnramTrace" | "rvUnramNorm"
                ) && self.token() == Token::LParen
                {
                    use crate::rv_unram::Unary;
                    let operation = match name.text.as_str() {
                        "rvUnramInv" => Unary::Inverse,
                        "rvUnramTeich" => Unary::Teichmuller,
                        "rvUnramTrace" => Unary::Trace,
                        _ => Unary::Norm,
                    };
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::UnramUnary(operation, Box::new(input)),
                        depth,
                    );
                }
                if matches!(
                    name.text.as_str(),
                    "rvLegendreParameter" | "rvLegendreUnit" | "rvLegendreFactor"
                ) {
                    use crate::legendre_closed::Operation;
                    let operation = match name.text.as_str() {
                        "rvLegendreParameter" => Operation::Parameter,
                        "rvLegendreUnit" => Operation::Unit,
                        _ => {
                            self.expect(Token::LBracket)?;
                            let i = self.number()?.0;
                            self.expect(Token::RBracket)?;
                            Operation::Factor(i)
                        }
                    };
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::LegendreClosed(operation, Box::new(input)),
                        depth,
                    );
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
                if name.text == "rvKummerHilbCharCoeff" && self.token() == Token::LBracket {
                    let query = self.kummer_hilbert_characteristic()?;
                    return expr(span, ExprForm::KummerHilbertCoefficient(Box::new(query)), 0);
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
                if name.text == "rvASTwoOmega" && self.eat(Token::LBracket)? {
                    let cover = self.label()?;
                    self.expect(Token::Comma)?;
                    let i = self.signed_wide(
                        "a signed Laurent exponent",
                        "Laurent differential exponent",
                    )?;
                    self.expect(Token::Comma)?;
                    let j = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::AsTwoOmega(cover, i, j), 0);
                }
                if name.text == "rvASTwoOrdinaryGen" && self.eat(Token::LBracket)? {
                    let cover = self.label()?;
                    self.expect(Token::Comma)?;
                    let index = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::AsTwoOrdinaryGen(cover, index), 0);
                }
                if matches!(
                    name.text.as_str(),
                    "rvASTwoCartierIter" | "rvASTwoTranslate"
                ) && self.eat(Token::LBracket)?
                {
                    let operation = if name.text == "rvASTwoTranslate" {
                        AsTwoDiffOp::Translate(
                            self.signed_wide("a signed translation", "Artin–Schreier translation")?,
                        )
                    } else {
                        AsTwoDiffOp::Iterate(self.number()?.0)
                    };
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::AsTwoDiffOperation(operation, Box::new(input)),
                        depth,
                    );
                }
                if name.text == "rvASTwoCartier" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::AsTwoDiffOperation(AsTwoDiffOp::Cartier, Box::new(input)),
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
                if matches!(name.text.as_str(), "rvASOddCoeff" | "rvASTwoCoeff")
                    && self.eat(Token::LBracket)?
                {
                    let cover = self.label()?;
                    self.expect(Token::Comma)?;
                    let (character, coefficient) = self.bidegree()?;
                    self.expect(Token::RBracket)?;
                    let form = if name.text == "rvASOddCoeff" {
                        ExprForm::AsOddCoefficient(cover, character, coefficient)
                    } else {
                        ExprForm::AsTwoCoefficient(cover, character, coefficient)
                    };
                    return expr(span, form, 0);
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
                if matches!(
                    name.text.as_str(),
                    "rmbarCoords" | "rmbarClass" | "rmbarLift"
                ) && self.token() == Token::LParen
                {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    let op = match name.text.as_str() {
                        "rmbarCoords" => RealExtOp::Coordinates,
                        "rmbarClass" => RealExtOp::Class,
                        _ => RealExtOp::Lift,
                    };
                    return expr(span, ExprForm::RealExtension(op, Box::new(input)), depth);
                }
                if matches!(
                    name.text.as_str(),
                    "rhoDivide" | "rhoEinf" | "rhoProductNull"
                ) && self.eat(Token::LBracket)?
                {
                    let (index, _) = self.number()?;
                    let op = match name.text.as_str() {
                        "rhoDivide" => RealExtOp::Divide(index),
                        "rhoEinf" => {
                            if self.eat(Token::Comma)? {
                                RealExtOp::Edge(Some(index), self.number()?.0)
                            } else {
                                RealExtOp::Edge(None, index)
                            }
                        }
                        _ => {
                            self.expect(Token::RBracket)?;
                            let inputs = self.real_extension_pair(nest, span, "rhoProductNull")?;
                            let depth = inputs[0].depth.max(inputs[1].depth) + 1;
                            return expr(span, ExprForm::RhoProductNull(index, inputs), depth);
                        }
                    };
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RealExtension(op, Box::new(input)), depth);
                }
                if name.text == "rhoReduce" && self.eat(Token::LBracket)? {
                    let (exponent, at) = self.number()?;
                    if exponent == 0 {
                        return Err(Error::new(at, "rho truncation exponent must be positive"));
                    }
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::RealRhoOperation(RealRhoOp::Reduce(exponent), Box::new(input)),
                        depth,
                    );
                }
                let real_rho = match name.text.as_str() {
                    "rhoTruncLift" => Some(RealRhoOp::Lift),
                    "rhoConnect" => Some(RealRhoOp::Connect),
                    "rhoConnectRep" => Some(RealRhoOp::ConnectRep),
                    "rhoBeta" => Some(RealRhoOp::Beta),
                    "rhoBetaRep" => Some(RealRhoOp::BetaRep),
                    "rhoBetaNull" => Some(RealRhoOp::BetaNull),
                    "rhoFrob" => Some(RealRhoOp::Frobenius),
                    "rhoFrobConnectNull" => Some(RealRhoOp::FrobeniusConnectNull),
                    _ => None,
                };
                if let Some(op) = real_rho.filter(|_| self.token() == Token::LParen) {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RealRhoOperation(op, Box::new(input)), depth);
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
                let operation = match (name.text.as_str(), self.token()) {
                    ("shMassey", Token::LParen) => Some(SheafHigherOp::Class),
                    ("shCurvature", Token::LParen) => Some(SheafHigherOp::Curvature),
                    ("shMasseySolve", Token::LParen) => Some(SheafHigherOp::Solve),
                    ("shInMassey", Token::LParen) => Some(SheafHigherOp::Contains),
                    ("shMasseyWitness" | "shMasseyFind", Token::LBracket) => {
                        self.advance()?;
                        let from = self.number()?.0;
                        self.expect(Token::Comma)?;
                        let to = self.number()?.0;
                        self.expect(Token::RBracket)?;
                        Some(SheafHigherOp::Witness {
                            from,
                            to,
                            targeted: name.text == "shMasseyFind",
                        })
                    }
                    _ => None,
                };
                if let Some(operation) = operation {
                    self.expect(Token::LParen)?;
                    let inputs = self.arguments(nest)?;
                    let mut depth = inputs
                        .iter()
                        .map(|value| value.depth + 1)
                        .max()
                        .unwrap_or(0);
                    let witnesses = if operation.explicit() {
                        self.witnesses(nest, &mut depth)?
                    } else {
                        Vec::new()
                    };
                    return expr(
                        span,
                        ExprForm::SheafHigher(operation, inputs, witnesses),
                        depth,
                    );
                }
                if matches!(
                    (name.text.as_str(), self.token()),
                    (
                        "shDefStart" | "shDefTerm" | "shDefGauge" | "shDefGaugeWitness",
                        Token::LBracket
                    ) | (
                        "shDefAppend" | "shDefCurvature" | "shDefObstruction" | "shDefLift",
                        Token::LParen
                    )
                ) {
                    return self.sheaf_deformation_expr(&name, span, nest);
                }
                if matches!(
                    name.text.as_str(),
                    "shKurA" | "shKurR" | "shKurKappa" | "shKurArc" | "shKurArcKappa"
                ) && self.token() == Token::LBracket
                {
                    return self.sheaf_kuranishi_expression(&name.text, span);
                }
                if matches!(
                    name.text.as_str(),
                    "shDiskId" | "shDiskGen" | "shDiskComp" | "shDiskCostalkLift"
                ) || crate::sheaf_disk::Bridge::detect(&name.text).is_some()
                {
                    return self.disk_expression(&name, span, nest);
                }
                if matches!(name.text.as_str(), "shParVec" | "shParRaw" | "shParGen")
                    || crate::sheaf_disk::parabolic::Unary::detect(&name.text).is_some()
                    || crate::sheaf_disk::parabolic::Transport::detect(&name.text).is_some()
                {
                    return self.parabolic_expression(&name, span, nest);
                }
                if matches!(
                    name.text.as_str(),
                    "shConvRaw"
                        | "shConvGen"
                        | "shConvProject"
                        | "shConvSection"
                        | "shConvRawT"
                        | "shConvT"
                        | "shConvMap"
                ) {
                    return self.convolution_expression(&name, span, nest);
                }
                if name.text == "shGlueApex" && self.eat(Token::LBracket)? {
                    let cone = self.label()?;
                    self.expect(Token::Comma)?;
                    let index = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::SheafConeApex(cone, index), 0);
                }
                if matches!(
                    name.text.as_str(),
                    "shGlueLift" | "shGlueTail" | "shConeTransfer" | "shConeUntransfer"
                ) && self.eat(Token::LBracket)?
                {
                    let family = self.label()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    let form = match name.text.as_str() {
                        "shGlueLift" => ExprForm::SheafConeTail(family, true, Box::new(input)),
                        "shGlueTail" => ExprForm::SheafConeTail(family, false, Box::new(input)),
                        "shConeTransfer" => {
                            ExprForm::SheafConeTransfer(family, false, Box::new(input))
                        }
                        _ => ExprForm::SheafConeTransfer(family, true, Box::new(input)),
                    };
                    return expr(span, form, depth);
                }
                if name.text == "shBCApply" && self.eat(Token::LBracket)? {
                    let comparison = self.label()?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::SheafBCApply(comparison, Box::new(input)),
                        depth,
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
                if name.text == "rmbarTerm" && self.eat(Token::LBracket)? {
                    let grade = self.motivic_grade()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RealBarTerm(grade, index, at), 0);
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
                if name.text == "rmrhoextgen" && self.eat(Token::LBracket)? {
                    let grade = self.real_rho_grade()?;
                    self.expect(Token::Comma)?;
                    let (index, at) = self.number()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RealRhoGenerator(grade, index, at), 0);
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
                if name.text == "rvHeckeKL" && self.eat(Token::LBracket)? {
                    let window = self.affine_window()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::HeckeBasis(true, window), 0);
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
                if name.text == "rvTLPath" && self.eat(Token::LBracket)? {
                    let (level, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let path = self.exponents()?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::TlPath(level, path), 0);
                }
                if matches!(name.text.as_str(), "rvTLSurgery" | "rvTLSurgeryTV")
                    && self.eat(Token::LBracket)?
                {
                    let kind = if name.text == "rvTLSurgery" {
                        crate::rv_tl::Surgery::Wrt
                    } else {
                        crate::rv_tl::Surgery::Doubled
                    };
                    let (level, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let (strands, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let word = self.signed_tuple("a signed Artin generator", "Artin generator")?;
                    self.expect(Token::Comma)?;
                    let framings = self.signed_tuple("a signed framing", "framing")?;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::TlSurgery(kind, level, strands, word, framings),
                        0,
                    );
                }
                if name.text == "rvTLFramed" && self.eat(Token::LBracket)? {
                    let (level, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let (strands, _) = self.number()?;
                    self.expect(Token::Comma)?;
                    let word = self.signed_tuple("a signed Artin generator", "Artin generator")?;
                    self.expect(Token::Comma)?;
                    let colors = self.exponents()?;
                    self.expect(Token::Comma)?;
                    let framings = self.signed_tuple("a signed framing", "framing")?;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::TlFramed(level, strands, word, colors, framings),
                        0,
                    );
                }
                if name.text == "rvTLColoredJones" && self.eat(Token::LBracket)? {
                    let (level, strands, color) = self.tl_grade()?;
                    self.expect(Token::Comma)?;
                    let word = self.signed_tuple("a signed Artin generator", "Artin generator")?;
                    self.expect(Token::RBracket)?;
                    return expr(
                        span,
                        ExprForm::TlColoredJones(level, strands, color, word),
                        0,
                    );
                }
                if (matches!(
                    name.text.as_str(),
                    "leeOrientation" | "leeDtRep" | "leeDtNull" | "leeDtDual"
                ) && self.token() == Token::LBracket)
                    || (matches!(name.text.as_str(), "leeDtDualD" | "leeDtPair")
                        && self.token() == Token::LParen)
                {
                    return self.lee_dt_expression(&name.text, span, nest + 1);
                }
                if matches!(name.text.as_str(), "leeRep" | "leeScale" | "leeFilter")
                    && self.token() == Token::LBracket
                    || matches!(
                        name.text.as_str(),
                        "leeLift" | "leePhi" | "leeD" | "leeNull"
                    ) && self.token() == Token::LParen
                {
                    return self.lee_expression(&name.text, span, nest + 1);
                }
                if matches!(
                    name.text.as_str(),
                    "khTerm"
                        | "khBasis"
                        | "khGen"
                        | "khNull"
                        | "khEulerAt"
                        | "khScale"
                        | "khABasis"
                        | "khAGen"
                        | "khANull"
                        | "khAnnular"
                        | "khToF2"
                        | "khPsi"
                        | "khKappaNull"
                        | "khKappaDual"
                        | "khHighest"
                        | "khMod"
                        | "khPageGen"
                        | "khPageClass"
                ) && self.token() == Token::LBracket
                    || matches!(
                        name.text.as_str(),
                        "khD"
                            | "khD0"
                            | "khDMinus2"
                            | "khForgetK"
                            | "khClass"
                            | "khRep"
                            | "khPair"
                            | "khE"
                            | "khF"
                            | "khH"
                            | "khVMinus"
                            | "khVZero"
                            | "khWPlus"
                            | "khLee0"
                            | "khLeePlus"
                            | "khPageRep"
                            | "khPageD"
                            | "khPageNext"
                            | "khPagePrev"
                            | "khPageBoundary"
                    ) && self.token() == Token::LParen
                {
                    return self.kh_expression(&name.text, span, nest + 1);
                }
                if matches!(name.text.as_str(), "rvTLTrace" | "rvTLJones")
                    && self.eat(Token::LBracket)?
                {
                    let form = if name.text == "rvTLTrace" {
                        let (level, strands, height) = self.tl_grade()?;
                        self.expect(Token::Comma)?;
                        let word =
                            self.signed_tuple("a signed Artin generator", "Artin generator")?;
                        ExprForm::TlTrace(level, strands, height, word)
                    } else {
                        let (level, _) = self.number()?;
                        self.expect(Token::Comma)?;
                        let (strands, _) = self.number()?;
                        self.expect(Token::Comma)?;
                        let word =
                            self.signed_tuple("a signed Artin generator", "Artin generator")?;
                        ExprForm::TlJones(level, strands, word)
                    };
                    self.expect(Token::RBracket)?;
                    return expr(span, form, 0);
                }
                if matches!(
                    name.text.as_str(),
                    "rvTLE" | "rvTLBraid" | "rvTLCoord" | "rvTLProject" | "rvCycloQGalois"
                ) && self.eat(Token::LBracket)?
                {
                    let operation = match name.text.as_str() {
                        "rvTLE" => TlOp::Generator(self.number()?.0),
                        "rvTLProject" => {
                            let first = self.number()?.0;
                            self.expect(Token::Comma)?;
                            let length = self.number()?.0;
                            TlOp::Project(first, length)
                        }
                        "rvTLBraid" => TlOp::Braid(
                            self.signed_tuple("a signed Artin generator", "Artin generator")?,
                        ),
                        "rvTLCoord" => {
                            let path = self.exponents()?;
                            self.expect(Token::RBracket)?;
                            let input = self.paren_expr(nest)?;
                            let depth = input.depth + 1;
                            return expr(
                                span,
                                ExprForm::TlCoordinate(path, Box::new(input)),
                                depth,
                            );
                        }
                        _ => {
                            let power =
                                self.signed_wide("a signed Galois power", "Galois power")?;
                            self.expect(Token::RBracket)?;
                            let input = self.paren_expr(nest)?;
                            let depth = input.depth + 1;
                            return expr(
                                span,
                                ExprForm::CycloRatioGalois(power, Box::new(input)),
                                depth,
                            );
                        }
                    };
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(
                        span,
                        ExprForm::TlOperation(operation, Box::new(input)),
                        depth,
                    );
                }
                if matches!(name.text.as_str(), "rvTLScale" | "rvCycloQMul")
                    && self.token() == Token::LParen
                {
                    let arguments = self.paren_pair(nest)?;
                    let depth = 1 + arguments[0].depth.max(arguments[1].depth);
                    let form = if name.text == "rvTLScale" {
                        ExprForm::TlScale(Box::new(arguments))
                    } else {
                        ExprForm::CycloRatioProduct(Box::new(arguments))
                    };
                    return expr(span, form, depth);
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
                    "rvWzwGalois" | "rvWzwCoeffGalois" if self.eat(Token::LBracket)? => {
                        let power = self.signed_wide("a signed Galois power", "Galois power")?;
                        self.expect(Token::RBracket)?;
                        if name.text == "rvWzwGalois" {
                            Some(WzwOp::Galois(power))
                        } else {
                            Some(WzwOp::CoeffGalois(power))
                        }
                    }
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
                if name.text == "rvCycloGalois" && self.eat(Token::LBracket)? {
                    let power = self.signed_wide("a signed Galois power", "Galois power")?;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::CycloGalois(power, Box::new(input)), depth);
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
                if matches!(name.text.as_str(), "rvCTMod" | "rvCTInitial")
                    && self.eat(Token::LBracket)?
                {
                    let period = self.label()?;
                    self.expect(Token::Comma)?;
                    let (prime, exponent) = self.cartier_level()?;
                    let form = if name.text == "rvCTMod" {
                        self.expect(Token::Comma)?;
                        ExprForm::RVCTQuery(
                            period,
                            prime,
                            exponent,
                            self.wide_unsigned("period index")?,
                        )
                    } else {
                        ExprForm::RVCTInitial(period, prime, exponent)
                    };
                    self.expect(Token::RBracket)?;
                    return expr(span, form, 0);
                }
                if name.text == "rvCTResidue" && self.eat(Token::LBracket)? {
                    let (prime, exponent) = self.cartier_level()?;
                    self.expect(Token::Comma)?;
                    let value = self.wide_unsigned("Cartier residue")?;
                    self.expect(Token::RBracket)?;
                    return expr(span, ExprForm::RVCTLiteral(prime, exponent, value), 0);
                }
                if name.text == "rvCTSection" && self.eat(Token::LBracket)? {
                    let digit = self.number()?.0;
                    self.expect(Token::RBracket)?;
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RVCTSection(digit, Box::new(input)), depth);
                }
                if name.text == "rvCTRead" && self.token() == Token::LParen {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::RVCTRead(Box::new(input)), depth);
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
                if matches!(
                    name.text.as_str(),
                    "rvWittSpectrum"
                        | "rvWittSpectrumLine"
                        | "rvWittSpectrumTensor"
                        | "rvWittSpectrumAdams"
                        | "rvWittSpectrumExterior"
                        | "rvWittSpectrumVerschiebung"
                        | "rvWittSpectrumNorm"
                        | "rvWittSpectrumCyclic"
                        | "rvWittSpectrumLie"
                        | "rvK3FrobHilbSpectrum"
                        | "rvK3SpechtSpectrum"
                        | "rvK3StratumSpectrum"
                ) {
                    return self.witt_spectrum_expression(&name.text, span, nest);
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
                if let Some(operation) = match name.text.as_str() {
                    "rvOreEulerFactor" if self.token() == Token::LParen => Some(OreOp::EulerFactor),
                    "rvOreEulerReduce" if self.token() == Token::LParen => Some(OreOp::EulerReduce),
                    "rvOreFactorial" if self.eat(Token::LBracket)? => {
                        let height = self.number()?.0;
                        self.expect(Token::RBracket)?;
                        Some(OreOp::Factorial(height))
                    }
                    _ => None,
                } {
                    let input = self.paren_expr(nest)?;
                    let depth = input.depth + 1;
                    return expr(span, ExprForm::OreUnary(operation, Box::new(input)), depth);
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

    fn wz_numerator(&mut self) -> Result<Vec<(usize, usize, i128)>> {
        self.expect(Token::LBrace)?;
        let mut terms = Vec::new();
        while self.token() != Token::RBrace {
            self.expect(Token::LParen)?;
            let (n, k) = self.bidegree()?;
            self.expect(Token::Comma)?;
            let c = self.signed_wide("a signed coefficient", "WZ numerator coefficient")?;
            self.expect(Token::RParen)?;
            terms.push((n, k, c));
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(terms)
    }

    fn wz_denominator(&mut self) -> Result<Vec<crate::rv_wz::Factor>> {
        self.expect(Token::LBrace)?;
        let mut factors = Vec::new();
        while self.token() != Token::RBrace {
            self.expect(Token::LParen)?;
            let n = self.signed_wide("a signed coefficient", "WZ denominator coefficient")?;
            self.expect(Token::Comma)?;
            let k = self.signed_wide("a signed coefficient", "WZ denominator coefficient")?;
            self.expect(Token::Comma)?;
            let constant =
                self.signed_wide("a signed coefficient", "WZ denominator coefficient")?;
            self.expect(Token::Comma)?;
            let multiplicity = self.number()?.0;
            self.expect(Token::RParen)?;
            factors.push(crate::rv_wz::Factor {
                n,
                k,
                constant,
                multiplicity,
            });
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(factors)
    }

    fn quantum_grade(&mut self) -> Result<(usize, usize)> {
        let (rank, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (ambient, _) = self.number()?;
        Ok((rank, ambient))
    }

    fn tl_grade(&mut self) -> Result<(usize, usize, usize)> {
        let (level, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (strands, _) = self.number()?;
        self.expect(Token::Comma)?;
        let (height, _) = self.number()?;
        Ok((level, strands, height))
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

    fn tropical_weights(&mut self) -> Result<Vec<i64>> {
        self.expect(Token::LBracket)?;
        let mut weights = Vec::new();
        while self.token() != Token::RBracket {
            weights.push(self.signed_number("tropical weight")?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBracket)?;
        Ok(weights)
    }

    fn rational_entry(
        &mut self,
        numerator_label: &str,
        denominator_label: &str,
    ) -> Result<(i128, i128)> {
        let pair = self.eat(Token::LParen)?;
        let numerator = self.signed_wide("an integer", numerator_label)?;
        let denominator = if pair {
            self.expect(Token::Comma)?;
            let value = self.signed_wide("a positive integer", denominator_label)?;
            self.expect(Token::RParen)?;
            value
        } else {
            1
        };
        Ok((numerator, denominator))
    }

    fn residue_rational_rows(&mut self) -> Result<Vec<Vec<(i128, i128)>>> {
        let mut rows = Vec::new();
        self.expect(Token::LBrace)?;
        while self.token() != Token::RBrace {
            self.expect(Token::LBracket)?;
            let mut row = Vec::new();
            while self.token() != Token::RBracket {
                row.push(self.rational_entry("residue row numerator", "residue row denominator")?);
                if !self.eat(Token::Comma)? {
                    break;
                }
            }
            self.expect(Token::RBracket)?;
            rows.push(row);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(rows)
    }

    fn tropical_rows(&mut self, nest: usize) -> Result<(Vec<Vec<ParsedExpr>>, usize)> {
        self.expect(Token::LBrace)?;
        let mut rows = Vec::new();
        let mut depth = 0;
        while self.token() != Token::RBrace {
            self.expect(Token::LBracket)?;
            let mut row = Vec::new();
            while self.token() != Token::RBracket {
                let entry = self.expr(nest + 1)?;
                depth = depth.max(entry.depth + 1);
                row.push(entry);
                if !self.eat(Token::Comma)? {
                    break;
                }
            }
            self.expect(Token::RBracket)?;
            rows.push(row);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok((rows, depth))
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

    fn paren_triple(&mut self, nest: usize) -> Result<[ParsedExpr; 3]> {
        self.expect(Token::LParen)?;
        let a = self.expr(nest + 1)?;
        self.expect(Token::Comma)?;
        let b = self.expr(nest + 1)?;
        self.expect(Token::Comma)?;
        let c = self.expr(nest + 1)?;
        self.expect(Token::RParen)?;
        Ok([a, b, c])
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

    fn real_extension_pair(
        &mut self,
        nest: usize,
        span: Span,
        name: &str,
    ) -> Result<Box<[ParsedExpr; 2]>> {
        self.expect(Token::LParen)?;
        let pair = self
            .arguments(nest)?
            .try_into()
            .map_err(|_| Error::new(span, format!("{name} requires two RMExt classes")))?;
        Ok(Box::new(pair))
    }

    fn mixed_h_generator(&mut self, span: Span) -> Result<ParsedExpr> {
        let name = self.label()?;
        self.expect(Token::Comma)?;
        let b = self.signed_degree()?;
        self.expect(Token::Comma)?;
        let (index, at) = self.number()?;
        self.expect(Token::RBracket)?;
        expr(span, ExprForm::MixedHGenerator(name, b, index, at), 0)
    }

    fn sheaf_deformation_expr(
        &mut self,
        name: &Label,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::sheaf_deform::Operation;
        let (form, depth) = match name.text.as_str() {
            "shDefStart" => {
                self.expect(Token::LBracket)?;
                let sheaf = self.label()?;
                self.expect(Token::RBracket)?;
                (ExprForm::SheafDefStart(sheaf), 0)
            }
            "shDefAppend" => {
                let args = self.paren_pair(nest)?;
                let depth = 1 + args[0].depth.max(args[1].depth);
                (ExprForm::SheafDefAppend(Box::new(args)), depth)
            }
            "shDefGauge" | "shDefGaugeWitness" => {
                self.expect(Token::LBracket)?;
                let exponent = self.number()?.0;
                self.expect(Token::RBracket)?;
                let args = self.paren_pair(nest)?;
                let depth = 1 + args[0].depth.max(args[1].depth);
                (
                    ExprForm::SheafDefGauge(
                        exponent,
                        name.text == "shDefGaugeWitness",
                        Box::new(args),
                    ),
                    depth,
                )
            }
            _ => {
                let operation = match name.text.as_str() {
                    "shDefTerm" => {
                        self.expect(Token::LBracket)?;
                        let order = self.number()?.0;
                        self.expect(Token::RBracket)?;
                        Operation::Term(order)
                    }
                    "shDefCurvature" => Operation::Curvature,
                    "shDefObstruction" => Operation::Obstruction,
                    "shDefLift" => Operation::Lift,
                    _ => unreachable!(),
                };
                let arg = self.paren_expr(nest)?;
                let depth = 1 + arg.depth;
                (ExprForm::SheafDefOperation(operation, Box::new(arg)), depth)
            }
        };
        expr(span, form, depth)
    }

    fn disk_label_index(&mut self) -> Result<(Label, usize)> {
        self.expect(Token::LBracket)?;
        let name = self.label()?;
        self.expect(Token::Comma)?;
        let index = self.number()?.0;
        self.expect(Token::RBracket)?;
        Ok((name, index))
    }

    fn convolution_command(
        &mut self,
        command: &str,
    ) -> Result<crate::sheaf_disk::convolution::Command> {
        use crate::sheaf_disk::convolution::{Command, Query};
        self.advance()?;
        Ok(match command {
            "shconv" => {
                let name = self.label()?;
                self.expect(Token::LBracket)?;
                let source = self.label()?;
                self.expect(Token::Comma)?;
                let lambda = self.signed_number("middle convolution scalar")?;
                self.expect(Token::RBracket)?;
                Command::Convolve(name, source, lambda)
            }
            "shconvmap" => Command::Map(self.label()?, self.disk_labels()?),
            "shconvrank" => {
                let [name] = self.disk_labels()?;
                Command::Rank(name)
            }
            "shconvcoords" => Command::Report(Box::new(self.paren_expr(0)?)),
            _ => {
                let [name] = self.disk_labels()?;
                let query = match command {
                    "shdiskinfinity" => Query::Infinity,
                    "shdiskrigidity" => Query::Rigidity,
                    "shdiskgroup" => Query::GroupOrder,
                    _ => unreachable!(),
                };
                Command::Sequence(name, query)
            }
        })
    }

    fn convolution_expression(
        &mut self,
        name: &Label,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::sheaf_disk::convolution::{Expression, Form};
        let (form, depth) = match name.text.as_str() {
            "shConvRaw" => {
                let [name] = self.disk_labels()?;
                (Expression::Raw(name, self.disk_vector_rows()?), 0)
            }
            "shConvGen" => {
                let (name, index) = self.disk_label_index()?;
                (Expression::Generator(name, index), 0)
            }
            "shConvProject" | "shConvSection" => {
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                let expression = if name.text == "shConvProject" {
                    Expression::Project(Box::new(input))
                } else {
                    Expression::Section(Box::new(input))
                };
                (expression, depth)
            }
            "shConvRawT" | "shConvT" => {
                self.expect(Token::LBracket)?;
                let index = self.signed_number("signed middle convolution generator")?;
                self.expect(Token::RBracket)?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                let form = if name.text == "shConvRawT" {
                    Form::Raw
                } else {
                    Form::Quotient
                };
                (Expression::Action(form, index, Box::new(input)), depth)
            }
            "shConvMap" => {
                let [map] = self.disk_labels()?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                (Expression::Map(map, Box::new(input)), depth)
            }
            _ => unreachable!(),
        };
        expr(span, ExprForm::Convolution(form), depth)
    }

    fn disk_vector_rows(&mut self) -> Result<Vec<Vec<i64>>> {
        self.expect(Token::LBrace)?;
        let mut rows = Vec::new();
        while self.token() != Token::RBrace {
            rows.push(self.sheaf_matrix()?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(rows)
    }

    fn parabolic_command(
        &mut self,
        command: &str,
    ) -> Result<crate::sheaf_disk::parabolic::Command> {
        use crate::sheaf_disk::parabolic::Command;
        self.advance()?;
        Ok(match command {
            "shdiskseq" => Command::Sequence(self.label()?, self.sheaf_chain()?),
            "shhurwitz" => {
                let name = self.label()?;
                let [source] = self.disk_labels()?;
                Command::Hurwitz(name, source, self.sheaf_matrix()?)
            }
            "shdiskpoint" => {
                let name = self.label()?;
                let (source, index) = self.disk_label_index()?;
                Command::Point(name, source, index)
            }
            "shdiskprod" => {
                let [source] = self.disk_labels()?;
                Command::Product(source)
            }
            "shdiskseqsame" | "shdiskmonosame" => {
                Command::Same(self.disk_labels()?, command == "shdiskseqsame")
            }
            "shdiskseqmap" => Command::Map(
                self.label()?,
                self.disk_labels()?,
                self.sheaf_cone_columns()?,
            ),
            "shdiskseqmapbraid" => Command::BraidMap(self.label()?, self.disk_labels()?),
            "shdiskseqmapat" => {
                let (map, index) = self.disk_label_index()?;
                Command::MapAt(map, index)
            }
            "shparrank" => {
                let [sequence] = self.disk_labels()?;
                Command::Rank(sequence)
            }
            "shparrows" | "shparvector" => {
                Command::Report(Box::new(self.paren_expr(0)?), command == "shparrows")
            }
            _ => unreachable!(),
        })
    }

    fn parabolic_expression(
        &mut self,
        name: &Label,
        span: Span,
        nest: usize,
    ) -> Result<ParsedExpr> {
        use crate::sheaf_disk::parabolic::{Expression, Transport, Unary};
        let (form, depth) = match name.text.as_str() {
            "shParVec" => {
                let [sequence] = self.disk_labels()?;
                (Expression::Nearby(sequence, self.sheaf_matrix()?), 0)
            }
            "shParRaw" => {
                let [sequence] = self.disk_labels()?;
                (Expression::Raw(sequence, self.disk_vector_rows()?), 0)
            }
            "shParGen" => {
                let (sequence, index) = self.disk_label_index()?;
                (Expression::Generator(sequence, index), 0)
            }
            operation => {
                let transport = Transport::detect(operation)
                    .map(|kind| self.disk_labels().map(|[name]| (kind, name)))
                    .transpose()?;
                let input = self.paren_expr(nest)?;
                let depth = 1 + input.depth;
                let form = if let Some((kind, name)) = transport {
                    Expression::Transport(kind, name, Box::new(input))
                } else {
                    Expression::Unary(
                        Unary::detect(operation).expect("known parabolic operator"),
                        Box::new(input),
                    )
                };
                (form, depth)
            }
        };
        expr(span, ExprForm::Parabolic(form), depth)
    }

    fn disk_command(&mut self, command: &str) -> Result<crate::sheaf_disk::Command> {
        use crate::sheaf_disk::{Command, Extension, Operator};
        self.advance()?;
        if command == "shdiskmapmat" {
            return Ok(Command::MapMatrix(Box::new(self.paren_expr(0)?)));
        }
        if command == "shdiskhomdim" {
            return Ok(Command::HomDimension(self.disk_labels()?));
        }
        let name = if matches!(command, "shdiskmat" | "shdiskranks") {
            None
        } else {
            Some(self.label()?)
        };
        Ok(match command {
            "shdisk" => {
                self.expect(Token::LBracket)?;
                let p = self.number()?.0;
                self.expect(Token::Comma)?;
                let n = self.number()?.0;
                self.expect(Token::Comma)?;
                let m = self.number()?.0;
                self.expect(Token::RBracket)?;
                Command::Declare(
                    name.unwrap(),
                    [p, n, m],
                    [self.sheaf_matrix()?, self.sheaf_matrix()?],
                )
            }
            "shdiskmap" => Command::Map(
                name.unwrap(),
                self.disk_labels()?,
                [self.sheaf_matrix()?, self.sheaf_matrix()?],
            ),
            "shdiskmat" => {
                let [disk, operator] = self.disk_labels()?;
                Command::Matrix(disk, Operator::from_label(&operator)?)
            }
            "shdiskranks" => {
                let [disk] = self.disk_labels()?;
                Command::Ranks(disk)
            }
            "shdiskrealize" => Command::Realize(name.unwrap(), Box::new(self.disk_labels()?)),
            "shdiskoperator" => {
                let [disk, operator, source, target] = self.disk_labels()?;
                Command::Operator(
                    name.unwrap(),
                    disk,
                    Operator::from_label(&operator)?,
                    [source, target],
                )
            }
            _ => {
                let [disk] = self.disk_labels()?;
                let name = name.unwrap();
                match command {
                    "shdiskdual" => Command::Dual(name, disk),
                    "shdiskstalk" | "shdiskcostalk" => {
                        Command::Complex(name, disk, command == "shdiskcostalk")
                    }
                    _ => {
                        let extension = match command {
                            "shdiskshriek" => Extension::Shriek,
                            "shdiskstar" => Extension::Star,
                            "shdiskmiddle" => Extension::Middle,
                            _ => unreachable!(),
                        };
                        Command::Extend(name, disk, extension)
                    }
                }
            }
        })
    }

    fn disk_labels<const N: usize>(&mut self) -> Result<[Label; N]> {
        self.expect(Token::LBracket)?;
        let mut names = Vec::with_capacity(N);
        for i in 0..N {
            if i != 0 {
                self.expect(Token::Comma)?;
            }
            names.push(self.label()?);
        }
        self.expect(Token::RBracket)?;
        Ok(names.try_into().ok().expect("fixed label count"))
    }

    fn disk_expression(&mut self, name: &Label, span: Span, nest: usize) -> Result<ParsedExpr> {
        use crate::sheaf_disk::{Bridge, Expression};
        let (form, depth) = match name.text.as_str() {
            "shDiskId" => {
                let [disk] = self.disk_labels()?;
                (Expression::Id(disk), 0)
            }
            "shDiskGen" => {
                self.expect(Token::LBracket)?;
                let from = self.label()?;
                self.expect(Token::Comma)?;
                let to = self.label()?;
                self.expect(Token::Comma)?;
                let index = self.number()?.0;
                self.expect(Token::RBracket)?;
                (Expression::Generator([from, to], index), 0)
            }
            "shDiskComp" => {
                let pair = self.paren_pair(nest)?;
                let depth = 1 + pair[0].depth.max(pair[1].depth);
                (Expression::Compose(Box::new(pair)), depth)
            }
            "shDiskCostalkLift" => {
                let [name] = self.disk_labels()?;
                let input = self.paren_expr(nest)?;
                let depth = input.depth + 1;
                (Expression::Lift(name, Box::new(input)), depth)
            }
            operation => {
                let bridge = Bridge::detect(operation).expect("known disk operator");
                let names = self.disk_labels()?;
                let input = if bridge.has_input() {
                    Some(Box::new(self.paren_expr(nest)?))
                } else {
                    None
                };
                let depth = input.as_ref().map_or(0, |input| input.depth + 1);
                (Expression::Bridge(bridge, names, input), depth)
            }
        };
        expr(span, ExprForm::Disk(form), depth)
    }

    fn sheaf_cone_columns(&mut self) -> Result<Vec<ParsedExpr>> {
        self.expect(Token::LBrace)?;
        let mut values = Vec::new();
        while self.token() != Token::RBrace {
            values.push(self.expr(0)?);
            if !self.eat(Token::Comma)? {
                break;
            }
        }
        self.expect(Token::RBrace)?;
        Ok(values)
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
