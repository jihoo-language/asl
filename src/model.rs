use std::cmp::Ordering;
use std::collections::HashMap;
use std::rc::Rc;

use crate::adams::Adams;
use crate::binary_field::BinaryField;
use crate::binary_linear;
use crate::bp_classes;
use crate::bp_height;
use crate::bp_mod_ext;
use crate::bp_page;
use crate::bp_rel;
use crate::bp_triple;
use crate::bp_vpage;
use crate::ce::CePage;
use crate::cobar;
use crate::cohomology::Cohomology;
use crate::comodule;
use crate::dependent_sets;
use crate::ec_morphism::Morphism as ECMap;
use crate::ec_poly::Scalar as ECScalar;
use crate::elliptic::{
    self, Curve as EllipticCurve, Isogeny, Kernel as ECKernel, Point as ECPoint,
};
use crate::elliptic_end_order::level as ec_level;
use crate::error::{Error, Result};
use crate::hw_algebra::Matrix as HWMatrix;
use crate::hw_closed;
use crate::hyperelliptic::{self, Residues};
use crate::hyperelliptic_chart::Chart;
use crate::hyperelliptic_jacobian::{self, Divisor};
use crate::hyperelliptic_kummer;
use crate::hyperelliptic_kummer_surface::enriques::points as enriques;
use crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::eigenpacket;
use crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::hecke;
use crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::period_polynomial as period;
use crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::residual_hecke;
use crate::hyperelliptic_kummer_surface::hilbert::elliptic_genus::{farey, inverse, maass};
use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::exterior as ks_exterior;
use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::exterior::nomizu;
use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::projector as ks_projector;
use crate::hyperelliptic_kummer_surface::qz;
use crate::hyperelliptic_kummer_surface::resolution::ResolvedPoint as KummerK3Point;
use crate::hyperelliptic_kummer_surface::{
    self, hilbert_square, mordell, nakajima, nikulin, pencil, unimodular, Point as KummerPoint,
};
use crate::hyperelliptic_kummer_surface::{clifford, cm_ring, rational_sparse, spinor};
use crate::jet_hom;
use crate::jet_page;
use crate::khovanov::{self, Form as KhForm};
use crate::legendre_cm;
use crate::map_cochain;
use crate::map_ext::{self, MapAdams};
use crate::map_may::MapMayPage;
use crate::may::MayPage;
use crate::milnor::MilnorBasis;
use crate::mixed_bars;
use crate::mixed_cone;
use crate::mixed_hom;
use crate::mixed_homology;
use crate::mixed_minimal;
use crate::mixed_tensor;
use crate::motivic;
use crate::odd_cobar::{self, Prime};
use crate::odd_ext::{self, OddExt};
use crate::odd_may::{Filtration, OddPage};
use crate::page_transition::Transition;
use crate::prime_linear::Coordinate;
use crate::profile::Profile;
use crate::q_koszul::{self, Grade, QExt};
use crate::q_module::{Module, Presentation};
use crate::rank_jets::Path;
use crate::real_bars;
use crate::real_ext;
use crate::real_extensions;
use crate::real_massey;
use crate::real_motivic;
use crate::real_rho;
use crate::rho_bock;
use crate::rho_localized;
use crate::rho_minimal;
use crate::rv_as_cover::Cover as RvAsCover;
use crate::rv_as_odd::Cover as RvAsOddCover;
use crate::rv_as_two::Cover as RvAsTwoCover;
use crate::rv_bar;
use crate::rv_bar_dual;
use crate::rv_bar_page;
use crate::rv_barcobar;
use crate::rv_crossed::{self, postnikov};
use crate::rv_crystal::{Crystal, Entries as CrystalEntries};
use crate::rv_curve::Curve as RvCurve;
use crate::rv_cyclic_homology;
use crate::rv_cyclotomic;
use crate::rv_fusion;
use crate::rv_hh;
use crate::rv_hh_homology;
use crate::rv_hh_mixed;
use crate::rv_hive::{self, Boundary as HiveBoundary};
use crate::rv_hom::{self, Hom};
use crate::rv_ideal::Ideal as RvIdeal;
use crate::rv_iwhecke;
use crate::rv_koszul::{self, Grade as RvGrade, Tor};
use crate::rv_laurent;
use crate::rv_lens;
use crate::rv_negative;
use crate::rv_ore;
use crate::rv_pcurv;
use crate::rv_period::{self, Grade as PeriodGrade, Level as PeriodLevel, Period};
use crate::rv_phodge::{Module as PHodgeModule, QuadraticSpec, Specification as PHodgeSpec};
use crate::rv_poly::{self, Minors, Polynomial};
use crate::rv_qgr;
use crate::rv_qp_quad::{Element as QuadElement, Field as QuadField};
use crate::rv_res_space::{self, Grade as ResGrade, Space as ResSpace};
use crate::rv_satake;
use crate::rv_sphhecke;
use crate::rv_tl;
use crate::rv_unram::{self, Ring as UnramRing};
use crate::rv_witt_cover::Cover as RvWittCover;
use crate::rv_witt_tower::Tower as RvWittTower;
use crate::rv_wzw;
use crate::set_family;
use crate::set_j;
use crate::set_logic;
use crate::set_topology;
use crate::sheaf::{self, Cochains, Form as SheafForm, Sheaf};
use crate::sheaf_cellular::{Cellular, Resolution, Total};
use crate::sheaf_deform;
use crate::sheaf_disk::{
    braid::{Sequence as DiskSequence, SequenceMap as DiskSequenceMap},
    convolution::{Convolution, ConvolutionMap, Form as ConvolutionForm},
    parabolic::{Form as ParabolicForm, Parabolic},
};
use crate::sheaf_ext::{Ends as SheafEnds, Hom as SheafHom};
use crate::sheaf_hyper;
use crate::sheaf_kuranishi;
use crate::sheaf_leray::Leray;
use crate::sheaf_page;
use crate::sheaf_support::{self, Supports};
use crate::stable_marks::{self, Reduction as StableCurve, Tree as StableTree};
use crate::steenrod::{self, Word};
use crate::syntax::Label;
use crate::tropical_matroid::{self, Grade as ValGrade, Matroid as ValMatroid};
use crate::v_bars;
use crate::v_koszul;
use crate::vbock::{self, BocksteinPage};
use crate::{bp, bp_polynomial};
use crate::{mixed, mixed_page};

const MAX_SPACE_POINTS: usize = 512;
const MAX_TYPES: usize = 4096;
const MAX_TYPE_DEPTH: usize = 64;
const MAX_ELEMENTS: usize = 4096;
const MAX_WORK: usize = 25_000_000;
const MAX_STORAGE: usize = 1_000_000;
pub(crate) const MAX_SHOWN_BYTES: usize = 65_536;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) struct TypeId(usize);

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum DiffFamily {
    ArtinSchreier,
    Witt,
}

impl DiffFamily {
    pub(crate) fn type_name(self) -> &'static str {
        match self {
            Self::ArtinSchreier => "ASDiff",
            Self::Witt => "WittDiff",
        }
    }

    pub(crate) fn operand(self) -> &'static str {
        match self {
            Self::ArtinSchreier => "an ASDiff differential",
            Self::Witt => "a WittDiff differential",
        }
    }

    pub(crate) fn scale_name(self) -> &'static str {
        match self {
            Self::ArtinSchreier => "rvASDiffScale",
            Self::Witt => "rvWittDiffScale",
        }
    }

    pub(crate) fn cartier_name(self) -> &'static str {
        match self {
            Self::ArtinSchreier => "rvASCartier",
            Self::Witt => "rvWittCartier",
        }
    }
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub(crate) enum Kind {
    Space(String),
    FiniteGroup(Rc<str>),
    PostnikovCochain(Rc<str>, usize),
    PostnikovCell(Rc<str>),
    SetNat(Rc<str>, Rc<str>),
    QModule(String),
    MixedHomology(mixed::Grade),
    MixedChain(mixed::Grade, bool),
    MixedMinimum(mixed_minimal::Grade, bool),
    MixedPage(mixed_page::Grade),
    Sheaf(SheafForm, Rc<str>, usize),
    SheafComplex(SheafForm, Rc<str>, usize),
    SheafDeformation(Rc<str>, usize),
    Convolution(ConvolutionForm, Rc<str>, Prime),
    Parabolic(ParabolicForm, Rc<str>, Prime),
    DiskHom(Rc<str>, Rc<str>, Prime),
    SheafHom(SheafForm, SheafEnds, usize),
    SheafPage(Rc<str>, sheaf_page::Grade),
    SheafLocal(Rc<str>, sheaf_support::Grade),
    Star(TypeId, usize),
    Product(TypeId, TypeId),
    Arrow(TypeId, TypeId),
    Cohomology(TypeId, usize),
    Operations(usize),
    Cobar(usize, usize),
    MCobar(motivic::Grade),
    RMCobar(motivic::Grade),
    RMRhoCobar(real_rho::Grade),
    RMRhoExt(real_rho::Grade),
    RhoBock(rho_bock::Grade),
    RhoLoc(motivic::Grade),
    RhoMin(motivic::Grade),
    MExt(motivic::Grade),
    RMExt(motivic::Grade),
    RMBar(motivic::Grade),
    RMTriple(motivic::Grade),
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
    CurveDiff(DiffFamily, Rc<str>, BinaryField),
    ASOddDiff(Rc<str>, Prime),
    ASTwoDiff(Rc<str>, Prime),
    ASOddDeRham(Rc<str>, Prime),
    ASOddTate(Rc<str>, Prime, i64),
    PrimeScalar(Prime),
    ASDeRham(Rc<str>, BinaryField),
    BinaryScalar(BinaryField),
    EOType(usize),
    RVLaurent(rv_hom::Grade),
    RVIwHecke(usize),
    RVSphHecke(usize),
    RVSatake(usize),
    RVHive(HiveBoundary),
    RVQGr(rv_qgr::Grade),
    RVQGrP(rv_pcurv::ModGrade),
    RVOre,
    RVWittSpectrum,
    RVSchurCharacter(u128),
    RVStableSchurCharacter,
    RVPartitionAlgebra(usize),
    RVPartitionLocal(usize),
    RVKaroubiMap(Rc<str>, Rc<str>),
    RVCTState(PeriodGrade),
    RVCTResidue(PeriodLevel),
    RVFusion(rv_fusion::Grade),
    RVCyclo(usize),
    ECCurve(Rc<str>),
    Jacobian(Rc<str>),
    Kummer(Rc<str>),
    KummerK3(Rc<str>),
    Enriques(Rc<str>),
    K3PencilBase(Rc<str>),
    K3PencilSection(Rc<str>),
    K3PencilKummer2(Rc<str>),
    K3PencilCoverPoint(Rc<str>),
    KummerDisc(Rc<str>),
    KummerWeil(Rc<str>),
    KummerDiv(Rc<str>),
    K3Lattice(Rc<str>),
    Nakajima(nakajima::Grade),
    HilbertCohomology(nakajima::Grade),
    HodgeTensor(nakajima::Grade),
    Verbitsky(nakajima::Grade),
    K3Clifford(Rc<str>, u8),
    OG6Clifford(Rc<str>, u8),
    KSExterior(Rc<str>, i128),
    KSNomizu(Rc<str>, i128, i128),
    KSEnd(Rc<str>),
    QZ,
    K3Spinor(Rc<str>, u8),
    Mukai(Rc<str>),
    OG10(Rc<str>),
    AbMukai(Rc<str>),
    OG6(Rc<str>),
    MukaiQuotient(Rc<str>),
    HilbertSquare(Rc<str>, i128),
    SiegelWeight,
    ModularHecke(u128),
    HeckeEigen(Box<str>),
    HeckeSymbol(Box<str>),
    MazurTate(period::mazur_tate::Grade),
    MazurAmice(period::amice::Grade),
    MazurOrdinary(Box<str>, u64),
    PeriodPolynomial(period::Grade),
    ECMap(Rc<str>, Rc<str>),
    ECPoint(Rc<str>),
    ECLevel(Rc<str>, bool),
    ECKernel(Rc<str>, u64),
    ECIso(Rc<str>),
    HWMatrix(Residues, usize),
    HWWMatrix(hw_closed::Grade),
    HWScalar(Residues),
    RVUnram(Rc<str>),
    RVCycloQ(usize),
    RVQpQuad(Rc<str>),
    RVResSpace(ResGrade),
    RVValMat(ValGrade),
    RVHyperArithmetic(crate::hyperelliptic_graph::Grade),
    RVHyperCover(crate::hyperelliptic_graph::Grade),
    RVStableTree(usize),
    RVStableMarks(u32, usize),
    Lee(Rc<str>, crate::lee::Grade),
    LeeDual(Rc<str>, i64),
    Khovanov(Rc<str>, khovanov::Grade, KhForm),
    CMOrder(legendre_cm::Grade),
    ECEndOrder(Rc<str>),
    CMClass(Rc<str>),
    RvTl(rv_tl::Grade),
    RVWzw(rv_wzw::Grade),
    JetHom(jet_hom::Grade),
    JetPage(jet_page::Grade),
    VKoszul(v_koszul::Grade),
    VExt(v_koszul::Grade),
    VBock(vbock::Grade),
    RVBarE(rv_bar_page::Grade),
    RvCyclicChain(rv_cyclic_homology::Grade),
    RvCyclicHomology(rv_cyclic_homology::Grade),
    RVNCChain(rv_negative::Grade),
    RvNegativeCyclic(rv_negative::Grade),
    RvHochschild(rv_hh::Grade),
    RVHHChain(rv_hh::Grade),
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
    BPVPage(bp_vpage::Grade),
    BPModExt(bp_mod_ext::Grade),
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
    pub(crate) fn real_grade(&self) -> Option<motivic::Grade> {
        match *self {
            Self::RMCobar(grade) | Self::RMExt(grade) => Some(grade),
            _ => None,
        }
    }

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
    GroupElement(usize),
    PostnikovCochain(Rc<postnikov::Cochain>),
    PostnikovCell(postnikov::Cell),
    Pair(Rc<(Value, Value)>),
    Map(Rc<[Value]>),
    SetNat(Rc<[usize]>),
    Class(Rc<[usize]>),
    Operation(Rc<[Word]>),
    Cobar(Rc<[cobar::Word]>),
    RMCobar(real_motivic::Cochain),
    RMTriple(Rc<real_massey::Coset>),
    RhoLoc(rho_localized::Cochain),
    RhoMin(rho_minimal::Support),
    PCobar(odd_cobar::Support),
    OddClass(Rc<[Coordinate]>),
    PrimeScalar(u32),
    SheafDeformation(Rc<sheaf_deform::Jet>),
    SheafC(sheaf::Support),
    MapCobar(map_cochain::Support),
    QKoszul(q_koszul::Support),
    RVHom(rv_hom::Support),
    CurveDiff(Rc<[binary_linear::Entry]>),
    ASDeRham(Rc<[binary_linear::Entry]>),
    BinaryScalar(usize),
    EOType(Rc<[usize]>),
    RVLaurent(rv_laurent::Support),
    RVIwHecke(rv_iwhecke::Support),
    RVSphHecke(rv_sphhecke::Support),
    RVSatake(rv_satake::Support),
    RVHive(rv_hive::Interior),
    RVQGr(rv_qgr::Support),
    RVQGrP(rv_pcurv::Support),
    RVOre(rv_ore::Support),
    RVWittSpectrum(Rc<crate::rv_witt_spectrum::Spectrum>),
    RVSchurCharacter(Rc<crate::rv_witt_spectrum::character::Character>),
    RVStableSchurCharacter(Rc<crate::rv_witt_spectrum::character::StableCharacter>),
    RVPartitionAlgebra(Rc<crate::partition_algebra::Algebra>),
    RVPartitionLocal(Rc<crate::partition_algebra::local::LocalAlgebra>),
    RVKaroubiMap(Rc<crate::partition_algebra::karoubi::Map>),
    PeriodState(rv_period::State),
    PeriodResidue(u128),
    RVFusion(rv_fusion::Support),
    RVCyclo(rv_cyclotomic::Polynomial),
    ECCurve(Rc<str>),
    Jacobian(Rc<Divisor>),
    Kummer(Rc<KummerPoint>),
    KummerK3(Rc<KummerK3Point>),
    Enriques(Rc<enriques::Orbit>),
    K3PencilBase(Rc<pencil::Parameter>),
    K3PencilSection(Rc<pencil::sections::Section>),
    K3PencilKummer2(Rc<pencil::sections::Descent>),
    K3PencilCoverPoint(Rc<pencil::sections::covering::Point>),
    KummerDisc(u8),
    KummerWeil(Rc<crate::hyperelliptic_kummer_surface::weil::Vector>),
    KummerDiv(Rc<crate::hyperelliptic_kummer_surface::polarized::Divisor>),
    K3Lattice(Rc<unimodular::Vector>),
    Nakajima(Rc<nakajima::State>),
    HilbertCohomology(Rc<nakajima::cohomology::Class>),
    HodgeTensor(Rc<hyperelliptic_kummer_surface::hodge_wick::Tensor>),
    K3Clifford(Rc<clifford::Element>),
    OG6Clifford(Rc<crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::spin::Element>),
    KSExterior(Rc<ks_exterior::Form>),
    KSNomizu(Rc<nomizu::Cochain>),
    KSEnd(Rc<ks_projector::End>),
    QZ(qz::Scalar),
    K3Spinor(Rc<spinor::Element>),
    Mukai(Rc<crate::hyperelliptic_kummer_surface::mukai::Vector>),
    OG10(Rc<crate::hyperelliptic_kummer_surface::mukai::moduli::og10::Class>),
    AbMukai(Rc<crate::hyperelliptic_kummer_surface::mukai::abelian::Vector>),
    OG6(Rc<crate::hyperelliptic_kummer_surface::mukai::abelian::og6::Class>),
    MukaiQuotient(Rc<crate::hyperelliptic_kummer_surface::mukai::Vector>),
    HilbertSquare(Rc<hilbert_square::Class>),
    SiegelWeight(crate::siegel::Weight),
    ModularHecke(Rc<hecke::Operator>),
    HeckeEigen(Rc<eigenpacket::Element>),
    HeckeSymbol(Rc<period::modular_symbol::Symbol>),
    MazurTate(Rc<period::mazur_tate::Distribution>),
    MazurAmice(Rc<period::amice::Amice>),
    MazurOrdinary(u64),
    PeriodPolynomial(Rc<period::Polynomial>),
    ECMap(Rc<ECMap>),
    ECPoint(Rc<ECPoint>),
    ECLevelVec(ec_level::Vector),
    ECLevelMat(ec_level::Matrix),
    ECKernel(Rc<ECKernel>),
    ECIso(Rc<str>),
    HWMatrix(Rc<HWMatrix<u64>>),
    HWWMatrix(Rc<hw_closed::Matrix>),
    HWScalar(u64),
    RVUnram(rv_unram::Coordinates),
    RVCycloQ(rv_lens::Rational),
    RVQpQuad(QuadElement),
    RVResSpace(Rc<ResSpace>),
    RVValMat(Rc<ValMatroid>),
    RVHyperArithmetic(Rc<crate::hyperelliptic_graph::arith::Arithmetic>),
    RVHyperCover(Rc<crate::hyperelliptic_graph::Cover>),
    RVStableTree(Rc<StableTree>),
    RVStableMarks(Rc<StableCurve>),
    LeeChain(crate::lee::Chain),
    KhChain(khovanov::Chain),
    KhClass(khovanov::Class),
    CMOrder(legendre_cm::Element),
    ECEndOrder(legendre_cm::Element),
    CMClass(cm_ring::Class),
    RvTl(rv_tl::Support),
    RVWzw(rv_wzw::Support),
    JetHom(jet_hom::Support),
    RVNCChain(rv_negative::Support),
    RVHHChain(rv_hh::Support),
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
    pub(crate) fn strict_less(&self, left: usize, right: usize) -> bool {
        left != right && self.order[left][right]
    }

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
    bp_vpage: Option<Rc<bp_vpage::Page>>,
    bp_mod_ext: Option<Rc<bp_mod_ext::Group>>,
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
    hive_witness: Option<rv_hive::Interior>,
    cardinality_bound: usize,
    heights: [Option<usize>; 2],
    inhabited: bool,
    depth: usize,
}

type MinorsKey = (TypeId, q_koszul::Indices, usize);

pub(crate) struct Model {
    types: Vec<TypeInfo>,
    finite_groups: HashMap<Rc<str>, Rc<rv_crossed::Group>>,
    crossed_modules: HashMap<Rc<str>, Rc<rv_crossed::Crossed>>,
    postnikov_systems: HashMap<Rc<str>, Rc<postnikov::Postnikov>>,
    set_families: HashMap<Rc<str>, Rc<set_family::Diagram>>,
    set_origins: HashMap<Rc<str>, dependent_sets::Origin>,
    set_recollements: HashMap<Rc<str>, Rc<crate::set_gluing::Recollement>>,
    set_monodromies: HashMap<Rc<str>, Rc<crate::set_monodromy::Monodromy>>,
    set_linearizations: HashMap<Rc<str>, Rc<str>>,
    set_logic: HashMap<Rc<str>, set_logic::Origin>,
    set_topologies: HashMap<Rc<str>, Rc<set_topology::Topology>>,
    set_reflections: HashMap<Rc<str>, Rc<set_topology::Reflection>>,
    set_reflection_certificates: set_topology::Certificates,
    set_sheafifications: HashMap<Rc<str>, Rc<set_j::Sheafification>>,
    set_sheaf_certificates: HashMap<Rc<str>, Option<Rc<set_j::Certificate>>>,
    set_elements: HashMap<TypeId, Rc<set_family::Elements>>,
    mixed_complexes: Vec<(String, Rc<mixed::Complex>)>,
    mixed_hh_origins: HashMap<mixed::Id, Rc<rv_hh_mixed::Origin>>,
    mixed_indices: HashMap<String, mixed::Id>,
    mixed_tensors: HashMap<mixed::Id, Rc<mixed_tensor::Tensor>>,
    mixed_cones: HashMap<mixed::Id, Rc<mixed_cone::Cone>>,
    mixed_homs: HashMap<mixed::Id, Rc<mixed_hom::Hom>>,
    mixed_minimals: Vec<(String, Rc<mixed_minimal::Minimal>)>,
    mixed_minimal_ids: HashMap<String, mixed_minimal::Id>,
    mixed_bars: HashMap<mixed::Id, Rc<mixed_bars::Bars>>,
    mixed_homology: HashMap<mixed::Grade, Rc<mixed_homology::Slice>>,
    mixed_pages: HashMap<TypeId, Rc<mixed_page::Page>>,
    sheaves: HashMap<String, Rc<Sheaf>>,
    sheaf_sites: HashMap<TypeId, Rc<sheaf::PosetData>>,
    sheaf_maps: HashMap<String, Rc<sheaf_hyper::Map>>,
    convolutions: HashMap<Rc<str>, Rc<Convolution>>,
    convolution_maps: HashMap<Rc<str>, Rc<ConvolutionMap>>,
    disk_sequences: HashMap<Rc<str>, Rc<DiskSequence>>,
    disk_sequence_maps: HashMap<Rc<str>, Rc<DiskSequenceMap>>,
    parabolic_groups: HashMap<Rc<str>, Rc<Parabolic>>,
    disk_quivers: HashMap<Rc<str>, Rc<crate::sheaf_disk::Disk>>,
    disk_homs: HashMap<(Rc<str>, Rc<str>), Rc<crate::sheaf_disk::Hom>>,
    disk_realizations: HashMap<Rc<str>, Rc<crate::sheaf_disk::Realization>>,
    cone_gluings: HashMap<String, Rc<crate::sheaf_cone::Cone>>,
    cone_homotopies: HashMap<String, Rc<crate::sheaf_cone::Homotopy>>,
    bc_squares: HashMap<String, Rc<crate::sheaf_base_change::Square>>,
    bc_comparisons: HashMap<String, Rc<crate::sheaf_base_change::Comparison>>,
    sheaf_kuranishi: HashMap<Rc<str>, Rc<sheaf_kuranishi::Chart>>,
    sheaf_homs: HashMap<SheafEnds, Rc<SheafHom>>,
    sheaf_complexes: HashMap<String, Total>,
    sheaf_resolutions: HashMap<String, Rc<Resolution>>,
    sheaf_supports: HashMap<String, Rc<Supports>>,
    leray_complexes: HashMap<String, Rc<Leray>>,
    serre_shifts: HashMap<String, usize>,
    type_ids: HashMap<Kind, TypeId>,
    milnor_bases: HashMap<usize, Rc<MilnorBasis>>,
    cyclotomics: HashMap<usize, rv_cyclotomic::Polynomial>,
    v_bars: HashMap<v_bars::Family, Rc<v_bars::Decomposition>>,
    kh_braids: HashMap<Rc<str>, Rc<khovanov::Braid>>,
    rv_ideals: HashMap<Rc<str>, Rc<RvIdeal>>,
    rv_curves: HashMap<String, Rc<RvCurve>>,
    rv_periods: HashMap<String, Rc<Period>>,
    rv_as_covers: HashMap<String, Rc<RvAsCover>>,
    rv_phodge_modules: HashMap<String, Rc<PHodgeModule>>,
    rv_qp_quadratics: HashMap<String, QuadField>,
    legendre_cm_families: HashMap<String, legendre_cm::Family>,
    rv_crystals: HashMap<String, Rc<Crystal>>,
    rv_as_odd_covers: HashMap<String, Rc<RvAsOddCover>>,
    rv_as_two_covers: HashMap<String, Rc<RvAsTwoCover>>,
    rv_witt_covers: HashMap<String, Rc<RvWittCover>>,
    elliptic_curves: HashMap<Rc<str>, Rc<EllipticCurve>>,
    elliptic_end_orders: HashMap<String, Rc<crate::elliptic_end_order::Certificate>>,
    ec_cm_cycles: HashMap<String, Rc<crate::elliptic_end_order::cycle::Cycle>>,
    ec_cm_actions: HashMap<String, Rc<crate::elliptic_end_order::realization::Realization>>,
    ec_torsion_levels: HashMap<String, Rc<ec_level::Level>>,
    elliptic_isogenies: HashMap<Rc<str>, Rc<Isogeny>>,
    elliptic_modular_polynomials: HashMap<u32, Rc<crate::elliptic_modular::Polynomial>>,
    hyperelliptic: HashMap<String, Rc<hyperelliptic::Curve>>,
    hw_closed_curves: HashMap<String, Rc<hw_closed::Curve>>,
    jacobian_charts: HashMap<String, Rc<Chart>>,
    jacobian_levels: HashMap<String, Rc<crate::jacobian_level::Basis>>,
    jacobian_power_levels: HashMap<String, Rc<crate::jacobian_level::Basis>>,
    enriques_quotients: HashMap<String, Rc<enriques::Quotient>>,
    k3_pencils: HashMap<String, Rc<pencil::Pencil>>,
    k3_cassels: HashMap<String, Rc<pencil::sections::covering::cassels::Pairing>>,
    k3_selmer_four: HashMap<String, Rc<pencil::sections::covering::four::Certificate>>,
    k3_pencil_covers: HashMap<String, Rc<pencil::sections::covering::Cover>>,
    k3_pencil_mw_images: HashMap<String, Rc<pencil::sections::covering::MordellImage>>,
    k3_mordell: HashMap<String, Rc<mordell::Mordell>>,
    k3_automorphisms: HashMap<String, Rc<hyperelliptic_kummer_surface::dynamics::Word>>,
    k3_periods: HashMap<String, Rc<hyperelliptic_kummer_surface::period::Period>>,
    cm_ring_orders: HashMap<String, cm_ring::Order>,
    karoubi_objects: HashMap<String, Rc<crate::partition_algebra::local::LocalAlgebra>>,
    deligne_square: Option<Rc<crate::partition_algebra::square::Square>>,
    deligne_towers: HashMap<usize, Rc<crate::partition_algebra::vacillating::Tower>>,
    hilbert_strata: HashMap<String, Rc<hyperelliptic_kummer_surface::stratum::Partition>>,
    specht_strata: HashMap<String, Rc<hyperelliptic_kummer_surface::specht::Decoration>>,
    cm_class_polynomials: HashMap<String, Rc<crate::cm_class_poly::Polynomial>>,
    cm_class_discriminants: HashMap<i128, Rc<crate::cm_class_poly::Polynomial>>,
    nakajima_gram: Option<Rc<nakajima::Gram>>,
    nakajima_inverse_gram: Option<Rc<nakajima::Gram>>,
    og6_nilpotents: HashMap<
        String,
        Rc<hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::nilpotent::Limit>,
    >,
    og6_cusps: HashMap<String, Rc<hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::Cusp>>,
    og6_lattices: HashMap<String, Rc<hyperelliptic_kummer_surface::mukai::abelian::og6::Lattice>>,
    mukai_moduli: HashMap<String, Rc<hyperelliptic_kummer_surface::mukai::moduli::Lattice>>,
    hilbert_square_form: Option<Rc<hilbert_square::Form>>,
    richelot_isogenies: HashMap<String, Rc<crate::richelot::Isogeny>>,
    richelot_descents: HashMap<String, Rc<crate::richelot::descent::Descent>>,
    bielliptic_quotients: HashMap<String, Rc<crate::bielliptic::Quotient>>,
    elliptic_gluings: HashMap<String, Rc<crate::elliptic_gluing::Gluing>>,
    hecke_eigenpackets: HashMap<String, Rc<eigenpacket::Packet>>,
    hecke_residuals: HashMap<String, Rc<residual_hecke::Specialization>>,
    hecke_eigenperiods: HashMap<String, Rc<period::eigenperiod::Eigenperiod>>,
    mazur_ordinary: HashMap<String, Rc<period::mazur_tate::Ordinary>>,
    period_bases: HashMap<period::Grade, Rc<period::Basis>>,
    integral_period_bases: HashMap<period::Grade, Rc<period::haberland::IntegralBasis>>,
    siegel_maass_forms: HashMap<String, Rc<maass::Form>>,
    siegel_chambers: HashMap<String, Rc<farey::Period>>,
    pub(crate) igusa_inverse_cache: inverse::Cache,
    rv_unram_rings: HashMap<String, Rc<UnramRing>>,
    hodge_wick_forms: HashMap<String, Rc<hyperelliptic_kummer_surface::hodge_wick::Form>>,
    rv_witt_towers: HashMap<String, Rc<RvWittTower>>,
    rv_tor: HashMap<TypeId, Rc<Tor>>,
    minor_bases: HashMap<MinorsKey, Rc<Minors>>,
    rv_bar_paths: Vec<(String, Rc<rv_bar::Path>)>,
    rv_cyclic_homology: HashMap<TypeId, Rc<rv_hh_homology::Slice<rv_negative::Term>>>,
    rv_negative_homology: HashMap<TypeId, Rc<rv_hh_homology::Slice<rv_negative::Term>>>,
    rv_hh_homology: HashMap<TypeId, Rc<rv_hh_homology::Slice<rv_hh::Word>>>,
    rv_path_ext: HashMap<TypeId, Rc<rv_bar_dual::Ext>>,
    rv_bar_pages: HashMap<TypeId, Rc<rv_bar_page::Page>>,
    rho_pages: HashMap<TypeId, Rc<rho_bock::Page>>,
    real_ext: HashMap<TypeId, Rc<real_ext::Group>>,
    real_bars: HashMap<(usize, i128), Rc<real_bars::Group>>,
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
            finite_groups: HashMap::new(),
            crossed_modules: HashMap::new(),
            postnikov_systems: HashMap::new(),
            set_families: HashMap::new(),
            set_origins: HashMap::new(),
            set_recollements: HashMap::new(),
            set_monodromies: HashMap::new(),
            set_linearizations: HashMap::new(),
            set_logic: HashMap::new(),
            set_topologies: HashMap::new(),
            set_reflections: HashMap::new(),
            set_reflection_certificates: HashMap::new(),
            set_sheafifications: HashMap::new(),
            set_sheaf_certificates: HashMap::new(),
            set_elements: HashMap::new(),
            mixed_complexes: Vec::new(),
            mixed_hh_origins: HashMap::new(),
            mixed_indices: HashMap::new(),
            mixed_tensors: HashMap::new(),
            mixed_cones: HashMap::new(),
            mixed_homs: HashMap::new(),
            mixed_minimals: Vec::new(),
            mixed_minimal_ids: HashMap::new(),
            mixed_bars: HashMap::new(),
            mixed_homology: HashMap::new(),
            mixed_pages: HashMap::new(),
            sheaves: HashMap::new(),
            sheaf_sites: HashMap::new(),
            sheaf_maps: HashMap::new(),
            convolutions: HashMap::new(),
            convolution_maps: HashMap::new(),
            disk_sequences: HashMap::new(),
            disk_sequence_maps: HashMap::new(),
            parabolic_groups: HashMap::new(),
            disk_quivers: HashMap::new(),
            disk_homs: HashMap::new(),
            disk_realizations: HashMap::new(),
            cone_gluings: HashMap::new(),
            cone_homotopies: HashMap::new(),
            bc_squares: HashMap::new(),
            bc_comparisons: HashMap::new(),
            sheaf_kuranishi: HashMap::new(),
            sheaf_homs: HashMap::new(),
            sheaf_complexes: HashMap::new(),
            sheaf_resolutions: HashMap::new(),
            sheaf_supports: HashMap::new(),
            leray_complexes: HashMap::new(),
            serre_shifts: HashMap::new(),
            type_ids: HashMap::new(),
            milnor_bases: HashMap::new(),
            cyclotomics: HashMap::new(),
            v_bars: HashMap::new(),
            kh_braids: HashMap::new(),
            rv_ideals: HashMap::new(),
            rv_curves: HashMap::new(),
            rv_periods: HashMap::new(),
            rv_as_covers: HashMap::new(),
            rv_phodge_modules: HashMap::new(),
            rv_qp_quadratics: HashMap::new(),
            legendre_cm_families: HashMap::new(),
            rv_crystals: HashMap::new(),
            rv_as_odd_covers: HashMap::new(),
            rv_as_two_covers: HashMap::new(),
            rv_witt_covers: HashMap::new(),
            elliptic_curves: HashMap::new(),
            elliptic_end_orders: HashMap::new(),
            ec_cm_cycles: HashMap::new(),
            ec_cm_actions: HashMap::new(),
            ec_torsion_levels: HashMap::new(),
            elliptic_isogenies: HashMap::new(),
            elliptic_modular_polynomials: HashMap::new(),
            hyperelliptic: HashMap::new(),
            hw_closed_curves: HashMap::new(),
            jacobian_charts: HashMap::new(),
            jacobian_levels: HashMap::new(),
            jacobian_power_levels: HashMap::new(),
            enriques_quotients: HashMap::new(),
            k3_pencils: HashMap::new(),
            k3_cassels: HashMap::new(),
            k3_selmer_four: HashMap::new(),
            k3_pencil_covers: HashMap::new(),
            k3_pencil_mw_images: HashMap::new(),
            k3_mordell: HashMap::new(),
            k3_automorphisms: HashMap::new(),
            k3_periods: HashMap::new(),
            cm_ring_orders: HashMap::new(),
            karoubi_objects: HashMap::new(),
            deligne_square: None,
            deligne_towers: HashMap::new(),
            hilbert_strata: HashMap::new(),
            specht_strata: HashMap::new(),
            cm_class_polynomials: HashMap::new(),
            cm_class_discriminants: HashMap::new(),
            nakajima_gram: None,
            nakajima_inverse_gram: None,
            mukai_moduli: HashMap::new(),
            og6_nilpotents: HashMap::new(),
            og6_cusps: HashMap::new(),
            og6_lattices: HashMap::new(),
            hilbert_square_form: None,
            richelot_isogenies: HashMap::new(),
            richelot_descents: HashMap::new(),
            bielliptic_quotients: HashMap::new(),
            elliptic_gluings: HashMap::new(),
            hecke_eigenpackets: HashMap::new(),
            hecke_residuals: HashMap::new(),
            hecke_eigenperiods: HashMap::new(),
            mazur_ordinary: HashMap::new(),
            period_bases: HashMap::new(),
            integral_period_bases: HashMap::new(),
            siegel_maass_forms: HashMap::new(),
            siegel_chambers: HashMap::new(),
            igusa_inverse_cache: inverse::Cache::default(),
            rv_unram_rings: HashMap::new(),
            hodge_wick_forms: HashMap::new(),
            rv_witt_towers: HashMap::new(),
            rv_tor: HashMap::new(),
            minor_bases: HashMap::new(),
            rv_bar_paths: Vec::new(),
            rv_cyclic_homology: HashMap::new(),
            rv_negative_homology: HashMap::new(),
            rv_hh_homology: HashMap::new(),
            rv_path_ext: HashMap::new(),
            rv_bar_pages: HashMap::new(),
            rho_pages: HashMap::new(),
            real_ext: HashMap::new(),
            real_bars: HashMap::new(),
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

    pub(crate) fn define_finite_group(
        &mut self,
        name: Label,
        group: rv_crossed::Group,
    ) -> Result<()> {
        if self.finite_groups.contains_key(name.text.as_str()) {
            return Err(Error::new(name.span, "finite group is already defined"));
        }
        self.allocate(group.cells() + name.text.len())?;
        self.finite_groups.insert(name.text.into(), Rc::new(group));
        Ok(())
    }

    pub(crate) fn finite_group(&self, name: &Label) -> Result<Rc<rv_crossed::Group>> {
        self.finite_groups
            .get(name.text.as_str())
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown finite group '{}'", name.text)))
    }

    pub(crate) fn finite_group_type(&mut self, name: &Label) -> Result<TypeId> {
        self.finite_group(name)?;
        self.inhabited_type(Kind::FiniteGroup(name.text.as_str().into()), 0)
    }

    pub(crate) fn define_crossed_module(
        &mut self,
        name: Label,
        crossed: rv_crossed::Crossed,
    ) -> Result<()> {
        if self.crossed_modules.contains_key(name.text.as_str()) {
            return Err(Error::new(name.span, "crossed module is already defined"));
        }
        self.allocate(crossed.cells() + name.text.len())?;
        self.crossed_modules
            .insert(name.text.into(), Rc::new(crossed));
        Ok(())
    }

    pub(crate) fn crossed_module(&self, name: &Label) -> Result<Rc<rv_crossed::Crossed>> {
        self.crossed_modules
            .get(name.text.as_str())
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown crossed module '{}'", name.text)))
    }

    pub(crate) fn define_postnikov(
        &mut self,
        name: Label,
        post: postnikov::Postnikov,
    ) -> Result<()> {
        if self.postnikov_systems.contains_key(name.text.as_str()) {
            return Err(Error::new(name.span, "Postnikov system is already defined"));
        }
        self.allocate(post.cells() + name.text.len())?;
        self.postnikov_systems
            .insert(name.text.into(), Rc::new(post));
        Ok(())
    }

    pub(crate) fn postnikov(&self, name: &Label) -> Result<Rc<postnikov::Postnikov>> {
        self.postnikov_systems
            .get(name.text.as_str())
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown Postnikov system '{}'", name.text),
                )
            })
    }

    pub(crate) fn postnikov_cochain_type(&mut self, name: &Label, degree: usize) -> Result<TypeId> {
        self.postnikov(name)?;
        if !matches!(degree, 2 | 3) {
            return Err(Error::new(
                name.span,
                "Postnikov cochain degree must be 2 or 3",
            ));
        }
        self.inhabited_type(Kind::PostnikovCochain(name.text.as_str().into(), degree), 0)
    }

    pub(crate) fn postnikov_cell_type(&mut self, name: &Label) -> Result<TypeId> {
        self.postnikov(name)?;
        self.inhabited_type(Kind::PostnikovCell(name.text.as_str().into()), 0)
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
            Kind::FiniteGroup(ref name) => self.finite_groups[name.as_ref()].order(),
            Kind::PostnikovCochain(ref name, degree) => {
                let post = &self.postnikov_systems[name.as_ref()];
                (post.prime.get() as usize).saturating_pow(post.size(degree) as u32)
            }
            Kind::PostnikovCell(ref name) => {
                let post = &self.postnikov_systems[name.as_ref()];
                post.order() * post.prime.get() as usize
            }
            Kind::QModule(_) => 1,
            Kind::SetNat(..) => {
                if inhabited {
                    usize::MAX
                } else {
                    0
                }
            }
            Kind::MixedMinimum(..) => usize::MAX,
            Kind::MixedChain(..) | Kind::MixedPage(..) | Kind::MixedHomology(..) => usize::MAX,
            Kind::Sheaf(_, ref name, k) if k >= self.sheaves[name.as_ref()].height => 1,
            Kind::SheafComplex(_, ref name, k)
                if k >= self.sheaf_complexes[name.as_ref()].height() =>
            {
                1
            }
            Kind::SheafPage(ref name, grade)
                if !grade.has_leading(self.sheaf_complexes[name.as_ref()].filtered()) =>
            {
                1
            }
            Kind::Convolution(form, ref name, prime) => {
                let dimension = self.convolutions[name].dimension(form);
                (prime.get() as usize).saturating_pow(u32::try_from(dimension).unwrap_or(u32::MAX))
            }
            Kind::Parabolic(form, ref name, prime) => {
                let dimension = self.parabolic_groups[name].dimension(form);
                (prime.get() as usize).saturating_pow(u32::try_from(dimension).unwrap_or(u32::MAX))
            }
            Kind::DiskHom(ref from, ref to, prime) => {
                let rank = self.disk_homs[&(from.clone(), to.clone())].rank();
                (prime.get() as usize).saturating_pow(u32::try_from(rank).unwrap_or(u32::MAX))
            }
            Kind::SheafDeformation(_, 0) => 1,
            Kind::SheafDeformation(..) => usize::MAX,
            Kind::SheafHom(_, ref ends, k) if k >= self.sheaf_homs[ends].height => 1,
            Kind::Sheaf(..)
            | Kind::SheafComplex(..)
            | Kind::SheafPage(..)
            | Kind::SheafHom(..)
            | Kind::SheafLocal(..) => usize::MAX,
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
            | Kind::BPPage(..)
            | Kind::BPVPage(..)
            | Kind::BPModExt(..) => usize::MAX,
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
            | Kind::MCobar(..)
            | Kind::RMCobar(..)
            | Kind::RMRhoCobar(..)
            | Kind::RhoLoc(..)
            | Kind::RhoMin(..)
            | Kind::Ext(..)
            | Kind::ExtAn(..)
            | Kind::MExt(..)
            | Kind::RMExt(..)
            | Kind::RMRhoExt(..)
            | Kind::RMBar(..)
            | Kind::RMTriple(..)
            | Kind::RhoBock(..)
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
            Kind::Lee(ref braid, grade) => {
                if self.kh_braids[braid].lee_possible(grade.i, grade.f) {
                    usize::MAX
                } else {
                    1
                }
            }
            Kind::LeeDual(ref braid, i) => {
                let braid = &self.kh_braids[braid];
                if braid.lee_possible(i, braid.lee_lower_bound()) {
                    usize::MAX
                } else {
                    1
                }
            }
            Kind::Khovanov(ref braid, grade, form) => {
                if self.kh_braids[braid].possible(grade, form) {
                    usize::MAX
                } else {
                    1
                }
            }
            Kind::RvTl(grade) if !grade.possible() => 1,
            Kind::QKoszul(..)
            | Kind::QExt(..)
            | Kind::JetHom(..)
            | Kind::JetPage(..)
            | Kind::RVHom(..)
            | Kind::RVLaurent(..)
            | Kind::RVIwHecke(..)
            | Kind::RVSphHecke(..)
            | Kind::RVSatake(..)
            | Kind::RVFusion(..)
            | Kind::RVCyclo(..)
            | Kind::RVCycloQ(..)
            | Kind::RVQpQuad(..)
            | Kind::CMOrder(..)
            | Kind::ECEndOrder(..)
            | Kind::RVWzw(..)
            | Kind::RvTl(..) => usize::MAX,
            Kind::RVResSpace(grade) => grade.count(),
            Kind::RVValMat(grade) if grade.rank == 0 || grade.rank == grade.ground => 1,
            Kind::RVValMat(..) => usize::MAX,
            Kind::RVStableTree(3) | Kind::RVStableMarks(_, 3) => 1,
            Kind::RVHyperArithmetic(..) | Kind::RVHyperCover(..) => usize::MAX,
            Kind::RVStableTree(..) | Kind::RVStableMarks(..) => usize::MAX,
            Kind::ECCurve(..) | Kind::ECIso(..) => 1,
            Kind::Jacobian(ref name) => {
                if self.hw_closed_curves[name.as_ref()].grade.genus == 0 {
                    1
                } else {
                    usize::MAX
                }
            }
            Kind::ECMap(..) | Kind::Kummer(..) | Kind::KummerK3(..) | Kind::KummerDiv(..) => {
                usize::MAX
            }
            Kind::Enriques(ref name) => self.enriques_quotients[name.as_ref()].count(),
            Kind::K3PencilBase(ref name) => self.k3_pencils[name.as_ref()].count(),
            Kind::K3PencilSection(..)
            | Kind::K3PencilKummer2(..)
            | Kind::K3PencilCoverPoint(..) => usize::MAX,
            Kind::KummerDisc(..) => usize::from(nikulin::ORDER),
            Kind::CMClass(ref name) => self.cm_ring_orders[name.as_ref()].cardinality_bound(),
            Kind::KummerWeil(..) => usize::MAX,
            Kind::PeriodPolynomial(grade) if grade.is_trivial() => 1,
            Kind::MazurOrdinary(_, modulus) => usize::try_from(modulus).unwrap_or(usize::MAX),
            Kind::MazurAmice(ref grade) => grade.cardinality(),
            Kind::K3Lattice(..)
            | Kind::Mukai(..)
            | Kind::MukaiQuotient(..)
            | Kind::OG10(..)
            | Kind::AbMukai(..)
            | Kind::OG6(..)
            | Kind::SiegelWeight
            | Kind::ModularHecke(..)
            | Kind::HeckeEigen(..)
            | Kind::HeckeSymbol(..)
            | Kind::MazurTate(..)
            | Kind::PeriodPolynomial(..) => usize::MAX,
            Kind::Nakajima(ref grade)
            | Kind::HilbertCohomology(ref grade)
            | Kind::HodgeTensor(ref grade)
            | Kind::Verbitsky(ref grade) => {
                if grade.supports_states() {
                    usize::MAX
                } else {
                    1
                }
            }
            Kind::HilbertSquare(_, degree) if hilbert_square::dimension(degree) == 0 => 1,
            Kind::HilbertSquare(..) => usize::MAX,
            Kind::K3Clifford(..) | Kind::OG6Clifford(..) | Kind::K3Spinor(..) => usize::MAX,
            Kind::KSExterior(_, 0..=64) => usize::MAX,
            Kind::KSExterior(..) => 1,
            Kind::KSNomizu(_, _, 0..=65) => usize::MAX,
            Kind::KSNomizu(..) => 1,
            Kind::KSEnd(..) => usize::MAX,
            Kind::QZ => usize::MAX,
            Kind::ECPoint(..) => usize::MAX,
            Kind::ECLevel(ref name, matrix) => (self.ec_torsion_levels[name.as_ref()].modulus
                as usize)
                .saturating_pow(if matrix { 4 } else { 2 }),
            Kind::ECKernel(_, ell) => {
                if inhabited {
                    ell as usize + 1
                } else {
                    0
                }
            }
            Kind::HWScalar(ring) => ring.count(1),
            Kind::HWMatrix(ring, size) => ring.count(size * size),
            Kind::HWWMatrix(ref grade) => {
                hw_closed::count(&self.rv_unram_rings[grade.ring.as_ref()], grade.genus)
            }
            Kind::RVUnram(ref name) => self.rv_unram_rings[name.as_ref()].count(),
            Kind::RVQGr(..)
            | Kind::RVQGrP(..)
            | Kind::RVOre
            | Kind::RVWittSpectrum
            | Kind::RVSchurCharacter(..)
            | Kind::RVStableSchurCharacter
            | Kind::RVPartitionAlgebra(..)
            | Kind::RVPartitionLocal(..)
            | Kind::RVKaroubiMap(..) => usize::MAX,
            Kind::RVCTState(..) => usize::MAX,
            Kind::RVCTResidue(level) => usize::try_from(level.modulus).unwrap_or(usize::MAX),
            Kind::RVHive(..) => {
                if inhabited {
                    usize::MAX
                } else {
                    0
                }
            }
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
            Kind::RvCyclicChain(..) | Kind::RvCyclicHomology(..) => usize::MAX,
            Kind::RVNCChain(..) | Kind::RvNegativeCyclic(..) => usize::MAX,
            Kind::RVBarCobar(..) | Kind::RVHHChain(..) | Kind::RvHochschild(..) => usize::MAX,
            Kind::RVBar(..) | Kind::RVBarE(..) | Kind::RVBarDual(..) | Kind::RVPathExt(..) => {
                usize::MAX
            }
            Kind::RVKoszul(..) | Kind::RVTor(..) => usize::MAX,
            Kind::RVPoly(ref selected, d) if selected.is_empty() && d != 0 => 1,
            Kind::RVPoly(ref selected, d) if selected.len() < 2 || d == 0 => 2,
            Kind::ASDeRham(ref name, field) => {
                let dimension = 2 * self.rv_as_covers[name.as_ref()].genus();
                field
                    .order()
                    .saturating_pow(u32::try_from(dimension).unwrap_or(u32::MAX))
            }
            Kind::BinaryScalar(field) => field.order(),
            Kind::EOType(genus) => support_space_size(genus),
            Kind::PrimeScalar(prime) => prime.get() as usize,
            Kind::ASOddDeRham(ref name, prime) | Kind::ASOddTate(ref name, prime, _) => {
                let dimension = 2 * self.rv_as_odd_covers[name.as_ref()].genus();
                (prime.get() as usize).saturating_pow(u32::try_from(dimension).unwrap_or(u32::MAX))
            }
            Kind::ASTwoDiff(ref name, prime) => {
                let genus = self.rv_as_two_covers[name.as_ref()].genus();
                (prime.get() as usize).saturating_pow(u32::try_from(genus).unwrap_or(u32::MAX))
            }
            Kind::ASOddDiff(ref name, prime) => {
                let genus = self.rv_as_odd_covers[name.as_ref()].genus();
                (prime.get() as usize).saturating_pow(u32::try_from(genus).unwrap_or(u32::MAX))
            }
            Kind::CurveDiff(family, ref name, field) => {
                let genus = self.curve_genus(family, name);
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
            bp_vpage: None,
            bp_mod_ext: None,
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
            hive_witness: None,
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
                | Kind::SetNat(..)
                | Kind::MixedChain(..)
                | Kind::MixedHomology(..)
                | Kind::MixedMinimum(..)
                | Kind::MixedPage(..)
                | Kind::Sheaf(..)
                | Kind::SheafComplex(..)
                | Kind::SheafDeformation(..)
                | Kind::Convolution(..)
                | Kind::Parabolic(..)
                | Kind::DiskHom(..)
                | Kind::SheafHom(..)
                | Kind::SheafPage(..)
                | Kind::SheafLocal(..)
                | Kind::Operations(..)
                | Kind::Cobar(..)
                | Kind::MCobar(..)
                | Kind::RMCobar(..)
                | Kind::RMRhoCobar(..)
                | Kind::RhoLoc(..)
                | Kind::RhoMin(..)
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
                | Kind::Ext(..)
                | Kind::ExtAn(..)
                | Kind::MExt(..)
                | Kind::RMExt(..)
                | Kind::RMRhoExt(..)
                | Kind::RMBar(..)
                | Kind::RMTriple(..)
                | Kind::RhoBock(..)
                | Kind::XExt(..)
                | Kind::MapCobar(..)
                | Kind::MapExt(..)
                | Kind::QKoszul(..)
                | Kind::JetHom(..)
                | Kind::JetPage(..)
                | Kind::QExt(..)
                | Kind::RVHom(..)
                | Kind::CurveDiff(..)
                | Kind::ASOddDiff(..)
                | Kind::ASTwoDiff(..)
                | Kind::ASOddDeRham(..)
                | Kind::ASOddTate(..)
                | Kind::PrimeScalar(..)
                | Kind::MazurOrdinary(..)
                | Kind::ASDeRham(..)
                | Kind::BinaryScalar(..)
                | Kind::EOType(..)
                | Kind::RVLaurent(..)
                | Kind::RVIwHecke(..)
                | Kind::RVSphHecke(..)
                | Kind::RVSatake(..)
                | Kind::RVHive(..)
                | Kind::RVQGr(..)
                | Kind::RVQGrP(..)
                | Kind::RVOre
                | Kind::RVWittSpectrum
                | Kind::RVSchurCharacter(..)
                | Kind::RVStableSchurCharacter
                | Kind::RVPartitionAlgebra(..)
                | Kind::RVPartitionLocal(..)
                | Kind::RVKaroubiMap(..)
                | Kind::RVCTState(..)
                | Kind::RVCTResidue(..)
                | Kind::RVFusion(..)
                | Kind::RVCyclo(..)
                | Kind::ECCurve(..)
                | Kind::Jacobian(..)
                | Kind::Kummer(..)
                | Kind::KummerK3(..)
                | Kind::Enriques(..)
                | Kind::K3PencilBase(..)
                | Kind::K3PencilSection(..)
                | Kind::K3PencilKummer2(..)
                | Kind::K3PencilCoverPoint(..)
                | Kind::KummerDisc(..)
                | Kind::KummerWeil(..)
                | Kind::KummerDiv(..)
                | Kind::K3Lattice(..)
                | Kind::Nakajima(..)
                | Kind::HilbertCohomology(..)
                | Kind::HodgeTensor(..)
                | Kind::Verbitsky(..)
                | Kind::K3Clifford(..)
                | Kind::OG6Clifford(..)
                | Kind::KSExterior(..)
                | Kind::KSNomizu(..)
                | Kind::FiniteGroup(..)
                | Kind::PostnikovCochain(..)
                | Kind::PostnikovCell(..)
                | Kind::KSEnd(..)
                | Kind::QZ
                | Kind::K3Spinor(..)
                | Kind::Mukai(..)
                | Kind::MukaiQuotient(..)
                | Kind::OG10(..)
                | Kind::AbMukai(..)
                | Kind::OG6(..)
                | Kind::HilbertSquare(..)
                | Kind::SiegelWeight
                | Kind::ModularHecke(..)
                | Kind::HeckeEigen(..)
                | Kind::HeckeSymbol(..)
                | Kind::MazurTate(..)
                | Kind::MazurAmice(..)
                | Kind::PeriodPolynomial(..)
                | Kind::ECMap(..)
                | Kind::ECPoint(..)
                | Kind::ECLevel(..)
                | Kind::ECKernel(..)
                | Kind::ECIso(..)
                | Kind::HWMatrix(..)
                | Kind::HWWMatrix(..)
                | Kind::HWScalar(..)
                | Kind::RVUnram(..)
                | Kind::RVCycloQ(..)
                | Kind::RVQpQuad(..)
                | Kind::RVResSpace(..)
                | Kind::RVValMat(..)
                | Kind::RVHyperArithmetic(..)
                | Kind::RVHyperCover(..)
                | Kind::RVStableTree(..)
                | Kind::RVStableMarks(..)
                | Kind::CMOrder(..)
                | Kind::ECEndOrder(..)
                | Kind::CMClass(..)
                | Kind::RVWzw(..)
                | Kind::Lee(..)
                | Kind::LeeDual(..)
                | Kind::Khovanov(..)
                | Kind::RvTl(..)
                | Kind::VKoszul(..)
                | Kind::VExt(..)
                | Kind::VBock(..)
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

    pub(crate) fn sheaf_site_data(&mut self, base: TypeId) -> Result<Rc<sheaf::PosetData>> {
        if let Some(data) = self.sheaf_sites.get(&base) {
            return Ok(data.clone());
        }
        let poset = self.poset(base).expect("named sheaf base poset");
        let data = Rc::new(sheaf::PosetData::build(self, &poset)?);
        self.sheaf_sites.insert(base, data.clone());
        Ok(data)
    }

    pub(crate) fn check_new_sheaf(&self, name: &Label) -> Result<()> {
        if self.sheaves.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("sheaf '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn check_new_mixed(&self, name: &Label) -> Result<()> {
        if self.mixed_indices.contains_key(&name.text) {
            return Err(Error::new(name.span, "mixed complex is already defined"));
        }
        Ok(())
    }

    pub(crate) fn define_mixed(&mut self, name: String, complex: mixed::Complex) -> mixed::Id {
        let id = mixed::Id(self.mixed_complexes.len());
        self.mixed_indices.insert(name.clone(), id);
        self.mixed_complexes.push((name, Rc::new(complex)));
        id
    }

    pub(crate) fn define_hh_mixed(
        &mut self,
        name: String,
        complex: mixed::Complex,
        origin: rv_hh_mixed::Origin,
    ) {
        let id = self.define_mixed(name, complex);
        self.mixed_hh_origins.insert(id, Rc::new(origin));
    }

    pub(crate) fn mixed_hh_origin(&self, id: mixed::Id) -> Result<Rc<rv_hh_mixed::Origin>> {
        self.mixed_hh_origins
            .get(&id)
            .cloned()
            .ok_or_else(|| Error::plain("operation requires a mixed complex imported by rvHHMixed"))
    }

    pub(crate) fn define_mixed_cone(
        &mut self,
        name: String,
        complex: mixed::Complex,
        cone: mixed_cone::Cone,
    ) {
        let id = self.define_mixed(name, complex);
        self.mixed_cones.insert(id, Rc::new(cone));
    }

    pub(crate) fn mixed_cone(&self, id: mixed::Id) -> Result<Rc<mixed_cone::Cone>> {
        self.mixed_cones
            .get(&id)
            .cloned()
            .ok_or_else(|| Error::plain("operation requires a mixed mapping cone"))
    }

    pub(crate) fn define_mixed_hom(
        &mut self,
        name: String,
        complex: mixed::Complex,
        hom: mixed_hom::Hom,
    ) {
        let id = self.define_mixed(name, complex);
        self.mixed_homs.insert(id, Rc::new(hom));
    }

    pub(crate) fn mixed_hom(&self, id: mixed::Id) -> Result<Rc<mixed_hom::Hom>> {
        self.mixed_homs
            .get(&id)
            .cloned()
            .ok_or_else(|| Error::plain("operation requires a mixed complex defined by mixhom"))
    }

    pub(crate) fn define_mixed_tensor(
        &mut self,
        name: String,
        complex: mixed::Complex,
        tensor: mixed_tensor::Tensor,
    ) {
        let id = self.define_mixed(name, complex);
        self.mixed_tensors.insert(id, Rc::new(tensor));
    }

    pub(crate) fn mixed_tensor(&self, id: mixed::Id) -> Result<Rc<mixed_tensor::Tensor>> {
        self.mixed_tensors
            .get(&id)
            .cloned()
            .ok_or_else(|| Error::plain("operation requires a mixed complex defined by mixtensor"))
    }

    pub(crate) fn mixed_id(&self, name: &Label) -> Result<mixed::Id> {
        self.mixed_indices
            .get(&name.text)
            .copied()
            .ok_or_else(|| Error::new(name.span, format!("unknown mixed complex '{}'", name.text)))
    }

    pub(crate) fn mixed_complex(&self, id: mixed::Id) -> Rc<mixed::Complex> {
        self.mixed_complexes[id.0].1.clone()
    }

    pub(crate) fn mixed_name(&self, id: mixed::Id) -> &str {
        &self.mixed_complexes[id.0].0
    }

    pub(crate) fn mixed_chain_type(&mut self, grade: mixed::Grade, total: bool) -> Result<TypeId> {
        self.inhabited_type(Kind::MixedChain(grade, total), 0)
    }

    pub(crate) fn mixed_homology_type(&mut self, grade: mixed::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::MixedHomology(grade), 0)
    }

    pub(crate) fn mixed_homology(
        &mut self,
        grade: mixed::Grade,
    ) -> Result<Rc<mixed_homology::Slice>> {
        if let Some(slice) = self.mixed_homology.get(&grade) {
            return Ok(slice.clone());
        }
        let slice = Rc::new(mixed_homology::Slice::new(self, grade)?);
        self.mixed_homology.insert(grade, slice.clone());
        Ok(slice)
    }

    pub(crate) fn mixed_bars(&mut self, id: mixed::Id) -> Result<Rc<mixed_bars::Bars>> {
        if let Some(bars) = self.mixed_bars.get(&id) {
            return Ok(bars.clone());
        }
        let complex = self.mixed_complex(id);
        let bars = Rc::new(mixed_bars::Bars::build(self, &complex)?);
        self.mixed_bars.insert(id, bars.clone());
        Ok(bars)
    }

    pub(crate) fn show_mixed_bars(&mut self, id: mixed::Id) -> Result<String> {
        let bars = self.mixed_bars(id)?;
        let mut out = Text::new(MAX_SHOWN_BYTES);
        out.push(&format!("mixbars[{}] = free [", self.mixed_name(id)))?;
        for (index, bar) in bars
            .bars()
            .iter()
            .filter(|bar| bar.torsion.is_none())
            .enumerate()
        {
            if index != 0 {
                out.push(", ")?;
            }
            out.push(&bar.birth.to_string())?;
        }
        out.push("], torsion [")?;
        for (index, (birth, torsion)) in bars
            .bars()
            .iter()
            .filter_map(|bar| bar.torsion.as_ref().map(|t| (bar.birth, t)))
            .enumerate()
        {
            if index != 0 {
                out.push(", ")?;
            }
            out.push(&format!("({}, {})", birth, torsion.length))?;
        }
        out.push("]")?;
        Ok(out.string)
    }

    pub(crate) fn check_new_mixed_minimal(&self, name: &Label) -> Result<()> {
        if self.mixed_minimal_ids.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                "minimal mixed model is already defined",
            ));
        }
        Ok(())
    }

    pub(crate) fn define_mixed_minimal(&mut self, name: String, data: mixed_minimal::Minimal) {
        let id = mixed_minimal::Id(self.mixed_minimals.len());
        self.mixed_minimal_ids.insert(name.clone(), id);
        self.mixed_minimals.push((name, Rc::new(data)));
    }

    pub(crate) fn mixed_minimal_id(&self, name: &Label) -> Result<mixed_minimal::Id> {
        self.mixed_minimal_ids
            .get(&name.text)
            .copied()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown minimal mixed model '{}'", name.text),
                )
            })
    }

    pub(crate) fn mixed_minimal_name(&self, id: mixed_minimal::Id) -> &str {
        &self.mixed_minimals[id.0].0
    }

    pub(crate) fn mixed_minimal(&self, id: mixed_minimal::Id) -> Rc<mixed_minimal::Minimal> {
        self.mixed_minimals[id.0].1.clone()
    }

    pub(crate) fn mixed_minimal_type(
        &mut self,
        grade: mixed_minimal::Grade,
        total: bool,
    ) -> Result<TypeId> {
        self.inhabited_type(Kind::MixedMinimum(grade, total), 0)
    }

    pub(crate) fn mixed_page_type(&mut self, grade: mixed_page::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::MixedPage(grade), 0)
    }

    pub(crate) fn mixed_page(&mut self, ty: TypeId) -> Result<Rc<mixed_page::Page>> {
        if let Some(page) = self.mixed_pages.get(&ty) {
            return Ok(page.clone());
        }
        let Kind::MixedPage(grade) = *self.kind(ty) else {
            unreachable!()
        };
        let page = Rc::new(mixed_page::Page::build(self, grade)?);
        self.mixed_pages.insert(ty, page.clone());
        Ok(page)
    }

    pub(crate) fn define_sheaf(&mut self, name: String, sheaf: Sheaf) {
        self.define_sheaf_shared(name, Rc::new(sheaf));
    }

    pub(crate) fn define_sheaf_shared(&mut self, name: String, sheaf: Rc<Sheaf>) {
        self.sheaves.insert(name, sheaf);
    }

    pub(crate) fn sheaf(&self, name: &Label) -> Result<Rc<Sheaf>> {
        self.sheaf_named(&name.text)
            .ok_or_else(|| Error::new(name.span, format!("unknown sheaf '{}'", name.text)))
    }

    pub(crate) fn sheaf_named(&self, name: &str) -> Option<Rc<Sheaf>> {
        self.sheaves.get(name).cloned()
    }

    pub(crate) fn sheaf_type(
        &mut self,
        form: SheafForm,
        name: &Label,
        degree: usize,
    ) -> Result<TypeId> {
        self.sheaf(name)?;
        self.known_sheaf_type(form, name.text.as_str().into(), degree)
    }

    pub(crate) fn known_sheaf_type(
        &mut self,
        form: SheafForm,
        name: Rc<str>,
        degree: usize,
    ) -> Result<TypeId> {
        self.inhabited_type(Kind::Sheaf(form, name, degree), 0)
    }

    pub(crate) fn check_new_sheaf_kuranishi(&self, name: &Label) -> Result<()> {
        if self.sheaf_kuranishi.contains_key(name.text.as_str()) {
            return Err(Error::new(
                name.span,
                format!("Kuranishi chart '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn define_sheaf_kuranishi(
        &mut self,
        name: &Label,
        chart: sheaf_kuranishi::Chart,
    ) -> Result<()> {
        self.check_new_sheaf_kuranishi(name)?;
        self.allocate(1)?;
        self.sheaf_kuranishi
            .insert(name.text.as_str().into(), Rc::new(chart));
        Ok(())
    }

    pub(crate) fn sheaf_kuranishi(&self, name: &Label) -> Result<Rc<sheaf_kuranishi::Chart>> {
        self.sheaf_kuranishi
            .get(name.text.as_str())
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown Kuranishi chart '{}'", name.text),
                )
            })
    }

    pub(crate) fn show_kuranishi_equation(&mut self, name: &Label, index: usize) -> Result<String> {
        let chart = self.sheaf_kuranishi(name)?;
        let mut out = Text::new(MAX_SHOWN_BYTES);
        out.push(&format!("shkuranishiEq[{}, {index}] = ", name.text))?;
        let mut first = true;
        for (mono, scalar) in chart.equation(self, index)? {
            self.spend(1 + mono.factors().len())?;
            if !first {
                out.push(" + ")?;
            }
            first = false;
            if scalar != 1 {
                out.push(&format!("{scalar}*"))?;
            }
            for (i, &(variable, power)) in mono.factors().iter().enumerate() {
                if i != 0 {
                    out.push("*")?;
                }
                out.push(&format!("x{variable}"))?;
                if power != 1 {
                    out.push(&format!("^{power}"))?;
                }
            }
        }
        if first {
            out.push("0")?;
        }
        Ok(out.string)
    }

    pub(crate) fn known_sheaf_deformation_type(
        &mut self,
        name: Rc<str>,
        order: usize,
    ) -> Result<TypeId> {
        let ends = SheafEnds::endomorphisms(name.clone());
        self.known_sheaf_hom_type(SheafForm::Cochain, ends, 1)?;
        self.inhabited_type(Kind::SheafDeformation(name, order), 0)
    }

    pub(crate) fn sheaf_hom_type(
        &mut self,
        form: SheafForm,
        source: &Label,
        target: &Label,
        degree: usize,
    ) -> Result<TypeId> {
        self.sheaf(source)?;
        self.sheaf(target)?;
        let ends = SheafEnds {
            source: source.text.as_str().into(),
            target: target.text.as_str().into(),
        };
        self.known_sheaf_hom_type(form, ends, degree)
    }

    pub(crate) fn known_sheaf_hom_type(
        &mut self,
        form: SheafForm,
        ends: SheafEnds,
        degree: usize,
    ) -> Result<TypeId> {
        self.ensure_sheaf_hom(&ends)?;
        self.inhabited_type(Kind::SheafHom(form, ends, degree), 0)
    }

    pub(crate) fn ensure_sheaf_hom(&mut self, ends: &SheafEnds) -> Result<Rc<SheafHom>> {
        if let Some(hom) = self.sheaf_homs.get(ends) {
            return Ok(hom.clone());
        }
        let source = self.sheaf_named(&ends.source).expect("registered sheaf");
        let target = self.sheaf_named(&ends.target).expect("registered sheaf");
        let data = Rc::new(SheafHom::build(self, source, target)?);
        self.sheaf_homs.insert(ends.clone(), data.clone());
        Ok(data)
    }

    pub(crate) fn sheaf_hom(&self, ends: &SheafEnds) -> Rc<SheafHom> {
        self.sheaf_homs[ends].clone()
    }

    pub(crate) fn check_new_sheaf_map(&self, name: &Label) -> Result<()> {
        if self.sheaf_maps.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("sheaf map '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn disk_named(&self, name: &str) -> Result<Rc<crate::sheaf_disk::Disk>> {
        self.disk_quivers
            .get(name)
            .cloned()
            .ok_or_else(|| Error::plain(format!("unknown disk quiver '{name}'")))
    }

    pub(crate) fn check_new_disk(&self, name: &Label) -> Result<()> {
        if self.disk_quivers.contains_key(name.text.as_str()) {
            return Err(Error::new(
                name.span,
                format!("disk quiver '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn define_disk(&mut self, name: String, disk: crate::sheaf_disk::Disk) {
        self.define_disk_shared(name, Rc::new(disk));
    }

    pub(crate) fn define_disk_shared(&mut self, name: String, disk: Rc<crate::sheaf_disk::Disk>) {
        self.disk_quivers.insert(name.into(), disk);
    }

    pub(crate) fn disk_sequence(&self, name: &str) -> Result<Rc<DiskSequence>> {
        self.disk_sequences
            .get(name)
            .cloned()
            .ok_or_else(|| Error::plain(format!("unknown disk sequence '{name}'")))
    }

    pub(crate) fn disk_sequence_map(&self, name: &str) -> Result<Rc<DiskSequenceMap>> {
        self.disk_sequence_maps
            .get(name)
            .cloned()
            .ok_or_else(|| Error::plain(format!("unknown disk sequence map '{name}'")))
    }

    pub(crate) fn check_new_disk_sequence(&self, name: &Label, map: bool) -> Result<()> {
        let exists = if map {
            self.disk_sequence_maps.contains_key(name.text.as_str())
        } else {
            self.disk_sequences.contains_key(name.text.as_str())
        };
        if exists {
            let item = if map {
                "disk sequence map"
            } else {
                "disk sequence"
            };
            return Err(Error::new(
                name.span,
                format!("{item} '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn define_disk_sequence(&mut self, name: String, sequence: DiskSequence) {
        self.disk_sequences.insert(name.into(), Rc::new(sequence));
    }

    pub(crate) fn define_disk_sequence_map(&mut self, name: String, map: DiskSequenceMap) {
        self.disk_sequence_maps.insert(name.into(), Rc::new(map));
    }

    pub(crate) fn define_convolution(&mut self, name: String, data: Convolution) {
        let name: Rc<str> = name.into();
        self.disk_sequences
            .insert(name.clone(), data.sequence.clone());
        self.convolutions.insert(name, Rc::new(data));
    }

    pub(crate) fn convolution(&self, name: &str) -> Result<Rc<Convolution>> {
        self.convolutions
            .get(name)
            .cloned()
            .ok_or_else(|| Error::plain(format!("unknown middle convolution '{name}'")))
    }

    pub(crate) fn define_convolution_map(&mut self, name: String, data: ConvolutionMap) {
        let name: Rc<str> = name.into();
        self.disk_sequence_maps
            .insert(name.clone(), data.descended.clone());
        self.convolution_maps.insert(name, Rc::new(data));
    }

    pub(crate) fn convolution_map(&self, name: &str) -> Result<Rc<ConvolutionMap>> {
        self.convolution_maps
            .get(name)
            .cloned()
            .ok_or_else(|| Error::plain(format!("unknown middle convolution map '{name}'")))
    }

    pub(crate) fn convolution_type(
        &mut self,
        name: Rc<str>,
        form: ConvolutionForm,
    ) -> Result<(TypeId, Rc<Convolution>)> {
        let data = self.convolution(&name)?;
        let ty = self.inhabited_type(Kind::Convolution(form, name, data.sequence.prime), 0)?;
        Ok((ty, data))
    }

    pub(crate) fn parabolic_type(
        &mut self,
        name: Rc<str>,
        form: ParabolicForm,
    ) -> Result<(TypeId, Rc<Parabolic>)> {
        let data = if let Some(group) = self.parabolic_groups.get(&name) {
            group.clone()
        } else {
            let sequence = self.disk_sequence(&name)?;
            let group = Rc::new(Parabolic::build(self, sequence)?);
            self.allocate(1)?;
            self.parabolic_groups.insert(name.clone(), group.clone());
            group
        };
        let ty = self.inhabited_type(Kind::Parabolic(form, name, data.sequence.prime), 0)?;
        Ok((ty, data))
    }

    pub(crate) fn disk_hom_type(
        &mut self,
        from: Rc<str>,
        to: Rc<str>,
    ) -> Result<(TypeId, Rc<crate::sheaf_disk::Hom>)> {
        let key = (from.clone(), to.clone());
        let hom = if let Some(hom) = self.disk_homs.get(&key) {
            hom.clone()
        } else {
            let a = self.disk_named(&from)?;
            let b = self.disk_named(&to)?;
            let hom = Rc::new(crate::sheaf_disk::Hom::build(self, a, b)?);
            self.allocate(1)?;
            self.disk_homs.insert(key, hom.clone());
            hom
        };
        let ty = self.inhabited_type(Kind::DiskHom(from, to, hom.source.prime), 0)?;
        Ok((ty, hom))
    }

    pub(crate) fn define_disk_realization(&mut self, value: crate::sheaf_disk::Realization) {
        self.disk_realizations
            .insert(value.name.clone(), Rc::new(value));
    }

    pub(crate) fn disk_realization(
        &self,
        name: &Label,
    ) -> Result<Rc<crate::sheaf_disk::Realization>> {
        self.disk_realizations
            .get(name.text.as_str())
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown realized disk '{}'", name.text)))
    }

    pub(crate) fn define_cone_gluing(&mut self, name: String, value: crate::sheaf_cone::Cone) {
        self.cone_gluings.insert(name, Rc::new(value));
    }

    pub(crate) fn cone_gluing(&self, name: &Label) -> Result<Rc<crate::sheaf_cone::Cone>> {
        self.cone_gluings
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown cone gluing '{}'", name.text)))
    }

    pub(crate) fn check_new_cone_homotopy(&self, name: &Label) -> Result<()> {
        if self.cone_homotopies.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("cone homotopy '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn define_cone_homotopy(
        &mut self,
        name: String,
        value: crate::sheaf_cone::Homotopy,
    ) {
        self.cone_homotopies.insert(name, Rc::new(value));
    }

    pub(crate) fn cone_homotopy(&self, name: &Label) -> Result<Rc<crate::sheaf_cone::Homotopy>> {
        self.cone_homotopies
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown cone homotopy '{}'", name.text)))
    }

    pub(crate) fn check_new_bc_square(&self, name: &Label) -> Result<()> {
        if self.bc_squares.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("base-change square '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn define_bc_square(
        &mut self,
        name: String,
        value: crate::sheaf_base_change::Square,
    ) {
        self.bc_squares.insert(name, Rc::new(value));
    }

    pub(crate) fn bc_square(&self, name: &Label) -> Result<Rc<crate::sheaf_base_change::Square>> {
        self.bc_squares.get(&name.text).cloned().ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown base-change square '{}'", name.text),
            )
        })
    }

    pub(crate) fn check_new_bc_comparison(&self, name: &Label) -> Result<()> {
        if self.bc_comparisons.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("base-change comparison '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn define_bc_comparison(
        &mut self,
        name: String,
        value: crate::sheaf_base_change::Comparison,
    ) {
        self.bc_comparisons.insert(name, Rc::new(value));
    }

    pub(crate) fn bc_comparison(
        &self,
        name: &Label,
    ) -> Result<Rc<crate::sheaf_base_change::Comparison>> {
        self.bc_comparisons.get(&name.text).cloned().ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown base-change comparison '{}'", name.text),
            )
        })
    }

    pub(crate) fn define_sheaf_map(&mut self, name: String, map: sheaf_hyper::Map) {
        self.define_sheaf_map_shared(name, Rc::new(map));
    }

    pub(crate) fn define_sheaf_map_shared(&mut self, name: String, map: Rc<sheaf_hyper::Map>) {
        self.sheaf_maps.insert(name, map);
    }

    pub(crate) fn sheaf_map(&self, name: &Label) -> Result<Rc<sheaf_hyper::Map>> {
        self.sheaf_maps
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown sheaf map '{}'", name.text)))
    }

    pub(crate) fn check_new_sheaf_complex(&self, name: &Label) -> Result<()> {
        if self.sheaf_complexes.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("sheaf complex '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn define_sheaf_complex(&mut self, name: String, complex: sheaf_hyper::Complex) {
        self.define_sheaf_complex_shared(name, Rc::new(complex));
    }

    pub(crate) fn define_sheaf_complex_shared(
        &mut self,
        name: String,
        complex: Rc<sheaf_hyper::Complex>,
    ) {
        self.sheaf_complexes.insert(name, Total::Roos(complex));
    }

    pub(crate) fn define_serre_shift(&mut self, name: String, shift: usize) {
        self.serre_shifts.insert(name, shift);
    }

    pub(crate) fn serre_shift(&self, name: &Label) -> Result<usize> {
        self.serre_shifts.get(&name.text).copied().ok_or_else(|| {
            Error::new(
                name.span,
                format!("'{}' is not a declared Serre complex", name.text),
            )
        })
    }

    pub(crate) fn define_leray_complex(&mut self, name: String, leray: Leray) {
        self.leray_complexes.insert(name, Rc::new(leray));
    }

    pub(crate) fn leray_complex(&self, name: &Label) -> Result<Rc<Leray>> {
        self.leray_complexes
            .get(&name.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("'{}' is not a declared Leray complex", name.text),
                )
            })
    }

    pub(crate) fn sheaf_complex(&self, name: &Label) -> Result<Rc<sheaf_hyper::Complex>> {
        match self.sheaf_total(name)? {
            Total::Roos(complex) => Ok(complex),
            Total::Cellular(_) => Err(Error::new(name.span, format!(
                "'{}' is a cellular total complex; this operation requires an ordinary Roos complex", name.text
            ))),
        }
    }

    pub(crate) fn sheaf_total(&self, name: &Label) -> Result<Total> {
        self.sheaf_total_named(&name.text)
            .ok_or_else(|| Error::new(name.span, format!("unknown sheaf complex '{}'", name.text)))
    }

    pub(crate) fn sheaf_total_named(&self, name: &str) -> Option<Total> {
        self.sheaf_complexes.get(name).cloned()
    }

    pub(crate) fn define_sheaf_cellular(&mut self, name: String, complex: Cellular) {
        self.sheaf_complexes
            .insert(name, Total::Cellular(Rc::new(complex)));
    }

    pub(crate) fn define_sheaf_resolution(&mut self, name: &Label, resolution: Resolution) {
        self.sheaf_resolutions
            .insert(name.text.clone(), Rc::new(resolution));
    }

    pub(crate) fn check_new_sheaf_resolution(&self, name: &Label) -> Result<()> {
        if self.sheaf_resolutions.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("sheaf resolution '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn sheaf_resolution(&self, name: &Label) -> Result<Rc<Resolution>> {
        self.sheaf_resolutions
            .get(&name.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown sheaf resolution '{}'", name.text),
                )
            })
    }

    pub(crate) fn sheaf_complex_type(
        &mut self,
        form: SheafForm,
        name: &Label,
        degree: usize,
    ) -> Result<TypeId> {
        self.sheaf_total(name)?;
        self.known_sheaf_complex_type(form, name.text.as_str().into(), degree)
    }

    pub(crate) fn known_sheaf_complex_type(
        &mut self,
        form: SheafForm,
        name: Rc<str>,
        degree: usize,
    ) -> Result<TypeId> {
        self.inhabited_type(Kind::SheafComplex(form, name, degree), 0)
    }

    pub(crate) fn define_sheaf_supports(&mut self, name: String, data: Supports) {
        self.sheaf_supports.insert(name, Rc::new(data));
    }

    pub(crate) fn sheaf_supports(&self, name: &Label) -> Result<Rc<Supports>> {
        self.sheaf_supports_named(&name.text).ok_or_else(|| {
            Error::new(
                name.span,
                format!(
                    "'{}' is not a declared extraordinary inverse image",
                    name.text
                ),
            )
        })
    }

    pub(crate) fn sheaf_supports_named(&self, name: &str) -> Option<Rc<Supports>> {
        self.sheaf_supports.get(name).cloned()
    }

    pub(crate) fn sheaf_local_type(
        &mut self,
        name: Rc<str>,
        grade: sheaf_support::Grade,
    ) -> Result<TypeId> {
        self.inhabited_type(Kind::SheafLocal(name, grade), 0)
    }

    pub(crate) fn sheaf_page_type(
        &mut self,
        name: &Label,
        grade: sheaf_page::Grade,
    ) -> Result<TypeId> {
        self.sheaf_total(name)?;
        self.known_sheaf_page_type(name.text.as_str().into(), grade)
    }

    pub(crate) fn known_sheaf_page_type(
        &mut self,
        name: Rc<str>,
        grade: sheaf_page::Grade,
    ) -> Result<TypeId> {
        self.inhabited_type(Kind::SheafPage(name, grade), 0)
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

    pub(crate) fn define_set_family(
        &mut self,
        name: &Label,
        diagram: set_family::Diagram,
    ) -> Result<()> {
        self.check_new_set_family(name)?;
        self.allocate(1)?;
        self.set_families
            .insert(Rc::from(name.text.as_str()), Rc::new(diagram));
        Ok(())
    }

    pub(crate) fn define_dependent_set(
        &mut self,
        name: &Label,
        diagram: set_family::Diagram,
        origin: dependent_sets::Origin,
    ) -> Result<()> {
        self.allocate(1)?;
        self.define_set_family(name, diagram)?;
        self.set_origins
            .insert(Rc::from(name.text.as_str()), origin);
        Ok(())
    }

    pub(crate) fn define_set_recollement(
        &mut self,
        name: &Label,
        data: crate::set_gluing::Recollement,
    ) -> Result<()> {
        if self.set_recollements.contains_key(name.text.as_str()) {
            return Err(Error::new(
                name.span,
                format!("set recollement '{}' is already defined", name.text),
            ));
        }
        self.allocate(1)?;
        self.set_recollements
            .insert(Rc::from(name.text.as_str()), Rc::new(data));
        Ok(())
    }

    pub(crate) fn set_recollement(
        &self,
        name: &Label,
    ) -> Result<Rc<crate::set_gluing::Recollement>> {
        self.set_recollements
            .get(name.text.as_str())
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown set recollement '{}'", name.text),
                )
            })
    }

    pub(crate) fn define_set_topology(
        &mut self,
        name: &Label,
        topology: set_topology::Topology,
    ) -> Result<()> {
        if self.set_topologies.contains_key(name.text.as_str()) {
            return Err(Error::new(
                name.span,
                format!("set topology '{}' is already defined", name.text),
            ));
        }
        self.allocate(1)?;
        self.set_topologies
            .insert(name.text.as_str().into(), Rc::new(topology));
        Ok(())
    }

    pub(crate) fn set_topology(&self, name: &str) -> Result<Rc<set_topology::Topology>> {
        self.set_topologies
            .get(name)
            .cloned()
            .ok_or_else(|| Error::plain(format!("unknown set topology '{name}'")))
    }

    pub(crate) fn define_set_reflection(
        &mut self,
        name: &Label,
        diagram: set_family::Diagram,
        reflection: set_topology::Reflection,
    ) -> Result<()> {
        self.allocate(1)?;
        self.define_set_family(name, diagram)?;
        self.set_reflections
            .insert(name.text.as_str().into(), Rc::new(reflection));
        Ok(())
    }

    pub(crate) fn set_reflection(&self, name: &str) -> Result<Rc<set_topology::Reflection>> {
        self.set_reflections
            .get(name)
            .cloned()
            .ok_or_else(|| Error::plain(format!("'{name}' is not a declared topology reflection")))
    }

    pub(crate) fn set_reflection_certificate(
        &mut self,
        topology: Rc<str>,
        family: Rc<str>,
    ) -> Result<Option<Rc<set_topology::Certificate>>> {
        let key = (topology, family);
        if let Some(certificate) = self.set_reflection_certificates.get(&key) {
            return Ok(certificate.clone());
        }
        let data = self.set_topology(&key.0)?;
        let diagram = self.set_family_named(&key.1).expect("registered family");
        let certificate = set_topology::Certificate::build(self, &data, &diagram)?.map(Rc::new);
        self.allocate(1)?;
        self.set_reflection_certificates
            .insert(key, certificate.clone());
        Ok(certificate)
    }

    pub(crate) fn define_set_sheafification(
        &mut self,
        name: &Label,
        diagram: set_family::Diagram,
        origin: set_j::Sheafification,
    ) -> Result<()> {
        self.allocate(1)?;
        self.define_set_family(name, diagram)?;
        self.set_sheafifications
            .insert(Rc::from(name.text.as_str()), Rc::new(origin));
        Ok(())
    }

    pub(crate) fn set_sheafification(&self, name: &str) -> Result<Rc<set_j::Sheafification>> {
        self.set_sheafifications.get(name).cloned().ok_or_else(|| {
            Error::plain(format!(
                "'{name}' is not a declared double-negation sheafification"
            ))
        })
    }

    pub(crate) fn set_sheaf_certificate(
        &mut self,
        name: Rc<str>,
    ) -> Result<Option<Rc<set_j::Certificate>>> {
        if let Some(certificate) = self.set_sheaf_certificates.get(&name) {
            return Ok(certificate.clone());
        }
        let family = self.set_family_named(&name).expect("registered family");
        let certificate = set_j::Certificate::build(self, &family)?.map(Rc::new);
        self.allocate(1)?;
        self.set_sheaf_certificates
            .insert(name, certificate.clone());
        Ok(certificate)
    }

    pub(crate) fn define_set_logic(
        &mut self,
        name: &Label,
        diagram: set_family::Diagram,
        origin: set_logic::Origin,
    ) -> Result<()> {
        self.allocate(1)?;
        self.define_set_family(name, diagram)?;
        self.set_logic.insert(Rc::from(name.text.as_str()), origin);
        Ok(())
    }

    pub(crate) fn set_omega(&self, name: &str) -> Result<Rc<set_logic::Omega>> {
        if let Some(set_logic::Origin::Omega(omega)) = self.set_logic.get(name) {
            return Ok(omega.clone());
        }
        Err(Error::plain(format!(
            "'{name}' is not a declared set truth classifier"
        )))
    }

    pub(crate) fn set_subobject(&self, name: &str) -> Result<Rc<set_logic::Subobject>> {
        if let Some(set_logic::Origin::Subobject(subobject)) = self.set_logic.get(name) {
            return Ok(subobject.clone());
        }
        Err(Error::plain(format!(
            "'{name}' is not a declared predicate comprehension"
        )))
    }

    pub(crate) fn set_origin(&self, name: &str) -> Option<dependent_sets::Origin> {
        self.set_origins.get(name).cloned()
    }

    pub(crate) fn check_new_set_family(&self, name: &Label) -> Result<()> {
        if self.set_families.contains_key(name.text.as_str()) {
            return Err(Error::new(
                name.span,
                format!("set family '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn set_family(&self, name: &Label) -> Result<Rc<set_family::Diagram>> {
        self.set_family_named(&name.text)
            .ok_or_else(|| Error::new(name.span, format!("unknown set family '{}'", name.text)))
    }

    pub(crate) fn set_family_named(&self, name: &str) -> Option<Rc<set_family::Diagram>> {
        self.set_families.get(name).cloned()
    }

    pub(crate) fn check_new_set_monodromy(&self, name: &Label) -> Result<()> {
        if self.set_monodromies.contains_key(name.text.as_str()) {
            return Err(Error::new(
                name.span,
                format!("set monodromy '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn define_set_monodromy(
        &mut self,
        name: &Label,
        action: crate::set_monodromy::Monodromy,
    ) -> Result<()> {
        self.allocate(1)?;
        self.set_monodromies
            .insert(Rc::from(name.text.as_str()), Rc::new(action));
        Ok(())
    }

    pub(crate) fn set_monodromy(
        &self,
        name: &Label,
    ) -> Result<Rc<crate::set_monodromy::Monodromy>> {
        self.set_monodromies
            .get(name.text.as_str())
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown set monodromy '{}'", name.text)))
    }

    pub(crate) fn define_set_linearization(
        &mut self,
        name: Label,
        family: &Label,
        sheaf: Sheaf,
    ) -> Result<()> {
        self.allocate(1)?;
        self.set_linearizations
            .insert(Rc::from(name.text.as_str()), Rc::from(family.text.as_str()));
        self.define_sheaf(name.text, sheaf);
        Ok(())
    }

    pub(crate) fn require_set_linearization(&self, name: &Label, family: &str) -> Result<()> {
        if self
            .set_linearizations
            .get(name.text.as_str())
            .is_none_or(|f| f.as_ref() != family)
        {
            return Err(Error::new(
                name.span,
                format!(
                    "'{}' must be a linearization of set family '{family}'",
                    name.text
                ),
            ));
        }
        Ok(())
    }

    pub(crate) fn set_elements(&self, name: &Label) -> Result<(TypeId, Rc<set_family::Elements>)> {
        self.named(&name.text)
            .and_then(|ty| self.set_elements.get(&ty).map(|data| (ty, data.clone())))
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown set elements space '{}'", name.text),
                )
            })
    }

    pub(crate) fn define_set_elements(&mut self, name: &Label, family: &Label) -> Result<()> {
        let diagram = self.set_family(family)?;
        let (ty, data) =
            set_family::Elements::build(self, name, Rc::from(family.text.as_str()), &diagram)?;
        self.set_elements.insert(ty, Rc::new(data));
        Ok(())
    }

    pub(crate) fn set_nat_type(&mut self, source: Rc<str>, target: Rc<str>) -> Result<TypeId> {
        self.register_set_nat_type(source, target, false)
    }

    pub(crate) fn register_set_nat_type(
        &mut self,
        source: Rc<str>,
        target: Rc<str>,
        validated_witness: bool,
    ) -> Result<TypeId> {
        let kind = Kind::SetNat(source.clone(), target.clone());
        if let Some(&ty) = self.type_ids.get(&kind) {
            return Ok(ty);
        }
        self.check_type_limit(0)?;
        let a = self
            .set_family_named(&source)
            .expect("registered source family");
        let b = self
            .set_family_named(&target)
            .expect("registered target family");
        if a.base != b.base {
            return Err(Error::plain(
                "SetNat requires two families on the same named site",
            ));
        }
        let inhabited = validated_witness || source == target || a.has_natural(self, &b)?;
        Ok(self.insert(kind, None, None, inhabited, 0))
    }

    pub(crate) fn check_named_space(&self, name: &Label, count: usize) -> Result<()> {
        if self.q_endpoint_named(&name.text).is_some() {
            return Err(Error::new(
                name.span,
                format!("space '{}' is already defined", name.text),
            ));
        }
        if count > MAX_SPACE_POINTS {
            return Err(Error::new(
                name.span,
                "space exceeds the limit of 512 points",
            ));
        }
        self.check_type_limit(0)
    }

    // Caller certifies the order and accounts for the constructed table.
    pub(crate) fn insert_derived_space(
        &mut self,
        name: &Label,
        names: Vec<String>,
        order: Vec<Vec<bool>>,
    ) -> TypeId {
        let inhabited = !names.is_empty();
        self.insert(
            Kind::Space(name.text.clone()),
            Some(Poset::new(names, order)),
            None,
            inhabited,
            0,
        )
    }

    pub(crate) fn define_space(&mut self, name: &Label, chains: &[Vec<Label>]) -> Result<()> {
        self.check_named_space(name, 0)?;

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
        self.insert_derived_space(name, names, order);
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

    pub(crate) fn mcobar_type(&mut self, grade: motivic::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::MCobar(grade), 0)
    }

    pub(crate) fn real_rho_cobar_type(&mut self, grade: real_rho::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RMRhoCobar(grade), 0)
    }

    pub(crate) fn real_rho_ext_type(&mut self, grade: real_rho::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RMRhoExt(grade), 0)
    }

    pub(crate) fn real_cobar_type(&mut self, grade: motivic::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RMCobar(grade), 0)
    }

    pub(crate) fn rho_min_type(&mut self, grade: motivic::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RhoMin(grade), 0)
    }

    pub(crate) fn rho_loc_type(&mut self, grade: motivic::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RhoLoc(grade), 0)
    }

    pub(crate) fn rho_page_type(&mut self, grade: rho_bock::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RhoBock(grade), 0)
    }

    pub(crate) fn rho_page(&mut self, ty: TypeId) -> Result<Rc<rho_bock::Page>> {
        if let Some(page) = self.rho_pages.get(&ty) {
            return Ok(page.clone());
        }
        let Kind::RhoBock(grade) = *self.kind(ty) else {
            unreachable!();
        };
        let page = Rc::new(rho_bock::Page::build(self, grade)?);
        self.rho_pages.insert(ty, page.clone());
        Ok(page)
    }

    pub(crate) fn rmbar_type(&mut self, grade: motivic::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RMBar(grade), 0)
    }

    pub(crate) fn rmext_type(&mut self, grade: motivic::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RMExt(grade), 0)
    }

    pub(crate) fn real_bars(&mut self, s: usize, k: i128) -> Result<Rc<real_bars::Group>> {
        if let Some(group) = self.real_bars.get(&(s, k)) {
            return Ok(group.clone());
        }
        let group = Rc::new(real_bars::Group::build(self, s, k)?);
        self.real_bars.insert((s, k), group.clone());
        Ok(group)
    }

    pub(crate) fn rm_triple_type(&mut self, grade: motivic::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RMTriple(grade), 0)
    }

    pub(crate) fn has_real_ext(&self, ty: TypeId) -> bool {
        self.real_ext.contains_key(&ty)
    }

    pub(crate) fn real_ext(&mut self, ty: TypeId) -> Result<Rc<real_ext::Group>> {
        if let Some(group) = self.real_ext.get(&ty) {
            return Ok(group.clone());
        }
        let group = Rc::new(match *self.kind(ty) {
            Kind::RMExt(grade) => {
                let bars = grade
                    .t
                    .checked_sub(i128::from(grade.w))
                    .and_then(|k| self.real_bars.get(&(grade.s, k)))
                    .cloned();
                real_ext::Group::build(self, grade, None, bars.as_deref())?
            }
            Kind::RMRhoExt(grade) => {
                real_ext::Group::build(self, grade.chain, Some(grade.exponent), None)?
            }
            _ => unreachable!(),
        });
        self.real_ext.insert(ty, group.clone());
        Ok(group)
    }

    pub(crate) fn mext_type(&mut self, grade: motivic::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::MExt(grade), 0)
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

    pub(crate) fn rv_cyclic_type(
        &mut self,
        grade: rv_cyclic_homology::Grade,
        homology: bool,
    ) -> Result<TypeId> {
        let kind = if homology {
            Kind::RvCyclicHomology(grade)
        } else {
            Kind::RvCyclicChain(grade)
        };
        self.inhabited_type(kind, 0)
    }

    pub(crate) fn rv_cyclic_homology(
        &mut self,
        ty: TypeId,
    ) -> Result<Rc<rv_hh_homology::Slice<rv_negative::Term>>> {
        if let Some(homology) = self.rv_cyclic_homology.get(&ty) {
            return Ok(homology.clone());
        }
        let Kind::RvCyclicHomology(grade) = self.kind(ty) else {
            unreachable!()
        };
        let grade = grade.clone();
        let homology = Rc::new(rv_cyclic_homology::build(self, &grade)?);
        self.rv_cyclic_homology.insert(ty, homology.clone());
        Ok(homology)
    }

    pub(crate) fn rv_negative_type(
        &mut self,
        grade: rv_negative::Grade,
        homology: bool,
    ) -> Result<TypeId> {
        let kind = if homology {
            Kind::RvNegativeCyclic(grade)
        } else {
            Kind::RVNCChain(grade)
        };
        self.inhabited_type(kind, 0)
    }

    pub(crate) fn rv_negative_homology(
        &mut self,
        ty: TypeId,
    ) -> Result<Rc<rv_hh_homology::Slice<rv_negative::Term>>> {
        if let Some(homology) = self.rv_negative_homology.get(&ty) {
            return Ok(homology.clone());
        }
        let Kind::RvNegativeCyclic(grade) = self.kind(ty) else {
            unreachable!()
        };
        let grade = grade.clone();
        let homology = Rc::new(rv_negative::build(self, &grade)?);
        self.rv_negative_homology.insert(ty, homology.clone());
        Ok(homology)
    }

    pub(crate) fn rv_hh_class_type(&mut self, grade: rv_hh::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RvHochschild(grade), 0)
    }

    pub(crate) fn rv_hh_homology(
        &mut self,
        ty: TypeId,
    ) -> Result<Rc<rv_hh_homology::Slice<rv_hh::Word>>> {
        if let Some(homology) = self.rv_hh_homology.get(&ty) {
            return Ok(homology.clone());
        }
        let Kind::RvHochschild(grade) = self.kind(ty) else {
            unreachable!()
        };
        let grade = grade.clone();
        let homology = Rc::new(rv_hh_homology::build(self, &grade)?);
        self.rv_hh_homology.insert(ty, homology.clone());
        Ok(homology)
    }

    pub(crate) fn rv_hh_type(&mut self, grade: rv_hh::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RVHHChain(grade), 0)
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

    pub(crate) fn wzw_type(&mut self, grade: rv_wzw::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RVWzw(grade), 0)
    }

    pub(crate) fn define_kh_braid(
        &mut self,
        name: &Label,
        strands: usize,
        word: &[i128],
    ) -> Result<()> {
        if self.kh_braids.contains_key(name.text.as_str()) {
            return Err(Error::new(name.span, "Khovanov braid already defined"));
        }
        let braid = khovanov::Braid::new(self, strands, word)?;
        self.kh_braids
            .insert(name.text.as_str().into(), Rc::new(braid));
        Ok(())
    }

    pub(crate) fn kh_braid(&self, name: &str) -> Result<Rc<khovanov::Braid>> {
        self.kh_braids
            .get(name)
            .cloned()
            .ok_or_else(|| Error::plain(format!("unknown Khovanov braid '{name}'")))
    }

    fn lee_elements(
        &mut self,
        braid: &khovanov::Braid,
        i: i64,
        f: Option<i64>,
    ) -> Result<Vec<Value>> {
        let f = f.unwrap_or_else(|| braid.lee_lower_bound());
        if braid.lee_possible(i, f)
            && khovanov::Context::new(braid)
                .lee_basis(self, i)?
                .first()
                .is_some_and(|(_, j)| *j >= f)
        {
            return Err(Error::plain(
                "cannot enumerate a nonzero rational Lee cochain group",
            ));
        }
        self.allocate(1)?;
        Ok(vec![Value::LeeChain(Rc::default())])
    }

    pub(crate) fn lee_type(&mut self, braid: Rc<str>, grade: crate::lee::Grade) -> Result<TypeId> {
        self.kh_braid(&braid)?;
        self.inhabited_type(Kind::Lee(braid, grade), 0)
    }

    pub(crate) fn lee_dual_type(&mut self, braid: Rc<str>, i: i64) -> Result<TypeId> {
        self.kh_braid(&braid)?;
        self.inhabited_type(Kind::LeeDual(braid, i), 0)
    }

    pub(crate) fn kh_type(
        &mut self,
        braid: Rc<str>,
        grade: khovanov::Grade,
        form: KhForm,
    ) -> Result<TypeId> {
        self.kh_braid(&braid)?;
        self.inhabited_type(Kind::Khovanov(braid, grade, form), 0)
    }

    pub(crate) fn tl_type(&mut self, grade: rv_tl::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RvTl(grade), 0)
    }

    pub(crate) fn cyclotomic_rational_type(&mut self, order: usize) -> Result<TypeId> {
        rv_cyclotomic::check_order(order)?;
        self.inhabited_type(Kind::RVCycloQ(order), 0)
    }

    pub(crate) fn cyclotomic_type(&mut self, order: usize) -> Result<TypeId> {
        rv_cyclotomic::check_order(order)?;
        self.inhabited_type(Kind::RVCyclo(order), 0)
    }

    pub(crate) fn cyclotomic(&mut self, order: usize) -> Result<rv_cyclotomic::Polynomial> {
        if let Some(modulus) = self.cyclotomics.get(&order) {
            return Ok(modulus.clone());
        }
        let modulus = rv_cyclotomic::modulus(self, order)?;
        self.allocate(2 * modulus.len())?;
        self.cyclotomics.insert(order, modulus.clone());
        Ok(modulus)
    }

    pub(crate) fn fusion_type(&mut self, grade: rv_fusion::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RVFusion(grade), 0)
    }

    pub(crate) fn ore_type(&mut self) -> Result<TypeId> {
        self.inhabited_type(Kind::RVOre, 0)
    }

    pub(crate) fn witt_spectrum_type(&mut self) -> Result<TypeId> {
        self.inhabited_type(Kind::RVWittSpectrum, 0)
    }

    pub(crate) fn schur_character_type(&mut self, weight: u128) -> Result<TypeId> {
        self.inhabited_type(Kind::RVSchurCharacter(weight), 0)
    }

    pub(crate) fn stable_schur_character_type(&mut self) -> Result<TypeId> {
        self.inhabited_type(Kind::RVStableSchurCharacter, 0)
    }

    pub(crate) fn partition_algebra_type(&mut self, rank: usize) -> Result<TypeId> {
        self.inhabited_type(Kind::RVPartitionAlgebra(rank), 0)
    }

    pub(crate) fn partition_local_type(&mut self, rank: usize) -> Result<TypeId> {
        self.inhabited_type(Kind::RVPartitionLocal(rank), 0)
    }

    pub(crate) fn karoubi_map_type(&mut self, source: &Label, target: &Label) -> Result<TypeId> {
        let a = self.karoubi_object(source)?;
        let b = self.karoubi_object(target)?;
        if a.rank() != b.rank() {
            return Err(Error::new(
                target.span,
                "Karoubi maps require retracts of the same partition rank",
            ));
        }
        self.inhabited_type(
            Kind::RVKaroubiMap(source.text.as_str().into(), target.text.as_str().into()),
            0,
        )
    }

    pub(crate) fn quantum_mod_type(&mut self, grade: rv_pcurv::ModGrade) -> Result<TypeId> {
        self.inhabited_type(Kind::RVQGrP(grade), 0)
    }

    pub(crate) fn quantum_grassmannian_type(&mut self, grade: rv_qgr::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::RVQGr(grade), 0)
    }

    pub(crate) fn hive_type(&mut self, boundary: HiveBoundary) -> Result<TypeId> {
        let kind = Kind::RVHive(boundary);
        if let Some(&ty) = self.type_ids.get(&kind) {
            return Ok(ty);
        }
        self.check_type_limit(0)?;
        let Kind::RVHive(boundary) = &kind else {
            unreachable!()
        };
        let witness = rv_hive::witness(self, boundary)?;
        self.allocate(3 * boundary.rank() + witness.as_ref().map_or(0, |point| point.len()))?;
        let ty = self.insert(kind, None, None, witness.is_some(), 0);
        self.types[ty.0].hive_witness = witness;
        Ok(ty)
    }

    pub(crate) fn hive_witness(&self, ty: TypeId) -> Result<Value> {
        self.types[ty.0]
            .hive_witness
            .clone()
            .map(Value::RVHive)
            .ok_or_else(|| Error::plain("hive boundary has no integral witness"))
    }

    pub(crate) fn satake_type(&mut self, rank: usize) -> Result<TypeId> {
        self.inhabited_type(Kind::RVSatake(rank), 0)
    }

    pub(crate) fn sph_hecke_type(&mut self, rank: usize) -> Result<TypeId> {
        self.inhabited_type(Kind::RVSphHecke(rank), 0)
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

    pub(crate) fn bp_mod_ext_type(&mut self, grade: bp_mod_ext::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::BPModExt(grade), 0)
    }

    pub(crate) fn bp_vpage_type(&mut self, grade: bp_vpage::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::BPVPage(grade), 0)
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

    pub(crate) fn ec_curve(&self, name: &str) -> Result<Rc<EllipticCurve>> {
        self.elliptic_curves
            .get(name)
            .cloned()
            .ok_or_else(|| Error::plain(format!("unknown elliptic curve '{name}'")))
    }

    pub(crate) fn ec_isogeny(&self, name: &str) -> Result<Rc<Isogeny>> {
        self.elliptic_isogenies
            .get(name)
            .cloned()
            .ok_or_else(|| Error::plain(format!("unknown cyclic isogeny '{name}'")))
    }

    pub(crate) fn define_ec_torsion_level(
        &mut self,
        name: Label,
        build: impl FnOnce(&mut Self) -> Result<ec_level::Level>,
    ) -> Result<()> {
        if self.ec_torsion_levels.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                "elliptic torsion level is already defined",
            ));
        }
        let level = build(self)?;
        self.ec_torsion_levels.insert(name.text, Rc::new(level));
        Ok(())
    }

    pub(crate) fn ec_torsion_level(&self, name: &str) -> Result<Rc<ec_level::Level>> {
        self.ec_torsion_levels
            .get(name)
            .cloned()
            .ok_or_else(|| Error::plain(format!("unknown elliptic torsion level '{name}'")))
    }

    pub(crate) fn ec_level_type(&mut self, name: &Label, matrix: bool) -> Result<TypeId> {
        self.ec_torsion_level(&name.text)
            .map_err(|error| error.at(name.span))?;
        self.inhabited_type(Kind::ECLevel(name.text.as_str().into(), matrix), 0)
    }

    pub(crate) fn define_ec_cm_action(
        &mut self,
        name: Label,
        build: impl FnOnce(&mut Self) -> Result<crate::elliptic_end_order::realization::Realization>,
    ) -> Result<()> {
        if self.ec_cm_actions.contains_key(&name.text) {
            return Err(Error::new(name.span, "CM action is already defined"));
        }
        let action = build(self)?;
        self.ec_cm_actions.insert(name.text, Rc::new(action));
        Ok(())
    }

    pub(crate) fn ec_cm_action(
        &self,
        name: &Label,
    ) -> Result<Rc<crate::elliptic_end_order::realization::Realization>> {
        self.ec_cm_actions
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown CM action '{}'", name.text)))
    }

    pub(crate) fn define_ec_cm_cycle(
        &mut self,
        name: Label,
        address: crate::elliptic_end_order::horizontal::Address,
    ) -> Result<()> {
        if self.ec_cm_cycles.contains_key(&name.text) {
            return Err(Error::new(name.span, "CM isogeny cycle is already defined"));
        }
        let cycle = crate::elliptic_end_order::cycle::Cycle::new(self, &name, address)?;
        self.ec_cm_cycles.insert(name.text, Rc::new(cycle));
        Ok(())
    }

    pub(crate) fn ec_cm_cycle(
        &self,
        name: &Label,
    ) -> Result<Rc<crate::elliptic_end_order::cycle::Cycle>> {
        self.ec_cm_cycles.get(&name.text).cloned().ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown CM isogeny cycle '{}'", name.text),
            )
        })
    }

    pub(crate) fn define_ec_end_order(
        &mut self,
        name: Label,
        certificate: crate::elliptic_end_order::Certificate,
    ) -> Result<()> {
        if self.elliptic_end_orders.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                "elliptic endomorphism-order certificate is already defined",
            ));
        }
        self.elliptic_end_orders
            .insert(name.text, Rc::new(certificate));
        Ok(())
    }

    pub(crate) fn ec_end_order(
        &self,
        name: &str,
    ) -> Result<Rc<crate::elliptic_end_order::Certificate>> {
        self.elliptic_end_orders.get(name).cloned().ok_or_else(|| {
            Error::plain(format!(
                "unknown elliptic endomorphism-order certificate '{name}'"
            ))
        })
    }

    pub(crate) fn ec_end_order_type(&mut self, name: &Label) -> Result<TypeId> {
        self.ec_end_order(&name.text)?;
        self.inhabited_type(Kind::ECEndOrder(name.text.as_str().into()), 0)
    }

    pub(crate) fn define_ec_curve(
        &mut self,
        name: Label,
        ring: Label,
        a: ECScalar,
        b: ECScalar,
    ) -> Result<()> {
        if self.elliptic_curves.contains_key(name.text.as_str()) {
            return Err(Error::new(name.span, "elliptic curve is already defined"));
        }
        let field = self.rv_unram(&ring)?;
        let curve = EllipticCurve::new(self, ring.text.into(), field, a, b)?;
        self.elliptic_curves
            .insert(name.text.into(), Rc::new(curve));
        Ok(())
    }

    pub(crate) fn define_ec_isogeny(
        &mut self,
        name: Label,
        source: Label,
        target: Label,
        kernel: Rc<ECKernel>,
    ) -> Result<()> {
        self.define_ec_isogeny_checked(name, source, target, kernel, |_, _, _| Ok(()))
    }

    pub(crate) fn define_ec_isogeny_checked(
        &mut self,
        name: Label,
        source: Label,
        target: Label,
        kernel: Rc<ECKernel>,
        check: impl FnOnce(&mut Model, &EllipticCurve, &EllipticCurve) -> Result<()>,
    ) -> Result<()> {
        if self.elliptic_isogenies.contains_key(name.text.as_str()) {
            return Err(Error::new(name.span, "cyclic isogeny is already defined"));
        }
        if self.elliptic_curves.contains_key(target.text.as_str()) {
            return Err(Error::new(
                target.span,
                "isogeny target name is already defined",
            ));
        }
        let curve = self.ec_curve(&source.text)?;
        let target: Rc<str> = target.text.into();
        let (isogeny, quotient) =
            Isogeny::new(self, &curve, source.text.into(), target.clone(), kernel)?;
        check(self, &curve, &quotient)?;
        self.elliptic_curves.insert(target, Rc::new(quotient));
        self.elliptic_isogenies
            .insert(name.text.into(), Rc::new(isogeny));
        Ok(())
    }

    pub(crate) fn ec_curve_type(&mut self, name: Rc<str>, point: bool) -> Result<TypeId> {
        self.ec_curve(&name)?;
        let kind = if point {
            Kind::ECPoint(name)
        } else {
            Kind::ECCurve(name)
        };
        self.inhabited_type(kind, 0)
    }

    pub(crate) fn ec_kernel_type(
        &mut self,
        name: Rc<str>,
        ell: u64,
        witness: bool,
    ) -> Result<TypeId> {
        let kind = Kind::ECKernel(name.clone(), ell);
        if let Some(&id) = self.type_ids.get(&kind) {
            return Ok(id);
        }
        let curve = self.ec_curve(&name)?;
        let inhabited = witness || curve.has_kernel(self, ell)?;
        Ok(self.insert(kind, None, None, inhabited, 0))
    }

    pub(crate) fn ec_map_type(&mut self, source: Rc<str>, target: Rc<str>) -> Result<TypeId> {
        if self.ec_curve(&source)?.ring_name != self.ec_curve(&target)?.ring_name {
            return Err(Error::plain(
                "elliptic morphisms require the same named residue field",
            ));
        }
        self.inhabited_type(Kind::ECMap(source, target), 0)
    }

    pub(crate) fn ec_push_type(
        &mut self,
        target: Rc<str>,
        ell: u64,
        input_kernel: TypeId,
    ) -> Result<TypeId> {
        let witness = self.types[input_kernel.0].inhabited;
        self.ec_kernel_type(target, ell, witness)
    }

    pub(crate) fn ec_isogeny_type(&mut self, name: Rc<str>) -> Result<TypeId> {
        self.ec_isogeny(&name)?;
        self.inhabited_type(Kind::ECIso(name), 0)
    }

    pub(crate) fn ec_coefficient_type(&mut self, curve: &EllipticCurve) -> Result<TypeId> {
        self.inhabited_type(Kind::RVUnram(curve.ring_name.clone()), 0)
    }

    pub(crate) fn show_ec_kernels(&mut self, name: &Label, ell: usize) -> Result<String> {
        let curve = self.ec_curve(&name.text)?;
        let kernels = elliptic::enumerate_kernels(self, &curve, ell)?;
        let mut out = Text::new(MAX_SHOWN_BYTES);
        out.push(&format!(
            "rvECIsogenies[{}, {ell}] = count={}, kernels=[",
            name.text,
            kernels.len()
        ))?;
        for (index, kernel) in kernels.iter().enumerate() {
            if index > 0 {
                out.push(", ")?;
            }
            self.write_ec_polynomial(&mut out, kernel)?;
        }
        out.push("]")?;
        Ok(out.string)
    }

    fn write_ec_polynomial(&mut self, out: &mut Text, polynomial: &[ECScalar]) -> Result<()> {
        out.push("[")?;
        for (i, s) in polynomial.iter().enumerate() {
            self.spend(1 + s.len())?;
            if i > 0 {
                out.push(", ")?;
            }
            out.push(&format!("{s:?}"))?;
        }
        out.push("]")
    }

    pub(crate) fn define_hyperelliptic(
        &mut self,
        name: Label,
        prime: usize,
        lift: Vec<i128>,
    ) -> Result<()> {
        self.require_fresh_hyperelliptic_name(&name)?;
        let curve = hyperelliptic::Curve::new(self, prime, lift)?;
        self.hyperelliptic.insert(name.text, Rc::new(curve));
        Ok(())
    }

    fn require_fresh_hyperelliptic_name(&self, name: &Label) -> Result<()> {
        if self.hyperelliptic.contains_key(&name.text)
            || self.hw_closed_curves.contains_key(&name.text)
        {
            return Err(Error::new(
                name.span,
                "hyperelliptic curve is already defined",
            ));
        }
        Ok(())
    }

    pub(crate) fn hyperelliptic(&self, name: &Label) -> Result<Rc<hyperelliptic::Curve>> {
        self.hyperelliptic.get(&name.text).cloned().ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown hyperelliptic curve '{}'", name.text),
            )
        })
    }

    pub(crate) fn hw_matrix_type(&mut self, ring: Residues, size: usize) -> Result<TypeId> {
        ring.matrix_size(size)?;
        self.inhabited_type(Kind::HWMatrix(ring, size), 0)
    }

    pub(crate) fn hw_scalar_type(&mut self, ring: Residues) -> Result<TypeId> {
        self.inhabited_type(Kind::HWScalar(ring), 0)
    }

    pub(crate) fn define_hww_curve(
        &mut self,
        name: Label,
        ring: &Label,
        coefficients: Vec<Vec<u64>>,
    ) -> Result<()> {
        self.require_fresh_hyperelliptic_name(&name)?;
        let coefficient_ring = self.rv_unram(ring)?;
        let curve = hw_closed::Curve::new(
            self,
            ring.text.as_str().into(),
            coefficient_ring,
            coefficients,
        )
        .map_err(|error| error.at(name.span))?;
        self.hw_closed_curves.insert(name.text, Rc::new(curve));
        Ok(())
    }

    pub(crate) fn hww_curve(&self, name: &Label) -> Result<Rc<hw_closed::Curve>> {
        self.hw_closed_curves
            .get(&name.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown unramified hyperelliptic curve '{}'", name.text),
                )
            })
    }

    pub(crate) fn define_jacobian_chart(
        &mut self,
        target: Label,
        source: &Label,
        branch: ECScalar,
    ) -> Result<()> {
        let curve = self.hww_curve(source)?;
        let (chart, polynomial) = Chart::build(self, source.text.as_str().into(), &curve, branch)
            .map_err(|error| error.at(source.span))?;
        let ring = Label {
            text: curve.grade.ring.to_string(),
            span: source.span,
        };
        let name = target.text.clone();
        self.define_hww_curve(target, &ring, polynomial)?;
        self.jacobian_charts.insert(name, Rc::new(chart));
        Ok(())
    }

    pub(crate) fn jacobian_chart(&self, name: &Label) -> Result<Rc<Chart>> {
        self.jacobian_charts
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown Jacobian chart '{}'", name.text)))
    }

    pub(crate) fn check_new_jacobian_level(&self, name: &Label) -> Result<()> {
        if self.jacobian_levels.contains_key(&name.text)
            || self.jacobian_power_levels.contains_key(&name.text)
        {
            return Err(Error::new(
                name.span,
                format!("Jacobian torsion basis '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn define_jacobian_level(
        &mut self,
        name: Label,
        basis: crate::jacobian_level::Basis,
    ) -> Result<()> {
        self.check_new_jacobian_level(&name)?;
        self.allocate(1 + name.text.len())?;
        self.jacobian_levels.insert(name.text, Rc::new(basis));
        Ok(())
    }

    pub(crate) fn define_jacobian_power_level(
        &mut self,
        name: Label,
        basis: crate::jacobian_level::Basis,
    ) -> Result<()> {
        self.check_new_jacobian_level(&name)?;
        self.allocate(1 + name.text.len())?;
        self.jacobian_power_levels.insert(name.text, Rc::new(basis));
        Ok(())
    }

    pub(crate) fn jacobian_power_level(
        &self,
        name: &Label,
    ) -> Result<Rc<crate::jacobian_level::Basis>> {
        if let Some(basis) = self.jacobian_power_levels.get(&name.text) {
            return Ok(Rc::clone(basis));
        }
        self.jacobian_level(name)
    }

    pub(crate) fn jacobian_level(&self, name: &Label) -> Result<Rc<crate::jacobian_level::Basis>> {
        self.jacobian_levels
            .get(&name.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown Jacobian torsion basis '{}'", name.text),
                )
            })
    }

    pub(crate) fn kummer_divisor_type(&mut self, basis: &Label) -> Result<TypeId> {
        crate::hyperelliptic_kummer_surface::polarized::basis(self, basis)?;
        self.inhabited_type(Kind::KummerDiv(basis.text.as_str().into()), 0)
    }

    pub(crate) fn define_k3_automorphism(
        &mut self,
        name: Label,
        word: hyperelliptic_kummer_surface::dynamics::Word,
    ) -> Result<()> {
        if self.k3_automorphisms.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("K3 automorphism word '{}' is already defined", name.text),
            ));
        }
        self.allocate(5 + name.text.len() + word.generators.len())?;
        self.k3_automorphisms.insert(name.text, Rc::new(word));
        Ok(())
    }

    pub(crate) fn define_k3_period(
        &mut self,
        name: Label,
        marking: Rc<str>,
        u: Rc<unimodular::Vector>,
        v: Rc<unimodular::Vector>,
    ) -> Result<()> {
        if self.k3_periods.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("K3 period '{}' is already defined", name.text),
            ));
        }
        let period = hyperelliptic_kummer_surface::period::Period::new(self, marking, u, v)?;
        self.k3_periods.insert(name.text, Rc::new(period));
        Ok(())
    }

    pub(crate) fn k3_period(
        &self,
        name: &Label,
    ) -> Result<Rc<hyperelliptic_kummer_surface::period::Period>> {
        self.k3_periods
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown K3 period '{}'", name.text)))
    }

    pub(crate) fn k3_automorphism(
        &self,
        name: &Label,
    ) -> Result<Rc<hyperelliptic_kummer_surface::dynamics::Word>> {
        self.k3_automorphisms
            .get(&name.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown K3 automorphism word '{}'", name.text),
                )
            })
    }

    pub(crate) fn define_karoubi_object(
        &mut self,
        name: Label,
        projector: Rc<crate::partition_algebra::local::LocalAlgebra>,
    ) -> Result<()> {
        if self.karoubi_objects.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Karoubi object '{}' is already defined", name.text),
            ));
        }
        crate::partition_algebra::karoubi::certify_object(self, &projector)
            .map_err(|error| error.at(name.span))?;
        self.allocate(4 + name.text.len())?;
        self.karoubi_objects.insert(name.text, projector);
        Ok(())
    }

    pub(crate) fn karoubi_object(
        &self,
        name: &Label,
    ) -> Result<Rc<crate::partition_algebra::local::LocalAlgebra>> {
        self.karoubi_objects
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown Karoubi object '{}'", name.text)))
    }

    pub(crate) fn deligne_square(
        &mut self,
    ) -> Result<Rc<crate::partition_algebra::square::Square>> {
        if let Some(square) = &self.deligne_square {
            return Ok(square.clone());
        }
        let square = Rc::new(crate::partition_algebra::square::Square::new(self)?);
        self.deligne_square = Some(square.clone());
        Ok(square)
    }

    pub(crate) fn deligne_tower(
        &mut self,
        rank: usize,
    ) -> Result<Rc<crate::partition_algebra::vacillating::Tower>> {
        if let Some(tower) = self.deligne_towers.get(&rank) {
            return Ok(tower.clone());
        }
        let tower = Rc::new(crate::partition_algebra::vacillating::Tower::new(
            self, rank,
        )?);
        self.deligne_towers.insert(rank, tower.clone());
        Ok(tower)
    }

    pub(crate) fn define_hilbert_stratum(
        &mut self,
        name: Label,
        partition: hyperelliptic_kummer_surface::stratum::Partition,
    ) -> Result<()> {
        if self.hilbert_strata.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Hilbert–Chow stratum '{}' is already defined", name.text),
            ));
        }
        self.allocate(8 + name.text.len() + 4 * partition.parts().len())?;
        self.hilbert_strata.insert(name.text, Rc::new(partition));
        Ok(())
    }

    pub(crate) fn hilbert_stratum(
        &self,
        name: &Label,
    ) -> Result<Rc<hyperelliptic_kummer_surface::stratum::Partition>> {
        self.hilbert_strata.get(&name.text).cloned().ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown Hilbert–Chow stratum '{}'", name.text),
            )
        })
    }

    pub(crate) fn define_specht_stratum(
        &mut self,
        name: Label,
        decoration: hyperelliptic_kummer_surface::specht::Decoration,
    ) -> Result<()> {
        if self.specht_strata.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Specht stratum '{}' is already defined", name.text),
            ));
        }
        self.allocate(name.text.len() + decoration.storage())?;
        self.specht_strata.insert(name.text, Rc::new(decoration));
        Ok(())
    }

    pub(crate) fn specht_stratum(
        &self,
        name: &Label,
    ) -> Result<Rc<hyperelliptic_kummer_surface::specht::Decoration>> {
        self.specht_strata
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown Specht stratum '{}'", name.text)))
    }

    pub(crate) fn show_jacobian_level_frobenius(
        &mut self,
        basis: &str,
        base: &str,
        power: usize,
        result: &crate::jacobian_level::Frobenius,
    ) -> Result<String> {
        let mut out = Text::new(MAX_SHOWN_BYTES);
        let multiplier = result.multiplier;
        out.push(&format!(
            "rvJacTorsionFrob[{basis}, {base}, {power}] = multiplier={multiplier}, matrix="
        ))?;
        self.write_jacobian_level_matrix(&mut out, &result.matrix)?;
        self.spend(result.characteristic.len())?;
        out.push(&format!(", characteristic={:?}", result.characteristic))?;
        Ok(out.string)
    }

    pub(crate) fn show_jacobian_power_frob(
        &mut self,
        basis: &str,
        base: &str,
        level: u64,
        power: usize,
        result: &crate::jacobian_level::Frobenius,
        inverse: &HWMatrix<u64>,
    ) -> Result<String> {
        let mut out = Text::new(MAX_SHOWN_BYTES);
        out.push(&format!(
            "rvJacPowerFrob[{basis}, {base}, {power}] = level={level}, multiplier={}, matrix=",
            result.multiplier
        ))?;
        self.write_jacobian_level_matrix(&mut out, &result.matrix)?;
        out.push(", inverse=")?;
        self.write_jacobian_level_matrix(&mut out, inverse)?;
        self.spend(result.characteristic.len())?;
        out.push(&format!(", characteristic={:?}", result.characteristic))?;
        Ok(out.string)
    }

    pub(crate) fn show_jacobian_power_transition(
        &mut self,
        source: &str,
        target: &str,
        root_exponent: u64,
        matrix: &HWMatrix<u64>,
        inverse: &HWMatrix<u64>,
    ) -> Result<String> {
        let mut out = Text::new(MAX_SHOWN_BYTES);
        out.push(&format!(
            "rvJacPowerTransition[{source}, {target}] = root_exponent={root_exponent}, matrix="
        ))?;
        self.write_jacobian_level_matrix(&mut out, matrix)?;
        out.push(", inverse=")?;
        self.write_jacobian_level_matrix(&mut out, inverse)?;
        Ok(out.string)
    }

    pub(crate) fn show_jacobian_richelot_matrix(
        &mut self,
        source: &str,
        target: &str,
        result: &crate::jacobian_level::richelot::Report,
    ) -> Result<String> {
        let mut out = Text::new(MAX_SHOWN_BYTES);
        out.push(&format!(
            "rvJacRichelotMatrix[{source}, {target}] = level={}, root exponent={}, multiplier={}, rank={}, column-action matrix=",
            result.prime.get(), result.root_exponent, result.multiplier, result.rank,
        ))?;
        self.write_jacobian_level_matrix(&mut out, &result.matrix)?;
        out.push(", column-action adjoint=")?;
        self.write_jacobian_level_matrix(&mut out, &result.adjoint)?;
        Ok(out.string)
    }

    fn write_jacobian_level_matrix(
        &mut self,
        out: &mut Text,
        matrix: &HWMatrix<u64>,
    ) -> Result<()> {
        out.push("[")?;
        for (index, row) in matrix.entries().chunks(matrix.size().max(1)).enumerate() {
            self.spend(1 + row.len())?;
            if index != 0 {
                out.push(", ")?;
            }
            out.push(&format!("{row:?}"))?;
        }
        out.push("]")
    }

    pub(crate) fn define_richelot(
        &mut self,
        target: Label,
        ring: &Label,
        coefficients: Vec<Vec<u64>>,
        isogeny: crate::richelot::Isogeny,
    ) -> Result<()> {
        self.define_hww_curve(target.clone(), ring, coefficients)?;
        self.richelot_isogenies
            .insert(target.text, Rc::new(isogeny));
        Ok(())
    }

    pub(crate) fn richelot(&self, target: &Label) -> Result<Rc<crate::richelot::Isogeny>> {
        self.richelot_isogenies
            .get(&target.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    target.span,
                    format!("curve '{}' has no stored Richelot splitting", target.text),
                )
            })
    }

    pub(crate) fn define_richelot_descent(
        &mut self,
        target: Label,
        ring: &Label,
        coefficients: Vec<Vec<u64>>,
        descent: crate::richelot::descent::Descent,
    ) -> Result<()> {
        self.define_hww_curve(target.clone(), ring, coefficients)?;
        self.richelot_descents.insert(target.text, Rc::new(descent));
        Ok(())
    }

    pub(crate) fn richelot_descent(
        &self,
        target: &Label,
    ) -> Result<Rc<crate::richelot::descent::Descent>> {
        self.richelot_descents
            .get(&target.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    target.span,
                    format!("curve '{}' has no stored Richelot descent", target.text),
                )
            })
    }

    pub(crate) fn define_bielliptic(
        &mut self,
        targets: [Label; 2],
        source: Rc<hw_closed::Curve>,
        matrix: [ECScalar; 4],
    ) -> Result<()> {
        let [plus, minus] = &targets;
        if plus.text == minus.text {
            return Err(Error::new(
                minus.span,
                "bielliptic target names must be distinct",
            ));
        }
        self.require_fresh_hyperelliptic_name(plus)?;
        self.require_fresh_hyperelliptic_name(minus)?;
        let quotients = crate::bielliptic::Quotient::build(self, source, matrix)?;
        for (name, quotient) in targets.into_iter().zip(quotients) {
            self.hw_closed_curves
                .insert(name.text.clone(), quotient.target.clone());
            self.bielliptic_quotients
                .insert(name.text, Rc::new(quotient));
        }
        Ok(())
    }

    pub(crate) fn bielliptic_quotient(
        &self,
        target: &Label,
    ) -> Result<Rc<crate::bielliptic::Quotient>> {
        self.bielliptic_quotients
            .get(&target.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    target.span,
                    format!("curve '{}' has no stored bielliptic quotient", target.text),
                )
            })
    }

    pub(crate) fn define_elliptic_gluing(
        &mut self,
        name: Label,
        factors: [crate::elliptic_gluing::Factor; 2],
        matrix: [ECScalar; 4],
    ) -> Result<()> {
        self.require_fresh_hyperelliptic_name(&name)?;
        let gluing =
            crate::elliptic_gluing::Gluing::build(self, name.text.as_str().into(), factors, matrix)
                .map_err(|error| error.at(name.span))?;
        self.hw_closed_curves
            .insert(name.text.clone(), gluing.curve.clone());
        self.elliptic_gluings.insert(name.text, Rc::new(gluing));
        Ok(())
    }

    pub(crate) fn elliptic_gluing(
        &self,
        name: &Label,
    ) -> Result<Rc<crate::elliptic_gluing::Gluing>> {
        self.elliptic_gluings
            .get(&name.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("curve '{}' has no stored elliptic gluing", name.text),
                )
            })
    }

    pub(crate) fn jacobian_curve(&self, name: &str) -> Result<Rc<hw_closed::Curve>> {
        let curve = self.residue_hyperelliptic_curve(name)?;
        hyperelliptic_jacobian::require_curve(&curve)?;
        Ok(curve)
    }

    fn residue_hyperelliptic_curve(&self, name: &str) -> Result<Rc<hw_closed::Curve>> {
        self.hw_closed_curves
            .get(name)
            .cloned()
            .ok_or_else(|| Error::plain(format!("unknown residue hyperelliptic curve '{name}'")))
    }

    pub(crate) fn kummer_curve(&self, name: &str) -> Result<Rc<hw_closed::Curve>> {
        let curve = self.residue_hyperelliptic_curve(name)?;
        hyperelliptic_kummer_surface::require_curve(&curve)?;
        Ok(curve)
    }

    pub(crate) fn show_jacobian_composition(
        &mut self,
        pair: &Divisor,
        quotient: &[ECScalar],
    ) -> Result<String> {
        let mut out = Text::new(MAX_SHOWN_BYTES);
        out.push("rvJacCompose = u=")?;
        self.write_ec_polynomial(&mut out, &pair.u)?;
        out.push(", v=")?;
        self.write_ec_polynomial(&mut out, &pair.v)?;
        out.push(", norm quotient=")?;
        self.write_ec_polynomial(&mut out, quotient)?;
        Ok(out.string)
    }

    pub(crate) fn jacobian_type(&mut self, name: Rc<str>) -> Result<TypeId> {
        self.jacobian_curve(&name)?;
        self.inhabited_type(Kind::Jacobian(name), 0)
    }

    pub(crate) fn kummer_type(&mut self, name: Rc<str>) -> Result<TypeId> {
        self.kummer_curve(&name)?;
        self.inhabited_type(Kind::Kummer(name), 0)
    }

    pub(crate) fn kummer_k3_type(&mut self, name: Rc<str>) -> Result<TypeId> {
        self.kummer_curve(&name)?;
        self.inhabited_type(Kind::KummerK3(name), 0)
    }

    pub(crate) fn define_k3_pencil(
        &mut self,
        name: Label,
        basis: Label,
        first: u8,
        second: u8,
    ) -> Result<()> {
        if self.k3_pencils.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("K3 pencil '{}' is already defined", name.text),
            ));
        }
        let pencil =
            pencil::Pencil::new(self, basis, first, second).map_err(|error| error.at(name.span))?;
        self.allocate(1 + name.text.len())?;
        self.k3_pencils.insert(name.text, Rc::new(pencil));
        Ok(())
    }

    pub(crate) fn k3_pencil(&self, name: &Label) -> Result<Rc<pencil::Pencil>> {
        self.k3_pencils
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown K3 pencil '{}'", name.text)))
    }

    pub(crate) fn k3_pencil_base_type(&mut self, name: &Label) -> Result<TypeId> {
        self.k3_pencil(name)?;
        self.inhabited_type(Kind::K3PencilBase(name.text.as_str().into()), 0)
    }

    pub(crate) fn k3_pencil_section_type(&mut self, name: &Label) -> Result<TypeId> {
        self.k3_pencil(name)?;
        self.inhabited_type(Kind::K3PencilSection(name.text.as_str().into()), 0)
    }

    pub(crate) fn k3_pencil_kummer2_type(&mut self, name: &Label) -> Result<TypeId> {
        self.k3_pencil(name)?;
        self.inhabited_type(Kind::K3PencilKummer2(name.text.as_str().into()), 0)
    }

    pub(crate) fn define_k3_cassels(
        &mut self,
        name: Label,
        cover: &Label,
        witnesses: Vec<pencil::sections::covering::cassels::AuxiliaryInput>,
    ) -> Result<()> {
        if self.k3_cassels.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Cassels witnesses '{}' are already defined", name.text),
            ));
        }
        let cover = self.k3_pencil_cover(cover)?;
        let pairing = pencil::sections::covering::cassels::Pairing::new(self, cover, witnesses)
            .map_err(|error| error.at(name.span))?;
        self.allocate(1 + name.text.len())?;
        self.k3_cassels.insert(name.text, Rc::new(pairing));
        Ok(())
    }

    pub(crate) fn k3_cassels(
        &self,
        name: &Label,
    ) -> Result<Rc<pencil::sections::covering::cassels::Pairing>> {
        self.k3_cassels.get(&name.text).cloned().ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown Cassels witnesses '{}'", name.text),
            )
        })
    }

    pub(crate) fn define_k3_selmer_four(
        &mut self,
        name: Label,
        image: &Label,
        witnesses: &[Label],
    ) -> Result<()> {
        if self.k3_selmer_four.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("four-Selmer certificate '{}' is already defined", name.text),
            ));
        }
        let image = self.k3_pencil_mw_image(image)?;
        let witnesses = witnesses
            .iter()
            .map(|witness| self.k3_cassels(witness))
            .collect::<Result<Vec<_>>>()?;
        let certificate =
            pencil::sections::covering::four::Certificate::new(self, image, &witnesses)
                .map_err(|error| error.at(name.span))?;
        self.allocate(1 + name.text.len())?;
        self.k3_selmer_four.insert(name.text, Rc::new(certificate));
        Ok(())
    }

    pub(crate) fn k3_selmer_four(
        &self,
        name: &Label,
    ) -> Result<Rc<pencil::sections::covering::four::Certificate>> {
        self.k3_selmer_four.get(&name.text).cloned().ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown four-Selmer certificate '{}'", name.text),
            )
        })
    }

    pub(crate) fn define_k3_pencil_cover(
        &mut self,
        name: Label,
        pencil_name: Label,
        descent: Rc<pencil::sections::Descent>,
    ) -> Result<()> {
        if self.k3_pencil_covers.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("K3 pencil covering '{}' is already defined", name.text),
            ));
        }
        let pencil = self.k3_pencil(&pencil_name)?;
        let cover = pencil::sections::covering::Cover::new(self, pencil_name, pencil, descent)
            .map_err(|error| error.at(name.span))?;
        self.allocate(1 + name.text.len())?;
        self.k3_pencil_covers.insert(name.text, Rc::new(cover));
        Ok(())
    }

    pub(crate) fn k3_pencil_cover(
        &self,
        name: &Label,
    ) -> Result<Rc<pencil::sections::covering::Cover>> {
        self.k3_pencil_covers
            .get(&name.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown K3 pencil covering '{}'", name.text),
                )
            })
    }

    pub(crate) fn k3_pencil_cover_point_type(&mut self, name: &Label) -> Result<TypeId> {
        self.k3_pencil_cover(name)?;
        self.inhabited_type(Kind::K3PencilCoverPoint(name.text.as_str().into()), 0)
    }

    pub(crate) fn define_k3_pencil_mw_image(
        &mut self,
        name: Label,
        pencil: &Label,
        sections: &[Rc<pencil::sections::Section>],
    ) -> Result<()> {
        if self.k3_pencil_mw_images.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!(
                    "K3 pencil Mordell-Weil certificate '{}' is already defined",
                    name.text
                ),
            ));
        }
        let named_pencil = self.k3_pencil(pencil)?;
        let image = pencil::sections::covering::MordellImage::new(
            self,
            pencil.clone(),
            named_pencil,
            sections,
        )
        .map_err(|error| error.at(name.span))?;
        self.allocate(1 + name.text.len())?;
        self.k3_pencil_mw_images.insert(name.text, Rc::new(image));
        Ok(())
    }

    pub(crate) fn k3_pencil_mw_image(
        &self,
        name: &Label,
    ) -> Result<Rc<pencil::sections::covering::MordellImage>> {
        self.k3_pencil_mw_images
            .get(&name.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown K3 pencil Mordell-Weil certificate '{}'", name.text),
                )
            })
    }

    pub(crate) fn define_k3_mordell(
        &mut self,
        name: Label,
        pencil: &Label,
        zero: u8,
    ) -> Result<()> {
        if self.k3_mordell.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("K3 Mordell–Weil marking '{}' is already defined", name.text),
            ));
        }
        let pencil = self.k3_pencil(pencil)?;
        let marked =
            mordell::Mordell::new(self, pencil, zero).map_err(|error| error.at(name.span))?;
        self.allocate(1 + name.text.len())?;
        self.k3_mordell.insert(name.text, Rc::new(marked));
        Ok(())
    }

    pub(crate) fn k3_mordell(&self, name: &Label) -> Result<Rc<mordell::Mordell>> {
        self.k3_mordell.get(&name.text).cloned().ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown K3 Mordell–Weil marking '{}'", name.text),
            )
        })
    }

    pub(crate) fn define_enriques(&mut self, declaration: enriques::Declare) -> Result<()> {
        let enriques::Declare(name, base, cover, basis, mask) = declaration;
        if self.enriques_quotients.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Enriques quotient '{}' is already defined", name.text),
            ));
        }
        let quotient = enriques::Quotient::new(self, base, cover, &basis, mask)?;
        self.allocate(1 + name.text.len())?;
        self.enriques_quotients.insert(name.text, Rc::new(quotient));
        Ok(())
    }

    pub(crate) fn enriques_quotient(&self, name: &str) -> Result<Rc<enriques::Quotient>> {
        self.enriques_quotients
            .get(name)
            .cloned()
            .ok_or_else(|| Error::plain(format!("unknown Enriques quotient '{name}'")))
    }

    pub(crate) fn enriques_type(&mut self, name: Rc<str>) -> Result<TypeId> {
        self.enriques_quotient(&name)?;
        self.inhabited_type(Kind::Enriques(name), 0)
    }

    pub(crate) fn kummer_disc_type(&mut self, name: &Label) -> Result<TypeId> {
        nikulin::basis(self, name)?;
        self.inhabited_type(Kind::KummerDisc(name.text.as_str().into()), 0)
    }

    pub(crate) fn kummer_weil_type(&mut self, name: &Label) -> Result<TypeId> {
        nikulin::basis_for(
            self,
            name,
            "Kummer Weil vectors require a certified genus-two level-2 basis",
        )?;
        self.inhabited_type(Kind::KummerWeil(name.text.as_str().into()), 0)
    }

    pub(crate) fn k3_lattice_type(&mut self, name: &Label) -> Result<TypeId> {
        nikulin::basis_for(
            self,
            name,
            "K3 lattices require a certified genus-two level-2 basis",
        )?;
        self.inhabited_type(Kind::K3Lattice(name.text.as_str().into()), 0)
    }

    pub(crate) fn og6_clifford_type(&mut self, name: &Label, parity: u8) -> Result<TypeId> {
        self.og6_cusp(name)?;
        self.inhabited_type(Kind::OG6Clifford(name.text.as_str().into(), parity), 0)
    }

    pub(crate) fn ks_exterior_type(&mut self, name: &Label, degree: i128) -> Result<TypeId> {
        self.og6_cusp(name)?;
        self.inhabited_type(Kind::KSExterior(name.text.as_str().into(), degree), 0)
    }

    pub(crate) fn ks_nomizu_type(&mut self, name: &Label, m: i128, degree: i128) -> Result<TypeId> {
        self.og6_nilpotent(name)?;
        self.inhabited_type(Kind::KSNomizu(name.text.as_str().into(), m, degree), 0)
    }

    pub(crate) fn ks_end_type(&mut self, name: &Label) -> Result<TypeId> {
        self.og6_cusp(name)?;
        self.inhabited_type(Kind::KSEnd(name.text.as_str().into()), 0)
    }

    pub(crate) fn qz_type(&mut self) -> Result<TypeId> {
        self.inhabited_type(Kind::QZ, 0)
    }

    pub(crate) fn k3_clifford_type(&mut self, basis: &Label, parity: u8) -> Result<TypeId> {
        nikulin::basis_for(
            self,
            basis,
            "K3 Clifford algebra requires a certified genus-two level-2 basis",
        )?;
        self.inhabited_type(Kind::K3Clifford(basis.text.as_str().into(), parity), 0)
    }

    pub(crate) fn k3_spinor_type(&mut self, basis: &Label, chirality: u8) -> Result<TypeId> {
        nikulin::basis_for(
            self,
            basis,
            "K3 spinors require a certified genus-two level-2 basis",
        )?;
        self.inhabited_type(Kind::K3Spinor(basis.text.as_str().into(), chirality), 0)
    }

    pub(crate) fn define_mukai_moduli(
        &mut self,
        name: Label,
        marking: Rc<str>,
        vector: Rc<hyperelliptic_kummer_surface::mukai::Vector>,
    ) -> Result<()> {
        if self.mukai_moduli.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Mukai moduli lattice '{}' is already defined", name.text),
            ));
        }
        let lattice =
            hyperelliptic_kummer_surface::mukai::moduli::Lattice::new(self, marking, vector)?;
        self.allocate(name.text.len() + lattice.cells())?;
        self.mukai_moduli.insert(name.text, Rc::new(lattice));
        Ok(())
    }

    pub(crate) fn mukai_moduli(
        &self,
        name: &Label,
    ) -> Result<Rc<hyperelliptic_kummer_surface::mukai::moduli::Lattice>> {
        self.mukai_moduli.get(&name.text).cloned().ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown Mukai moduli lattice '{}'", name.text),
            )
        })
    }

    pub(crate) fn abelian_mukai_type(&mut self, name: &Label) -> Result<TypeId> {
        nikulin::basis_for(
            self,
            name,
            "abelian Mukai lattices require a certified genus-two level-2 basis",
        )?;
        self.inhabited_type(Kind::AbMukai(name.text.as_str().into()), 0)
    }

    pub(crate) fn define_og6(
        &mut self,
        name: Label,
        marking: Rc<str>,
        vector: Rc<hyperelliptic_kummer_surface::mukai::abelian::Vector>,
    ) -> Result<()> {
        if self.og6_lattices.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("OG6 lattice '{}' is already defined", name.text),
            ));
        }
        let lattice =
            hyperelliptic_kummer_surface::mukai::abelian::og6::Lattice::new(self, marking, vector)?;
        self.allocate(name.text.len() + lattice.cells())?;
        self.og6_lattices.insert(name.text, Rc::new(lattice));
        Ok(())
    }

    pub(crate) fn define_og6_nilpotent(
        &mut self,
        name: Label,
        limit: hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::nilpotent::Limit,
    ) -> Result<()> {
        if self.og6_nilpotents.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("OG6 nilpotent '{}' is already defined", name.text),
            ));
        }
        self.allocate(name.text.len() + limit.cells())?;
        self.og6_nilpotents.insert(name.text, Rc::new(limit));
        Ok(())
    }

    pub(crate) fn og6_nilpotent(
        &self,
        name: &Label,
    ) -> Result<Rc<hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::nilpotent::Limit>> {
        self.og6_nilpotents
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown OG6 nilpotent '{}'", name.text)))
    }

    pub(crate) fn define_og6_cusp(
        &mut self,
        name: Label,
        lattice_name: &Label,
        inputs: [Rc<hyperelliptic_kummer_surface::mukai::abelian::og6::Class>; 3],
    ) -> Result<()> {
        if self.og6_cusps.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("OG6 cusp '{}' is already defined", name.text),
            ));
        }
        let lattice = self.og6_lattice(lattice_name)?;
        let cusp = hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::Cusp::new(
            self,
            name.text.as_str().into(),
            lattice_name.text.as_str().into(),
            &lattice,
            inputs,
        )?;
        self.allocate(name.text.len() + cusp.cells())?;
        self.og6_cusps.insert(name.text, Rc::new(cusp));
        Ok(())
    }

    pub(crate) fn og6_cusp(
        &self,
        name: &Label,
    ) -> Result<Rc<hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::Cusp>> {
        self.og6_cusps
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown OG6 cusp '{}'", name.text)))
    }

    pub(crate) fn og6_lattice(
        &self,
        name: &Label,
    ) -> Result<Rc<hyperelliptic_kummer_surface::mukai::abelian::og6::Lattice>> {
        self.og6_lattices
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown OG6 lattice '{}'", name.text)))
    }

    pub(crate) fn og6_type(&mut self, name: &Label) -> Result<TypeId> {
        self.og6_lattice(name)?;
        self.inhabited_type(Kind::OG6(name.text.as_str().into()), 0)
    }

    pub(crate) fn og10_type(&mut self, name: &Label) -> Result<TypeId> {
        self.mukai_moduli(name)?.require_og10()?;
        self.inhabited_type(Kind::OG10(name.text.as_str().into()), 0)
    }

    pub(crate) fn mukai_quotient_type(&mut self, name: &Label) -> Result<TypeId> {
        self.mukai_moduli(name)?
            .require_isotropic("RVMukaiQuotient")?;
        self.inhabited_type(Kind::MukaiQuotient(name.text.as_str().into()), 0)
    }

    pub(crate) fn mukai_type(&mut self, name: &Label) -> Result<TypeId> {
        nikulin::basis_for(
            self,
            name,
            "Mukai lattices require a certified genus-two level-2 basis",
        )?;
        self.inhabited_type(Kind::Mukai(name.text.as_str().into()), 0)
    }

    pub(crate) fn nakajima_type(
        &mut self,
        name: &Label,
        length: i128,
        degree: i128,
    ) -> Result<TypeId> {
        nikulin::basis_for(
            self,
            name,
            "Nakajima states require a certified genus-two level-2 basis",
        )?;
        let grade = nakajima::Grade {
            basis: name.text.as_str().into(),
            length,
            degree,
        };
        self.inhabited_type(Kind::Nakajima(grade), 0)
    }

    pub(crate) fn hilbert_cohomology_type(
        &mut self,
        name: &Label,
        length: i128,
        degree: i128,
    ) -> Result<TypeId> {
        nikulin::basis_for(
            self,
            name,
            "Hilbert cohomology requires a certified genus-two level-2 basis",
        )?;
        let grade = nakajima::Grade {
            basis: name.text.as_str().into(),
            length,
            degree,
        };
        self.inhabited_type(Kind::HilbertCohomology(grade), 0)
    }

    pub(crate) fn register_hodge_wick(
        &mut self,
        name: Label,
        multiplicity: usize,
        polynomial: Vec<i128>,
        gram: Vec<Vec<hyperelliptic_kummer_surface::hodge_wick::Entry>>,
    ) -> Result<()> {
        if self.hodge_wick_forms.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Hodge Wick marking '{}' is already defined", name.text),
            ));
        }
        let form = hyperelliptic_kummer_surface::hodge_wick::Form::new(
            self,
            multiplicity,
            polynomial,
            gram,
        )
        .map_err(|error| error.at(name.span))?;
        self.hodge_wick_forms.insert(name.text, Rc::new(form));
        Ok(())
    }

    pub(crate) fn hodge_wick(
        &self,
        name: &str,
    ) -> Result<Rc<hyperelliptic_kummer_surface::hodge_wick::Form>> {
        self.hodge_wick_forms
            .get(name)
            .cloned()
            .ok_or_else(|| Error::plain(format!("unknown Hodge Wick marking '{name}'")))
    }

    pub(crate) fn verbitsky_type(&mut self, grade: nakajima::Grade) -> Result<TypeId> {
        self.hodge_wick(&grade.basis)?.validate_grade(&grade)?;
        self.inhabited_type(Kind::Verbitsky(grade), 0)
    }

    pub(crate) fn hodge_tensor_type(&mut self, grade: nakajima::Grade) -> Result<TypeId> {
        self.hodge_wick(&grade.basis)?.validate_grade(&grade)?;
        self.inhabited_type(Kind::HodgeTensor(grade), 0)
    }

    pub(crate) fn nakajima_gram(&mut self) -> Result<Rc<nakajima::Gram>> {
        if let Some(gram) = &self.nakajima_gram {
            return Ok(gram.clone());
        }
        let gram = Rc::new(nakajima::Gram::new(self)?);
        self.nakajima_gram = Some(gram.clone());
        Ok(gram)
    }

    pub(crate) fn nakajima_inverse_gram(&mut self) -> Result<Rc<nakajima::Gram>> {
        if let Some(gram) = &self.nakajima_inverse_gram {
            return Ok(gram.clone());
        }
        let gram = self.nakajima_gram()?;
        let inverse = Rc::new(gram.inverse(self)?);
        self.nakajima_inverse_gram = Some(inverse.clone());
        Ok(inverse)
    }

    pub(crate) fn hilbert_square_type(&mut self, name: &Label, degree: i128) -> Result<TypeId> {
        nikulin::basis_for(
            self,
            name,
            "Hilbert-square cohomology requires a certified genus-two level-2 basis",
        )?;
        self.inhabited_type(Kind::HilbertSquare(name.text.as_str().into(), degree), 0)
    }

    pub(crate) fn hilbert_square_form(&mut self) -> Result<Rc<hilbert_square::Form>> {
        if let Some(form) = &self.hilbert_square_form {
            return Ok(form.clone());
        }
        let form = Rc::new(hilbert_square::Form::new(self)?);
        self.hilbert_square_form = Some(form.clone());
        Ok(form)
    }

    pub(crate) fn siegel_weight_type(&mut self) -> Result<TypeId> {
        self.inhabited_type(Kind::SiegelWeight, 0)
    }

    pub(crate) fn modular_hecke_type(&mut self, weight: u128) -> Result<TypeId> {
        hecke::validate_weight(weight)?;
        self.inhabited_type(Kind::ModularHecke(weight), 0)
    }

    pub(crate) fn register_hecke_eigenpacket(
        &mut self,
        name: Label,
        weight: u128,
        generator: Rc<hecke::Operator>,
    ) -> Result<()> {
        if self.hecke_eigenpackets.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Hecke eigenpacket '{}' is already defined", name.text),
            ));
        }
        let packet = eigenpacket::Packet::new(self, weight, generator)
            .map_err(|error| error.at(name.span))?;
        self.allocate(name.text.len() + 1)?;
        self.hecke_eigenpackets.insert(name.text, Rc::new(packet));
        Ok(())
    }

    pub(crate) fn hecke_eigenpacket(&self, name: &Label) -> Result<Rc<eigenpacket::Packet>> {
        self.hecke_eigenpackets
            .get(&name.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown Hecke eigenpacket '{}'", name.text),
                )
            })
    }

    pub(crate) fn hecke_eigen_type(&mut self, name: &Label) -> Result<TypeId> {
        self.hecke_eigenpacket(name)?;
        self.inhabited_type(Kind::HeckeEigen(name.text.as_str().into()), 0)
    }

    pub(crate) fn register_hecke_residual(
        &mut self,
        name: Label,
        source: Label,
        prime: usize,
        root: i128,
    ) -> Result<()> {
        if self.hecke_residuals.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!(
                    "residual Hecke specialization '{}' is already defined",
                    name.text
                ),
            ));
        }
        let packet = self.hecke_eigenpacket(&source)?;
        let cells = name.text.len() + source.text.len() + 8;
        let specialization = residual_hecke::Specialization::new(self, source, packet, prime, root)
            .map_err(|error| error.at(name.span))?;
        self.allocate(cells)?;
        self.hecke_residuals
            .insert(name.text, Rc::new(specialization));
        Ok(())
    }

    pub(crate) fn hecke_residual(
        &self,
        name: &Label,
    ) -> Result<Rc<residual_hecke::Specialization>> {
        self.hecke_residuals
            .get(&name.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown residual Hecke specialization '{}'", name.text),
                )
            })
    }

    pub(crate) fn register_hecke_eigenperiod(
        &mut self,
        declaration: period::eigenperiod::Declare,
    ) -> Result<()> {
        let period::eigenperiod::Declare(name, source, sign, pivot) = declaration;
        if self.hecke_eigenperiods.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Hecke eigenperiod '{}' is already defined", name.text),
            ));
        }
        let packet = self.hecke_eigenpacket(&source)?;
        let cells = name.text.len() + source.text.len() + 8;
        let eigenperiod = period::eigenperiod::Eigenperiod::new(self, source, packet, sign, pivot)
            .map_err(|error| error.at(name.span))?;
        self.allocate(cells)?;
        self.hecke_eigenperiods
            .insert(name.text, Rc::new(eigenperiod));
        Ok(())
    }

    pub(crate) fn hecke_eigenperiod(
        &self,
        name: &Label,
    ) -> Result<Rc<period::eigenperiod::Eigenperiod>> {
        self.hecke_eigenperiods
            .get(&name.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown Hecke eigenperiod '{}'", name.text),
                )
            })
    }

    pub(crate) fn mazur_tate_type(&mut self, grade: period::mazur_tate::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::MazurTate(grade), 0)
    }

    pub(crate) fn mazur_amice_type(&mut self, grade: period::amice::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::MazurAmice(grade), 0)
    }

    pub(crate) fn mazur_ordinary_type(&mut self, name: &Label) -> Result<TypeId> {
        let modulus = self.mazur_ordinary(name)?.modulus;
        self.inhabited_type(Kind::MazurOrdinary(name.text.as_str().into(), modulus), 0)
    }

    pub(crate) fn register_mazur_ordinary(
        &mut self,
        name: Label,
        period: Label,
        embedding: period::mazur_tate::Embedding,
        precision: usize,
    ) -> Result<()> {
        if self.mazur_ordinary.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!(
                    "ordinary Mazur–Tate measure '{}' is already defined",
                    name.text
                ),
            ));
        }
        let measure = period::mazur_tate::Ordinary::new(self, &period, embedding, precision)
            .map_err(|error| error.at(name.span))?;
        self.allocate(name.text.len() + 20)?;
        self.mazur_ordinary.insert(name.text, Rc::new(measure));
        Ok(())
    }

    pub(crate) fn mazur_ordinary(&self, name: &Label) -> Result<Rc<period::mazur_tate::Ordinary>> {
        self.mazur_ordinary.get(&name.text).cloned().ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown ordinary Mazur–Tate measure '{}'", name.text),
            )
        })
    }

    pub(crate) fn hecke_symbol_type(&mut self, name: &Label) -> Result<TypeId> {
        self.hecke_eigenperiod(name)?;
        self.inhabited_type(Kind::HeckeSymbol(name.text.as_str().into()), 0)
    }

    pub(crate) fn period_type(&mut self, grade: period::Grade) -> Result<TypeId> {
        self.inhabited_type(Kind::PeriodPolynomial(grade), 0)
    }

    pub(crate) fn period_basis(&mut self, grade: period::Grade) -> Result<Rc<period::Basis>> {
        if let Some(basis) = self.period_bases.get(&grade) {
            return Ok(basis.clone());
        }
        let basis = Rc::new(period::Basis::new(self, grade)?);
        self.period_bases.insert(grade, basis.clone());
        Ok(basis)
    }

    pub(crate) fn integral_period_basis(
        &mut self,
        grade: period::Grade,
    ) -> Result<Rc<period::haberland::IntegralBasis>> {
        if let Some(basis) = self.integral_period_bases.get(&grade) {
            return Ok(basis.clone());
        }
        let basis = Rc::new(period::haberland::IntegralBasis::new(self, grade)?);
        self.integral_period_bases.insert(grade, basis.clone());
        Ok(basis)
    }

    pub(crate) fn register_siegel_chamber(
        &mut self,
        name: Label,
        coordinates: [i128; 3],
    ) -> Result<()> {
        if self.siegel_chambers.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Siegel chamber '{}' is already defined", name.text),
            ));
        }
        let period = farey::Period::new(self, coordinates).map_err(|error| error.at(name.span))?;
        self.allocate(name.text.len() + 10)?;
        self.siegel_chambers.insert(name.text, Rc::new(period));
        Ok(())
    }

    pub(crate) fn siegel_chamber(&self, name: &Label) -> Result<Rc<farey::Period>> {
        self.siegel_chambers
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown Siegel chamber '{}'", name.text)))
    }

    pub(crate) fn register_siegel_maass(
        &mut self,
        name: Label,
        weight: u128,
        terms: Vec<(maass::Monomial, i128)>,
    ) -> Result<()> {
        if self.siegel_maass_forms.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Maaß lift '{}' is already defined", name.text),
            ));
        }
        let form = maass::Form::new(self, weight, terms).map_err(|error| error.at(name.span))?;
        self.allocate(name.text.len() + 1)?;
        self.siegel_maass_forms.insert(name.text, Rc::new(form));
        Ok(())
    }

    pub(crate) fn siegel_maass(&self, name: &Label) -> Result<Rc<maass::Form>> {
        self.siegel_maass_forms
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown Maaß lift '{}'", name.text)))
    }

    pub(crate) fn show_jacobian_kummer(
        &mut self,
        support: &hyperelliptic_kummer::Support,
    ) -> Result<String> {
        let mut out = Text::new(MAX_SHOWN_BYTES);
        out.push("rvJacKummer = B=")?;
        self.write_ec_polynomial(&mut out, &support.polynomial)?;
        out.push(", nonsquare orbits by degree=[")?;
        for (i, (degree, count)) in support.counts.iter().enumerate() {
            self.spend(1)?;
            if i > 0 {
                out.push(", ")?;
            }
            out.push(&format!("({degree}, {count})"))?;
        }
        out.push("]")?;
        Ok(out.string)
    }

    pub(crate) fn hww_matrix_type(&mut self, grade: hw_closed::Grade) -> Result<TypeId> {
        let ring = self
            .rv_unram_named(&grade.ring)
            .expect("registered unramified ring");
        hw_closed::matrix_size(&ring, grade.genus)?;
        self.inhabited_type(Kind::HWWMatrix(grade), 0)
    }

    pub(crate) fn register_rv_unram(
        &mut self,
        name: Label,
        p: usize,
        precision: usize,
        polynomial: Vec<i128>,
    ) -> Result<()> {
        if self.rv_unram_rings.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("unramified ring '{}' is already defined", name.text),
            ));
        }
        let ring =
            UnramRing::new(self, p, precision, polynomial).map_err(|error| error.at(name.span))?;
        self.rv_unram_rings.insert(name.text, Rc::new(ring));
        Ok(())
    }

    pub(crate) fn rv_unram_named(&self, name: &str) -> Option<Rc<UnramRing>> {
        self.rv_unram_rings.get(name).cloned()
    }

    pub(crate) fn rv_unram(&self, name: &Label) -> Result<Rc<UnramRing>> {
        self.rv_unram_named(&name.text).ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown unramified ring '{}'", name.text),
            )
        })
    }

    pub(crate) fn rv_unram_type(&mut self, name: &Label) -> Result<TypeId> {
        self.rv_unram(name)?;
        self.inhabited_type(Kind::RVUnram(name.text.as_str().into()), 0)
    }

    pub(crate) fn rv_witt_tower(&self, name: &Label) -> Result<Rc<RvWittTower>> {
        self.rv_witt_towers
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown Witt tower '{}'", name.text)))
    }

    pub(crate) fn check_new_rv_witt_tower(&self, name: &Label) -> Result<()> {
        if self.rv_witt_towers.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Witt tower '{}' is already defined", name.text),
            ));
        }
        Ok(())
    }

    pub(crate) fn define_rv_witt_tower(&mut self, name: String, tower: RvWittTower) {
        self.rv_witt_towers.insert(name, Rc::new(tower));
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

    pub(crate) fn rv_witt_cover_named(&self, name: &str) -> Option<Rc<RvWittCover>> {
        self.rv_witt_covers.get(name).cloned()
    }

    fn curve_genus(&self, family: DiffFamily, name: &str) -> usize {
        match family {
            DiffFamily::ArtinSchreier => self.rv_as_covers[name].genus(),
            DiffFamily::Witt => self.rv_witt_covers[name].genus(),
        }
    }

    pub(crate) fn eo_type(&mut self, genus: usize) -> Result<TypeId> {
        self.inhabited_type(Kind::EOType(genus), 0)
    }

    pub(crate) fn binary_scalar_type(&mut self, field: BinaryField) -> Result<TypeId> {
        self.inhabited_type(Kind::BinaryScalar(field), 0)
    }

    pub(crate) fn rv_as_derham_type(&mut self, name: &Label, modulus: usize) -> Result<TypeId> {
        self.rv_as_cover(name)?;
        let field = BinaryField::new(self, modulus)?;
        self.rv_as_derham_field(name.text.as_str().into(), field)
            .map_err(|error| error.at(name.span))
    }

    pub(crate) fn rv_as_derham_field(
        &mut self,
        name: Rc<str>,
        field: BinaryField,
    ) -> Result<TypeId> {
        self.rv_as_covers[name.as_ref()].polynomial_chart()?;
        self.inhabited_type(Kind::ASDeRham(name, field), 0)
    }

    pub(crate) fn rv_curve_diff_type(
        &mut self,
        family: DiffFamily,
        name: &Label,
        modulus: usize,
    ) -> Result<TypeId> {
        match family {
            DiffFamily::ArtinSchreier => {
                self.rv_as_cover(name)?;
            }
            DiffFamily::Witt => {
                self.rv_witt_cover(name)?;
            }
        }
        let field = BinaryField::new(self, modulus)?;
        self.inhabited_type(Kind::CurveDiff(family, name.text.as_str().into(), field), 0)
    }

    pub(crate) fn prime_scalar_type(&mut self, prime: Prime) -> Result<TypeId> {
        self.inhabited_type(Kind::PrimeScalar(prime), 0)
    }

    pub(crate) fn rv_as_odd_derham_type(&mut self, name: &Label) -> Result<TypeId> {
        let cover = self.rv_as_odd_cover(name)?;
        self.rv_as_odd_derham_named(name.text.as_str().into(), &cover)
    }

    pub(crate) fn rv_as_odd_derham_named(
        &mut self,
        name: Rc<str>,
        cover: &RvAsOddCover,
    ) -> Result<TypeId> {
        if cover.genus().checked_mul(2).is_none() {
            return Err(Error::plain(
                "odd Artin–Schreier de Rham dimension is out of range",
            ));
        }
        self.inhabited_type(Kind::ASOddDeRham(name, cover.prime()), 0)
    }

    pub(crate) fn rv_as_odd_tate_type(&mut self, name: &Label, degree: i64) -> Result<TypeId> {
        let cover = self.rv_as_odd_cover(name)?;
        self.rv_as_odd_tate_named(name.text.as_str().into(), &cover, degree)
    }

    pub(crate) fn rv_as_odd_tate_named(
        &mut self,
        name: Rc<str>,
        cover: &RvAsOddCover,
        degree: i64,
    ) -> Result<TypeId> {
        if cover.genus().checked_mul(2).is_none() {
            return Err(Error::plain(
                "odd Artin–Schreier de Rham dimension is out of range",
            ));
        }
        self.inhabited_type(Kind::ASOddTate(name, cover.prime(), degree), 0)
    }

    pub(crate) fn rv_as_odd_diff_type(&mut self, name: &Label) -> Result<TypeId> {
        let cover = self.rv_as_odd_cover(name)?;
        self.inhabited_type(Kind::ASOddDiff(name.text.as_str().into(), cover.prime()), 0)
    }

    pub(crate) fn rv_as_two_cover_named(&self, name: &str) -> Option<Rc<RvAsTwoCover>> {
        self.rv_as_two_covers.get(name).cloned()
    }

    pub(crate) fn rv_as_two_diff_type(&mut self, name: &Label) -> Result<TypeId> {
        let cover = self.rv_as_two_cover(name)?;
        self.inhabited_type(Kind::ASTwoDiff(name.text.as_str().into(), cover.prime()), 0)
    }

    pub(crate) fn rv_as_two_cover(&self, name: &Label) -> Result<Rc<RvAsTwoCover>> {
        self.rv_as_two_cover_named(&name.text).ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown two-pole Artin–Schreier cover '{}'", name.text),
            )
        })
    }

    pub(crate) fn register_rv_as_two_cover(
        &mut self,
        name: Label,
        prime: usize,
        terms: Vec<(i128, i128)>,
    ) -> Result<()> {
        if self.rv_as_two_covers.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!(
                    "two-pole Artin–Schreier cover '{}' is already defined",
                    name.text
                ),
            ));
        }
        let cover = RvAsTwoCover::build(self, prime, terms).map_err(|error| error.at(name.span))?;
        self.rv_as_two_covers.insert(name.text, Rc::new(cover));
        Ok(())
    }

    pub(crate) fn rv_as_odd_cover_named(&self, name: &str) -> Option<Rc<RvAsOddCover>> {
        self.rv_as_odd_covers.get(name).cloned()
    }

    pub(crate) fn rv_as_odd_cover(&self, name: &Label) -> Result<Rc<RvAsOddCover>> {
        self.rv_as_odd_covers
            .get(&name.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown odd Artin–Schreier cover '{}'", name.text),
                )
            })
    }

    pub(crate) fn check_new_rv_as_odd_cover(&self, name: &Label) -> Result<()> {
        if self.rv_as_odd_covers.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!(
                    "odd Artin–Schreier cover '{}' is already defined",
                    name.text
                ),
            ));
        }
        Ok(())
    }

    pub(crate) fn define_rv_as_odd_cover(&mut self, name: String, cover: RvAsOddCover) {
        self.rv_as_odd_covers.insert(name, Rc::new(cover));
    }

    pub(crate) fn rv_crystal(&self, name: &Label) -> Result<Rc<Crystal>> {
        self.rv_crystals
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown F-crystal '{}'", name.text)))
    }

    pub(crate) fn register_rv_crystal(
        &mut self,
        name: Label,
        prime: usize,
        genus: usize,
        entries: CrystalEntries,
    ) -> Result<()> {
        if self.rv_crystals.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("F-crystal '{}' is already defined", name.text),
            ));
        }
        let crystal =
            Crystal::build(self, prime, genus, entries).map_err(|error| error.at(name.span))?;
        self.rv_crystals.insert(name.text, Rc::new(crystal));
        Ok(())
    }

    pub(crate) fn define_rv_period(
        &mut self,
        name: Label,
        variables: usize,
        terms: Vec<rv_period::Literal>,
    ) -> Result<()> {
        if self.rv_periods.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("Laurent period '{}' is already defined", name.text),
            ));
        }
        let period = Period::build(self, variables, terms)?;
        self.rv_periods.insert(name.text, Rc::new(period));
        Ok(())
    }

    pub(crate) fn rv_period(&self, name: &Label) -> Result<Rc<Period>> {
        self.rv_periods
            .get(&name.text)
            .cloned()
            .ok_or_else(|| Error::new(name.span, format!("unknown Laurent period '{}'", name.text)))
    }

    pub(crate) fn rv_period_for(&self, grade: &PeriodGrade) -> Rc<Period> {
        self.rv_periods[grade.period.as_ref()].clone()
    }

    pub(crate) fn rv_ct_state_grade(
        &mut self,
        name: Label,
        p: usize,
        s: usize,
    ) -> Result<PeriodGrade> {
        self.rv_period(&name)?;
        let level = PeriodLevel::new(self, p, s)?;
        Ok(PeriodGrade {
            period: name.text.into(),
            level,
        })
    }

    pub(crate) fn rv_ct_state_type(&mut self, grade: PeriodGrade) -> Result<TypeId> {
        self.inhabited_type(Kind::RVCTState(grade), 0)
    }

    pub(crate) fn rv_ct_residue_type(&mut self, level: PeriodLevel) -> Result<TypeId> {
        self.inhabited_type(Kind::RVCTResidue(level), 0)
    }

    pub(crate) fn register_qp_quad(
        &mut self,
        name: Label,
        p: usize,
        d: i128,
        r: i128,
    ) -> Result<()> {
        if self.rv_qp_quadratics.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("quadratic field '{}' is already defined", name.text),
            ));
        }
        let field = QuadField::new(self, p, d, r).map_err(|error| error.at(name.span))?;
        self.allocate(5)?;
        self.rv_qp_quadratics.insert(name.text, field);
        Ok(())
    }

    pub(crate) fn qp_quad_named(&self, name: &str) -> Option<QuadField> {
        self.rv_qp_quadratics.get(name).copied()
    }

    pub(crate) fn qp_quad(&self, name: &Label) -> Result<QuadField> {
        self.qp_quad_named(&name.text).ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown quadratic field '{}'", name.text),
            )
        })
    }

    pub(crate) fn rv_res_space_type(&mut self, grade: ResGrade) -> Result<TypeId> {
        self.inhabited_type(Kind::RVResSpace(grade), 0)
    }

    pub(crate) fn rv_hyper_arith_type(
        &mut self,
        grade: crate::hyperelliptic_graph::Grade,
    ) -> Result<TypeId> {
        self.inhabited_type(Kind::RVHyperArithmetic(grade), 0)
    }

    pub(crate) fn rv_hyper_cover_type(
        &mut self,
        grade: crate::hyperelliptic_graph::Grade,
    ) -> Result<TypeId> {
        self.inhabited_type(Kind::RVHyperCover(grade), 0)
    }

    pub(crate) fn rv_stable_tree_type(&mut self, marks: usize) -> Result<TypeId> {
        stable_marks::check_count(marks)?;
        self.inhabited_type(Kind::RVStableTree(marks), 0)
    }

    pub(crate) fn rv_stable_curve_type(&mut self, prime: u32, marks: usize) -> Result<TypeId> {
        stable_marks::check_count(marks)?;
        self.inhabited_type(Kind::RVStableMarks(prime, marks), 0)
    }

    pub(crate) fn rv_valuated_matroid_type(&mut self, grade: ValGrade) -> Result<TypeId> {
        self.inhabited_type(Kind::RVValMat(grade), 0)
    }

    pub(crate) fn qp_quad_type(&mut self, name: &Label) -> Result<TypeId> {
        self.qp_quad(name)?;
        self.inhabited_type(Kind::RVQpQuad(name.text.as_str().into()), 0)
    }

    pub(crate) fn register_legendre_cm(
        &mut self,
        name: Label,
        family: legendre_cm::Family,
    ) -> Result<()> {
        self.allocate(8)?;
        self.register_qp_quad(
            name.clone(),
            family.prime,
            family.discriminant,
            family.trace,
        )?;
        self.legendre_cm_families.insert(name.text, family);
        Ok(())
    }

    pub(crate) fn legendre_cm_named(&self, name: &str) -> Option<legendre_cm::Family> {
        self.legendre_cm_families.get(name).copied()
    }

    pub(crate) fn legendre_cm(&self, name: &Label) -> Result<legendre_cm::Family> {
        self.legendre_cm_named(&name.text).ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown Legendre CM field '{}'", name.text),
            )
        })
    }

    pub(crate) fn cm_order_type(&mut self, name: &Label, exponent: usize) -> Result<TypeId> {
        self.legendre_cm(name)?;
        legendre_cm::require_positive(exponent)?;
        let grade = legendre_cm::Grade {
            family: name.text.as_str().into(),
            exponent,
        };
        self.inhabited_type(Kind::CMOrder(grade), 0)
    }

    pub(crate) fn elliptic_modular_polynomial(
        &mut self,
        level: Prime,
    ) -> Result<Rc<crate::elliptic_modular::Polynomial>> {
        if let Some(polynomial) = self.elliptic_modular_polynomials.get(&level.get()) {
            return Ok(Rc::clone(polynomial));
        }
        let polynomial = Rc::new(crate::elliptic_modular::Polynomial::new(self, level)?);
        self.elliptic_modular_polynomials
            .insert(level.get(), Rc::clone(&polynomial));
        Ok(polynomial)
    }

    pub(crate) fn cm_class_polynomial(
        &self,
        name: &Label,
    ) -> Result<Rc<crate::cm_class_poly::Polynomial>> {
        self.cm_class_polynomials
            .get(&name.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown CM class polynomial '{}'", name.text),
                )
            })
    }

    pub(crate) fn cm_class_polynomial_at(
        &self,
        discriminant: i128,
    ) -> Option<Rc<crate::cm_class_poly::Polynomial>> {
        self.cm_class_discriminants.get(&discriminant).cloned()
    }

    pub(crate) fn define_cm_class_polynomial(
        &mut self,
        name: Label,
        polynomial: crate::cm_class_poly::Polynomial,
    ) -> Result<()> {
        if self.cm_class_polynomials.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("CM class polynomial '{}' is already defined", name.text),
            ));
        }
        let polynomial = Rc::new(polynomial);
        if let Some(existing) = self.cm_class_discriminants.get(&polynomial.discriminant) {
            if !polynomial.same(existing) {
                return Err(Error::new(
                    name.span,
                    "CM class polynomial lifts to the same discriminant disagree",
                ));
            }
            self.cm_class_polynomials
                .insert(name.text, Rc::clone(existing));
        } else {
            self.cm_class_discriminants
                .insert(polynomial.discriminant, Rc::clone(&polynomial));
            self.cm_class_polynomials.insert(name.text, polynomial);
        }
        Ok(())
    }

    pub(crate) fn define_cm_ring_order(&mut self, name: Label, discriminant: i128) -> Result<()> {
        if self.cm_ring_orders.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("CM ring order '{}' is already defined", name.text),
            ));
        }
        let order = cm_ring::Order::new(discriminant)?;
        self.allocate(5 + name.text.len())?;
        self.cm_ring_orders.insert(name.text, order);
        Ok(())
    }

    pub(crate) fn cm_ring_order(&self, name: &Label) -> Result<cm_ring::Order> {
        self.cm_ring_orders
            .get(&name.text)
            .copied()
            .ok_or_else(|| Error::new(name.span, format!("unknown CM ring order '{}'", name.text)))
    }

    pub(crate) fn cm_class_type(&mut self, name: &Label) -> Result<TypeId> {
        self.cm_ring_order(name)?;
        self.inhabited_type(Kind::CMClass(name.text.as_str().into()), 0)
    }

    pub(crate) fn register_rv_phodge_quad(
        &mut self,
        name: Label,
        field: QuadField,
        spec: QuadraticSpec<QuadElement>,
    ) -> Result<()> {
        if self.rv_phodge_modules.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("filtered module '{}' is already defined", name.text),
            ));
        }
        let module = PHodgeModule::build_quadratic(self, field, spec)
            .map_err(|error| error.at(name.span))?;
        self.rv_phodge_modules.insert(name.text, Rc::new(module));
        Ok(())
    }

    pub(crate) fn rv_phodge(&self, name: &Label) -> Result<Rc<PHodgeModule>> {
        self.rv_phodge_modules
            .get(&name.text)
            .cloned()
            .ok_or_else(|| {
                Error::new(
                    name.span,
                    format!("unknown filtered module '{}'", name.text),
                )
            })
    }

    pub(crate) fn register_rv_phodge(&mut self, name: Label, spec: PHodgeSpec) -> Result<()> {
        if self.rv_phodge_modules.contains_key(&name.text) {
            return Err(Error::new(
                name.span,
                format!("filtered module '{}' is already defined", name.text),
            ));
        }
        let module = PHodgeModule::build(self, spec).map_err(|error| error.at(name.span))?;
        self.rv_phodge_modules.insert(name.text, Rc::new(module));
        Ok(())
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

    pub(crate) fn bp_mod_ext(&mut self, ty: TypeId) -> Result<Rc<bp_mod_ext::Group>> {
        if let Some(group) = &self.types[ty.0].bp_mod_ext {
            return Ok(group.clone());
        }
        let Kind::BPModExt(grade) = *self.kind(ty) else {
            unreachable!()
        };
        let group = Rc::new(bp_mod_ext::Group::build(self, grade)?);
        self.types[ty.0].bp_mod_ext = Some(group.clone());
        Ok(group)
    }

    pub(crate) fn bp_vpage(&mut self, ty: TypeId) -> Result<Rc<bp_vpage::Page>> {
        if let Some(data) = &self.types[ty.0].bp_vpage {
            return Ok(data.clone());
        }
        let Kind::BPVPage(grade) = *self.kind(ty) else {
            unreachable!()
        };
        let data = Rc::new(bp_vpage::Page::build(self, grade)?);
        self.types[ty.0].bp_vpage = Some(data.clone());
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
        let data = match *self.kind(ty) {
            Kind::MExt(grade) => Adams::build_motivic(self, grade)?,
            _ => {
                let (profile, s, t) = self.kind(ty).ext_grade().unwrap();
                Adams::build(self, profile, s, t)?
            }
        };
        let data = Rc::new(data);
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
            | Kind::MixedChain(..)
            | Kind::MixedHomology(..)
            | Kind::MixedMinimum(..)
            | Kind::MixedPage(..)
            | Kind::Sheaf(..)
            | Kind::SheafComplex(..)
            | Kind::SheafDeformation(..)
            | Kind::Convolution(..)
            | Kind::Parabolic(..)
            | Kind::DiskHom(..)
            | Kind::SheafHom(..)
            | Kind::SheafPage(..)
            | Kind::SheafLocal(..)
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
            | Kind::BPVPage(..)
            | Kind::BPModExt(..)
            | Kind::CobarAn(..)
            | Kind::MCobar(..)
            | Kind::RMCobar(..)
            | Kind::RMRhoCobar(..)
            | Kind::RhoLoc(..)
            | Kind::RhoMin(..)
            | Kind::Ext(..)
            | Kind::ExtAn(..)
            | Kind::MExt(..)
            | Kind::RMExt(..)
            | Kind::RMRhoExt(..)
            | Kind::RMBar(..)
            | Kind::RMTriple(..)
            | Kind::RhoBock(..)
            | Kind::XExt(..)
            | Kind::MapCobar(..)
            | Kind::MapExt(..)
            | Kind::QKoszul(..)
            | Kind::JetHom(..)
            | Kind::JetPage(..)
            | Kind::QExt(..)
            | Kind::RVHom(..)
            | Kind::CurveDiff(..)
            | Kind::ASOddDiff(..)
            | Kind::ASTwoDiff(..)
            | Kind::ASOddDeRham(..)
            | Kind::ASOddTate(..)
            | Kind::PrimeScalar(..)
            | Kind::MazurOrdinary(..)
            | Kind::ASDeRham(..)
            | Kind::BinaryScalar(..)
            | Kind::EOType(..)
            | Kind::RVLaurent(..)
            | Kind::RVIwHecke(..)
            | Kind::RVSphHecke(..)
            | Kind::RVSatake(..)
            | Kind::RVHive(..)
            | Kind::RVQGr(..)
            | Kind::RVQGrP(..)
            | Kind::RVOre
            | Kind::RVWittSpectrum
            | Kind::RVSchurCharacter(..)
            | Kind::RVStableSchurCharacter
            | Kind::RVPartitionAlgebra(..)
            | Kind::RVPartitionLocal(..)
            | Kind::RVKaroubiMap(..)
            | Kind::RVCTState(..)
            | Kind::RVCTResidue(..)
            | Kind::RVFusion(..)
            | Kind::RVCyclo(..)
            | Kind::ECCurve(..)
            | Kind::Jacobian(..)
            | Kind::Kummer(..)
            | Kind::KummerK3(..)
            | Kind::Enriques(..)
            | Kind::K3PencilBase(..)
            | Kind::K3PencilSection(..)
            | Kind::K3PencilKummer2(..)
            | Kind::K3PencilCoverPoint(..)
            | Kind::KummerDisc(..)
            | Kind::KummerWeil(..)
            | Kind::KummerDiv(..)
            | Kind::K3Lattice(..)
            | Kind::Nakajima(..)
            | Kind::HilbertCohomology(..)
            | Kind::HodgeTensor(..)
            | Kind::Verbitsky(..)
            | Kind::K3Clifford(..)
            | Kind::OG6Clifford(..)
            | Kind::KSExterior(..)
            | Kind::KSNomizu(..)
            | Kind::FiniteGroup(..)
            | Kind::PostnikovCochain(..)
            | Kind::PostnikovCell(..)
            | Kind::KSEnd(..)
            | Kind::QZ
            | Kind::K3Spinor(..)
            | Kind::Mukai(..)
            | Kind::MukaiQuotient(..)
            | Kind::OG10(..)
            | Kind::AbMukai(..)
            | Kind::OG6(..)
            | Kind::HilbertSquare(..)
            | Kind::SiegelWeight
            | Kind::ModularHecke(..)
            | Kind::HeckeEigen(..)
            | Kind::HeckeSymbol(..)
            | Kind::MazurTate(..)
            | Kind::MazurAmice(..)
            | Kind::PeriodPolynomial(..)
            | Kind::ECMap(..)
            | Kind::ECPoint(..)
            | Kind::ECLevel(..)
            | Kind::ECEndOrder(..)
            | Kind::ECKernel(..)
            | Kind::ECIso(..)
            | Kind::HWMatrix(..)
            | Kind::HWWMatrix(..)
            | Kind::HWScalar(..)
            | Kind::RVUnram(..)
            | Kind::RVCycloQ(..)
            | Kind::RVResSpace(..)
            | Kind::RVValMat(..)
            | Kind::RVHyperArithmetic(..)
            | Kind::RVHyperCover(..)
            | Kind::RVStableTree(..)
            | Kind::RVStableMarks(..)
            | Kind::RVQpQuad(..)
            | Kind::CMOrder(..)
            | Kind::RVWzw(..)
            | Kind::Lee(..)
            | Kind::CMClass(..)
            | Kind::LeeDual(..)
            | Kind::Khovanov(..)
            | Kind::RvTl(..)
            | Kind::VKoszul(..)
            | Kind::VExt(..)
            | Kind::VBock(..)
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
            | Kind::SetNat(..) => 1,
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
            Kind::FiniteGroup(name) => out.push(&format!("RVGroup({name})")),
            Kind::PostnikovCochain(name, degree) => {
                out.push(&format!("RVPostnikovCochain({name},{degree})"))
            }
            Kind::PostnikovCell(name) => out.push(&format!("RVPostnikovCell({name})")),
            Kind::SetNat(source, target) => out.push(&format!("SetNat({source}, {target})")),
            Kind::MixedHomology(grade) => out.push(&format!(
                "MixedH({}, {})",
                self.mixed_name(grade.complex),
                grade.b
            )),
            Kind::MixedChain(grade, total) => {
                let name = if *total { "MixedU" } else { "MixedChain" };
                out.push(&format!(
                    "{name}({}, {})",
                    self.mixed_name(grade.complex),
                    grade.b
                ))
            }
            Kind::MixedMinimum(grade, total) => {
                let name = if *total { "MixedMinU" } else { "MixedMin" };
                out.push(&format!(
                    "{name}({}, {})",
                    self.mixed_minimal_name(grade.minimal),
                    grade.b
                ))
            }
            Kind::MixedPage(grade) => out.push(&format!(
                "MixedE({}, {}, {}, {})",
                self.mixed_name(grade.chain.complex),
                grade.r,
                grade.chain.b,
                grade.f
            )),
            Kind::Sheaf(form, name, degree) => {
                out.push(&format!("{}({name}, {degree})", form.name()))
            }
            Kind::SheafDeformation(name, order) => out.push(&format!("ShDef({name}, {order})")),
            Kind::Convolution(form, name, _) => out.push(&format!("{}({name})", form.name())),
            Kind::Parabolic(form, name, _) => out.push(&format!("{}({name})", form.name())),
            Kind::DiskHom(from, to, _) => out.push(&format!("ShDiskHom({from}, {to})")),
            Kind::SheafHom(form, ends, degree) => {
                let name = match form {
                    SheafForm::Cochain => "ShHomC",
                    SheafForm::Class => "ShExt",
                };
                out.push(&format!(
                    "{name}({}, {}, {degree})",
                    ends.source, ends.target
                ))
            }
            Kind::SheafComplex(form, name, degree) => {
                out.push(&format!("{}({name}, {degree})", form.total_name()))
            }
            Kind::SheafLocal(name, grade) => {
                let local = self
                    .sheaf_supports_named(name)
                    .expect("registered support data");
                out.push(&format!(
                    "ShLocal({name}, {}, {}, {})",
                    local.point_name(grade.point),
                    grade.mode.name(),
                    grade.degree
                ))
            }
            Kind::SheafPage(name, grade) => out.push(&format!(
                "ShE({name}, {}, {}, {})",
                grade.r, grade.p, grade.q
            )),
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
            &Kind::MCobar(motivic::Grade { s, t, w }) => {
                out.push(&format!("MCobar({s}, {t}, {w})"))
            }
            &Kind::RhoMin(motivic::Grade { s, t, w }) => {
                out.push(&format!("RhoMin({s}, {t}, {w})"))
            }
            &Kind::RhoLoc(motivic::Grade { s, t, w }) => {
                out.push(&format!("RhoLoc({s}, {t}, {w})"))
            }
            &Kind::RMCobar(motivic::Grade { s, t, w }) => {
                out.push(&format!("RMCobar({s}, {t}, {w})"))
            }
            &Kind::MExt(motivic::Grade { s, t, w }) => out.push(&format!("MExt({s}, {t}, {w})")),
            &Kind::RMRhoCobar(real_rho::Grade {
                exponent,
                chain: motivic::Grade { s, t, w },
            }) => out.push(&format!("RMRhoCobar({exponent}, {s}, {t}, {w})")),
            &Kind::RMRhoExt(real_rho::Grade {
                exponent,
                chain: motivic::Grade { s, t, w },
            }) => out.push(&format!("RMRhoExt({exponent}, {s}, {t}, {w})")),
            &Kind::RMExt(motivic::Grade { s, t, w }) => out.push(&format!("RMExt({s}, {t}, {w})")),
            &Kind::RMBar(motivic::Grade { s, t, w }) => out.push(&format!("RMBar({s}, {t}, {w})")),
            &Kind::RMTriple(motivic::Grade { s, t, w }) => {
                out.push(&format!("RMTriple({s}, {t}, {w})"))
            }
            &Kind::RhoBock(rho_bock::Grade {
                chain: motivic::Grade { s, t, w },
                r,
                f,
            }) => out.push(&format!("RhoBock({r}, {s}, {t}, {w}, {f})")),
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
            Kind::RvCyclicChain(grade) | Kind::RvCyclicHomology(grade) => {
                let name = grade
                    .mode
                    .type_name(matches!(self.kind(ty), Kind::RvCyclicHomology(_)));
                out.push(&format!(
                    "{name}({}, {}, {:?})",
                    grade.ordinary.ideal, grade.ordinary.b, grade.ordinary.weight
                ))
            }
            Kind::RVNCChain(grade) | Kind::RvNegativeCyclic(grade) => {
                let name = if matches!(self.kind(ty), Kind::RVNCChain(_)) {
                    "RVNCChain"
                } else {
                    "RVNC"
                };
                out.push(&format!(
                    "{name}({}, {}, {}, {:?})",
                    grade.ordinary.ideal, grade.length, grade.ordinary.b, grade.ordinary.weight
                ))
            }
            Kind::RvHochschild(grade) => out.push(&format!(
                "RVHH({}, {}, {:?})",
                grade.ideal, grade.b, grade.weight
            )),
            Kind::RVHHChain(grade) => out.push(&format!(
                "RVHHChain({}, {}, {:?})",
                grade.ideal, grade.b, grade.weight
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
            Kind::ASDeRham(name, field) => {
                out.push(&format!("ASDeRham({name}, {})", field.modulus()))
            }
            Kind::BinaryScalar(field) => out.push(&format!("BinaryScalar({})", field.modulus())),
            Kind::EOType(genus) => out.push(&format!("EOType({genus})")),
            Kind::ASOddDiff(name, _) => out.push(&format!("ASOddDiff({name})")),
            Kind::ASTwoDiff(name, _) => out.push(&format!("ASTwoDiff({name})")),
            Kind::ASOddDeRham(name, _) => out.push(&format!("ASOddDeRham({name})")),
            Kind::ASOddTate(name, _, degree) => out.push(&format!("ASOddTate({name}, {degree})")),
            Kind::PrimeScalar(prime) => out.push(&format!("PrimeScalar({})", prime.get())),
            Kind::CurveDiff(family, name, field) => out.push(&format!(
                "{}({name}, {})",
                family.type_name(),
                field.modulus()
            )),
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
            Kind::BPVPage(grade) => out.push(&grade.to_string()),
            Kind::BPModExt(grade) => out.push(&grade.to_string()),
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
            Kind::RVSphHecke(rank) => out.push(&format!("RVSphHecke({rank})")),
            Kind::RVSatake(rank) => out.push(&format!("RVSatake({rank})")),
            Kind::RVFusion(grade) => {
                out.push(&format!("RVFusion({}, {})", grade.rank, grade.level))
            }
            Kind::RVWzw(grade) => out.push(&format!(
                "RVWzw({}, {})",
                grade.fusion.rank, grade.fusion.level
            )),
            Kind::Lee(name, grade) => out.push(&format!("LeeC({name}, {}, {})", grade.i, grade.f)),
            Kind::LeeDual(name, i) => out.push(&format!("LeeDual({name}, {i})")),
            Kind::Khovanov(name, grade, form) => {
                match form {
                    KhForm::Modular(prime, k) => {
                        return out.push(&format!(
                            "KhPC({name}, {}, {}, {}, {k})",
                            prime.get(),
                            grade.i,
                            grade.j
                        ));
                    }
                    KhForm::Page(spec) => {
                        return out.push(&format!(
                            "KhPage({name}, {}, {}, {}, {}, {})",
                            spec.prime.get(),
                            spec.r,
                            grade.i,
                            grade.j,
                            spec.k
                        ));
                    }
                    _ => {}
                }
                out.push(&format!("{}({name}, {}, {}", form.name(), grade.i, grade.j))?;
                if let Some(k) = form.annular().or_else(|| form.cutoff()) {
                    out.push(&format!(", {k}"))?;
                }
                out.push(")")
            }
            Kind::RvTl(grade) => out.push(&format!(
                "RVTL({}, {}, {})",
                grade.level, grade.strands, grade.height
            )),
            Kind::RVResSpace(grade) => out.push(&format!(
                "RVResSpace({}, {})",
                grade.prime.get(),
                grade.width,
            )),
            Kind::RVHyperArithmetic(grade) => out.push(&format!(
                "RVHyperArithmetic({}, {})",
                grade.prime.get(),
                grade.marks,
            )),
            Kind::RVHyperCover(grade) => out.push(&format!(
                "RVHyperCover({}, {})",
                grade.prime.get(),
                grade.marks
            )),
            Kind::RVStableTree(marks) => out.push(&format!("RVStableTree({marks})")),
            Kind::RVStableMarks(prime, marks) => {
                out.push(&format!("RVStableMarks({prime}, {marks})"))
            }
            Kind::RVValMat(grade) => out.push(&format!(
                "RVValMat({}, {}, {})",
                grade.prime, grade.ground, grade.rank
            )),
            Kind::RVQpQuad(name) => out.push(&format!("RVQpQuad({name})")),
            Kind::CMOrder(grade) => {
                out.push(&format!("RVCMOrder({}, {})", grade.family, grade.exponent))
            }
            Kind::ECEndOrder(name) => out.push(&format!("RVECEnd({name})")),
            Kind::CMClass(name) => out.push(&format!("RVCMClass({name})")),
            Kind::ECCurve(name) => out.push(&format!("RVECCurve({name})")),
            Kind::Jacobian(name) => out.push(&format!("RVJac({name})")),
            Kind::Kummer(name) => out.push(&format!("RVKummer({name})")),
            Kind::KummerK3(name) => out.push(&format!("RVKummerK3({name})")),
            Kind::Enriques(name) => out.push(&format!("RVEnriques({name})")),
            Kind::K3PencilBase(name) => out.push(&format!("RVK3PencilBase({name})")),
            Kind::K3PencilSection(name) => out.push(&format!("RVK3PencilSection({name})")),
            Kind::K3PencilKummer2(name) => out.push(&format!("RVK3PencilKummer2({name})")),
            Kind::K3PencilCoverPoint(name) => out.push(&format!("RVK3PencilCoverPoint({name})")),
            Kind::KummerDisc(name) => out.push(&format!("RVKummerDisc({name})")),
            Kind::KummerWeil(name) => out.push(&format!("RVKummerWeil({name})")),
            Kind::KummerDiv(name) => out.push(&format!("RVKummerDiv({name})")),
            Kind::K3Lattice(name) => out.push(&format!("RVK3Lattice({name})")),
            Kind::Verbitsky(grade) => out.push(&format!(
                "RVVerbitsky({},{},{})",
                grade.basis, grade.length, grade.degree
            )),
            Kind::HodgeTensor(grade) => out.push(&format!(
                "RVHodgeTensor({},{},{})",
                grade.basis, grade.length, grade.degree
            )),
            Kind::HilbertCohomology(grade) => out.push(&format!(
                "RVHilbCoh({},{},{})",
                grade.basis, grade.length, grade.degree
            )),
            Kind::Nakajima(grade) => out.push(&format!(
                "RVNakajima({},{},{})",
                grade.basis, grade.length, grade.degree
            )),
            Kind::K3Clifford(name, parity) => out.push(&format!("RVK3Clifford({name},{parity})")),
            Kind::OG6Clifford(name, parity) => out.push(&format!("RVOG6Clifford({name},{parity})")),
            Kind::KSExterior(name, degree) => out.push(&format!("RVKSExterior({name},{degree})")),
            Kind::KSNomizu(name, m, degree) => {
                out.push(&format!("RVKSNomizu({name},{m},{degree})"))
            }
            Kind::KSEnd(name) => out.push(&format!("RVKSEnd({name})")),
            Kind::QZ => out.push("RVQZ"),
            Kind::K3Spinor(name, chirality) => out.push(&format!("RVK3Spinor({name},{chirality})")),
            Kind::Mukai(name) => out.push(&format!("RVMukai({name})")),
            Kind::OG10(name) => out.push(&format!("RVOG10({name})")),
            Kind::AbMukai(name) => out.push(&format!("RVAbMukai({name})")),
            Kind::OG6(name) => out.push(&format!("RVOG6({name})")),
            Kind::MukaiQuotient(name) => out.push(&format!("RVMukaiQuotient({name})")),
            Kind::HilbertSquare(name, degree) => {
                out.push(&format!("RVHilbSquare({name}, {degree})"))
            }
            Kind::SiegelWeight => out.push("RVSiegelWeight"),
            Kind::ModularHecke(weight) => out.push(&format!("RVHecke({weight})")),
            Kind::HeckeEigen(name) => out.push(&format!("RVHeckeEigen({name})")),
            Kind::HeckeSymbol(name) => out.push(&format!("RVHeckeSymbol({name})")),
            Kind::MazurTate(grade) => out.push(&format!(
                "RVMazurTate({},{},{},{})",
                grade.period,
                grade.prime.get(),
                grade.level,
                grade.index,
            )),
            Kind::MazurAmice(grade) => out.push(&format!(
                "RVMazurAmice({},{},{},{})",
                grade.ordinary, grade.level, grade.index, grade.character
            )),
            Kind::MazurOrdinary(name, _) => out.push(&format!("RVMazurOrdinary({name})")),
            Kind::PeriodPolynomial(grade) => {
                out.push(&format!("RVPeriod({},{})", grade.weight, grade.sign.text()))
            }
            Kind::ECMap(a, b) => out.push(&format!("RVECMap({a}, {b})")),
            Kind::ECPoint(name) => out.push(&format!("RVECPoint({name})")),
            Kind::ECLevel(name, matrix) => {
                let label = if *matrix {
                    "RVECLevelMat"
                } else {
                    "RVECLevelVec"
                };
                out.push(&format!("{label}({name})"))
            }
            Kind::ECIso(name) => out.push(&format!("RVECIso({name})")),
            Kind::ECKernel(name, ell) => out.push(&format!("RVECKernel({name}, {ell})")),
            Kind::HWMatrix(ring, size) => out.push(&format!(
                "RVHWMatrix({}, {size}, {})",
                ring.prime.get(),
                ring.precision
            )),
            Kind::HWWMatrix(grade) => {
                out.push(&format!("RVHWWMatrix({}, {})", grade.ring, grade.genus))
            }
            Kind::HWScalar(ring) => out.push(&format!(
                "RVHWScalar({}, {})",
                ring.prime.get(),
                ring.precision
            )),
            Kind::RVUnram(name) => out.push(&format!("RVUnram({name})")),
            Kind::RVCycloQ(order) => out.push(&format!("RVCycloQ({order})")),
            Kind::RVCyclo(order) => out.push(&format!("RVCyclo({order})")),
            Kind::RVCTState(grade) => out.push(&format!(
                "RVCTState({}, {}, {})",
                grade.period, grade.level.prime, grade.level.precision
            )),
            Kind::RVCTResidue(level) => out.push(&format!(
                "RVCTResidue({}, {})",
                level.prime, level.precision
            )),
            Kind::RVOre => out.push("RVOre"),
            Kind::RVWittSpectrum => out.push("RVWittSpectrum"),
            Kind::RVSchurCharacter(weight) => out.push(&format!("RVSchurCharacter({weight})")),
            Kind::RVStableSchurCharacter => out.push("RVStableSchurCharacter"),
            Kind::RVPartitionAlgebra(rank) => out.push(&format!("RVPartitionAlgebra({rank})")),
            Kind::RVPartitionLocal(rank) => out.push(&format!("RVPartitionLocal({rank})")),
            Kind::RVKaroubiMap(a, b) => out.push(&format!("RVKaroubiMap({a}, {b})")),
            Kind::RVQGrP(grade) => out.push(&format!(
                "RVQGrP({}, {}, {})",
                grade.prime.get(),
                grade.grassmannian.rank,
                grade.grassmannian.ambient
            )),
            Kind::RVQGr(grade) => out.push(&format!("RVQGr({}, {})", grade.rank, grade.ambient)),
            Kind::RVHive(HiveBoundary([lambda, mu, nu])) => {
                out.push(&format!("RVHive({lambda:?}, {mu:?}, {nu:?})"))
            }
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
            (Value::BinaryScalar(left), Value::BinaryScalar(right)) => Ok(left == right),
            (Value::EOType(left), Value::EOType(right)) => {
                self.spend(left.len().min(right.len()))?;
                Ok(left == right)
            }
            (Value::PrimeScalar(left), Value::PrimeScalar(right)) => Ok(left == right),
            (Value::GroupElement(a), Value::GroupElement(b)) => Ok(a == b),
            (Value::PostnikovCell(a), Value::PostnikovCell(b)) => Ok(a == b),
            (Value::PostnikovCochain(a), Value::PostnikovCochain(b)) => {
                self.spend(a.0.len().min(b.0.len()))?;
                Ok(a == b)
            }
            (Value::SetNat(left), Value::SetNat(right)) => {
                self.spend(left.len().min(right.len()))?;
                Ok(left == right)
            }
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
            (Value::RhoMin(left), Value::RhoMin(right)) => {
                self.spend(left.len().min(right.len()))?;
                Ok(Rc::ptr_eq(left, right) || left == right)
            }
            (Value::RhoLoc(left), Value::RhoLoc(right)) => Ok(Rc::ptr_eq(left, right)
                || rho_localized::compare(self, left, right)? == Ordering::Equal),
            (Value::RMCobar(left), Value::RMCobar(right)) => Ok(Rc::ptr_eq(left, right)
                || real_motivic::compare(self, left, right)? == Ordering::Equal),
            (Value::SheafDeformation(left), Value::SheafDeformation(right)) => {
                Ok(Rc::ptr_eq(left, right) || left.compare(self, right)?.is_eq())
            }
            (Value::SheafC(left), Value::SheafC(right)) => {
                Ok(Rc::ptr_eq(left, right) || sheaf::compare(self, left, right)?.is_eq())
            }
            (Value::OddClass(left), Value::OddClass(right)) => {
                Ok(Rc::ptr_eq(left, right)
                    || odd_ext::compare(self, left, right)? == Ordering::Equal)
            }
            (Value::PCobar(left), Value::PCobar(right)) => Ok(Rc::ptr_eq(left, right)
                || odd_cobar::compare(self, left, right)? == Ordering::Equal),
            (Value::RMTriple(left), Value::RMTriple(right)) => {
                Ok(Rc::ptr_eq(left, right) || left.compare(self, right)?.is_eq())
            }
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
            (Value::RVNCChain(left), Value::RVNCChain(right)) => Ok(Rc::ptr_eq(left, right)
                || rv_negative::compare(self, left, right)? == Ordering::Equal),
            (Value::RVHHChain(left), Value::RVHHChain(right)) => Ok(
                Rc::ptr_eq(left, right) || rv_hh::compare(self, left, right)? == Ordering::Equal
            ),
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
            (Value::RVHive(left), Value::RVHive(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(true);
                }
                self.spend(left.len().min(right.len()))?;
                Ok(left == right)
            }
            (Value::RVFusion(left), Value::RVFusion(right)) => {
                Ok(rv_fusion::compare(self, left, right)?.is_eq())
            }
            (Value::RVWzw(left), Value::RVWzw(right)) => {
                Ok(rv_wzw::compare(self, left, right)?.is_eq())
            }
            (Value::LeeChain(left), Value::LeeChain(right)) => {
                self.spend(left.len().min(right.len()))?;
                Ok(left == right)
            }
            (Value::KhChain(left), Value::KhChain(right)) => {
                Ok(khovanov::compare_chains(self, left, right)?.is_eq())
            }
            (Value::KhClass(left), Value::KhClass(right)) => {
                Ok(khovanov::compare_classes(self, left, right)?.is_eq())
            }
            (Value::RvTl(left), Value::RvTl(right)) => {
                Ok(rv_tl::compare(self, left, right)?.is_eq())
            }
            (Value::RVQpQuad(left), Value::RVQpQuad(right)) => Ok(left == right),
            (Value::RVHyperArithmetic(left), Value::RVHyperArithmetic(right)) => {
                Ok(Rc::ptr_eq(left, right) || left.compare(self, right)?.is_eq())
            }
            (Value::RVHyperCover(left), Value::RVHyperCover(right)) => {
                Ok(Rc::ptr_eq(left, right) || left.compare(self, right)?.is_eq())
            }
            (Value::RVStableTree(left), Value::RVStableTree(right)) => {
                Ok(Rc::ptr_eq(left, right) || left.compare(self, right)?.is_eq())
            }
            (Value::RVStableMarks(left), Value::RVStableMarks(right)) => {
                Ok(Rc::ptr_eq(left, right) || left.compare(self, right)?.is_eq())
            }
            (Value::RVResSpace(left), Value::RVResSpace(right)) => {
                Ok(Rc::ptr_eq(left, right) || left.compare(self, right)?.is_eq())
            }
            (Value::RVValMat(left), Value::RVValMat(right)) => {
                Ok(Rc::ptr_eq(left, right) || left.compare(self, right)?.is_eq())
            }
            (Value::CMOrder(left), Value::CMOrder(right)) => Ok(left == right),
            (Value::ECEndOrder(left), Value::ECEndOrder(right)) => Ok(left == right),
            (Value::ECLevelVec(a), Value::ECLevelVec(b)) => Ok(a == b),
            (Value::ECLevelMat(a), Value::ECLevelMat(b)) => Ok(a == b),
            (Value::CMClass(left), Value::CMClass(right)) => Ok(left == right),
            (Value::HWMatrix(a), Value::HWMatrix(b)) => {
                self.spend(a.entries().len().min(b.entries().len()))?;
                Ok(a == b)
            }
            (Value::HWWMatrix(a), Value::HWWMatrix(b)) => {
                Ok(hw_closed::compare(self, a, b)?.is_eq())
            }
            (Value::HWScalar(a), Value::HWScalar(b)) => Ok(a == b),
            (Value::ECCurve(a), Value::ECCurve(b)) | (Value::ECIso(a), Value::ECIso(b)) => {
                Ok(a == b)
            }
            (Value::ECPoint(a), Value::ECPoint(b)) => {
                self.spend(1)?;
                Ok(a == b)
            }
            (Value::ECKernel(a), Value::ECKernel(b)) => {
                self.spend(a.polynomial.len().min(b.polynomial.len()) * a.polynomial[0].len())?;
                Ok(a.polynomial == b.polynomial)
            }
            (Value::Jacobian(a), Value::Jacobian(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::Kummer(a), Value::Kummer(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::KummerK3(a), Value::KummerK3(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::Enriques(a), Value::Enriques(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::K3PencilBase(a), Value::K3PencilBase(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::K3PencilCoverPoint(a), Value::K3PencilCoverPoint(b)) => {
                Ok(a.compare(self, b)?.is_eq())
            }
            (Value::K3PencilSection(a), Value::K3PencilSection(b)) => {
                Ok(a.compare(self, b)?.is_eq())
            }
            (Value::K3PencilKummer2(a), Value::K3PencilKummer2(b)) => {
                Ok(a.compare(self, b)?.is_eq())
            }
            (Value::KummerDisc(a), Value::KummerDisc(b)) => Ok(a == b),
            (Value::KummerWeil(a), Value::KummerWeil(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::KummerDiv(a), Value::KummerDiv(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::K3Lattice(a), Value::K3Lattice(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::RVWittSpectrum(a), Value::RVWittSpectrum(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::RVPartitionLocal(a), Value::RVPartitionLocal(b)) => {
                Ok(a.compare(self, b)?.is_eq())
            }
            (Value::RVKaroubiMap(a), Value::RVKaroubiMap(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::RVPartitionAlgebra(a), Value::RVPartitionAlgebra(b)) => {
                Ok(a.compare(self, b)?.is_eq())
            }
            (Value::RVStableSchurCharacter(a), Value::RVStableSchurCharacter(b)) => {
                Ok(a.compare(self, b)?.is_eq())
            }
            (Value::RVSchurCharacter(a), Value::RVSchurCharacter(b)) => {
                Ok(a.compare(self, b)?.is_eq())
            }
            (Value::Nakajima(a), Value::Nakajima(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::HodgeTensor(a), Value::HodgeTensor(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::HilbertCohomology(a), Value::HilbertCohomology(b)) => {
                Ok(a.compare(self, b)?.is_eq())
            }
            (Value::K3Clifford(a), Value::K3Clifford(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::OG6Clifford(a), Value::OG6Clifford(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::KSExterior(a), Value::KSExterior(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::KSNomizu(a), Value::KSNomizu(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::KSEnd(a), Value::KSEnd(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::QZ(a), Value::QZ(b)) => Ok(a.compare(self, *b)?.is_eq()),
            (Value::K3Spinor(a), Value::K3Spinor(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::OG10(a), Value::OG10(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::AbMukai(a), Value::AbMukai(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::OG6(a), Value::OG6(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::Mukai(a), Value::Mukai(b))
            | (Value::MukaiQuotient(a), Value::MukaiQuotient(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::HilbertSquare(a), Value::HilbertSquare(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::SiegelWeight(a), Value::SiegelWeight(b)) => Ok(a == b),
            (Value::ModularHecke(a), Value::ModularHecke(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::HeckeEigen(a), Value::HeckeEigen(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::HeckeSymbol(a), Value::HeckeSymbol(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::MazurTate(a), Value::MazurTate(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::MazurAmice(a), Value::MazurAmice(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::MazurOrdinary(a), Value::MazurOrdinary(b)) => Ok(a == b),
            (Value::PeriodPolynomial(a), Value::PeriodPolynomial(b)) => {
                Ok(a.compare(self, b)?.is_eq())
            }
            (Value::ECMap(a), Value::ECMap(b)) => Ok(a.compare(self, b)?.is_eq()),
            (Value::RVUnram(left), Value::RVUnram(right)) => {
                Ok(rv_unram::compare(self, left, right)?.is_eq())
            }
            (Value::RVCycloQ(left), Value::RVCycloQ(right)) => {
                Ok(rv_lens::compare(self, left, right)?.is_eq())
            }
            (Value::RVCyclo(left), Value::RVCyclo(right)) => {
                Ok(rv_cyclotomic::compare(self, left, right)?.is_eq())
            }
            (Value::PeriodState(left), Value::PeriodState(right)) => {
                Ok(Rc::ptr_eq(left, right) || rv_period::compare(self, left, right)?.is_eq())
            }
            (Value::PeriodResidue(left), Value::PeriodResidue(right)) => Ok(left == right),
            (Value::RVOre(left), Value::RVOre(right)) => {
                Ok(rv_ore::compare(self, left, right)?.is_eq())
            }
            (Value::RVQGrP(left), Value::RVQGrP(right)) => {
                Ok(rv_pcurv::compare(self, left, right)?.is_eq())
            }
            (Value::RVQGr(left), Value::RVQGr(right)) => {
                Ok(rv_qgr::compare(self, left, right)?.is_eq())
            }
            (Value::RVSatake(left), Value::RVSatake(right)) => Ok(Rc::ptr_eq(left, right)
                || rv_satake::compare(self, left, right)? == Ordering::Equal),
            (Value::RVSphHecke(left), Value::RVSphHecke(right)) => Ok(Rc::ptr_eq(left, right)
                || rv_sphhecke::compare(self, left, right)? == Ordering::Equal),
            (Value::RVIwHecke(left), Value::RVIwHecke(right)) => Ok(Rc::ptr_eq(left, right)
                || rv_iwhecke::compare(self, left, right)? == Ordering::Equal),
            (Value::RVLaurent(left), Value::RVLaurent(right)) => Ok(Rc::ptr_eq(left, right)
                || rv_laurent::compare(self, left, right)? == Ordering::Equal),
            (Value::RVHom(left), Value::RVHom(right))
            | (Value::CurveDiff(left), Value::CurveDiff(right))
            | (Value::ASDeRham(left), Value::ASDeRham(right))
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
            (Value::GroupElement(a), Value::GroupElement(b)) => Ok(a.cmp(b)),
            (Value::PostnikovCell(a), Value::PostnikovCell(b)) => Ok(a.cmp(b)),
            (Value::PostnikovCochain(a), Value::PostnikovCochain(b)) => {
                self.spend(a.0.len().min(b.0.len()))?;
                Ok(a.cmp(b))
            }
            (Value::PrimeScalar(left), Value::PrimeScalar(right)) => Ok(left.cmp(right)),
            (Value::SetNat(left), Value::SetNat(right)) => {
                self.spend(left.len().min(right.len()))?;
                Ok(left.cmp(right))
            }
            (Value::BinaryScalar(left), Value::BinaryScalar(right)) => Ok(left.cmp(right)),
            (Value::EOType(left), Value::EOType(right)) => {
                self.spend(left.len().min(right.len()))?;
                Ok(left.cmp(right))
            }
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
            (Value::RhoMin(left), Value::RhoMin(right)) => {
                self.spend(left.len().min(right.len()))?;
                Ok(left.cmp(right))
            }
            (Value::RhoLoc(left), Value::RhoLoc(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    rho_localized::compare(self, left, right)
                }
            }
            (Value::RMCobar(left), Value::RMCobar(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    real_motivic::compare(self, left, right)
                }
            }
            (Value::SheafDeformation(left), Value::SheafDeformation(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    left.compare(self, right)
                }
            }
            (Value::SheafC(left), Value::SheafC(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    sheaf::compare(self, left, right)
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
            (Value::RMTriple(left), Value::RMTriple(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    left.compare(self, right)
                }
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
            (Value::RVNCChain(left), Value::RVNCChain(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    rv_negative::compare(self, left, right)
                }
            }
            (Value::RVHHChain(left), Value::RVHHChain(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    rv_hh::compare(self, left, right)
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
            (Value::RVHive(left), Value::RVHive(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(Ordering::Equal);
                }
                self.spend(left.len().min(right.len()))?;
                Ok(left.cmp(right))
            }
            (Value::RVFusion(left), Value::RVFusion(right)) => {
                rv_fusion::compare(self, left, right)
            }
            (Value::RVWzw(left), Value::RVWzw(right)) => rv_wzw::compare(self, left, right),
            (Value::LeeChain(left), Value::LeeChain(right)) => {
                self.spend(left.len().min(right.len()))?;
                Ok(left.cmp(right))
            }
            (Value::KhChain(left), Value::KhChain(right)) => {
                khovanov::compare_chains(self, left, right)
            }
            (Value::KhClass(left), Value::KhClass(right)) => {
                khovanov::compare_classes(self, left, right)
            }
            (Value::RvTl(left), Value::RvTl(right)) => rv_tl::compare(self, left, right),
            (Value::RVQpQuad(left), Value::RVQpQuad(right)) => Ok(left.cmp(right)),
            (Value::RVHyperArithmetic(left), Value::RVHyperArithmetic(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    left.compare(self, right)
                }
            }
            (Value::RVHyperCover(left), Value::RVHyperCover(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    left.compare(self, right)
                }
            }
            (Value::RVStableTree(left), Value::RVStableTree(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    left.compare(self, right)
                }
            }
            (Value::RVStableMarks(left), Value::RVStableMarks(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    left.compare(self, right)
                }
            }
            (Value::RVResSpace(left), Value::RVResSpace(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    left.compare(self, right)
                }
            }
            (Value::RVValMat(left), Value::RVValMat(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    left.compare(self, right)
                }
            }
            (Value::CMOrder(left), Value::CMOrder(right)) => Ok(left.cmp(right)),
            (Value::ECEndOrder(left), Value::ECEndOrder(right)) => Ok(left.cmp(right)),
            (Value::ECLevelVec(a), Value::ECLevelVec(b)) => Ok(a.cmp(b)),
            (Value::ECLevelMat(a), Value::ECLevelMat(b)) => Ok(a.cmp(b)),
            (Value::CMClass(left), Value::CMClass(right)) => Ok(left.cmp(right)),
            (Value::HWMatrix(a), Value::HWMatrix(b)) => {
                self.spend(a.entries().len().min(b.entries().len()))?;
                Ok(a.entries().cmp(b.entries()))
            }
            (Value::HWWMatrix(a), Value::HWWMatrix(b)) => hw_closed::compare(self, a, b),
            (Value::HWScalar(a), Value::HWScalar(b)) => Ok(a.cmp(b)),
            (Value::ECCurve(a), Value::ECCurve(b)) | (Value::ECIso(a), Value::ECIso(b)) => {
                Ok(a.cmp(b))
            }
            (Value::ECPoint(a), Value::ECPoint(b)) => {
                self.spend(1)?;
                Ok(a.cmp(b))
            }
            (Value::ECKernel(a), Value::ECKernel(b)) => {
                self.spend(a.polynomial.len().min(b.polynomial.len()) * a.polynomial[0].len())?;
                Ok(a.polynomial.cmp(&b.polynomial))
            }
            (Value::Jacobian(a), Value::Jacobian(b)) => a.compare(self, b),
            (Value::Kummer(a), Value::Kummer(b)) => a.compare(self, b),
            (Value::KummerK3(a), Value::KummerK3(b)) => a.compare(self, b),
            (Value::Enriques(a), Value::Enriques(b)) => a.compare(self, b),
            (Value::K3PencilBase(a), Value::K3PencilBase(b)) => a.compare(self, b),
            (Value::K3PencilSection(a), Value::K3PencilSection(b)) => a.compare(self, b),
            (Value::K3PencilKummer2(a), Value::K3PencilKummer2(b)) => a.compare(self, b),
            (Value::K3PencilCoverPoint(a), Value::K3PencilCoverPoint(b)) => a.compare(self, b),
            (Value::KummerDisc(a), Value::KummerDisc(b)) => Ok(a.cmp(b)),
            (Value::KummerWeil(a), Value::KummerWeil(b)) => a.compare(self, b),
            (Value::KummerDiv(a), Value::KummerDiv(b)) => a.compare(self, b),
            (Value::K3Lattice(a), Value::K3Lattice(b)) => a.compare(self, b),
            (Value::Nakajima(a), Value::Nakajima(b)) => a.compare(self, b),
            (Value::HilbertCohomology(a), Value::HilbertCohomology(b)) => a.compare(self, b),
            (Value::HodgeTensor(a), Value::HodgeTensor(b)) => a.compare(self, b),
            (Value::K3Clifford(a), Value::K3Clifford(b)) => a.compare(self, b),
            (Value::OG6Clifford(a), Value::OG6Clifford(b)) => a.compare(self, b),
            (Value::KSExterior(a), Value::KSExterior(b)) => a.compare(self, b),
            (Value::KSNomizu(a), Value::KSNomizu(b)) => a.compare(self, b),
            (Value::KSEnd(a), Value::KSEnd(b)) => a.compare(self, b),
            (Value::QZ(a), Value::QZ(b)) => a.compare(self, *b),
            (Value::K3Spinor(a), Value::K3Spinor(b)) => a.compare(self, b),
            (Value::OG10(a), Value::OG10(b)) => a.compare(self, b),
            (Value::AbMukai(a), Value::AbMukai(b)) => a.compare(self, b),
            (Value::OG6(a), Value::OG6(b)) => a.compare(self, b),
            (Value::Mukai(a), Value::Mukai(b))
            | (Value::MukaiQuotient(a), Value::MukaiQuotient(b)) => a.compare(self, b),
            (Value::HilbertSquare(a), Value::HilbertSquare(b)) => a.compare(self, b),
            (Value::SiegelWeight(a), Value::SiegelWeight(b)) => Ok(a.cmp(b)),
            (Value::ModularHecke(a), Value::ModularHecke(b)) => a.compare(self, b),
            (Value::HeckeEigen(a), Value::HeckeEigen(b)) => a.compare(self, b),
            (Value::HeckeSymbol(a), Value::HeckeSymbol(b)) => a.compare(self, b),
            (Value::MazurTate(a), Value::MazurTate(b)) => a.compare(self, b),
            (Value::MazurAmice(a), Value::MazurAmice(b)) => a.compare(self, b),
            (Value::MazurOrdinary(a), Value::MazurOrdinary(b)) => Ok(a.cmp(b)),
            (Value::PeriodPolynomial(a), Value::PeriodPolynomial(b)) => a.compare(self, b),
            (Value::ECMap(a), Value::ECMap(b)) => a.compare(self, b),
            (Value::RVUnram(left), Value::RVUnram(right)) => rv_unram::compare(self, left, right),
            (Value::RVCycloQ(left), Value::RVCycloQ(right)) => rv_lens::compare(self, left, right),
            (Value::RVCyclo(left), Value::RVCyclo(right)) => {
                rv_cyclotomic::compare(self, left, right)
            }
            (Value::PeriodState(left), Value::PeriodState(right)) => {
                if Rc::ptr_eq(left, right) {
                    Ok(Ordering::Equal)
                } else {
                    rv_period::compare(self, left, right)
                }
            }
            (Value::PeriodResidue(left), Value::PeriodResidue(right)) => Ok(left.cmp(right)),
            (Value::RVOre(left), Value::RVOre(right)) => rv_ore::compare(self, left, right),
            (Value::RVWittSpectrum(a), Value::RVWittSpectrum(b)) => a.compare(self, b),
            (Value::RVSchurCharacter(a), Value::RVSchurCharacter(b)) => a.compare(self, b),
            (Value::RVPartitionAlgebra(a), Value::RVPartitionAlgebra(b)) => a.compare(self, b),
            (Value::RVPartitionLocal(a), Value::RVPartitionLocal(b)) => a.compare(self, b),
            (Value::RVKaroubiMap(a), Value::RVKaroubiMap(b)) => a.compare(self, b),
            (Value::RVStableSchurCharacter(a), Value::RVStableSchurCharacter(b)) => {
                a.compare(self, b)
            }
            (Value::RVQGrP(left), Value::RVQGrP(right)) => rv_pcurv::compare(self, left, right),
            (Value::RVQGr(left), Value::RVQGr(right)) => rv_qgr::compare(self, left, right),
            (Value::RVSatake(left), Value::RVSatake(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(Ordering::Equal);
                }
                rv_satake::compare(self, left, right)
            }
            (Value::RVSphHecke(left), Value::RVSphHecke(right)) => {
                if Rc::ptr_eq(left, right) {
                    return Ok(Ordering::Equal);
                }
                rv_sphhecke::compare(self, left, right)
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
            | (Value::CurveDiff(left), Value::CurveDiff(right))
            | (Value::ASDeRham(left), Value::ASDeRham(right))
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
                Kind::FiniteGroup(ref name) => {
                    let count = self.finite_groups[name.as_ref()].order();
                    self.allocate(count)?;
                    (0..count).map(Value::GroupElement).collect()
                }
                Kind::PostnikovCell(ref name) => {
                    let post = self.postnikov_systems[name.as_ref()].clone();
                    let count = post.order() * post.prime.get() as usize;
                    self.check_element_limit(count)?;
                    self.allocate(2 * count)?;
                    (0..post.order())
                        .flat_map(|object| {
                            (0..post.prime.get()).map(move |coordinate| {
                                Value::PostnikovCell(postnikov::Cell { object, coordinate })
                            })
                        })
                        .collect()
                }
                Kind::PostnikovCochain(ref name, degree) => {
                    let post = self.postnikov_systems[name.as_ref()].clone();
                    let dimension = post.size(degree);
                    let prime = post.prime.get() as usize;
                    let count = self.types[ty.0].cardinality_bound;
                    self.check_element_limit(count)?;
                    self.spend(dimension * count)?;
                    (0..count)
                        .map(|mut index| {
                            let mut values = Vec::new();
                            for j in 0..dimension {
                                let residue = (index % prime) as u32;
                                index /= prime;
                                if residue != 0 {
                                    values.push((j, residue));
                                }
                            }
                            postnikov::Cochain(values).persist(self)
                        })
                        .collect::<Result<Vec<_>>>()?
                }
                Kind::QModule(_) => {
                    return Err(Error::plain(
                        "coefficient modules are allowed only as Koszul endpoints",
                    ))
                }
                Kind::SetNat(ref source, ref target) => {
                    let source = self.set_family_named(source).expect("registered source");
                    let target = self.set_family_named(target).expect("registered target");
                    source.natural_elements(self, &target)?
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
                Kind::RMBar(grade) => {
                    let bars = real_extensions::bars(self, grade)?;
                    let living = bars.living(self, grade.t)?;
                    self.enumerate_supports(0, living.len(), |support| {
                        Value::Class(support.iter().map(|&i| living[i]).collect())
                    })?
                }
                Kind::RMExt(..) | Kind::RMRhoExt(..) => {
                    let rank = self.real_ext(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::RhoBock(..) => {
                    let rank = self.rho_page(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::Ext(..) | Kind::ExtAn(..) | Kind::MExt(..) => {
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
                Kind::RVSphHecke(..) => {
                    return Err(Error::plain("cannot enumerate a spherical Hecke algebra"))
                }
                Kind::RVWzw(..) => {
                    return Err(Error::plain(
                        "cannot enumerate a cyclotomic torus state module",
                    ))
                }
                Kind::Lee(ref name, grade) => {
                    let braid = self.kh_braid(name)?;
                    self.lee_elements(&braid, grade.i, Some(grade.f))?
                }
                Kind::LeeDual(ref name, i) => {
                    let braid = self.kh_braid(name)?;
                    self.lee_elements(&braid, i, None)?
                }
                Kind::Khovanov(ref name, grade, form) => {
                    let braid = self.kh_braid(name)?;
                    if let KhForm::Modular(prime, k) = form {
                        braid.modular_elements(self, grade, prime, k, MAX_ELEMENTS)?
                    } else if let KhForm::Page(spec) = form {
                        braid.page_elements(self, grade, spec, MAX_ELEMENTS)?
                    } else if let Some(k) = form.cutoff() {
                        braid.binary_elements(self, grade, k, MAX_ELEMENTS)?
                    } else if form.is_class() {
                        braid
                            .group(self, grade, form.annular())?
                            .elements(self, MAX_ELEMENTS)?
                    } else {
                        braid.chain_elements(self, grade, form.annular())?
                    }
                }
                Kind::RvTl(grade) => {
                    if grade.possible() {
                        return Err(Error::plain(
                            "cannot enumerate a cyclotomic Temperley–Lieb path sector",
                        ));
                    }
                    vec![rv_tl::zero()]
                }
                Kind::RVResSpace(grade) => {
                    let count = grade.count();
                    self.check_element_limit(count)?;
                    rv_res_space::enumerate(self, grade, count)?
                }
                Kind::RVStableTree(3) => vec![Value::RVStableTree(StableTree::three(self)?)],
                Kind::RVStableMarks(prime, 3) => {
                    vec![Value::RVStableMarks(StableCurve::three(self, prime)?)]
                }
                Kind::RVHyperArithmetic(..) => {
                    return Err(Error::plain(
                        "cannot enumerate arithmetic hyperelliptic covers",
                    ));
                }
                Kind::RVHyperCover(..) => {
                    return Err(Error::plain("cannot enumerate hyperelliptic metric graphs"));
                }
                Kind::RVStableTree(..) | Kind::RVStableMarks(..) => {
                    return Err(Error::plain(
                        "cannot enumerate the space of stable marked reductions",
                    ));
                }
                Kind::RVValMat(grade) if grade.rank == 0 || grade.rank == grade.ground => {
                    vec![Value::RVValMat(tropical_matroid::Matroid::build(
                        self,
                        grade,
                        |_, _| Ok(Some(0)),
                    )?)]
                }
                Kind::RVValMat(..) => {
                    return Err(Error::plain("cannot enumerate a tropical Grassmannian"))
                }
                Kind::RVQpQuad(..) => {
                    return Err(Error::plain("cannot enumerate a quadratic field"));
                }
                Kind::CMOrder(..) => {
                    return Err(Error::plain("cannot enumerate an integral CM order"));
                }
                Kind::ECEndOrder(..) => {
                    return Err(Error::plain(
                        "cannot enumerate an integral elliptic endomorphism order",
                    ));
                }
                Kind::CMClass(ref name) => self.cm_ring_orders[name.as_ref()]
                    .elements(self, MAX_ELEMENTS)?
                    .into_iter()
                    .map(Value::CMClass)
                    .collect(),
                Kind::HWScalar(ring) => {
                    let size = ring.count(1);
                    self.check_element_limit(size)?;
                    self.spend(size)?;
                    self.allocate(size)?;
                    (0..size).map(|c| Value::HWScalar(c as u64)).collect()
                }
                Kind::ECCurve(ref name) => vec![Value::ECCurve(name.clone())],
                Kind::ECIso(ref name) => vec![Value::ECIso(name.clone())],
                Kind::Jacobian(ref name) => {
                    let curve = self.jacobian_curve(name)?;
                    hyperelliptic_jacobian::elements(self, &curve, MAX_ELEMENTS)?
                }
                Kind::Kummer(ref name) => {
                    let curve = self.kummer_curve(name)?;
                    hyperelliptic_kummer_surface::elements(self, &curve, MAX_ELEMENTS)?
                }
                Kind::KummerK3(ref name) => {
                    let curve = self.kummer_curve(name)?;
                    hyperelliptic_kummer_surface::resolution::elements(self, &curve, MAX_ELEMENTS)?
                }
                Kind::Enriques(ref name) => {
                    self.enriques_quotient(name)?.elements(self, MAX_ELEMENTS)?
                }
                Kind::K3PencilSection(..) | Kind::K3PencilKummer2(..) => {
                    return Err(Error::plain(
                        "cannot enumerate a generic K3 pencil section or square-class type",
                    ));
                }
                Kind::K3PencilCoverPoint(..) => {
                    return Err(Error::plain(
                        "cannot enumerate a generic K3 pencil covering",
                    ));
                }
                Kind::K3PencilBase(ref name) => self.k3_pencils[name.as_ref()]
                    .clone()
                    .elements(self, MAX_ELEMENTS)?,
                Kind::KummerDisc(..) => {
                    let count = usize::from(nikulin::ORDER);
                    self.spend(count)?;
                    self.allocate(count)?;
                    (0..nikulin::ORDER).map(Value::KummerDisc).collect()
                }
                Kind::KummerWeil(..) => {
                    return Err(Error::plain(
                        "cannot enumerate a rational Kummer Weil representation",
                    ))
                }
                Kind::KummerDiv(..) => {
                    return Err(Error::plain(
                        "cannot enumerate an integral marked Kummer divisor lattice",
                    ));
                }
                Kind::K3Lattice(..) => {
                    return Err(Error::plain("cannot enumerate an integral K3 lattice"))
                }
                Kind::Mukai(..) => return Err(Error::plain("cannot enumerate a Mukai lattice")),
                Kind::AbMukai(..) => {
                    return Err(Error::plain("cannot enumerate an abelian Mukai lattice"))
                }
                Kind::OG6(..) => return Err(Error::plain("cannot enumerate an OG6 lattice")),
                Kind::OG10(..) => return Err(Error::plain("cannot enumerate an OG10 lattice")),
                Kind::MukaiQuotient(..) => {
                    return Err(Error::plain("cannot enumerate a Mukai quotient lattice"))
                }
                Kind::HodgeTensor(ref grade) | Kind::Verbitsky(ref grade) => {
                    if grade.supports_states() {
                        return Err(Error::plain(
                            "cannot enumerate a rational Hodge Wick component",
                        ));
                    }
                    vec![hyperelliptic_kummer_surface::hodge_wick::Tensor::default()
                        .persist(self)?]
                }
                Kind::HilbertCohomology(ref grade) => {
                    if grade.supports_states() {
                        return Err(Error::plain(
                            "cannot enumerate a rational Hilbert cohomology component",
                        ));
                    }
                    vec![nakajima::cohomology::Class::zero().persist(self)?]
                }
                Kind::Nakajima(ref grade) => {
                    if grade.supports_states() {
                        return Err(Error::plain(
                            "cannot enumerate an integral Nakajima component",
                        ));
                    }
                    vec![nakajima::State::zero().persist(self)?]
                }
                Kind::HilbertSquare(_, degree) if hilbert_square::dimension(degree) == 0 => {
                    vec![hilbert_square::Class::zero(degree).persist(self)?]
                }
                Kind::HilbertSquare(..) => {
                    return Err(Error::plain(
                        "cannot enumerate a rational Hilbert-square component",
                    ));
                }
                Kind::PeriodPolynomial(grade) if grade.is_trivial() => {
                    vec![period::Polynomial::zero(grade).persist(self)?]
                }
                Kind::PeriodPolynomial(..) => {
                    return Err(Error::plain(
                        "cannot enumerate rational Manin period polynomials",
                    ));
                }
                Kind::MazurTate(..) => {
                    return Err(Error::plain(
                        "cannot enumerate an exact Mazur–Tate group ring",
                    ));
                }
                Kind::MazurAmice(..) => {
                    return Err(Error::plain(
                        "cannot enumerate an Amice cyclotomic group ring",
                    ));
                }
                Kind::MazurOrdinary(_, modulus) => {
                    let size = usize::try_from(modulus).map_err(|_| {
                        Error::plain(
                            "ordinary Mazur–Tate residue ring exceeds the enumeration limit",
                        )
                    })?;
                    if size > MAX_ELEMENTS {
                        return Err(Error::plain(
                            "ordinary Mazur–Tate residue ring exceeds the enumeration limit",
                        ));
                    }
                    self.spend(size)?;
                    self.allocate(size)?;
                    (0..modulus).map(Value::MazurOrdinary).collect()
                }
                Kind::HeckeSymbol(..) => {
                    return Err(Error::plain(
                        "cannot enumerate a rational Hecke modular-symbol space",
                    ));
                }
                Kind::HeckeEigen(..) => {
                    return Err(Error::plain(
                        "cannot enumerate a rational Hecke eigenpacket algebra",
                    ));
                }
                Kind::ModularHecke(..) => {
                    return Err(Error::plain("cannot enumerate the level-one Hecke algebra"));
                }
                Kind::SiegelWeight => {
                    return Err(Error::plain("cannot enumerate the Siegel weight lattice"))
                }
                Kind::ECMap(..) => {
                    return Err(Error::plain(
                        "an elliptic rational morphism type cannot be enumerated",
                    ));
                }
                Kind::OG6Clifford(..) => {
                    return Err(Error::plain(
                        "cannot enumerate an OG6 integral Clifford component",
                    ))
                }
                Kind::KSExterior(_, 0..=64) => {
                    return Err(Error::plain(
                        "cannot enumerate an integral Kuga–Satake exterior component",
                    ));
                }
                Kind::KSExterior(..) => vec![ks_exterior::Form::zero().persist(self)?],
                Kind::KSNomizu(_, _, 0..=65) => {
                    return Err(Error::plain(
                        "cannot enumerate a rational Kuga–Satake Nomizu cochain space",
                    ))
                }
                Kind::KSNomizu(..) => vec![nomizu::Cochain::zero().persist(self)?],
                Kind::KSEnd(..) => {
                    return Err(Error::plain(
                        "cannot enumerate a rational Kuga–Satake endomorphism space",
                    ));
                }
                Kind::QZ => return Err(Error::plain("cannot enumerate the Q/Z scalar group")),
                Kind::K3Clifford(..) => {
                    return Err(Error::plain(
                        "cannot enumerate the rational K3 Clifford algebra",
                    ))
                }
                Kind::K3Spinor(..) => {
                    return Err(Error::plain(
                        "cannot enumerate a rational K3 half-spin module",
                    ))
                }
                Kind::ECPoint(ref name) => self.ec_curve(name)?.elements(self, MAX_ELEMENTS)?,
                Kind::ECLevel(ref name, matrix) => {
                    let n = self.ec_torsion_levels[name.as_ref()].modulus;
                    let count = (n as usize).saturating_pow(if matrix { 4 } else { 2 });
                    self.check_element_limit(count)?;
                    self.spend(count)?;
                    self.allocate(count * if matrix { 4 } else { 2 })?;
                    (0..count)
                        .map(|mut code| {
                            let mut next = || {
                                let value = code as u64 % n;
                                code /= n as usize;
                                value
                            };
                            let a = [next(), next()];
                            if matrix {
                                Value::ECLevelMat(ec_level::Matrix([a, [next(), next()]]))
                            } else {
                                Value::ECLevelVec(a)
                            }
                        })
                        .collect()
                }
                Kind::ECKernel(ref name, ell) => {
                    let curve = self.ec_curve(name)?;
                    let mut kernels = elliptic::enumerate_kernels(self, &curve, ell as usize)?;
                    self.check_element_limit(kernels.len())?;
                    self.spend(kernels.len() * (1 + kernels.len().max(1).ilog2() as usize))?;
                    kernels.sort_unstable();
                    kernels
                        .into_iter()
                        .map(|k| {
                            let k = ECKernel::from_certified(self, &curve, ell, k)?;
                            Ok(Value::ECKernel(Rc::new(k)))
                        })
                        .collect::<Result<_>>()?
                }
                Kind::HWMatrix(ring, size) => ring.matrix_elements(self, size, MAX_ELEMENTS)?,
                Kind::HWWMatrix(ref grade) => {
                    let ring = self
                        .rv_unram_named(&grade.ring)
                        .expect("registered unramified ring");
                    let count = hw_closed::count(&ring, grade.genus);
                    self.check_element_limit(count)?;
                    hw_closed::enumerate(self, &ring, grade.genus, count)?
                }
                Kind::RVUnram(ref name) => {
                    let ring = self
                        .rv_unram_named(name)
                        .expect("registered unramified ring");
                    let count = ring.count();
                    self.check_element_limit(count)?;
                    rv_unram::enumerate(self, &ring, count)?
                }
                Kind::RVCycloQ(..) => {
                    return Err(Error::plain("cannot enumerate a rational cyclotomic field"))
                }
                Kind::RVCyclo(..) => {
                    return Err(Error::plain("cannot enumerate an integral cyclotomic ring"));
                }
                Kind::RVFusion(..) => {
                    return Err(Error::plain("cannot enumerate a Verlinde fusion ring"));
                }
                Kind::RVCTState(ref grade) => {
                    self.rv_period_for(grade)
                        .elements(self, grade.level, MAX_ELEMENTS)?
                }
                Kind::RVCTResidue(level) => {
                    let size = usize::try_from(level.modulus)
                        .ok()
                        .filter(|&n| n <= MAX_ELEMENTS)
                        .ok_or_else(|| {
                            Error::plain(
                                "Cartier residues exceed the finite type enumeration limit",
                            )
                        })?;
                    self.spend(size)?;
                    self.allocate(size)?;
                    (0..level.modulus).map(Value::PeriodResidue).collect()
                }
                Kind::RVOre => return Err(Error::plain("cannot enumerate an Ore operator ring")),
                Kind::RVWittSpectrum => {
                    return Err(Error::plain("cannot enumerate a Witt spectrum ring"))
                }
                Kind::RVPartitionAlgebra(..)
                | Kind::RVPartitionLocal(..)
                | Kind::RVKaroubiMap(..) => {
                    return Err(Error::plain(
                        "cannot enumerate a partition algebra or its Karoubi maps",
                    ))
                }
                Kind::RVSchurCharacter(..) | Kind::RVStableSchurCharacter => {
                    return Err(Error::plain("cannot enumerate a Schur character group"))
                }
                Kind::RVQGrP(..) => {
                    return Err(Error::plain(
                        "cannot enumerate a polynomial Schubert module",
                    ))
                }
                Kind::RVQGr(..) => {
                    return Err(Error::plain("cannot enumerate a quantum Grassmannian ring"))
                }
                Kind::RVHive(ref boundary) => {
                    let boundary = boundary.clone();
                    rv_hive::elements(self, &boundary, MAX_ELEMENTS)?
                }
                Kind::RVSatake(..) => {
                    return Err(Error::plain(
                        "cannot enumerate a symmetric Laurent polynomial ring",
                    ))
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
                Kind::ASDeRham(ref name, field) => {
                    let dimension = 2 * self.rv_as_covers[name.as_ref()].genus();
                    binary_linear::elements(self, field, dimension, MAX_ELEMENTS, Value::ASDeRham)?
                }
                Kind::EOType(genus) => {
                    let size = support_space_size(genus);
                    self.check_element_limit(size)?;
                    self.spend(size * (genus + 1))?;
                    self.allocate(size * (genus + 1))?;
                    (0..size)
                        .map(|mask| {
                            let mut rank = 0;
                            let sequence = (0..genus)
                                .rev()
                                .map(|position| {
                                    rank += (mask >> position) & 1;
                                    rank
                                })
                                .collect();
                            Value::EOType(sequence)
                        })
                        .collect()
                }
                Kind::BinaryScalar(field) => {
                    let size = field.order();
                    self.check_element_limit(size)?;
                    self.spend(size)?;
                    self.allocate(size)?;
                    (0..size).map(Value::BinaryScalar).collect()
                }
                Kind::PrimeScalar(prime) => {
                    let size = prime.get() as usize;
                    self.check_element_limit(size)?;
                    self.spend(size)?;
                    self.allocate(size)?;
                    (0..prime.get()).map(Value::PrimeScalar).collect()
                }
                Kind::ASOddTate(ref name, prime, degree) => {
                    let cover = self.rv_as_odd_covers[name.as_ref()].clone();
                    let rank = cover.tate(self)?.group(self, &cover, degree)?.rank();
                    odd_ext::elements(self, prime, rank, MAX_ELEMENTS)?
                }
                Kind::ASOddDeRham(ref name, prime) => {
                    let dimension = 2 * self.rv_as_odd_covers[name.as_ref()].genus();
                    odd_ext::elements(self, prime, dimension, MAX_ELEMENTS)?
                }
                Kind::ASTwoDiff(ref name, prime) => {
                    let genus = self.rv_as_two_covers[name.as_ref()].genus();
                    odd_ext::elements(self, prime, genus, MAX_ELEMENTS)?
                }
                Kind::ASOddDiff(ref name, prime) => {
                    let genus = self.rv_as_odd_covers[name.as_ref()].genus();
                    odd_ext::elements(self, prime, genus, MAX_ELEMENTS)?
                }
                Kind::CurveDiff(family, ref name, field) => {
                    let genus = self.curve_genus(family, name);
                    binary_linear::elements(self, field, genus, MAX_ELEMENTS, Value::CurveDiff)?
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
                Kind::RMTriple(..) => {
                    return Err(Error::plain(
                        "RMTriple types do not enumerate ambient cosets",
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
                Kind::RvCyclicChain(ref grade) => {
                    let grade = grade.clone();
                    let basis =
                        rv_cyclic_homology::basis(self, &grade, MAX_ELEMENTS.ilog2() as usize)?;
                    self.allocate(basis.iter().map(rv_negative::term_cells).sum())?;
                    self.enumerate_supports(0, basis.len(), |class| {
                        Value::RVNCChain(class.iter().map(|&i| basis[i].clone()).collect())
                    })?
                }
                Kind::RvCyclicHomology(..) => {
                    let rank = self.rv_cyclic_homology(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::RVNCChain(ref grade) => {
                    let grade = grade.clone();
                    let basis = rv_negative::basis(
                        self,
                        &grade,
                        grade.ordinary.b.into(),
                        MAX_ELEMENTS.ilog2() as usize,
                    )?;
                    self.allocate(basis.iter().map(rv_negative::term_cells).sum())?;
                    self.enumerate_supports(0, basis.len(), |class| {
                        Value::RVNCChain(class.iter().map(|&i| basis[i].clone()).collect())
                    })?
                }
                Kind::RvNegativeCyclic(..) => {
                    let rank = self.rv_negative_homology(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::RvHochschild(..) => {
                    let rank = self.rv_hh_homology(ty)?.rank();
                    self.enumerate_supports(0, rank, Value::Class)?
                }
                Kind::RVHHChain(ref grade) => {
                    let grade = grade.clone();
                    let basis = rv_hh::basis(self, &grade, MAX_ELEMENTS.ilog2() as usize)?;
                    self.allocate(basis.iter().map(|word| rv_hh::word_cells(word)).sum())?;
                    self.enumerate_supports(0, basis.len(), |class| {
                        Value::RVHHChain(class.iter().map(|&item| basis[item].clone()).collect())
                    })?
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
                Kind::MixedChain(grade, total) => {
                    self.mixed_complex(grade.complex)
                        .elements(self, grade, total, MAX_ELEMENTS)?
                }
                Kind::MixedMinimum(grade, total) => self.mixed_minimal(grade.minimal).elements(
                    self,
                    grade.b,
                    total,
                    MAX_ELEMENTS,
                )?,
                Kind::MixedHomology(grade) => {
                    let rank = self.mixed_homology(grade)?.rank();
                    let prime = self.mixed_complex(grade.complex).prime;
                    odd_ext::elements(self, prime, rank, MAX_ELEMENTS)?
                }
                Kind::MixedPage(grade) => {
                    let rank = self.mixed_page(ty)?.rank();
                    odd_ext::elements(
                        self,
                        self.mixed_complex(grade.chain.complex).prime,
                        rank,
                        MAX_ELEMENTS,
                    )?
                }
                Kind::Sheaf(form, ref name, degree) => {
                    let data = self.sheaf_named(name).expect("registered sheaf");
                    match form {
                        SheafForm::Cochain => data.elements(self, degree, MAX_ELEMENTS)?,
                        SheafForm::Class => {
                            let rank = data.rank(self, degree)?;
                            odd_ext::elements(self, data.prime, rank, MAX_ELEMENTS)?
                        }
                    }
                }
                Kind::SheafDeformation(ref name, order) => {
                    let hom = self.sheaf_hom(&SheafEnds::endomorphisms(name.clone()));
                    sheaf_deform::elements(self, &hom, order, MAX_ELEMENTS)?
                }
                Kind::Convolution(form, ref name, prime) => {
                    let dimension = self.convolutions[name].dimension(form);
                    odd_ext::elements(self, prime, dimension, MAX_ELEMENTS)?
                }
                Kind::Parabolic(form, ref name, prime) => {
                    let dimension = self.parabolic_groups[name].dimension(form);
                    odd_ext::elements(self, prime, dimension, MAX_ELEMENTS)?
                }
                Kind::DiskHom(ref from, ref to, prime) => {
                    let rank = self.disk_homs[&(from.clone(), to.clone())].rank();
                    odd_ext::elements(self, prime, rank, MAX_ELEMENTS)?
                }
                Kind::SheafHom(form, ref ends, degree) => {
                    let data = self.sheaf_hom(ends);
                    match form {
                        SheafForm::Cochain => data.elements(self, degree, MAX_ELEMENTS)?,
                        SheafForm::Class => {
                            let rank = data.cohomology.rank(self, data.as_ref(), degree)?;
                            odd_ext::elements(self, data.prime(), rank, MAX_ELEMENTS)?
                        }
                    }
                }
                Kind::SheafComplex(form, ref name, degree) => {
                    let data = self
                        .sheaf_total_named(name)
                        .expect("registered sheaf complex");
                    match form {
                        SheafForm::Cochain => {
                            data.filtered().elements(self, degree, MAX_ELEMENTS)?
                        }
                        SheafForm::Class => {
                            let rank = data.rank(self, degree)?;
                            odd_ext::elements(self, data.filtered().prime(), rank, MAX_ELEMENTS)?
                        }
                    }
                }
                Kind::SheafLocal(ref name, grade) => {
                    let data = self
                        .sheaf_supports_named(name)
                        .expect("registered support data");
                    let rank = data.rank(self, grade)?;
                    odd_ext::elements(self, data.prime(), rank, MAX_ELEMENTS)?
                }
                Kind::SheafPage(ref name, grade) => {
                    let data = self
                        .sheaf_total_named(name)
                        .expect("registered sheaf complex");
                    let rank = data.filtered().page(self, grade)?.rank();
                    odd_ext::elements(self, data.filtered().prime(), rank, MAX_ELEMENTS)?
                }
                Kind::Operations(degree) => steenrod::elements(self, degree, MAX_ELEMENTS)?,
                Kind::RhoMin(grade) => {
                    let basis = rho_minimal::basis(self, grade, MAX_ELEMENTS.ilog2() as usize)?;
                    self.enumerate_supports(0, basis.len(), |support| {
                        Value::RhoMin(support.iter().map(|&index| basis[index]).collect())
                    })?
                }
                Kind::RhoLoc(grade) => {
                    let basis = rho_localized::basis(self, grade, MAX_ELEMENTS.ilog2() as usize)?;
                    self.enumerate_supports(0, basis.len(), |support| {
                        Value::RhoLoc(support.iter().map(|&index| basis[index].clone()).collect())
                    })?
                }
                Kind::RMRhoCobar(grade) => {
                    let basis = real_rho::basis(self, grade, MAX_ELEMENTS.ilog2() as usize)?;
                    self.enumerate_supports(0, basis.len(), |support| {
                        Value::RMCobar(support.iter().map(|&i| basis[i].clone()).collect())
                    })?
                }
                Kind::RMCobar(grade) => {
                    let basis = real_motivic::basis(self, grade, MAX_ELEMENTS.ilog2() as usize)?;
                    self.enumerate_supports(0, basis.len(), |support| {
                        Value::RMCobar(support.iter().map(|&index| basis[index].clone()).collect())
                    })?
                }
                Kind::MCobar(grade) => {
                    let basis = motivic::basis(self, grade, MAX_ELEMENTS.ilog2() as usize)?;
                    self.enumerate_supports(0, basis.len(), |support| {
                        Value::Cobar(support.iter().map(|&index| basis[index].clone()).collect())
                    })?
                }
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
                Kind::BPModExt(grade) => {
                    let rank = self.bp_mod_ext(ty)?.rank();
                    odd_ext::elements(self, grade.prime(), rank, MAX_ELEMENTS)?
                }
                Kind::BPVPage(grade) => {
                    let rank = self.bp_vpage(ty)?.rank();
                    odd_ext::elements(self, grade.prime(), rank, MAX_ELEMENTS)?
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

    pub(crate) fn show_satake_coefficient(
        &self,
        coefficient: &[rv_satake::Term],
    ) -> Result<String> {
        let mut out = Text::new(MAX_SHOWN_BYTES);
        out.push("rvSatakeCoeff = {")?;
        for (index, term) in coefficient.iter().enumerate() {
            if index != 0 {
                out.push(", ")?;
            }
            out.push(&format!("{}: {}", term.power, term.coefficient))?;
        }
        out.push("}")?;
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

    fn write_hecke<'a>(
        &mut self,
        out: &mut Text,
        ty: TypeId,
        basis: &str,
        terms: impl Iterator<Item = (&'a [i128], &'a [(usize, i128)])>,
    ) -> Result<()> {
        self.write_type(out, ty, 0)?;
        out.push("{")?;
        let mut first = true;
        for (index, polynomial) in terms {
            for &(degree, coefficient) in polynomial {
                self.spend(1 + index.len())?;
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
                out.push(&format!("{basis}{index:?}"))?;
            }
        }
        if first {
            out.push("0")?;
        }
        out.push("}")
    }

    fn write_schubert<C: std::fmt::Display>(
        &mut self,
        out: &mut Text,
        ty: TypeId,
        terms: &[(rv_qgr::Schubert, C)],
    ) -> Result<()> {
        self.write_type(out, ty, 0)?;
        out.push("{")?;
        for (index, (term, coefficient)) in terms.iter().enumerate() {
            self.spend(1 + term.shape.len())?;
            if index != 0 {
                out.push(", ")?;
            }
            out.push(&format!(
                "({}, {:?}): {coefficient}",
                term.quantum, term.shape
            ))?;
        }
        out.push("}")
    }

    fn write_cube_terms<T: std::fmt::Display>(
        &mut self,
        out: &mut Text,
        ty: TypeId,
        terms: &[(khovanov::Cell, T)],
    ) -> Result<()> {
        self.write_type(out, ty, 0)?;
        out.push("{")?;
        if terms.is_empty() {
            out.push("0")?;
        }
        for (i, (cell, coefficient)) in terms.iter().enumerate() {
            self.spend(1 + cell.xs.len() + cell.state.count_ones() as usize)?;
            if i != 0 {
                out.push(", ")?;
            }
            let mut state = cell.state;
            let mut ones = Vec::new();
            while state != 0 {
                ones.push(state.trailing_zeros());
                state &= state - 1;
            }
            out.push(&format!("({ones:?}, {:?}): {coefficient}", cell.xs))?;
        }
        out.push("}")
    }

    fn write_sparse_rational<D: rational_sparse::Domain, K: Copy + Ord + std::fmt::Display>(
        &mut self,
        out: &mut Text,
        ty: TypeId,
        element: &rational_sparse::Coefficients<D, K>,
    ) -> Result<()> {
        self.write_type(out, ty, 0)?;
        out.push("{[")?;
        for (i, &(index, numerator)) in element.terms().iter().enumerate() {
            self.spend(1)?;
            if i != 0 {
                out.push(", ")?;
            }
            out.push(&format!(
                "({index}, {})",
                rational_sparse::signed_text(numerator)
            ))?;
        }
        out.push(&format!("] / {}}}", element.denominator_text()))
    }

    fn write_nakajima_terms(
        &mut self,
        out: &mut Text,
        state: &nakajima::State,
        denominator: u128,
    ) -> Result<()> {
        out.push("{")?;
        for (i, (monomial, coefficient)) in state.terms().enumerate() {
            self.spend(1 + monomial.len())?;
            if i > 0 {
                out.push(", ")?;
            }
            if denominator == 1 {
                out.push(&format!("{coefficient}*["))?;
            } else {
                out.push(&format!("{coefficient}/{denominator}*["))?;
            }
            for (j, factor) in monomial.iter().enumerate() {
                if j > 0 {
                    out.push(", ")?;
                }
                out.push(&format!(
                    "({},{},{})",
                    factor.mode.length, factor.mode.color, factor.exponent
                ))?;
            }
            out.push("]")?;
        }
        out.push("}")
    }

    fn write_hecke_element(out: &mut Text, element: &eigenpacket::Element) -> Result<()> {
        out.push("[")?;
        for (i, coefficient) in element.coordinates().enumerate() {
            if i != 0 {
                out.push(", ")?;
            }
            out.push(&coefficient)?;
        }
        out.push("]/")?;
        out.push(&element.denominator())
    }

    fn write_value(&mut self, out: &mut Text, ty: TypeId, value: &Value) -> Result<()> {
        match value {
            Value::GroupElement(i) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{{i}}}"))
            }
            Value::PostnikovCell(cell) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{{}, {}}}", cell.object, cell.coordinate))
            }
            Value::PostnikovCochain(cochain) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                self.spend(cochain.0.len())?;
                out.push(&format!("{:?}", cochain.0))?;
                out.push("}")
            }
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
            Value::SheafDeformation(jet) => {
                let Kind::SheafDeformation(name, _) = self.kind(ty) else {
                    unreachable!()
                };
                let hom = self.type_ids[&Kind::SheafHom(
                    SheafForm::Cochain,
                    SheafEnds::endomorphisms(name.clone()),
                    1,
                )];
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                for (i, terms) in jet.coefficients().iter().enumerate() {
                    if i != 0 {
                        out.push(", ")?;
                    }
                    out.push(&format!("t^{}=", i + 1))?;
                    self.write_value(out, hom, &Value::SheafC(terms.clone()))?;
                }
                out.push("}")
            }
            Value::SheafC(terms) => {
                let (poset, total, hom) = match self.kind(ty) {
                    Kind::SheafHom(_, ends, _) => {
                        let data = self.sheaf_hom(ends);
                        (data.source.poset.clone(), None, Some(data))
                    }
                    Kind::Sheaf(_, name, _) => (
                        self.sheaf_named(name)
                            .expect("registered sheaf")
                            .poset
                            .clone(),
                        None,
                        None,
                    ),
                    Kind::SheafComplex(_, name, _) => {
                        let total = self
                            .sheaf_total_named(name)
                            .expect("registered sheaf complex");
                        (total.poset().clone(), Some(total), None)
                    }
                    _ => unreachable!(),
                };
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                for (index, (generator, coefficient)) in terms.iter().enumerate() {
                    if index != 0 {
                        out.push(" + ")?;
                    }
                    if *coefficient != 1 {
                        out.push(&format!("{coefficient}*"))?;
                    }
                    if total.is_some() {
                        out.push(&format!("q{}:", generator.vertical))?;
                    }
                    match &generator.place {
                        sheaf::Place::Chain(chain) => {
                            out.push("[")?;
                            for (i, &point) in chain.iter().enumerate() {
                                if i != 0 {
                                    out.push("<")?;
                                }
                                out.push(&poset.names[point])?;
                            }
                            out.push("]")?;
                        }
                        sheaf::Place::Cell { ordinal, .. } => {
                            out.push(total.as_ref().unwrap().cell_label(*ordinal))?;
                        }
                    }
                    if let Some(hom) = &hom {
                        let (column, row) = hom.decode(generator);
                        out.push(&format!("#{column}->{row}"))?;
                    } else {
                        out.push(&format!("#{}", generator.index))?;
                    }
                }
                out.push("}")
            }
            Value::PrimeScalar(coefficient) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{{coefficient}}}"))
            }
            Value::OddClass(support) => {
                let mixed = if let Kind::MixedChain(grade, total) = *self.kind(ty) {
                    Some((self.mixed_complex(grade.complex), grade.b, total))
                } else {
                    None
                };
                let two_chart = match self.kind(ty) {
                    Kind::ASTwoDiff(name, _) if !support.is_empty() => {
                        let cover = self
                            .rv_as_two_cover_named(name)
                            .expect("registered two-pole cover");
                        self.spend(support.len())?;
                        Some(cover.cartier(self)?)
                    }
                    _ => None,
                };
                let chart = match self.kind(ty) {
                    Kind::ASOddDiff(name, _) | Kind::ASOddDeRham(name, _)
                        if !support.is_empty() =>
                    {
                        let de_rham = matches!(self.kind(ty), Kind::ASOddDeRham(..));
                        let cover = self
                            .rv_as_odd_cover_named(name)
                            .expect("registered odd Artin–Schreier cover");
                        self.spend(support.len())?;
                        Some((cover.cartier(self)?, cover.genus(), de_rham))
                    }
                    _ => None,
                };
                let minimal = if let Kind::MixedMinimum(grade, true) = *self.kind(ty) {
                    Some((self.mixed_minimal(grade.minimal), grade.b))
                } else {
                    None
                };
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
                    if let Some((complex, b, total)) = &mixed {
                        if *total {
                            let exponent = complex.filtration(generator, (*b).into()).unwrap();
                            out.push(&format!("u^{exponent}*"))?;
                        }
                        out.push(complex.name(generator))?;
                    } else {
                        if let Some(chart) = &two_chart {
                            let (i, j) = chart.exponent(generator);
                            out.push(&format!("omega({i},{j})"))?;
                        } else if let Some((chart, genus, de_rham)) = &chart {
                            let (label, index) = if *de_rham && generator >= *genus {
                                ("eta", generator - genus)
                            } else if *de_rham {
                                ("alpha", generator)
                            } else {
                                ("omega", generator)
                            };
                            let (b, r) = chart.exponent(index);
                            out.push(&format!("{label}({b},{r})"))?;
                        } else {
                            if let Some((data, b)) = &minimal {
                                let exponent = data.filtration(generator, *b).unwrap();
                                out.push(&format!("u^{exponent}*"))?;
                            }
                            out.push(&format!("g{generator}"))?;
                        }
                    }
                }
                out.push("}")
            }
            Value::Class(support) => {
                let bars = match *self.kind(ty) {
                    Kind::RMBar(grade) if !support.is_empty() => {
                        self.spend(support.len())?;
                        Some((real_extensions::bars(self, grade)?, grade.t))
                    }
                    _ => None,
                };
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if support.is_empty() {
                    out.push("0")?;
                }
                for (index, generator) in support.iter().enumerate() {
                    if index != 0 {
                        out.push(" + ")?;
                    }
                    if let Some((bars, t)) = &bars {
                        let exponent = bars.get(*generator)?.exponent_at(*t)?;
                        out.push(&format!("rho^{exponent}*b{generator}"))?;
                    } else {
                        out.push(&format!("g{generator}"))?;
                    }
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
                    if let Kind::MCobar(grade) = *self.kind(ty) {
                        self.spend(cobar::word_cells(word))?;
                        write_motivic_word(out, grade, word)?;
                    } else {
                        write_cobar_word(out, word)?;
                    }
                }
                out.push("}")
            }
            Value::RhoMin(terms) => {
                self.spend(terms.len())?;
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                for (index, (f, i)) in terms.iter().enumerate() {
                    if index > 0 {
                        out.push(" + ")?;
                    }
                    out.push(&format!("rho^{f}*g{i}"))?;
                }
                out.push("}")
            }
            Value::RhoLoc(terms) => {
                let Kind::RhoLoc(grade) = *self.kind(ty) else {
                    unreachable!();
                };
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                for (index, term) in terms.iter().rev().enumerate() {
                    if index != 0 {
                        out.push(" + ")?;
                    }
                    self.spend(term.cells())?;
                    out.push(&format!("rholocword[{}, {}]{{", grade.t, term.vertex()))?;
                    for (j, (u, xi)) in term.factors().enumerate() {
                        if j != 0 {
                            out.push(", ")?;
                        }
                        out.push(&format!("({u}, "))?;
                        write_number_list(out, xi)?;
                        out.push(")")?;
                    }
                    out.push("}")?;
                }
                out.push("}")
            }
            Value::RMCobar(terms) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                for (index, term) in terms.iter().rev().enumerate() {
                    if index != 0 {
                        out.push(" + ")?;
                    }
                    self.spend(term.cells())?;
                    out.push(&format!("rmcword[{}, {}]{{", term.rho, term.tau))?;
                    write_motivic_factors(out, &term.word)?;
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
            Value::RMTriple(coset) => {
                let Kind::RMTriple(grade) = *self.kind(ty) else {
                    unreachable!()
                };
                self.write_type(out, ty, 0)?;
                out.push("{base = ")?;
                let ext_type = self.rmext_type(grade)?;
                self.write_value(out, ext_type, &coset.base())?;
                out.push(&format!(", dimension = {}}}", coset.dimension()))
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
            Value::RVNCChain(terms) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                for (index, term) in terms.iter().enumerate() {
                    if index != 0 {
                        out.push(" + ")?;
                    }
                    out.push(&format!("u^{} * ", term.power))?;
                    write_hh_word(out, &term.word)?;
                }
                out.push("}")
            }
            Value::RVHHChain(words) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if words.is_empty() {
                    out.push("0")?;
                }
                for (index, word) in words.iter().enumerate() {
                    if index != 0 {
                        out.push(" + ")?;
                    }
                    write_hh_word(out, word)?;
                }
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
            Value::RVIwHecke(terms) => self.write_hecke(
                out,
                ty,
                "T",
                terms
                    .iter()
                    .map(|term| (term.window.as_ref(), term.coefficient.as_ref())),
            ),
            Value::RVSphHecke(terms) => self.write_hecke(
                out,
                ty,
                "S",
                terms
                    .iter()
                    .map(|term| (term.coweight.as_ref(), term.coefficient.as_ref())),
            ),
            Value::RVWzw(terms) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                for (index, (weight, coefficient)) in terms.iter().enumerate() {
                    self.spend(1 + weight.len() + coefficient.len())?;
                    if index != 0 {
                        out.push(", ")?;
                    }
                    out.push(&format!("{weight:?}: {coefficient:?}"))?;
                }
                out.push("}")
            }
            Value::KhChain(terms) => self.write_cube_terms(out, ty, terms),
            Value::LeeChain(terms) => self.write_cube_terms(out, ty, terms),
            Value::KhClass(terms) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                for (i, (generator, coefficient)) in terms.iter().enumerate() {
                    self.spend(1)?;
                    if i != 0 {
                        out.push(", ")?;
                    }
                    out.push(&format!("g{generator}: {coefficient}"))?;
                }
                out.push("}")
            }
            Value::RvTl(terms) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                for (index, (path, coefficient)) in terms.iter().enumerate() {
                    self.spend(1 + path.len() + coefficient.numerator.len())?;
                    if index != 0 {
                        out.push(", ")?;
                    }
                    out.push(&format!(
                        "{path:?}: {:?} / {}",
                        coefficient.numerator, coefficient.denominator
                    ))?;
                }
                out.push("}")
            }
            Value::RVResSpace(value) => {
                let grade = value.grade;
                self.spend(value.rows.len() * (1 + grade.width))?;
                out.push(&format!(
                    "RVResSpace({}, {}){{",
                    grade.prime.get(),
                    grade.width
                ))?;
                for (i, row) in value.rows.iter().enumerate() {
                    if i != 0 {
                        out.push(", ")?;
                    }
                    out.push(&format!("{row:?}"))?;
                }
                out.push("}")
            }
            Value::RVHyperArithmetic(value) => out.push(&value.summary(self)?),
            Value::RVHyperCover(value) => out.push(&value.summary(self)?),
            Value::RVStableTree(value) => out.push(&value.summary(self)?),
            Value::RVStableMarks(value) => out.push(&value.summary(self)?),
            Value::RVValMat(value) => {
                let grade = value.grade;
                out.push(&format!(
                    "RVValMat({}, {}, {}){{",
                    grade.prime, grade.ground, grade.rank
                ))?;
                for (i, (basis, weight)) in value.entries().iter().enumerate() {
                    self.spend(1 + basis.len())?;
                    if i != 0 {
                        out.push(", ")?;
                    }
                    out.push(&format!("{basis:?}: {weight}"))?;
                }
                out.push("}")
            }
            Value::RVQpQuad(value) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{{}, {}, {}}}", value.a, value.b, value.c))
            }
            Value::CMOrder(value) | Value::ECEndOrder(value) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{{}, {}}}", value.constant, value.generator))
            }
            Value::ECCurve(name) => {
                self.write_type(out, ty, 0)?;
                let curve = self.ec_curve(name)?;
                out.push(&format!("{{a={:?}, b={:?}}}", curve.a, curve.b))
            }
            Value::Jacobian(divisor) => {
                self.write_type(out, ty, 0)?;
                out.push("{u=")?;
                self.write_ec_polynomial(out, &divisor.u)?;
                out.push(", v=")?;
                self.write_ec_polynomial(out, &divisor.v)?;
                out.push("}")
            }
            Value::Kummer(point) => {
                self.write_type(out, ty, 0)?;
                self.write_ec_polynomial(out, point.coordinates())
            }
            Value::KummerK3(point) => {
                self.write_type(out, ty, 0)?;
                out.push("{base=")?;
                self.write_ec_polynomial(out, point.base().coordinates())?;
                if let Some(direction) = point.direction() {
                    out.push(", direction=")?;
                    self.write_ec_polynomial(out, direction)?;
                }
                out.push("}")
            }
            Value::K3PencilSection(section) => {
                self.write_type(out, ty, 0)?;
                match section.as_ref() {
                    pencil::sections::Section::Infinity => out.push("{O}"),
                    pencil::sections::Section::Affine(x, y) => {
                        out.push("{x=")?;
                        self.write_ec_polynomial(out, &x.numerator)?;
                        out.push("/")?;
                        self.write_ec_polynomial(out, &x.denominator)?;
                        out.push(", y=")?;
                        self.write_ec_polynomial(out, &y.numerator)?;
                        out.push("/")?;
                        self.write_ec_polynomial(out, &y.denominator)?;
                        out.push("}")
                    }
                }
            }
            Value::K3PencilCoverPoint(point) => {
                self.write_type(out, ty, 0)?;
                out.push("[")?;
                for (i, coordinate) in point.coordinates().iter().enumerate() {
                    if i != 0 {
                        out.push(", ")?;
                    }
                    self.write_ec_polynomial(out, coordinate)?;
                }
                out.push("]")
            }
            Value::K3PencilKummer2(descent) => {
                self.write_type(out, ty, 0)?;
                out.push("[")?;
                for (i, (nonsquare, support)) in descent.coordinates().enumerate() {
                    if i != 0 {
                        out.push(", ")?;
                    }
                    out.push(if nonsquare {
                        "(nonsquare, ["
                    } else {
                        "(square, ["
                    })?;
                    for (j, factor) in support.iter().enumerate() {
                        if j != 0 {
                            out.push(", ")?;
                        }
                        self.write_ec_polynomial(out, factor)?;
                    }
                    out.push("])")?;
                }
                out.push("]")
            }
            Value::K3PencilBase(parameter) => {
                self.write_type(out, ty, 0)?;
                self.write_ec_polynomial(out, parameter.coordinates())
            }
            Value::Enriques(orbit) => {
                self.write_type(out, ty, 0)?;
                out.push(if orbit.split {
                    "{split, base="
                } else {
                    "{nonsplit, base="
                })?;
                self.write_ec_polynomial(out, orbit.representative.base().coordinates())?;
                if let Some(direction) = orbit.representative.direction() {
                    out.push(", direction=")?;
                    self.write_ec_polynomial(out, direction)?;
                }
                out.push("}")
            }
            Value::KummerDisc(mask) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{{mask}}}"))
            }
            Value::KummerWeil(vector) => {
                self.write_type(out, ty, 0)?;
                let fraction = vector.fraction();
                self.spend(fraction.numerator.len())?;
                out.push(&format!(
                    "{{{:?} / {}}}",
                    fraction.numerator, fraction.denominator
                ))
            }
            Value::KummerDiv(divisor) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                for (i, coefficient) in divisor.coefficients().iter().enumerate() {
                    self.spend(1)?;
                    if i != 0 {
                        out.push(", ")?;
                    }
                    out.push(&coefficient.to_string())?;
                }
                out.push("}")
            }
            Value::HilbertSquare(class) => {
                self.write_type(out, ty, 0)?;
                out.push("{coordinates=[")?;
                self.spend(class.coordinates().len())?;
                for (index, coordinate) in class.coordinates().iter().enumerate() {
                    if index != 0 {
                        out.push(", ")?;
                    }
                    out.push(&coordinate.to_string())?;
                }
                out.push("]}")
            }
            Value::K3Lattice(vector) => {
                let (w, b) = vector.coordinates();
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{w={w:?}, b={b:?}}}"))
            }
            Value::K3Clifford(element) => self.write_sparse_rational(out, ty, element.as_ref()),
            Value::OG6Clifford(element) => self.write_sparse_rational(out, ty, element.as_ref()),
            Value::KSExterior(form) => self.write_sparse_rational(out, ty, form.as_ref()),
            Value::KSNomizu(form) => self.write_sparse_rational(out, ty, form.as_ref()),
            Value::KSEnd(end) => self.write_sparse_rational(out, ty, end.as_ref()),
            Value::QZ(scalar) => out.push(&format!("RVQZ({})", scalar.text())),
            Value::K3Spinor(element) => self.write_sparse_rational(out, ty, element.as_ref()),
            Value::AbMukai(vector) => {
                let (r, wedge, s) = vector.coordinates();
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{r={r}, wedge={wedge:?}, s={s}}}"))
            }
            Value::OG6(class) => {
                let (vector, k) = class.coordinates();
                let (r, wedge, s) = vector.coordinates();
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{r={r}, wedge={wedge:?}, s={s}, k={k}}}"))
            }
            Value::CMClass(class) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{:?}", class.0))
            }
            Value::OG10(class) => {
                let (vector, k) = class.coordinates();
                let (w, b) = vector.c.coordinates();
                self.write_type(out, ty, 0)?;
                out.push(&format!(
                    "{{r={}, w={w:?}, b={b:?}, s={}, k={k}}}",
                    vector.r, vector.s
                ))
            }
            Value::Mukai(vector) | Value::MukaiQuotient(vector) => {
                let (w, b) = vector.c.coordinates();
                self.write_type(out, ty, 0)?;
                out.push(&format!(
                    "{{r={}, w={w:?}, b={b:?}, s={}}}",
                    vector.r, vector.s
                ))
            }
            Value::Nakajima(state) => {
                self.write_type(out, ty, 0)?;
                self.write_nakajima_terms(out, state, 1)
            }
            Value::HodgeTensor(tensor) => {
                let (Kind::HodgeTensor(grade) | Kind::Verbitsky(grade)) = self.kind(ty) else {
                    unreachable!()
                };
                let form = self.hodge_wick(&grade.basis)?;
                self.write_type(out, ty, 0)?;
                out.push(&tensor.text(self, &form)?)
            }
            Value::HilbertCohomology(class) => {
                self.write_type(out, ty, 0)?;
                self.write_nakajima_terms(out, class.numerator(), class.denominator())
            }
            Value::PeriodPolynomial(polynomial) => {
                self.write_type(out, ty, 0)?;
                out.push("{[")?;
                for (i, coefficient) in polynomial.coordinates().iter().enumerate() {
                    if i != 0 {
                        out.push(", ")?;
                    }
                    out.push(&coefficient.to_string())?;
                }
                out.push("]}")
            }
            Value::HeckeEigen(element) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                Self::write_hecke_element(out, element)?;
                out.push("}")
            }
            Value::MazurTate(distribution) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                for (i, (residue, coefficient)) in distribution.terms().enumerate() {
                    if i != 0 {
                        out.push(", ")?;
                    }
                    out.push(&format!("{residue}: "))?;
                    Self::write_hecke_element(out, coefficient)?;
                }
                out.push("}")
            }
            Value::MazurAmice(transform) => {
                self.write_type(out, ty, 0)?;
                out.push("{[")?;
                for (i, coefficient) in transform.coefficients().iter().enumerate() {
                    if i != 0 {
                        out.push(", ")?;
                    }
                    out.push(&coefficient.to_string())?;
                }
                out.push("]}")
            }
            Value::MazurOrdinary(residue) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{{residue}}}"))
            }
            Value::HeckeSymbol(symbol) => {
                self.write_type(out, ty, 0)?;
                out.push("{[")?;
                for (i, coefficient) in symbol.elements().iter().enumerate() {
                    if i != 0 {
                        out.push(", ")?;
                    }
                    Self::write_hecke_element(out, coefficient)?;
                }
                out.push("]}")
            }
            Value::ModularHecke(operator) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                for (i, (index, coefficient)) in operator.terms().enumerate() {
                    if i != 0 {
                        out.push(", ")?;
                    }
                    out.push(&format!("T_{index}: {}", hecke::decimal(coefficient)))?;
                }
                out.push("}")
            }
            Value::SiegelWeight(weight) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{{weight}}}"))
            }
            Value::ECMap(map) => {
                self.write_type(out, ty, 0)?;
                let Some(coords) = &map.coordinates else {
                    return out.push("{infinity map}");
                };
                out.push(&format!("{{degree={}, xnum=", map.degree()))?;
                self.write_ec_polynomial(out, &coords.x.numerator)?;
                out.push(", xden=")?;
                self.write_ec_polynomial(out, &coords.x.denominator)?;
                out.push(", ynum=")?;
                self.write_ec_polynomial(out, &coords.y_multiplier.numerator)?;
                out.push(", yden=")?;
                self.write_ec_polynomial(out, &coords.y_multiplier.denominator)?;
                if let Some(kernel) = map.kernel() {
                    out.push(", kernel=")?;
                    self.write_ec_polynomial(out, &kernel.polynomial)?;
                }
                out.push("}")
            }
            Value::ECLevelVec(vector) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{{vector:?}}}"))
            }
            Value::ECLevelMat(matrix) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{columns={:?}}}", matrix.0))
            }
            Value::ECPoint(point) => {
                self.write_type(out, ty, 0)?;
                match point.as_ref() {
                    ECPoint::Infinity => out.push("{infinity}"),
                    ECPoint::Affine(x, y) => out.push(&format!("{{x={x:?}, y={y:?}}}")),
                }
            }
            Value::ECKernel(kernel) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                self.write_ec_polynomial(out, &kernel.polynomial)?;
                out.push("}")
            }
            Value::ECIso(name) => {
                self.write_type(out, ty, 0)?;
                let map = self.ec_isogeny(name)?;
                out.push(&format!(
                    "{{source={}, target={}, degree={}, xnum=",
                    map.source, map.target, map.kernel.ell
                ))?;
                self.write_ec_polynomial(out, &map.numerator)?;
                out.push(", xden=")?;
                self.write_ec_polynomial(out, &map.denominator)?;
                out.push(", dnum=")?;
                self.write_ec_polynomial(out, &map.derivative_numerator)?;
                out.push(", dden=")?;
                self.write_ec_polynomial(out, &map.derivative_denominator)?;
                out.push("}")
            }
            Value::HWMatrix(matrix) => {
                self.write_type(out, ty, 0)?;
                self.spend(matrix.entries().len())?;
                out.push("{[")?;
                let n = matrix.size();
                for i in 0..n {
                    if i != 0 {
                        out.push(", ")?;
                    }
                    out.push(&format!("{:?}", &matrix.entries()[i * n..(i + 1) * n]))?;
                }
                out.push("]}")
            }
            Value::HWWMatrix(matrix) => {
                self.write_type(out, ty, 0)?;
                out.push("{[")?;
                for (i, row) in matrix.entries().chunks(matrix.size().max(1)).enumerate() {
                    if i != 0 {
                        out.push(", ")?;
                    }
                    out.push("[")?;
                    for (j, coefficient) in row.iter().enumerate() {
                        self.spend(coefficient.len())?;
                        if j != 0 {
                            out.push(", ")?;
                        }
                        out.push(&format!("{coefficient:?}"))?;
                    }
                    out.push("]")?;
                }
                out.push("]}")
            }
            Value::HWScalar(value) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{{value}}}"))
            }
            Value::RVUnram(values) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{values:?}"))
            }
            Value::RVCycloQ(ratio) => {
                self.spend(ratio.numerator.len())?;
                self.write_type(out, ty, 0)?;
                out.push("{[")?;
                for (index, coefficient) in ratio.numerator.iter().enumerate() {
                    if index != 0 {
                        out.push(", ")?;
                    }
                    out.push(&coefficient.to_string())?;
                }
                out.push(&format!("] / {}}}", ratio.denominator))
            }
            Value::RVCyclo(coefficients) => {
                self.spend(coefficients.len())?;
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                for (index, coefficient) in coefficients.iter().enumerate() {
                    if index != 0 {
                        out.push(", ")?;
                    }
                    out.push(&coefficient.to_string())?;
                }
                out.push("}")
            }
            Value::RVFusion(terms) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                for (index, (weight, coefficient)) in terms.iter().enumerate() {
                    self.spend(1 + weight.len())?;
                    if index != 0 {
                        out.push(", ")?;
                    }
                    out.push(&format!("{weight:?}: {coefficient}"))?;
                }
                out.push("}")
            }
            Value::PeriodResidue(value) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{{value}}}"))
            }
            Value::PeriodState(terms) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                for (i, (term, coefficient)) in terms.iter().enumerate() {
                    if i != 0 {
                        out.push(" + ")?;
                    }
                    out.push(&format!(
                        "{coefficient}*z^{}{:?}",
                        term.power, term.exponents
                    ))?;
                }
                out.push("}")
            }
            Value::RVWittSpectrum(value) => out.push(&value.render()),
            Value::RVSchurCharacter(value) => out.push(&value.render()),
            Value::RVStableSchurCharacter(value) => out.push(&value.render()),
            Value::RVPartitionAlgebra(value) => out.push(&value.render()),
            Value::RVPartitionLocal(value) => out.push(&value.render()),
            Value::RVKaroubiMap(value) => out.push(&value.render()),
            Value::RVOre(terms) => {
                out.push("RVOre{")?;
                self.spend(terms.len())?;
                for (index, (powers, coefficient)) in terms.iter().enumerate() {
                    if index != 0 {
                        out.push(", ")?;
                    }
                    out.push(&format!("({}, {}): {coefficient}", powers.q, powers.theta))?;
                }
                out.push("}")
            }
            Value::RVQGr(terms) => self.write_schubert(out, ty, terms),
            Value::RVQGrP(terms) => self.write_schubert(out, ty, terms),
            Value::RVHive(interior) => {
                self.spend(interior.len())?;
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                for (index, vertex) in interior.iter().enumerate() {
                    if index != 0 {
                        out.push(", ")?;
                    }
                    out.push(&vertex.to_string())?;
                }
                out.push("}")
            }
            Value::RVSatake(terms) => {
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                if terms.is_empty() {
                    out.push("0")?;
                }
                for (index, term) in terms.iter().enumerate() {
                    self.spend(1 + term.exponents.len())?;
                    if index != 0 {
                        out.push(if term.coefficient < 0 { " - " } else { " + " })?;
                    } else if term.coefficient < 0 {
                        out.push("-")?;
                    }
                    let magnitude = term.coefficient.unsigned_abs();
                    let mut has_factor = magnitude != 1;
                    if has_factor {
                        out.push(&magnitude.to_string())?;
                    }
                    if term.power != 0 {
                        if has_factor {
                            out.push("*")?;
                        }
                        if term.power == 1 {
                            out.push("v")?;
                        } else {
                            out.push(&format!("v^{}", term.power))?;
                        }
                        has_factor = true;
                    }
                    for (index, &power) in term.exponents.iter().enumerate() {
                        if power != 0 {
                            if has_factor {
                                out.push("*")?;
                            }
                            out.push(&format!("X{}", index + 1))?;
                            if power != 1 {
                                out.push(&format!("^{power}"))?;
                            }
                            has_factor = true;
                        }
                    }
                    if !has_factor {
                        out.push("1")?;
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
            Value::ASDeRham(support) => {
                let Kind::ASDeRham(name, _) = self.kind(ty) else {
                    unreachable!()
                };
                let genus = self.rv_as_covers[name.as_ref()].genus();
                self.spend(support.len())?;
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
                    if coordinate < genus {
                        out.push(&format!("alpha{coordinate}"))?;
                    } else {
                        out.push(&format!("eta{}", coordinate - genus + 1))?;
                    }
                }
                out.push("}")
            }
            Value::EOType(sequence) => {
                self.spend(sequence.len())?;
                self.write_type(out, ty, 0)?;
                out.push("{")?;
                for (position, entry) in sequence.iter().enumerate() {
                    if position != 0 {
                        out.push(", ")?;
                    }
                    out.push(&entry.to_string())?;
                }
                out.push("}")
            }
            Value::BinaryScalar(coefficient) => {
                self.write_type(out, ty, 0)?;
                out.push(&format!("{{{coefficient}}}"))
            }
            Value::RVHom(support) | Value::JetClass(support) | Value::CurveDiff(support) => {
                let symbol = if matches!(value, Value::CurveDiff(_)) {
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
            Value::SetNat(rows) => {
                let Kind::SetNat(source, _) = self.kind(ty) else {
                    unreachable!()
                };
                let diagram = self.set_family_named(source).expect("registered source");
                out.push("{")?;
                for (point, name) in diagram.site.names.iter().enumerate() {
                    if point != 0 {
                        out.push(", ")?;
                    }
                    out.push(name)?;
                    out.push(": [")?;
                    let entries = diagram.components(point, rows);
                    self.spend(entries.len() + 1)?;
                    for (index, image) in entries.iter().enumerate() {
                        if index != 0 {
                            out.push(", ")?;
                        }
                        out.push(&image.to_string())?;
                    }
                    out.push("]")?;
                }
                out.push("}")
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

fn write_hh_word(out: &mut Text, word: &rv_hh::Word) -> Result<()> {
    write_hh_factor(out, &word[0])?;
    out.push("[")?;
    for (index, factor) in word[1..].iter().enumerate() {
        if index != 0 {
            out.push(" | ")?;
        }
        write_hh_factor(out, factor)?;
    }
    out.push("]")
}

fn write_hh_factor(out: &mut Text, factor: &rv_hh::Factor) -> Result<()> {
    out.push(&format!(
        "({}, {:?}):g{}",
        factor.s, factor.weight, factor.generator
    ))
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

fn write_motivic_word(out: &mut Text, grade: motivic::Grade, word: &cobar::Word) -> Result<()> {
    out.push(&format!("mcword[{}]{{", grade.tau_power(word)))?;
    write_motivic_factors(out, word)
}

fn write_motivic_factors(out: &mut Text, word: &cobar::Word) -> Result<()> {
    for (index, factor) in word.iter().enumerate() {
        if index > 0 {
            out.push(", ")?;
        }
        let mut xi: Vec<_> = factor.iter().map(|power| power / 2).collect();
        let mut tau: Vec<_> = factor.iter().map(|power| power % 2).collect();
        xi.truncate(
            xi.iter()
                .rposition(|&power| power != 0)
                .map_or(0, |i| i + 1),
        );
        tau.truncate(
            tau.iter()
                .rposition(|&power| power != 0)
                .map_or(0, |i| i + 1),
        );
        out.push("(")?;
        write_number_list(out, &xi)?;
        out.push(", ")?;
        write_number_list(out, &tau)?;
        out.push(")")?;
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
