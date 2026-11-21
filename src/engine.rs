use std::collections::{HashMap, HashSet};
use std::rc::Rc;

use crate::adams;
use crate::binary_field::{BinaryField, Embedding};
use crate::bp_height;
use crate::bp_page;
use crate::bp_triple;
use crate::bracket;
use crate::ce;
use crate::cobar;
use crate::cobar_squares;
use crate::cohomology;
use crate::comodule::{self, Generator};
use crate::cyclic::{self, Layout};
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
use crate::model::{Kind, Model, Poset, SuspensionPoint, TypeId, Value};
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
use crate::restriction_transfer;
use crate::rv_arithmetic;
use crate::rv_as_cover::Cover as RvAsCover;
use crate::rv_bar;
use crate::rv_bar_page;
use crate::rv_curve::{self, Curve as RvCurve};
use crate::rv_descent;
use crate::rv_field;
use crate::rv_gauss;
use crate::rv_hochster;
use crate::rv_hodge;
use crate::rv_hom;
use crate::rv_ideal::Ideal as RvIdeal;
use crate::rv_iwahori;
use crate::rv_koszul;
use crate::rv_laurent;
use crate::rv_newton;
use crate::rv_poly;
use crate::rv_rees;
use crate::rv_scheme;
use crate::rv_shintani;
use crate::rv_t_bracket;
use crate::rv_t_system;
use crate::rv_t_transfer;
use crate::rv_tree;
use crate::steenrod;
use crate::syntax::BarFamilySpec;
use crate::syntax::{
    BinaryOp, BpExtOp, BpRelOp, ClassForm, CyclicForm, Entries, ExprForm, JetHomSpec, JetOp,
    JetPageSpec, JetPathArg, Label, MapMaySpec, ParsedExpr, ParsedType, QSpec, RvIdealOp,
    RvSchemeOp, Statement, StatementForm, SuspensionForm, TypeForm, VBockSpec, Witness,
};
use crate::syntax::{RvBarPageSpec, RvBarSpec, RvDegree, RvGradeSpec};
use crate::transfer;
use crate::v_bars;
use crate::v_koszul;
use crate::vbock;
use crate::{bp, bp_classes, bp_ext, bp_polynomial::Polynomial as BpPolynomial, bp_rel};

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

enum TypedForm {
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
    BpExtOperation(BpExtOp, Box<TypedExpr>),
    BpMassey(bp::Massey, Vec<TypedExpr>),
    BpTriple(bp_triple::Layout, Box<[TypedExpr; 3]>),
    BpTripleBase(Box<TypedExpr>),
    BpRaise(Box<TypedExpr>),
    BpPageClass(Box<TypedExpr>),
    BpPageOperation(bp_page::Operation, Box<TypedExpr>),
    RvBarWord(Vec<usize>, Vec<usize>, Vec<TypedExpr>),
    RvBarLength(usize, Box<TypedExpr>),
    RvMComponent(Box<TypedExpr>),
    RvKoszulTerm(usize, Rc<RvIdeal>, Box<TypedExpr>),
    RvHasse(Vec<usize>, Box<TypedExpr>),
    RvVanish(BinaryField, Vec<usize>, Box<TypedExpr>),
    RvInMinors(TypeId, q_koszul::Indices, usize, Box<TypedExpr>),
    RvIdeal(RvIdealOp, Rc<RvIdeal>, Box<TypedExpr>),
    CobarCup(Profile, usize, Box<TypedExpr>, Box<TypedExpr>),
    ExtSquare(usize, Box<TypedExpr>),
    Nakamura(usize, Box<TypedExpr>),
    MayLiftTo(usize, usize, Box<TypedExpr>),
    Coact(usize, usize, Box<TypedExpr>),
    MapExtInduced(Box<TypedExpr>),
    MapWord(Generator, Generator, Box<TypedExpr>),
    Equator,
    Suspend(Box<TypedExpr>),
    Shift(Box<TypedExpr>),
    PScale(Prime, usize, Box<TypedExpr>),
    JetOperation(JetOp, Box<TypedExpr>),
    Bracket(Vec<TypedExpr>),
    Transfer(ScalarField, transfer::Operation, Vec<TypedExpr>),
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

    fn run_statement(&mut self, statement: StatementForm) -> Result<Option<String>> {
        let line = match statement {
            StatementForm::QModule(name, selected, presentation) => {
                self.model.define_q_module(&name, selected, presentation)?;
                return Ok(None);
            }
            StatementForm::Space(name, chains) => {
                self.model.define_space(&name, &chains)?;
                return Ok(None);
            }
            StatementForm::Let(name, expr) => {
                if self.bindings.contains_key(&name.text) {
                    return Err(Error::new(
                        name.span,
                        format!("value '{}' is already defined", name.text),
                    ));
                }
                let expr = self.checker().expr(expr)?;
                let value = self.eval(&expr)?;
                self.model.allocate(1)?;
                self.bindings.insert(
                    name.text,
                    Binding {
                        index: self.values.len(),
                        ty: expr.ty,
                    },
                );
                self.values.push(value);
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
            StatementForm::RvAsCover(name, variables, expressions) => {
                self.model.check_new_rv_as_cover(&name)?;
                let selected = q_koszul::indices(&mut self.model, variables)?;
                let forms = self.rv_generator_values(expressions, &selected, "rvascover")?;
                let cover = RvAsCover::build(&mut self.model, selected.len(), forms)?;
                self.model.define_rv_as_cover(name.text, cover);
                return Ok(None);
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
                        if a.is_empty() || b.is_empty() {
                            return Ok(odd_ext::zero());
                        }
                        let source = self.model.bp_page(left.ty)?;
                        let other = self.model.bp_page(right.ty)?;
                        let target = self.model.bp_page(expr.ty)?;
                        let product = source.cup(&mut self.model, &other, &target, &a, &b)?;
                        odd_ext::persist(&mut self.model, product)
                    }
                    (BinaryOp::BpMultiply, Value::BP(a), Value::BP(b)) => {
                        let result = a.multiply(&mut self.model, &b)?;
                        self.persist_bp(expr.ty, result)
                    }
                    (BinaryOp::Add, Value::RVPoly(a), Value::RVPoly(b)) => {
                        let polynomial = rv_poly::add(&mut self.model, &a, &b)?;
                        Ok(Value::RVPoly(rv_poly::persist(
                            &mut self.model,
                            polynomial,
                        )?))
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
                    (BinaryOp::Add, Value::Cobar(a), Value::Cobar(b)) => {
                        cobar::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::OddClass(a), Value::OddClass(b)) => {
                        let p = match *self.model.kind(expr.ty) {
                            Kind::PExt(p, ..)
                            | Kind::PExtAn(p, ..)
                            | Kind::PMay(p, ..)
                            | Kind::PMayAn(p, ..)
                            | Kind::PCeAn(p, ..) => p,
                            Kind::BPPage(grade) => Prime::validated(grade.shape.prime),
                            _ => unreachable!(),
                        };
                        odd_ext::add(&mut self.model, p, &a, &b)
                    }
                    (BinaryOp::Add, Value::PCobar(a), Value::PCobar(b)) => {
                        let (p, _, _, _) = self.model.kind(expr.ty).pcobar_grade().unwrap();
                        odd_cobar::add(&mut self.model, p, &a, &b)
                    }
                    (BinaryOp::Add, Value::RVLaurent(a), Value::RVLaurent(b)) => {
                        rv_laurent::add(&mut self.model, &a, &b)
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
                        adams::yoneda(&mut self.model, left.ty, &a, right.ty, &b, expr.ty)
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
            TypedForm::Cobar(op, input) => match (op, self.eval(input)?) {
                (CobarOp::Differential, Value::Cobar(words)) => {
                    let (profile, _, _) = self.model.kind(input.ty).cobar_grade().unwrap();
                    cobar::differential(&mut self.model, profile, &words)
                }
                (CobarOp::Differential, Value::PCobar(terms)) => {
                    let (p, profile, _, _) = self.model.kind(input.ty).pcobar_grade().unwrap();
                    odd_cobar::differential(&mut self.model, p, profile, &terms)
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
                    if class.is_empty() {
                        return Ok(cobar::zero());
                    }
                    let adams = self.model.adams(input.ty)?;
                    let words = adams.representative(&mut self.model, &class)?;
                    self.model.allocate(words.len())?;
                    Ok(Value::Cobar(words.into()))
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
                    let (profile, s, t) = self.model.kind(input.ty).cobar_grade().unwrap();
                    let ty = self.model.profile_ext_type(profile, s, t)?;
                    let adams = self.model.adams(ty)?;
                    let primitive = adams
                        .primitive(&mut self.model, &words)?
                        .ok_or_else(|| Error::plain("extPrimitive requires a cobar boundary"))?;
                    self.model.allocate(primitive.len())?;
                    Ok(Value::Cobar(primitive.into()))
                }
                _ => unreachable!(),
            },
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

impl Checker<'_> {
    fn ty(&mut self, ty: ParsedType) -> Result<TypeId> {
        self.model.spend(1)?;
        match ty.kind {
            TypeForm::Name(name) => self
                .model
                .named(&name.text)
                .ok_or_else(|| Error::new(name.span, format!("unknown space '{}'", name.text))),
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
            TypeForm::BpPage(mode, spec) => {
                let grade = self.bp_page_grade(mode, spec)?;
                self.model.bp_page_type(grade)
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
            TypeForm::BP(spec, depth) => {
                let mut shape = bp::Shape::new(self.model, spec.0, spec.1, spec.2)?;
                if let Some(depth) = depth {
                    shape = shape.quotient(depth)?;
                }
                self.model.bp_type(shape)
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

    fn expr_inner(&mut self, expr: ExprForm, span: Span) -> Result<(TypeId, TypedForm)> {
        self.model.spend(1)?;
        match expr {
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
                    Kind::Cobar(..) | Kind::CobarAn(..) => cobar::zero(),
                    Kind::PCobar(..) | Kind::PCobarAn(..) => odd_cobar::zero(),
                    Kind::PExt(..) | Kind::PExtAn(..) | Kind::PMay(..) | Kind::PMayAn(..) | Kind::PCeAn(..) | Kind::BPPage(..) => odd_ext::zero(),
                    Kind::MapCobar(..) => Value::MapCobar([].into()),
                    Kind::JetHom(..) => jet_hom::zero(),
                    Kind::JetPage(..) => jet_page::zero(),
                    Kind::QKoszul(..) | Kind::VKoszul(..) => Value::QKoszul([].into()),
                    Kind::RVBarE(..) => cohomology::zero(),
                    Kind::RVBar(..) => Value::RVBar([].into()),
                    Kind::RVKoszul(..) => Value::RVKoszul([].into()),
                    Kind::RVTor(..) => cohomology::zero(),
                    Kind::RVPoly(..) => Value::RVPoly([].into()),
                    Kind::BP(..) => bp::zero(),
                    Kind::BPRel(..) => bp_rel::zero(),
                    Kind::BPExt(..) | Kind::BPRelExt(..) => bp_classes::zero(),
                    Kind::QExt(..) | Kind::VExt(..) | Kind::VBock(..) => cohomology::zero(),
                    Kind::RVHom(..) => rv_hom::zero(),
                    Kind::RVLaurent(..) => rv_laurent::zero(),
                    Kind::Cohomology(..) | Kind::Stable(..) | Kind::Ext(..) | Kind::ExtAn(..) | Kind::XExt(..) | Kind::MapExt(..) | Kind::May(..) | Kind::MayAn(..) | Kind::MapMay(..) | Kind::CeAn(..) => {
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
            ExprForm::PCobarWord(p, factors) => {
                let prime = Prime::new(self.model, p)?;
                let ((s, t), value) = odd_cobar::literal(self.model, prime, factors)?;
                let ty = self.model.pcobar_type(prime, s, t)?;
                Ok((ty, TypedForm::Constant(value)))
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
                    | Kind::PCeAn(p, ..) => p,
                    Kind::BPPage(grade) => Prime::validated(grade.shape.prime),
                    _ => return Err(self.invalid_operand(
                        &value,
                        "pscale",
                        "a PCobar cochain, PExt class, PMay class, PCEAn class, or BP page class",
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
                            Kind::BPExt(..) | Kind::BPRel(..) | Kind::BPRelExt(..)
                        ) {
                            return Err(self.invalid_operand(
                                &input,
                                "bpScale",
                                "a BPExt class, BPRel cochain, or BPRelExt class",
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
                if matches!(operation, bp::Operation::Negate)
                    && matches!(kind, Kind::BPExt(..) | Kind::BPRel(..) | Kind::BPRelExt(..))
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
        operation.require_arity(arity)?;
        let (field, terms, degrees) = self.scalar_inputs(operation.name(), inputs)?;
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
        let (field, terms, degrees) = self.scalar_inputs("resHigher", inputs)?;
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

    fn scalar_inputs(
        &mut self,
        operation: &str,
        inputs: Vec<ParsedExpr>,
    ) -> Result<(ScalarField, Vec<TypedExpr>, Vec<transfer::Bidegree>)> {
        let arity = inputs.len();
        let mut terms = Vec::with_capacity(arity);
        let mut degrees = Vec::with_capacity(arity);
        let mut source = None;
        for input in inputs {
            let term = self.expr(input)?;
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
                    .ok_or_else(|| self.invalid_operand(&term, operation, expected))?;
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
                        return Err(self.invalid_operand(&term, operation, "a PExt class"));
                    }
                    (ScalarField::ModTwo(_), ScalarField::Odd(..)) => {
                        return Err(self.invalid_operand(&term, operation, "an Ext class"));
                    }
                }
            }
            source = Some(field);
            degrees.push((s, t));
            terms.push(term);
        }
        let field = source
            .ok_or_else(|| Error::plain(format!("{operation} requires at least one Ext class")))?;
        Ok((field, terms, degrees))
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

    fn bp_relative_shape(&mut self, spec: crate::syntax::BpRelSpec) -> Result<bp_rel::Shape> {
        bp_rel::Shape::new(self.model, spec.0, spec.1, spec.2, spec.3)
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

    fn cobar_result(
        &mut self,
        op: &CobarOp,
        argument: &TypedExpr,
        name: &str,
        span: Span,
    ) -> Result<TypeId> {
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
                        | Kind::PCobar(..)
                        | Kind::PCobarAn(..)
                        | Kind::PExt(..)
                        | Kind::PExtAn(..)
                        | Kind::PMay(..)
                        | Kind::PMayAn(..)
                        | Kind::PCeAn(..)
                        | Kind::BPPage(..)
                        | Kind::CobarAn(..)
                        | Kind::MapCobar(..)
                        | Kind::QKoszul(..)
                        | Kind::JetHom(..)
                        | Kind::RVBarE(..)
                        | Kind::RVBar(..)
                        | Kind::RVKoszul(..)
                        | Kind::RVTor(..)
                        | Kind::JetPage(..)
                        | Kind::RVPoly(..)
                        | Kind::BP(..)
                        | Kind::BPRel(..)
                        | Kind::BPExt(..)
                        | Kind::BPRelExt(..)
                        | Kind::QExt(..)
                        | Kind::RVHom(..)
                        | Kind::RVLaurent(..)
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
            BinaryOp::BpCup => {
                if matches!(self.model.kind(left.ty), Kind::BPExt(..))
                    || matches!(self.model.kind(right.ty), Kind::BPExt(..))
                {
                    let left = self.bp_ext_operand(left, "bpCup")?;
                    let right = self.bp_ext_operand(right, "bpCup")?;
                    return self.model.bp_ext_type(left.cup(right)?);
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
