use std::collections::{HashMap, HashSet};
use std::rc::Rc;

use crate::adams;
use crate::bracket;
use crate::ce;
use crate::cobar;
use crate::cobar_squares;
use crate::cohomology;
use crate::comodule::{self, Generator};
use crate::cyclic::{self, Layout};
use crate::error::{Error, Result, Span};
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
use crate::odd_transfer;
use crate::page_transition::{self, Direction as PageDirection, PageFamily};
use crate::profile::Profile;
use crate::q_koszul::{self, Grade};
use crate::q_module;
use crate::restriction_transfer;
use crate::steenrod;
use crate::syntax::{
    BinaryOp, ClassForm, CyclicForm, Entries, ExprForm, MapMaySpec, ParsedExpr, ParsedType, QSpec,
    Statement, StatementForm, SuspensionForm, TypeForm, VBockSpec, Witness,
};
use crate::transfer;
use crate::v_koszul;
use crate::vbock;

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
    Odd(Prime),
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
    VShift(i64, Box<TypedExpr>),
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
    Bracket(Vec<TypedExpr>),
    Transfer(ScalarField, transfer::Operation, Vec<TypedExpr>),
    HigherRestriction(Profile, Profile, Vec<TypedExpr>),
    MapTransfer(map_transfer::Operation, Vec<TypedExpr>),
    MapSystem(map_system::Operation, usize, Vec<DefiningEntry>),
    OddSystem(odd_massey::Operation, usize, Vec<DefiningEntry>),
    MapBracket {
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
            StatementForm::PExtDimension(p, s, t) => {
                let prime = Prime::new(&mut self.model, p)?;
                let ty = self.model.odd_ext_type(prime, s, t)?;
                let rank = self.model.odd_ext(ty)?.rank();
                format!("pextdim[{p}, {s}, {t}] = {rank}")
            }
            StatementForm::PExtAnDimension(p, n, s, t) => {
                let prime = Prime::new(&mut self.model, p)?;
                let ty = self
                    .model
                    .profile_odd_ext_type(prime, Profile::An(n), s, t)?;
                let rank = self.model.odd_ext(ty)?.rank();
                format!("pextandim[{p}, {n}, {s}, {t}] = {rank}")
            }
            StatementForm::PMayDimension(p, r, s, t, f) => {
                let prime = Prime::new(&mut self.model, p)?;
                let ty = self.model.pmay_type(prime, r, s, t, f)?;
                let rank = self.model.pmay(ty)?.rank();
                format!("pmaydim[{p}, {r}, {s}, {t}, {f}] = {rank}")
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
                    (BinaryOp::Add, Value::Cobar(a), Value::Cobar(b)) => {
                        cobar::add(&mut self.model, &a, &b)
                    }
                    (BinaryOp::Add, Value::OddClass(a), Value::OddClass(b)) => {
                        let (Kind::PExt(p, ..) | Kind::PExtAn(p, ..) | Kind::PMay(p, ..)) =
                            *self.model.kind(expr.ty)
                        else {
                            unreachable!()
                        };
                        odd_ext::add(&mut self.model, p, &a, &b)
                    }
                    (BinaryOp::Add, Value::PCobar(a), Value::PCobar(b)) => {
                        let (p, _, _, _) = self.model.kind(expr.ty).pcobar_grade().unwrap();
                        odd_cobar::add(&mut self.model, p, &a, &b)
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
                    (BinaryOp::MayMul, Value::OddClass(a), Value::OddClass(b)) => {
                        odd_may::multiply(&mut self.model, left.ty, &a, right.ty, &b, expr.ty)
                    }
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
                (CobarOp::MayClass, Value::PCobar(terms)) => {
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
                (CobarOp::MayDifferential, Value::OddClass(class)) => {
                    odd_may::differential(&mut self.model, input.ty, &class, expr.ty)
                }
                (
                    CobarOp::PageTransition(PageFamily::May, PageDirection::Next),
                    Value::OddClass(class),
                ) => odd_may::next(&mut self.model, input.ty, &class, expr.ty),
                (
                    CobarOp::PageTransition(PageFamily::May, PageDirection::Previous),
                    Value::OddClass(class),
                ) => odd_may::previous(&mut self.model, input.ty, &class, expr.ty),
                (CobarOp::MayRepresentative, Value::OddClass(class)) => {
                    if class.is_empty() {
                        return Ok(odd_cobar::zero());
                    }
                    let page = self.model.pmay(input.ty)?;
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
            TypedForm::QRestrict(input) => {
                let value = self.eval(input)?;
                q_koszul::restrict(&mut self.model, input.ty, value, expr.ty)
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
                        odd_ext::restrict(&mut self.model, input.ty, &class, expr.ty)
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
            TypedForm::Transfer(ScalarField::Odd(prime), operation, inputs) => {
                let inputs = inputs
                    .iter()
                    .map(|term| {
                        let Kind::PExt(_, s, t) = *self.model.kind(term.ty) else {
                            unreachable!()
                        };
                        let Value::OddClass(coordinates) = self.eval(term)? else {
                            unreachable!()
                        };
                        Ok(((s, t), coordinates))
                    })
                    .collect::<Result<Vec<odd_transfer::Input>>>()?;
                odd_transfer::evaluate(&mut self.model, *prime, &inputs, *operation)
            }
            TypedForm::HigherRestriction(source, target, inputs) => {
                let inputs = self.scalar_arguments(inputs)?;
                restriction_transfer::evaluate(&mut self.model, *source, *target, &inputs)
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
            TypedForm::MapBracket {
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
                if let Some(query) = query {
                    let (_, coordinates) = eval_class(query)?;
                    map_bracket::member(&mut self.model, *target, &inputs, &coordinates)
                } else {
                    map_bracket::predicate(&mut self.model, *target, &inputs)
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
                    Kind::PExt(..) | Kind::PExtAn(..) | Kind::PMay(..) => odd_ext::zero(),
                    Kind::MapCobar(..) => Value::MapCobar([].into()),
                    Kind::QKoszul(..) | Kind::VKoszul(..) => Value::QKoszul([].into()),
                    Kind::QExt(..) | Kind::VExt(..) | Kind::VBock(..) => cohomology::zero(),
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
                    | Kind::PMay(p, ..) => p,
                    _ => {
                        return Err(self.invalid_operand(
                            &value,
                            "pscale",
                            "a PCobar cochain, PExt class, or PMay class",
                        ))
                    }
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
                let rank = self.model.pmay(ty)?.rank();
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
                let Kind::CobarAn(n, s, t) = *self.model.kind(value.ty) else {
                    return Err(self.invalid_operand(&value, "ceClass", "a CobarAn cochain"));
                };
                let ty = self.model.ce_type(n, r, s, t, f)?;
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
            ExprForm::MapBracket(operation, inputs) => self.map_bracket(operation, inputs),
            ExprForm::MapTransfer(operation, inputs) => self.map_transfer(operation, inputs),
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
                        if let Some(operation) = map_system::Operation::from_name(&name.text) {
                            return self.map_system(operation, vec![*argument], Vec::new(), span);
                        }
                        if let Some(operation) = odd_massey::Operation::from_name(&name.text) {
                            return self.odd_system(operation, vec![*argument], Vec::new(), span);
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
        let range = self.model.named("S").unwrap();
        let (ty, query) = if let Some(query) = query {
            let query = self.expr(query)?;
            self.require_type(target, &query, operation.name())?;
            (range, Some(Box::new(query)))
        } else {
            (self.model.arrow(target, range)?, None)
        };
        Ok((
            ty,
            TypedForm::MapBracket {
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
        let (field, terms, degrees) = self.scalar_inputs(operation.name(), inputs, true)?;
        let (s, t) = operation.bidegree(degrees, arity)?;
        let ty = match (field, operation) {
            (ScalarField::ModTwo(profile), transfer::Operation::Product) => {
                self.model.profile_ext_type(profile, s, t)?
            }
            (ScalarField::ModTwo(profile), transfer::Operation::Lift) => {
                self.model.profile_cobar_type(profile, s, t)?
            }
            (ScalarField::Odd(prime), transfer::Operation::Product) => {
                self.model.odd_ext_type(prime, s, t)?
            }
            (ScalarField::Odd(prime), transfer::Operation::Lift) => {
                self.model.pcobar_type(prime, s, t)?
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
        let (ScalarField::ModTwo(source), terms, degrees) =
            self.scalar_inputs("resHigher", inputs, false)?
        else {
            unreachable!()
        };
        let target = Profile::An(n);
        if !source.can_restrict_to(target) {
            return Err(Error::new(
                span,
                "resHigher cannot restrict from a smaller Steenrod algebra",
            ));
        }
        let (s, t) = restriction_transfer::bidegree(degrees, arity)?;
        let ty = self.model.profile_ext_type(target, s, t)?;
        Ok((ty, TypedForm::HigherRestriction(source, target, terms)))
    }

    fn scalar_inputs(
        &mut self,
        operation: &str,
        inputs: Vec<ParsedExpr>,
        allow_odd: bool,
    ) -> Result<(ScalarField, Vec<TypedExpr>, Vec<transfer::Bidegree>)> {
        let arity = inputs.len();
        let mut terms = Vec::with_capacity(arity);
        let mut degrees = Vec::with_capacity(arity);
        let mut source = None;
        for input in inputs {
            let term = self.expr(input)?;
            let kind = self.model.kind(term.ty);
            let (field, s, t) = match *kind {
                Kind::PExt(prime, s, t) if allow_odd => (ScalarField::Odd(prime), s, t),
                _ => {
                    let expected = if matches!(source, Some(ScalarField::Odd(_))) {
                        "a PExt class"
                    } else {
                        "an Ext class"
                    };
                    let (profile, s, t) = kind
                        .ext_grade()
                        .ok_or_else(|| self.invalid_operand(&term, operation, expected))?;
                    (ScalarField::ModTwo(profile), s, t)
                }
            };
            if let Some(previous) = source {
                match (previous, field) {
                    (ScalarField::ModTwo(p), ScalarField::ModTwo(q)) => {
                        require_profile(p, q, operation, term.span)?;
                    }
                    (ScalarField::Odd(p), ScalarField::Odd(q)) if p == q => {}
                    (ScalarField::Odd(_), ScalarField::Odd(_)) => {
                        return Err(Error::new(
                            term.span,
                            format!("{operation} requires PExt classes at the same prime"),
                        ));
                    }
                    (ScalarField::Odd(_), ScalarField::ModTwo(_)) => {
                        return Err(self.invalid_operand(&term, operation, "a PExt class"));
                    }
                    (ScalarField::ModTwo(_), ScalarField::Odd(_)) => {
                        return Err(self.invalid_operand(&term, operation, "an Ext class"));
                    }
                }
            }
            source = Some(field);
            degrees.push((s, t));
            terms.push(term);
        }
        let profile = source
            .ok_or_else(|| Error::plain(format!("{operation} requires at least one Ext class")))?;
        Ok((profile, terms, degrees))
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
                let (n, _, s, t, _) = self.ce_degree(argument, name)?;
                return self.model.profile_cobar_type(Profile::An(n), s, t);
            }
            CobarOp::CeDifferential => {
                let (n, r, s, t, f) = self.ce_degree(argument, name)?;
                let s = checked_degree(s, 1, name, span)?;
                let f = checked_degree(f, r, name, span)?;
                return self.model.ce_type(n, r, s, t, f);
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
                if let Kind::PMay(p, _, s, t, _) = *self.model.kind(argument.ty) {
                    return self.model.pcobar_type(p, s, t);
                }
                if let Kind::MapMay(source, target, _, s, t, _) = *self.model.kind(argument.ty) {
                    return self.model.map_cobar_type(source, target, s, t);
                }
                let (profile, _, s, t, _) = self.may_degree(argument, name)?;
                return self.model.profile_cobar_type(profile, s, t);
            }
            CobarOp::MayDifferential => {
                if let Kind::PMay(p, r, s, t, f) = *self.model.kind(argument.ty) {
                    let s = checked_degree(s, 1, name, span)?;
                    let f = checked_degree(f, r, name, span)?;
                    return self.model.pmay_type(p, r, s, t, f);
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
                        | Kind::CobarAn(..)
                        | Kind::MapCobar(..)
                        | Kind::QKoszul(..)
                        | Kind::QExt(..)
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
                if let Kind::PMay(p, r, s, t, f) = *self.model.kind(left.ty) {
                    let Kind::PMay(q, page, u, v, g) = *self.model.kind(right.ty) else {
                        return Err(self.invalid_operand(right, "mayMul", "a PMay class"));
                    };
                    if p != q {
                        return Err(Error::new(right.span, "mayMul requires the same prime"));
                    }
                    if r != page {
                        return Err(Error::new(
                            right.span,
                            "mayMul requires classes on the same page",
                        ));
                    }
                    let s = checked_degree(s, u, "mayMul", span)?;
                    let t = checked_degree(t, v, "mayMul", span)?;
                    let f = checked_degree(f, g, "mayMul", span)?;
                    return self.model.pmay_type(p, r, s, t, f);
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
