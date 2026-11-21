use std::collections::HashMap;
use std::rc::Rc;

use crate::cohomology;
use crate::error::{Error, Result, Span};
use crate::model::{Kind, Model, Poset, SuspensionPoint, TypeId, Value};
use crate::observations::{betti, homotopy};
use crate::steenrod;
use crate::syntax::{
    BinaryOp, ClassForm, Entries, ExprForm, ParsedExpr, ParsedType, Statement, StatementForm,
    SuspensionForm, TypeForm,
};

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

enum TypedForm {
    Global(usize),
    Local(usize),
    Point(usize),
    Constant(Value),
    Binary(BinaryOp, Box<TypedExpr>, Box<TypedExpr>),
    Equator,
    Suspend(Box<TypedExpr>),
    Shift(Box<TypedExpr>),
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
                    (BinaryOp::Mul, Value::Operation(a), Value::Operation(b)) => {
                        steenrod::multiply(&mut self.model, &a, &b)
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
                        if index >= rank {
                            return Err(Error::new(
                                at,
                                format!(
                                    "generator index {index} is out of range for {} of rank {rank}",
                                    self.model.type_name(ty),
                                ),
                            ));
                        }
                        self.model.allocate(1)?;
                        Value::Class([index].into())
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
                    Kind::Cohomology(..) | Kind::Stable(..) => cohomology::zero(),
                    _ => {
                        return Err(Error::new(
                            span,
                            format!(
                                "zero requires a cohomology or Steenrod algebra type, found {}",
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
                if matches!(&function.kind, ExprForm::Name(name) if name.text == "susp")
                    && !self.bindings.contains_key("susp")
                    && !self.locals.iter().any(|(name, _)| name == "susp")
                {
                    let argument = self.expr(*argument)?;
                    let (base, degree) = self.class_type(&argument, "susp")?;
                    let degree = checked_degree(degree, 1, "susp", span)?;
                    let suspension = self.model.suspension_type(base)?;
                    let ty = self.model.cohomology_type(suspension, degree)?;
                    return Ok((ty, TypedForm::Suspend(Box::new(argument))));
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

    fn operation_degree(&self, value: &TypedExpr, op: &str) -> Result<usize> {
        match *self.model.kind(value.ty) {
            Kind::Operations(degree) => Ok(degree),
            _ => Err(self.invalid_operand(value, op, "a Steenrod operation")),
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
                    Kind::Cohomology(..) | Kind::Stable(..) | Kind::Operations(..)
                ) {
                    return Err(self.invalid_operand(
                        left,
                        "add",
                        "a cohomology class or Steenrod operation",
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
