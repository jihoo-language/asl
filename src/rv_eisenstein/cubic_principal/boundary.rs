use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};

use super::{add, arithmetic, divide, integer, multiply, neg, power, q, subtract, Rat, Q};

type Children = BTreeMap<u8, Node>;

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(super) enum Node {
    Constant(Rat),
    Branch(Children),
}

const ZERO: Node = Node::Constant(Rat::ZERO);

pub(super) fn arity(depth: usize) -> u8 {
    (Q + i128::from(depth == 0)) as u8
}

impl Node {
    pub(super) fn collapse(depth: usize, mut children: Children) -> Self {
        children.retain(|_, v| *v != ZERO);
        if children.is_empty() {
            return ZERO;
        }
        if children.len() == usize::from(arity(depth)) {
            if let Some(Self::Constant(first)) = children.values().next() {
                if children.values().all(|v| v == &Self::Constant(*first)) {
                    return Self::Constant(*first);
                }
            }
        }
        Self::Branch(children)
    }

    fn child(&self, index: u8) -> &Self {
        match self {
            Self::Constant(_) => self,
            Self::Branch(children) => children.get(&index).unwrap_or(&ZERO),
        }
    }

    fn depth(&self) -> usize {
        match self {
            Self::Constant(_) => 0,
            Self::Branch(children) => 1 + children.values().map(Self::depth).max().unwrap_or(0),
        }
    }

    fn cells(&self) -> usize {
        match self {
            Self::Constant(_) => 5,
            Self::Branch(children) => 7 + children.values().map(|c| c.cells() + 2).sum::<usize>(),
        }
    }

    fn leaves(&self) -> usize {
        match self {
            Self::Constant(x) => usize::from(x.numerator != 0),
            Self::Branch(x) => x.values().map(Self::leaves).sum(),
        }
    }

    fn lookup(&self, key: &[u8]) -> &Self {
        let mut node = self;
        for &d in key {
            match node {
                Self::Constant(_) => break,
                Self::Branch(children) => node = children.get(&d).unwrap_or(&ZERO),
            }
        }
        node
    }

    fn value(&self, key: &[u8]) -> Result<Rat> {
        match self.lookup(key) {
            Self::Constant(c) => Ok(*c),
            Self::Branch(_) => Err(Error::plain(
                "boundary function varies on the requested cylinder",
            )),
        }
    }

    fn mean(&self, model: &mut Model, depth: usize) -> Result<Rat> {
        model.spend(1)?;
        let Self::Branch(children) = self else {
            let Self::Constant(c) = self else {
                unreachable!()
            };
            return Ok(*c);
        };
        let mut total = Rat::ZERO;
        for child in children.values() {
            total = add(total, child.mean(model, depth + 1)?)?;
        }
        divide(total, integer(i128::from(arity(depth))))
    }

    fn conditional(&self, model: &mut Model, level: usize, depth: usize) -> Result<Self> {
        if depth == level {
            return Ok(Self::Constant(self.mean(model, depth)?));
        }
        let Self::Branch(children) = self else {
            return Ok(self.clone());
        };
        let mut out = Children::new();
        for (&index, child) in children {
            out.insert(index, child.conditional(model, level, depth + 1)?);
        }
        Ok(Self::collapse(depth, out))
    }

    fn visit_leaves(
        &self,
        path: &mut Vec<u8>,
        visit: &mut impl FnMut(&[u8], Rat) -> Result<()>,
    ) -> Result<()> {
        match self {
            Self::Constant(c) if c.numerator != 0 => visit(path, *c),
            Self::Constant(_) => Ok(()),
            Self::Branch(children) => {
                for (&i, child) in children {
                    path.push(i);
                    child.visit_leaves(path, visit)?;
                    path.pop();
                }
                Ok(())
            }
        }
    }

    fn pointwise(
        &self,
        model: &mut Model,
        rhs: &Self,
        depth: usize,
        operation: Binary,
    ) -> Result<Self> {
        model.spend(1)?;
        if let (Self::Constant(a), Self::Constant(b)) = (self, rhs) {
            return Ok(Self::Constant(match operation {
                Binary::Add => add(*a, *b)?,
                Binary::Multiply => multiply(*a, *b)?,
            }));
        }
        match operation {
            Binary::Add if *self == ZERO => return Ok(rhs.clone()),
            Binary::Add if *rhs == ZERO => return Ok(self.clone()),
            Binary::Multiply if *self == ZERO || *rhs == ZERO => return Ok(ZERO),
            _ => (),
        }
        let mut output = Children::new();
        for index in 0..arity(depth) {
            let v = self
                .child(index)
                .pointwise(model, rhs.child(index), depth + 1, operation)?;
            output.insert(index, v);
        }
        Ok(Self::collapse(depth, output))
    }

    pub(super) fn add_indicator(
        &mut self,
        model: &mut Model,
        key: &[u8],
        depth: usize,
        c: Rat,
    ) -> Result<()> {
        model.spend(1)?;
        if c.numerator == 0 {
            return Ok(());
        }
        let Some((&head, tail)) = key.split_first() else {
            *self = self.pointwise(model, &Self::Constant(c), depth, Binary::Add)?;
            return Ok(());
        };
        let mut children = match std::mem::replace(self, ZERO) {
            Self::Branch(children) => children,
            Self::Constant(v) if v.numerator == 0 => Children::new(),
            Self::Constant(v) => (0..arity(depth)).map(|k| (k, Self::Constant(v))).collect(),
        };
        children
            .entry(head)
            .or_insert(ZERO)
            .add_indicator(model, tail, depth + 1, c)?;
        *self = Self::collapse(depth, children);
        Ok(())
    }
}

#[derive(Clone, Copy)]
enum Binary {
    Add,
    Multiply,
}

#[derive(Clone)]
pub(crate) struct Function {
    pub(super) node: Node,
    presentation: usize,
}

impl PartialEq for Function {
    fn eq(&self, other: &Self) -> bool {
        self.node == other.node
    }
}
impl Eq for Function {}
impl PartialOrd for Function {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Function {
    fn cmp(&self, other: &Self) -> Ordering {
        self.node.cmp(&other.node)
    }
}

impl Function {
    pub(crate) fn scalar(value: Rat) -> Self {
        Self::from_node(Node::Constant(value), 0)
    }

    pub(super) fn from_node(node: Node, presentation: usize) -> Self {
        Self {
            presentation: presentation.max(node.depth()),
            node,
        }
    }

    pub(super) fn presentation(&self) -> usize {
        self.presentation
    }

    pub(crate) fn characteristic(model: &mut Model, path: &[u8]) -> Result<Self> {
        let mut f = ZERO;
        f.add_indicator(model, path, 0, integer(1))?;
        Ok(Self::from_node(f, path.len()))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(1 + self.node.cells())?;
        Ok(Value::CubicBoundary(Rc::new(self)))
    }

    pub(crate) fn description(&self) -> String {
        format!(
            "level={}, minimal_level={}, nonzero_cylinders={}",
            self.presentation,
            self.node.depth(),
            self.node.leaves()
        )
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let node = self.node.pointwise(model, &rhs.node, 0, Binary::Add)?;
        Ok(Self::from_node(
            node,
            self.presentation.max(rhs.presentation),
        ))
    }

    pub(super) fn product(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let node = self.node.pointwise(model, &rhs.node, 0, Binary::Multiply)?;
        Ok(Self::from_node(
            node,
            self.presentation.max(rhs.presentation),
        ))
    }

    pub(crate) fn scale(&self, model: &mut Model, scalar: Rat) -> Result<Self> {
        let node = self
            .node
            .pointwise(model, &Node::Constant(scalar), 0, Binary::Multiply)?;
        Ok(Self::from_node(node, self.presentation))
    }

    pub(crate) fn refine(&self, level: usize) -> Result<Self> {
        if level < self.presentation {
            return Err(Error::plain(
                "boundary refinement requires a level at least the current presentation",
            ));
        }
        Ok(Self::from_node(self.node.clone(), level))
    }

    pub(crate) fn coarsen(&self, level: usize) -> Result<Self> {
        if level > self.presentation {
            return Err(Error::plain(
                "boundary coarsening requires a level at most the current presentation",
            ));
        }
        if level < self.node.depth() {
            return Err(Error::plain(
                "boundary function is not constant on all cylinders at that level",
            ));
        }
        Ok(Self::from_node(self.node.clone(), level))
    }

    pub(crate) fn average(&self, model: &mut Model, level: usize) -> Result<Self> {
        let node = self.node.conditional(model, level, 0)?;
        Ok(Self::from_node(node, level))
    }

    pub(crate) fn difference(&self, model: &mut Model, level: usize) -> Result<Self> {
        let scaled = self.average(model, level)?;
        if level == 0 {
            return Ok(scaled);
        }
        let earlier = self.average(model, level - 1)?.scale(model, integer(-1))?;
        scaled.add(model, &earlier)
    }

    pub(super) fn visit(&self, visit: &mut impl FnMut(&[u8], Rat) -> Result<()>) -> Result<()> {
        self.node.visit_leaves(&mut Vec::new(), visit)
    }

    pub(crate) fn value(&self, key: &[u8]) -> Result<Rat> {
        self.node.value(key)
    }

    pub(crate) fn integrate(&self, model: &mut Model, key: &[u8]) -> Result<Rat> {
        multiply(
            mass(key.len())?,
            self.node.lookup(key).mean(model, key.len())?,
        )
    }

    pub(crate) fn poisson(&self, model: &mut Model, path: &[u8], t: Rat) -> Result<Rat> {
        let mut result = Rat::ZERO;
        self.visit(&mut |leaf, c| {
            model.spend(1)?;
            result = add(result, multiply(c, poisson_measure(leaf, path, t)?)?)?;
            Ok(())
        })?;
        Ok(result)
    }
}

pub(crate) fn mass(depth: usize) -> Result<Rat> {
    if depth == 0 {
        return Ok(integer(1));
    }
    let denominator = multiply(integer(Q + 1), power(q(), as_exponent(depth - 1)?)?)?;
    divide(integer(1), denominator)
}

fn as_exponent(depth: usize) -> Result<i64> {
    i64::try_from(depth).map_err(|_| arithmetic())
}

fn poisson_measure(leaf: &[u8], path: &[u8], t: Rat) -> Result<Rat> {
    let d = as_exponent(path.len())?;
    let k = leaf.iter().zip(path).take_while(|(a, b)| a == b).count();
    if k < leaf.len() || k == path.len() {
        return multiply(mass(leaf.len())?, power(t, 2 * as_exponent(k)? - d)?);
    }
    let mut result = multiply(mass(path.len())?, power(t, d)?)?;
    for k in leaf.len()..path.len() {
        let side = subtract(mass(k)?, mass(k + 1)?)?;
        result = add(result, multiply(side, power(t, 2 * as_exponent(k)? - d)?)?)?;
    }
    Ok(result)
}

pub(super) fn complementary(
    node: &mut Node,
    model: &mut Model,
    complement: bool,
    key: &[u8],
    c: Rat,
) -> Result<()> {
    if complement {
        node.add_indicator(model, &[], 0, c)?;
        node.add_indicator(model, key, 0, neg(c)?)
    } else {
        node.add_indicator(model, key, 0, c)
    }
}
