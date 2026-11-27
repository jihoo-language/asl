use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::Model;

use super::boundary::{arity, mass, Function, Node};
use super::{add, divide, integer, multiply, power, q, subtract, Rat, Q};

pub(crate) struct Parameter {
    z2: Rat,
    s: Rat,
    normalizer: Rat,
}

impl Parameter {
    pub(crate) fn new(z: Rat) -> Result<Self> {
        let z2 = multiply(z, z)?;
        if z2 == integer(1) {
            return Err(Error::plain(
                "normalized Knapp–Stein operator has a pole at z=1 or z=-1",
            ));
        }
        let s = divide(integer(Q * Q), z2)?;
        let normalizer = divide(
            subtract(integer(Q * Q), s)?,
            multiply(integer(Q + 1), subtract(q(), s)?)?,
        )?;
        Ok(Self { z2, s, normalizer })
    }

    pub(crate) fn eigenvalue(&self, r: usize) -> Result<Rat> {
        if r == 0 {
            return Ok(integer(1));
        }
        let numerator = subtract(integer(Q * Q), self.z2)?;
        if numerator.numerator == 0 {
            return Ok(Rat::ZERO);
        }
        let denominator = multiply(q(), subtract(self.z2, integer(1))?)?;
        let degree = i64::try_from(r - 1).map_err(|_| super::arithmetic())?;
        multiply(
            divide(numerator, denominator)?,
            power(divide(q(), self.z2)?, degree)?,
        )
    }

    fn tail(&self, r: usize) -> Result<Rat> {
        if r == 0 {
            return Ok(self.normalizer);
        }
        let ratio = divide(self.s, q())?;
        let denominator = subtract(integer(1), ratio)?;
        let initial = divide(integer(Q - 1), integer(Q + 1))?;
        let degree = i64::try_from(r).map_err(|_| super::arithmetic())?;
        divide(multiply(initial, power(ratio, degree)?)?, denominator)
    }

    pub(crate) fn transform(&self, model: &mut Model, f: &Function) -> Result<Function> {
        if let Node::Constant(_) = &f.node {
            return Ok(f.clone());
        }
        let data = Averages::new(model, &f.node, 0)?;
        let weights = Weights::new(self, data.depth)?;
        let by_kernel = data.kernel(model, &weights, 0, Rat::ZERO)?;
        let by_spectrum = data.spectrum(model, &weights, 0, data.mean)?;
        if by_kernel != by_spectrum {
            return Err(Error::plain(
                "normalized Knapp–Stein kernel disagrees with the martingale spectrum",
            ));
        }
        Ok(Function::from_node(by_kernel, f.presentation()))
    }

    pub(crate) fn cylinder_coefficient(
        &self,
        model: &mut Model,
        source: &[u8],
        target: &[u8],
    ) -> Result<Rat> {
        let common = source
            .iter()
            .zip(target)
            .take_while(|(a, b)| a == b)
            .count();
        let direct = if common == source.len() {
            Some(self.tail(source.len())?)
        } else if common < target.len() {
            let exponent = i64::try_from(common).map_err(|_| super::arithmetic())?;
            Some(multiply(mass(source.len())?, power(self.s, exponent)?)?)
        } else {
            None
        };
        let indicator = Function::characteristic(model, source)?;
        let image = self.transform(model, &indicator)?;
        let checked = image.value(target)?;
        if let Some(direct) = direct {
            if divide(direct, self.normalizer)? != checked {
                return Err(Error::plain(
                    "normalized Knapp–Stein cylinder coefficient is inconsistent",
                ));
            }
        }
        Ok(checked)
    }
}

struct Weights {
    prefix: Vec<Rat>,
    tail: Vec<Rat>,
    alpha: Vec<Rat>,
}

impl Weights {
    fn new(parameter: &Parameter, depth: usize) -> Result<Self> {
        let mut prefix = Vec::with_capacity(depth + 1);
        let mut tail = Vec::with_capacity(depth + 1);
        let mut alpha = Vec::with_capacity(depth + 1);
        let mut weighted_power = integer(1);
        for j in 0..=depth {
            prefix.push(divide(weighted_power, parameter.normalizer)?);
            tail.push(divide(parameter.tail(j)?, parameter.normalizer)?);
            alpha.push(parameter.eigenvalue(j)?);
            if j < depth {
                let ratio = divide(parameter.s, integer(i128::from(arity(j))))?;
                weighted_power = multiply(weighted_power, ratio)?;
            }
        }
        Ok(Self {
            prefix,
            tail,
            alpha,
        })
    }
}

struct Averages {
    mean: Rat,
    depth: usize,
    children: BTreeMap<u8, Averages>,
}

static ZERO: Averages = Averages {
    mean: Rat::ZERO,
    depth: 0,
    children: BTreeMap::new(),
};

impl Averages {
    fn new(model: &mut Model, node: &Node, level: usize) -> Result<Self> {
        model.spend(1)?;
        match node {
            Node::Constant(mean) => Ok(Self {
                mean: *mean,
                depth: 0,
                children: BTreeMap::new(),
            }),
            Node::Branch(slots) => {
                let mut total = Rat::ZERO;
                let mut depth = 0;
                let mut children = BTreeMap::new();
                for (&i, child) in slots {
                    let child = Self::new(model, child, level + 1)?;
                    total = add(total, child.mean)?;
                    depth = depth.max(child.depth + 1);
                    children.insert(i, child);
                }
                let mean = divide(total, integer(i128::from(arity(level))))?;
                Ok(Self {
                    mean,
                    depth,
                    children,
                })
            }
        }
    }

    fn kernel(
        &self,
        model: &mut Model,
        weights: &Weights,
        level: usize,
        offset: Rat,
    ) -> Result<Node> {
        model.spend(1)?;
        if self.children.is_empty() {
            return Ok(Node::Constant(add(
                offset,
                multiply(self.mean, weights.tail[level])?,
            )?));
        }
        let branches = arity(level);
        let mut result = BTreeMap::new();
        for i in 0..branches {
            let node = self.children.get(&i).unwrap_or(&ZERO);
            let side = subtract(self.mean, divide(node.mean, integer(i128::from(branches)))?)?;
            let contribution = multiply(weights.prefix[level], side)?;
            result.insert(
                i,
                node.kernel(model, weights, level + 1, add(offset, contribution)?)?,
            );
        }
        Ok(Node::collapse(level, result))
    }

    fn spectrum(
        &self,
        model: &mut Model,
        weights: &Weights,
        level: usize,
        sum: Rat,
    ) -> Result<Node> {
        model.spend(1)?;
        if self.children.is_empty() {
            return Ok(Node::Constant(sum));
        }
        let mut result = BTreeMap::new();
        for i in 0..arity(level) {
            let node = self.children.get(&i).unwrap_or(&ZERO);
            let increment = multiply(weights.alpha[level + 1], subtract(node.mean, self.mean)?)?;
            result.insert(
                i,
                node.spectrum(model, weights, level + 1, add(sum, increment)?)?,
            );
        }
        Ok(Node::collapse(level, result))
    }
}
