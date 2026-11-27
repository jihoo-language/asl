use std::collections::BTreeSet;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

use super::cubic_division::Extension;
use super::cubic_principal::{transport_cylinders, Function};
use super::cubic_tree::{Operator, Tree, Vertex};

mod ihara;

const MAX_WORD: usize = 20_000;

pub(crate) struct Declaration(
    pub(crate) Label,
    pub(crate) Label,
    pub(crate) Vec<(usize, ParsedExpr)>,
);

pub(crate) enum Expression {
    Literal(Label, Vec<i128>),
    Direction(Label, i128),
    CoreVertex(Label, usize),
    Unary(Unary, Box<ParsedExpr>),
    Binary(Binary, Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Inverse,
    Power(i64),
    Cyclic,
    Prefix,
    Axis,
    Primitive,
}

#[derive(Clone, Copy)]
pub(crate) enum Binary {
    Multiply,
    Conjugate,
    Apply,
    Boundary,
}

pub(crate) enum Report {
    Data(Label, ReportKind),
    Word(WordReport, Box<ParsedExpr>),
    PingPong(Label, i128, Box<ParsedExpr>),
}

pub(crate) enum ReportKind {
    Marking,
    Core,
    Nonbacktracking(usize, usize),
    Weighted(i128, i128),
    Ihara,
    Coefficient(usize),
    Trace(usize),
    Primitive(usize),
    Enumerate(usize),
}

#[derive(Clone, Copy)]
pub(crate) enum WordReport {
    Lengths,
    Path,
    Axis,
}

type Letter = i8;

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Word(Vec<Letter>);

impl Word {
    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(2 + self.0.len())?;
        Ok(Value::CubicSchottky(Rc::new(self)))
    }

    pub(crate) fn description(&self) -> String {
        format!("{:?}", self.0)
    }

    pub(crate) fn inverse(&self, model: &mut Model) -> Result<Self> {
        model.spend(self.0.len())?;
        Ok(Self(self.0.iter().rev().map(|&a| -a).collect()))
    }

    pub(crate) fn multiply(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        model.spend(self.0.len() + rhs.0.len())?;
        let mut result = self.0.clone();
        for &a in &rhs.0 {
            push_reduced(&mut result, a);
        }
        check_word_size(result.len())?;
        Ok(Self(result))
    }

    fn cyclic_parts(&self) -> (&[Letter], &[Letter]) {
        let mut peel = 0;
        while 2 * peel + 1 < self.0.len() && self.0[peel] == -self.0[self.0.len() - peel - 1] {
            peel += 1;
        }
        (&self.0[..peel], &self.0[peel..self.0.len() - peel])
    }

    pub(crate) fn cyclic(&self, model: &mut Model) -> Result<Self> {
        model.spend(self.0.len())?;
        Ok(Self(self.cyclic_parts().1.to_vec()))
    }

    pub(crate) fn prefix(&self, model: &mut Model) -> Result<Self> {
        model.spend(self.0.len())?;
        Ok(Self(self.cyclic_parts().0.to_vec()))
    }

    pub(crate) fn power(&self, model: &mut Model, k: i64) -> Result<Self> {
        if k == 0 || self.0.is_empty() {
            return Ok(Self(Vec::new()));
        }
        let (prefix, core) = self.cyclic_parts();
        let repeat = usize::try_from(k.unsigned_abs()).map_err(|_| word_limit())?;
        let count = core
            .len()
            .checked_mul(repeat)
            .and_then(|x| x.checked_add(2 * prefix.len()))
            .ok_or_else(word_limit)?;
        check_word_size(count)?;
        model.spend(count)?;
        let mut output = Vec::with_capacity(count);
        output.extend_from_slice(prefix);
        for _ in 0..repeat {
            if k < 0 {
                output.extend(core.iter().rev().map(|&a| -a));
            } else {
                output.extend_from_slice(core);
            }
        }
        output.extend(prefix.iter().rev().map(|&a| -a));
        Ok(Self(output))
    }

    pub(crate) fn primitive(&self, model: &mut Model) -> Result<bool> {
        model.spend(self.0.len())?;
        let core = self.cyclic_parts().1;
        Ok(!core.is_empty() && period(core) == core.len())
    }

    pub(crate) fn conjugate(&self, model: &mut Model, rhs: &Self) -> Result<bool> {
        model.spend(self.0.len() + rhs.0.len())?;
        let lhs = self.cyclic_parts().1;
        let rhs = rhs.cyclic_parts().1;
        if lhs.len() != rhs.len() {
            return Ok(false);
        }
        let a = least_rotation(lhs);
        let b = least_rotation(rhs);
        Ok(lhs[a..]
            .iter()
            .chain(&lhs[..a])
            .eq(rhs[b..].iter().chain(&rhs[..b])))
    }
}

fn push_reduced(stack: &mut Vec<Letter>, a: Letter) {
    if stack.last() == Some(&-a) {
        stack.pop();
    } else {
        stack.push(a);
    }
}

fn check_word_size(length: usize) -> Result<()> {
    if length > MAX_WORD {
        Err(word_limit())
    } else {
        Ok(())
    }
}

fn word_limit() -> Error {
    Error::plain("materialized Schottky words support at most 20000 letters")
}

fn period(word: &[Letter]) -> usize {
    if word.is_empty() {
        return 0;
    }
    let mut failure = vec![0usize; word.len()];
    for i in 1..word.len() {
        let mut j = failure[i - 1];
        while j > 0 && word[i] != word[j] {
            j = failure[j - 1];
        }
        failure[i] = j + usize::from(word[i] == word[j]);
    }
    let p = word.len() - failure[word.len() - 1];
    if word.len().is_multiple_of(p) {
        p
    } else {
        word.len()
    }
}

fn least_rotation(word: &[Letter]) -> usize {
    let n = word.len();
    let (mut a, mut b, mut k) = (0, 1, 0);
    while a < n && b < n && k < n {
        match word[(a + k) % n].cmp(&word[(b + k) % n]) {
            std::cmp::Ordering::Equal => k += 1,
            std::cmp::Ordering::Greater => {
                a = (a + k + 1).max(b + 1);
                k = 0;
            }
            std::cmp::Ordering::Less => {
                b = (b + k + 1).max(a + 1);
                k = 0;
            }
        }
    }
    a.min(b)
}

struct Generator {
    length: usize,
    matrix: Operator,
    inverse: Operator,
    directions: [Vertex; 2],
    indices: [u8; 2],
}

pub(crate) struct Group {
    pub(crate) division_name: Box<str>,
    pub(crate) division: Rc<Extension>,
    generators: Vec<Generator>,
}

impl Group {
    pub(crate) fn new(
        model: &mut Model,
        division_name: Box<str>,
        division: Rc<Extension>,
        inputs: Vec<(usize, Operator)>,
    ) -> Result<Self> {
        if inputs.len() > 14 {
            return Err(Error::plain(
                "a 28-valent root permits at most fourteen signed Schottky direction pairs",
            ));
        }
        let mut generators = Vec::with_capacity(inputs.len());
        let tree = Tree::new(model, &division)?;
        let standard = [tree.at(model, 1, 27)?, tree.at(model, 1, 0)?];
        let mut used = BTreeSet::new();
        for (length, matrix) in inputs {
            if length == 0 {
                return Err(Error::plain(
                    "Schottky translation lengths must be positive",
                ));
            }
            if !matrix.is_integral() {
                return Err(Error::plain(
                    "Schottky markings require integral GL2(O_D) matrix operators",
                ));
            }
            let inverse = tree.inverse(model, &matrix)?;
            let directions = [
                tree.apply(model, &matrix, &standard[0])?,
                tree.apply(model, &matrix, &standard[1])?,
            ];
            let indices = [
                tree.address(model, &directions[0])?[0],
                tree.address(model, &directions[1])?[0],
            ];
            if !used.insert(indices[0]) || !used.insert(indices[1]) {
                return Err(Error::plain(
                    "Schottky ping-pong requires all signed first directions to be distinct",
                ));
            }
            generators.push(Generator {
                length,
                matrix,
                inverse,
                directions,
                indices,
            });
        }
        model.allocate(8 + 500 * generators.len())?;
        Ok(Self {
            division_name,
            division,
            generators,
        })
    }

    fn letter(&self, a: i128) -> Result<Letter> {
        if a == 0 || a.unsigned_abs() > self.generators.len() as u128 {
            return Err(Error::plain(
                "Schottky word letter is outside the signed generator range",
            ));
        }
        Ok(a as Letter)
    }

    fn marking(&self, name: &str) -> String {
        let lengths: Vec<_> = self.lengths().collect();
        let directions: Vec<_> = self.generators.iter().map(|g| g.indices).collect();
        format!("rvCubicSchottkyData[{name}] = division={}, rank={}, lengths={lengths:?}, directions={directions:?}, ping_pong=certified",
            self.division_name,self.generators.len())
    }

    pub(crate) fn literal(&self, model: &mut Model, inputs: &[i128]) -> Result<Word> {
        model.spend(inputs.len())?;
        let mut result = Vec::new();
        for &a in inputs {
            push_reduced(&mut result, self.letter(a)?);
            check_word_size(result.len())?;
        }
        Ok(Word(result))
    }

    fn generator(&self, a: Letter) -> &Generator {
        &self.generators[usize::from(a.unsigned_abs()) - 1]
    }

    fn lengths(&self) -> impl ExactSizeIterator<Item = usize> + '_ {
        self.generators.iter().map(|g| g.length)
    }

    pub(crate) fn direction(&self, a: i128) -> Result<Vertex> {
        let a = self.letter(a)?;
        Ok(self.generator(a).directions[usize::from(a < 0)].clone())
    }

    fn displacement(&self, word: &[Letter]) -> u128 {
        word.iter().map(|&a| self.generator(a).length as u128).sum()
    }

    fn radius(&self, word: &[Letter]) -> Result<usize> {
        let v = self.displacement(word);
        if v > self.division.precision as u128 {
            return Err(Error::plain(
                "Schottky root path requires more lattice precision",
            ));
        }
        Ok(v as usize)
    }

    fn apply_letter(
        &self,
        model: &mut Model,
        tree: &Tree<'_>,
        a: Letter,
        v: &Vertex,
    ) -> Result<Vertex> {
        let generator = self.generator(a);
        let mut image = tree.apply(model, &generator.inverse, v)?;
        if generator.length.saturating_sub(image.level()) > self.division.precision {
            return Err(Error::plain(
                "Schottky translation requires more lattice precision",
            ));
        }
        for _ in 0..generator.length {
            model.spend(1)?;
            image = if a > 0 {
                tree.translate(model, &image)?
            } else {
                tree.translate_back(model, &image)?
            };
        }
        tree.apply(model, &generator.matrix, &image)
    }

    fn action(
        &self,
        model: &mut Model,
        tree: &Tree<'_>,
        word: &[Letter],
        v: &Vertex,
    ) -> Result<Vertex> {
        let mut image = v.clone();
        for &a in word.iter().rev() {
            image = self.apply_letter(model, tree, a, &image)?;
        }
        Ok(image)
    }

    pub(crate) fn apply(&self, model: &mut Model, word: &Word, v: &Vertex) -> Result<Vertex> {
        let tree = Tree::new(model, &self.division)?;
        self.action(model, &tree, &word.0, v)
    }

    pub(crate) fn boundary(
        &self,
        model: &mut Model,
        word: &Word,
        f: &Function,
    ) -> Result<Function> {
        let tree = Tree::new(model, &self.division)?;
        transport_cylinders(model, &tree, f, |model, v| {
            self.action(model, &tree, &word.0, v)
        })
    }

    pub(crate) fn axis(&self, model: &mut Model, word: &Word) -> Result<Vertex> {
        let tree = Tree::new(model, &self.division)?;
        let prefix = word.cyclic_parts().0;
        self.radius(prefix)?;
        self.action(model, &tree, prefix, &tree.root())
    }

    pub(crate) fn ping_pong(&self, model: &mut Model, a: i128, v: &Vertex) -> Result<String> {
        let tree = Tree::new(model, &self.division)?;
        let a = self.letter(a)?;
        let generator = self.generator(a);
        if v.level() == 0 {
            return Err(Error::plain(
                "Schottky ping-pong needs a nonroot ray vertex",
            ));
        }
        let first = tree.ancestor(v, 1)?;
        let sign = usize::from(a < 0);
        if first == generator.directions[1 - sign] {
            return Err(Error::plain(
                "Schottky ping-pong input lies in the excluded repelling cylinder",
            ));
        }
        let image = self.apply_letter(model, &tree, a, v)?;
        if image.level() == 0 || tree.ancestor(&image, 1)? != generator.directions[sign] {
            return Err(inconsistent());
        }
        Ok(format!(
            "rvCubicSchottkyPingPong = letter={a}, image_radius={}, attracting_direction={}, verified=true",
            image.level(),generator.indices[sign],
        ))
    }

    pub(crate) fn word_report(
        &self,
        model: &mut Model,
        word: &Word,
        kind: WordReport,
    ) -> Result<String> {
        model.spend(word.0.len())?;
        let (prefix, core) = word.cyclic_parts();
        let rooted = self.displacement(&word.0);
        let translation = self.displacement(core);
        match kind {
            WordReport::Lengths => {
                let power = if core.is_empty() {
                    0
                } else {
                    core.len() / period(core)
                };
                Ok(format!("rvCubicSchottkyLengths = letters={}, displacement={rooted}, translation={translation}, prefix={prefix:?}, core={core:?}, power={power}",
                    word.0.len()))
            }
            WordReport::Axis => {
                let tree = Tree::new(model, &self.division)?;
                let v = self.axis(model, word)?;
                let endpoint = self.action(model, &tree, &word.0, &v)?;
                let distance = tree.distance(&v, &endpoint)?;
                if distance as u128 != translation {
                    return Err(inconsistent());
                }
                Ok(format!("rvCubicSchottkyAxisCheck = vertex_radius={}, image_radius={}, displacement={distance}, verified=true",
                    v.level(),endpoint.level()))
            }
            WordReport::Path => {
                let tree = Tree::new(model, &self.division)?;
                let radius = self.radius(&word.0)?;
                let root = tree.root();
                let end = self.action(model, &tree, &word.0, &root)?;
                if end.level() != radius {
                    return Err(inconsistent());
                }
                let mut junctions = vec![0];
                for i in 1..=word.0.len() {
                    let prefix = &word.0[..i];
                    let distance = self.radius(prefix)?;
                    let junction = self.action(model, &tree, prefix, &root)?;
                    if junction != tree.ancestor(&end, distance)? {
                        return Err(inconsistent());
                    }
                    junctions.push(distance);
                }
                let address = tree.address(model, &end)?;
                Ok(format!("rvCubicSchottkyPath = radius={radius}, junctions={junctions:?}, address={address:?}, verified=true"))
            }
        }
    }

    pub(crate) fn core_vertex(&self, model: &mut Model, index: usize) -> Result<Vertex> {
        let tree = Tree::new(model, &self.division)?;
        let mut remaining = index;
        if remaining == 0 {
            return Ok(tree.root());
        }
        for generator in &self.generators {
            if remaining < generator.length {
                Tree::require_radius(self.division.precision, remaining)?;
                let mut ray = vec![0; remaining];
                ray[0] = 27;
                let endpoint = tree.vertex_at_address(model, &ray)?;
                return tree.apply(model, &generator.matrix, &endpoint);
            }
            remaining -= generator.length - 1;
        }
        Err(Error::plain("Schottky core vertex index is out of range"))
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str, kind: ReportKind) -> Result<String> {
        if matches!(kind, ReportKind::Marking) {
            return Ok(self.marking(name));
        }
        if let ReportKind::Weighted(a, b) = kind {
            let a = self.letter(a)?;
            let b = self.letter(b)?;
            let coefficient = if a == -b {
                "0".to_owned()
            } else {
                format!("u^{}", self.generator(b).length)
            };
            return Ok(format!(
                "rvCubicSchottkyWeighted[{name},{a},{b}] = {coefficient}"
            ));
        }
        let core = ihara::Core::new(model, self.lengths().collect())?;
        core.report(model, name, kind)
    }
}

fn inconsistent() -> Error {
    Error::plain("inconsistent Schottky ping-pong, marked axis or geodesic certificate")
}
