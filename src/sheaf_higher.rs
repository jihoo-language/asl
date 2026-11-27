use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result, Span};
use crate::model::{Kind, Model, TypeId, Value};
use crate::odd_cobar::Prime;
use crate::odd_ext;
use crate::prime_linear::{Basis, Coordinate};
use crate::sheaf::{self, Cochain, Cochains, Form, Layer, Term};
use crate::sheaf_ext::{self, Ends, Hom};

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Class,
    Curvature,
    Solve,
    Contains,
    Witness {
        from: usize,
        to: usize,
        targeted: bool,
    },
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Class => "shMassey",
            Self::Curvature => "shCurvature",
            Self::Solve => "shMasseySolve",
            Self::Contains => "shInMassey",
            Self::Witness {
                targeted: false, ..
            } => "shMasseyWitness",
            Self::Witness { targeted: true, .. } => "shMasseyFind",
        }
    }

    pub(crate) fn explicit(self) -> bool {
        matches!(self, Self::Class | Self::Curvature)
    }

    pub(crate) fn targeted(self) -> bool {
        matches!(self, Self::Contains | Self::Witness { targeted: true, .. })
    }
}

pub(crate) struct Signature {
    subjects: Vec<Rc<str>>,
    locations: Vec<Span>,
}

impl Signature {
    pub(crate) fn new(model: &Model, arguments: &[(TypeId, Span)]) -> Result<Self> {
        if arguments.len() < 2 {
            return Err(Error::plain(
                "a sheaf defining system requires at least two degree-one inputs",
            ));
        }
        let side = arguments.len() + 1;
        sheaf::check_cells(12usize.saturating_mul(side).saturating_mul(side))?;
        let mut subjects = Vec::with_capacity(side);
        let mut locations = Vec::with_capacity(arguments.len());
        for &(ty, span) in arguments {
            let Kind::SheafHom(_, ends, 1) = model.kind(ty) else {
                return Err(Error::new(
                    span,
                    "sheaf defining system inputs must have degree one in ShExt or ShHomC",
                ));
            };
            if let Some(previous) = subjects.last() {
                if previous != &ends.source {
                    return Err(Error::new(span, format!(
                        "sheaf defining system requires composable inputs, expected source {previous}, found {}", ends.source
                    )));
                }
            } else {
                subjects.push(ends.source.clone());
            }
            subjects.push(ends.target.clone());
            locations.push(span);
        }
        Ok(Self {
            subjects,
            locations,
        })
    }

    fn arity(&self) -> usize {
        self.subjects.len() - 1
    }

    fn ends(&self, from: usize, to: usize) -> Ends {
        Ends {
            source: self.subjects[from].clone(),
            target: self.subjects[to].clone(),
        }
    }

    pub(crate) fn top_type(&self, model: &mut Model, form: Form) -> Result<TypeId> {
        model.known_sheaf_hom_type(form, self.ends(0, self.arity()), 2)
    }

    pub(crate) fn witness_type(&self, model: &mut Model, from: usize, to: usize) -> Result<TypeId> {
        if to > self.arity()
            || !to
                .checked_sub(from)
                .is_some_and(|span| (2..self.arity()).contains(&span))
        {
            return Err(Error::plain(format!(
                "invalid proper interval [{from},{to}] in sheaf defining system"
            )));
        }
        model.known_sheaf_hom_type(Form::Cochain, self.ends(from, to), 1)
    }
}

pub(crate) struct Witness {
    pub(crate) from: usize,
    pub(crate) to: usize,
    pub(crate) span: Span,
    pub(crate) value: Value,
}

pub(crate) fn evaluate(
    model: &mut Model,
    signature: &Signature,
    operation: Operation,
    arguments: Vec<Value>,
    witnesses: Vec<Witness>,
    expected: Option<&[Coordinate]>,
) -> Result<Value> {
    let mut system = System::new(model, signature, arguments)?;
    if operation.explicit() {
        for witness in &witnesses {
            let Value::SheafC(value) = &witness.value else {
                unreachable!()
            };
            system.write(witness.from, witness.to, value.to_vec())?;
        }
        for (from, to) in system.intervals() {
            let defect = system.curvature(model, from, to)?;
            let actual = system.entry(from, to);
            if actual.is_empty() && defect.is_empty() {
                continue;
            }
            let hom = system.hom(model, from, to)?;
            let differential = hom.differential_support(model, actual)?;
            let defect = sum(model, system.prime, &differential, &defect)?;
            if !defect.is_empty() {
                let error = Error::plain(format!(
                    "sheaf defining system equation fails on interval [{from},{to}]"
                ));
                if let Some(entry) = witnesses
                    .iter()
                    .find(|entry| (entry.from, entry.to) == (from, to))
                {
                    return Err(error.at(entry.span));
                }
                return Err(error);
            }
        }
        let top = system.curvature(model, 0, signature.arity())?;
        return match operation {
            Operation::Curvature => sheaf::persist(model, top),
            _ => {
                let coordinates = system.class(model, &top)?;
                odd_ext::persist(model, coordinates)
            }
        };
    }
    let result = system.solve(model, expected)?;
    if matches!(operation, Operation::Contains) {
        return Ok(Value::Point(usize::from(result.is_some())));
    }
    let Some(class) = result else {
        let message = if expected.is_some() {
            "no sheaf defining system realizes the requested class"
        } else {
            "no sheaf defining system exists"
        };
        return Err(Error::plain(message));
    };
    match operation {
        Operation::Witness { from, to, .. } => sheaf::persist(model, system.take(from, to)),
        _ => odd_ext::persist(model, class),
    }
}

struct System<'a> {
    signature: &'a Signature,
    entries: Vec<Vec<Term>>,
    prime: Prime,
    cells: usize,
}

impl<'a> System<'a> {
    fn new(model: &mut Model, signature: &'a Signature, arguments: Vec<Value>) -> Result<Self> {
        let side = signature.subjects.len();
        let first = model.ensure_sheaf_hom(&signature.ends(0, 1))?;
        let prime = first.prime();
        let cells = 12 * side * side;
        model.spend(cells)?;
        let mut system = Self {
            signature,
            entries: vec![Vec::new(); side * side],
            prime,
            cells,
        };
        for (i, value) in arguments.into_iter().enumerate() {
            let hom = system.hom(model, i, i + 1)?;
            let terms = match value {
                Value::SheafC(raw) => {
                    if !hom.differential_support(model, &raw)?.is_empty() {
                        return Err(Error::new(
                            signature.locations[i],
                            "adjacent sheaf defining-system input must be closed",
                        ));
                    }
                    raw.to_vec()
                }
                Value::OddClass(class) => {
                    hom.cohomology
                        .representative_support(model, hom.as_ref(), 1, &class)?
                }
                _ => unreachable!(),
            };
            system.write(i, i + 1, terms)?;
        }
        Ok(system)
    }

    fn at(&self, from: usize, to: usize) -> usize {
        from * self.signature.subjects.len() + to
    }

    fn entry(&self, from: usize, to: usize) -> &[Term] {
        &self.entries[self.at(from, to)]
    }

    fn write(&mut self, from: usize, to: usize, terms: Vec<Term>) -> Result<()> {
        let index = self.at(from, to);
        let next = self.cells - term_cells(&self.entries[index]) + term_cells(&terms);
        sheaf::check_cells(next)?;
        self.entries[index] = terms;
        self.cells = next;
        Ok(())
    }

    fn take(&mut self, from: usize, to: usize) -> Vec<Term> {
        let index = self.at(from, to);
        std::mem::take(&mut self.entries[index])
    }

    fn intervals(&self) -> Vec<(usize, usize)> {
        let arity = self.signature.arity();
        (2..arity)
            .flat_map(|width| (0..=arity - width).map(move |from| (from, from + width)))
            .collect()
    }

    fn hom(&self, model: &mut Model, from: usize, to: usize) -> Result<Rc<Hom>> {
        model.ensure_sheaf_hom(&self.signature.ends(from, to))
    }

    fn curvature(&self, model: &mut Model, from: usize, to: usize) -> Result<Vec<Term>> {
        let mut result = Cochain::default();
        model.spend(to - from)?;
        for mid in from + 1..to {
            let left = self.entry(from, mid);
            let right = self.entry(mid, to);
            if left.is_empty() || right.is_empty() {
                continue;
            }
            let inner = self.hom(model, from, mid)?;
            let outer = self.hom(model, mid, to)?;
            let product = sheaf_ext::compose(model, (&outer, 1, right), (&inner, 1, left))?;
            for (term, coefficient) in product {
                result.add(model, self.prime, term, coefficient)?;
            }
        }
        Ok(result.into_terms())
    }

    fn class(&self, model: &mut Model, curvature: &[Term]) -> Result<Vec<Coordinate>> {
        if curvature.is_empty() {
            return Ok(Vec::new());
        }
        let hom = self.hom(model, 0, self.signature.arity())?;
        hom.cohomology
            .class_coordinates(model, hom.as_ref(), 2, curvature)
    }

    fn solve(
        &mut self,
        model: &mut Model,
        expected: Option<&[Coordinate]>,
    ) -> Result<Option<Vec<Coordinate>>> {
        let intervals = self.intervals();
        let mut frames: Vec<Frame> = Vec::new();
        let mut frame_cells = 0;
        loop {
            model.spend(1)?;
            let active = frames.len();
            if let Some(&(from, to)) = intervals.get(active) {
                let defect = self.curvature(model, from, to)?;
                let base = if defect.is_empty() {
                    Some(Vec::new())
                } else {
                    let complex = self.hom(model, from, to)?;
                    let primitive = Primitive::new(model, complex.as_ref())?;
                    sheaf::check_cells(self.cells + frame_cells + primitive.cells())?;
                    primitive.solve(model, self.prime, &defect)?
                };
                if let Some(base) = base {
                    self.write(from, to, base.clone())?;
                    frame_cells += term_cells(&base);
                    frames.push(Frame { base, digits: None });
                    sheaf::check_cells(self.cells + frame_cells)?;
                    continue;
                }
            } else {
                let top = self.curvature(model, 0, self.signature.arity())?;
                let class = self.class(model, &top)?;
                if expected.is_none_or(|value| class == value) {
                    return Ok(Some(class));
                }
            }

            // Traverse only H^1 choices. Degree-zero changes of basis in
            // this interval add coboundaries and change only larger intervals.
            loop {
                let Some(mut frame) = frames.pop() else {
                    return Ok(None);
                };
                frame_cells -= frame.cells();
                let (from, to) = intervals[frames.len()];
                let hom = self.hom(model, from, to)?;
                if let Some(next) = frame.advance(model, &hom)? {
                    self.write(from, to, next)?;
                    frame_cells += frame.cells();
                    sheaf::check_cells(self.cells + frame_cells)?;
                    frames.push(frame);
                    break;
                }
                self.write(from, to, Vec::new())?;
            }
        }
    }
}

struct Frame {
    base: Vec<Term>,
    digits: Option<Vec<u32>>,
}

impl Frame {
    fn cells(&self) -> usize {
        term_cells(&self.base) + self.digits.as_ref().map_or(0, Vec::len)
    }

    fn advance(&mut self, model: &mut Model, hom: &Hom) -> Result<Option<Vec<Term>>> {
        if self.digits.is_none() {
            let rank = hom.cohomology.rank(model, hom, 1)?;
            self.digits = Some(vec![0; rank]);
        }
        let digits = self.digits.as_mut().unwrap();
        let prime = hom.prime();
        for digit in digits.iter_mut() {
            model.spend(1)?;
            *digit += 1;
            if *digit != prime.get() {
                model.spend(digits.len())?;
                let value: Vec<_> = digits
                    .iter()
                    .copied()
                    .enumerate()
                    .filter(|&(_, c)| c != 0)
                    .collect();
                let rep = hom
                    .cohomology
                    .representative_support(model, hom, 1, &value)?;
                return Ok(Some(sum(model, prime, &self.base, &rep)?));
            }
            *digit = 0;
        }
        Ok(None)
    }
}

struct Primitive {
    source: Layer,
    target: Layer,
    images: Basis,
}

impl Primitive {
    fn new(model: &mut Model, hom: &Hom) -> Result<Self> {
        let source = hom.layer(model, 1, MAX_PIVOT_CELLS)?;
        let target = hom.layer(model, 2, MAX_PIVOT_CELLS)?;
        let mut images = Basis::default();
        for index in 0..source.size {
            let column = hom.indexed_column(model, &source.generator(index), &target)?;
            images.insert_with_witness(model, hom.prime(), column, vec![(index, 1)])?;
            sheaf::check_cells(source.cells() + target.cells() + images.cells())?;
        }
        Ok(Self {
            source,
            target,
            images,
        })
    }

    fn cells(&self) -> usize {
        self.source.cells() + self.target.cells() + self.images.cells()
    }

    fn solve(&self, model: &mut Model, prime: Prime, defect: &[Term]) -> Result<Option<Vec<Term>>> {
        let raw = defect
            .iter()
            .map(|(term, coefficient)| {
                Ok((self.target.index(model, term)?, prime.negate(*coefficient)))
            })
            .collect::<Result<Vec<_>>>()?;
        let (rest, primitive) = self.images.reduce_with_witness(model, prime, raw)?;
        Ok(rest.is_empty().then(|| {
            primitive
                .into_iter()
                .map(|(index, c)| (self.source.generator(index), c))
                .collect()
        }))
    }
}

fn term_cells(terms: &[Term]) -> usize {
    terms.iter().map(|(term, _)| term.cells()).sum()
}

fn sum(model: &mut Model, prime: Prime, left: &[Term], right: &[Term]) -> Result<Vec<Term>> {
    let mut output = Cochain::default();
    for (term, scalar) in left.iter().chain(right) {
        output.add(model, prime, term.clone(), *scalar)?;
    }
    Ok(output.into_terms())
}
