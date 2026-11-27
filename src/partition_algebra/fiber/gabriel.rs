use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::lee_linear::{Echelon, Rational};
use crate::model::Model;

use super::loewy::{clear_denominators, integral_row};
use super::{add_multiple, ratio, Coordinates, Diagram, Element, Fiber, Mode, Parameter};

mod polynomial;

const ONE: Rational = Rational::integer(1);
const NEGATIVE_ONE: Rational = Rational::integer(-1);

#[derive(Clone, Copy)]
pub(crate) enum Generator {
    Primitive(u128),
    Vertex(u128),
    Arrow(u128),
    Square([u128; 3]),
}

pub(crate) struct Gabriel {
    parameter: Parameter,
    primitive: Vec<Coordinates>,
    groups: Vec<Vec<usize>>,
    cartan: Vec<Vec<usize>>,
    radical: Vec<Vec<usize>>,
    squares: Vec<Vec<Span>>,
    arrows: Vec<Arrow>,
    quadratic: Vec<Quadratic>,
    relations: Vec<Coordinates>,
    cubic: Vec<[usize; 3]>,
    minimal_cubic: Vec<[usize; 3]>,
}

struct Arrow {
    source: usize,
    target: usize,
    lift: Coordinates,
}

struct Quadratic {
    arrows: [usize; 2],
    coordinates: Coordinates,
}

#[derive(Default)]
struct Span {
    basis: Vec<Coordinates>,
    echelon: Echelon,
}

impl Span {
    fn insert(&mut self, model: &mut Model, value: Coordinates) -> Result<bool> {
        let (row, denominator) = clear_denominators(&value)?;
        if self
            .echelon
            .insert(model, row, vec![(self.basis.len(), denominator)])?
            .is_none()
        {
            self.basis.push(value);
            Ok(true)
        } else {
            Ok(false)
        }
    }

    fn coordinates(&self, model: &mut Model, value: &Coordinates) -> Result<Coordinates> {
        let (row, denominator) = clear_denominators(value)?;
        let relation = self
            .echelon
            .relation(model, row, vec![(self.basis.len(), denominator)])?
            .ok_or_else(|| {
                Error::plain("partition path product escaped its radical-square corner")
            })?;
        let mut coordinates = normalized_relation(relation)?;
        coordinates.pop_last();
        for value in coordinates.values_mut() {
            *value = super::negate(*value)?;
        }
        Ok(coordinates)
    }
}

fn normalized_relation(row: crate::lee_linear::Vector) -> Result<Coordinates> {
    let (_, leading) = *row
        .last()
        .ok_or_else(|| Error::plain("partition linear dependence has no leading term"))?;
    row.into_iter()
        .map(|(i, c)| Ok((i, ratio(c, leading)?)))
        .collect()
}

impl Gabriel {
    pub(super) fn new(model: &mut Model, fiber: &Fiber) -> Result<Self> {
        let parameter = fiber.parameter;
        if parameter.rank > 2 {
            return Err(Error::plain(
                "partition primitive decompositions and Gabriel quivers are limited to rank 2",
            ));
        }
        let size = fiber.size();
        let diagrams: Vec<_> = (0..size).map(basis).collect();
        let identity = diagrams[fiber.table.index[&Diagram::identity(parameter.rank)]].clone();
        let mut primitive = Vec::new();
        split(model, fiber, &diagrams, identity.clone(), &mut primitive)?;
        certify_primitives(model, fiber, &identity, &primitive)?;
        let groups = simple_groups(model, fiber, &diagrams, &primitive)?;
        let simple: Vec<_> = groups.iter().map(Vec::len).collect();
        let semisimple = size - fiber.radical(model)?.rows.len();
        if simple.iter().map(|&d| d * d).sum::<usize>() != semisimple {
            return Err(Error::plain(
                "partition simple dimensions do not exhaust the semisimple quotient",
            ));
        }
        let reps: Vec<_> = groups.iter().map(|group| &primitive[group[0]]).collect();
        let radical_kernel = fiber.radical(model)?;
        let powers = fiber.powers(model)?;
        if !powers.basis(3).is_empty() {
            return Err(Error::plain(
                "partition Gabriel presentation requires radical cube zero",
            ));
        }
        let mut cartan = Vec::new();
        let mut radical = Vec::new();
        let mut squares = Vec::new();
        let mut arrows = Vec::new();
        for (target, &e) in reps.iter().enumerate() {
            let mut carts = Vec::new();
            let mut rads = Vec::new();
            let mut sqs = Vec::new();
            for (source, &f) in reps.iter().enumerate() {
                let whole = corner(model, fiber, e, f, diagrams.iter())?;
                let square = corner(model, fiber, e, f, powers.basis(2).iter())?;
                let mut lifted = Span::default();
                for term in &square.basis {
                    lifted.insert(model, term.clone())?;
                }
                for x in radical_kernel.rows.values() {
                    let value = sandwich(model, fiber, e, x, f)?;
                    if lifted.insert(model, value.clone())? {
                        arrows.push(Arrow {
                            source,
                            target,
                            lift: value,
                        });
                    }
                }
                if whole.basis.len().checked_sub(lifted.basis.len())
                    != Some(usize::from(source == target))
                {
                    return Err(Error::plain(
                        "partition basic corner did not have the expected semisimple dimension",
                    ));
                }
                carts.push(whole.basis.len());
                rads.push(lifted.basis.len());
                sqs.push(square);
            }
            cartan.push(carts);
            radical.push(rads);
            squares.push(sqs);
        }
        let weighted: usize = cartan
            .iter()
            .enumerate()
            .map(|(j, row)| {
                row.iter()
                    .enumerate()
                    .map(|(i, &c)| c * simple[i] * simple[j])
                    .sum::<usize>()
            })
            .sum();
        if weighted != size {
            return Err(Error::plain(
                "partition Cartan matrix does not reconstruct the diagram dimension",
            ));
        }
        let paths = path_presentation(model, fiber, &arrows, &squares)?;
        let basic: usize = cartan.iter().flatten().sum();
        if groups.len() + arrows.len() + paths.quadratic.len() - paths.relations.len() != basic {
            return Err(Error::plain(
                "partition path presentation does not reconstruct the basic dimension",
            ));
        }
        let cells = primitive
            .iter()
            .chain(arrows.iter().map(|arrow| &arrow.lift))
            .chain(squares.iter().flatten().flat_map(|span| &span.basis))
            .chain(paths.quadratic.iter().map(|p| &p.coordinates))
            .chain(paths.relations.iter())
            .map(|v| 4 + 10 * v.len())
            .sum::<usize>();
        model.allocate(
            cells + 6 * (paths.cubic.len() + paths.minimal_cubic.len()) + 8 * groups.len().pow(2),
        )?;
        Ok(Self {
            parameter,
            primitive,
            groups,
            cartan,
            radical,
            squares,
            arrows,
            quadratic: paths.quadratic,
            relations: paths.relations,
            cubic: paths.cubic,
            minimal_cubic: paths.minimal_cubic,
        })
    }

    pub(crate) fn data(&self) -> String {
        let simples: Vec<_> = self.groups.iter().map(Vec::len).collect();
        let basic: usize = self.cartan.iter().flatten().sum();
        let projectives: Vec<_> = (0..simples.len())
            .map(|i| {
                self.cartan
                    .iter()
                    .zip(&simples)
                    .map(|(row, d)| row[i] * d)
                    .sum::<usize>()
            })
            .collect();
        let mut layers = Vec::new();
        let mut ext = Vec::new();
        for (j, row) in self.cartan.iter().enumerate() {
            let mut projective_layers = Vec::new();
            let mut ext_row = Vec::new();
            for (i, &total) in row.iter().enumerate() {
                let radical = self.radical[j][i];
                let square = self.squares[j][i].basis.len();
                projective_layers.push([total - radical, radical - square, square]);
                ext_row.push(radical - square);
            }
            layers.push(projective_layers);
            ext.push(ext_row);
        }
        let mut ext2 = vec![vec![0; simples.len()]; simples.len()];
        for relation in &self.relations {
            let (&i, _) = relation.last_key_value().unwrap();
            let [a, b] = self.quadratic[i].arrows;
            ext2[self.arrows[b].target][self.arrows[a].source] += 1;
        }
        for &[a, _, c] in &self.minimal_cubic {
            ext2[self.arrows[c].target][self.arrows[a].source] += 1;
        }
        let arrows: Vec<_> = self.arrows.iter().map(|a| (a.source, a.target)).collect();
        format!("rvPartitionGabriel[{},{}] = simples={simples:?}, primitives={:?}, projectives={projectives:?}, basicDimension={basic}, cartan={:?}, layers={layers:?}, ext1={ext:?}, ext2={ext2:?}, arrows={arrows:?}", self.parameter.rank, self.parameter.n, self.groups, self.cartan)
    }

    pub(crate) fn relations(&self) -> String {
        let paths: Vec<_> = self.quadratic.iter().map(|q| q.arrows).collect();
        let coordinates = self
            .quadratic
            .iter()
            .map(|q| report_coordinates(&q.coordinates))
            .collect::<Vec<_>>()
            .join(", ");
        let relations = self
            .relations
            .iter()
            .map(report_coordinates)
            .collect::<Vec<_>>()
            .join(", ");
        format!("rvPartitionGabrielRelations[{},{}] = quadraticPaths={paths:?}, products=[{coordinates}], quadraticRelations=[{relations}], zeroCubicPaths={:?}, minimalCubicRelations={:?}", self.parameter.rank, self.parameter.n, self.cubic, self.minimal_cubic)
    }

    pub(super) fn generator(&self, fiber: &Fiber, generator: Generator) -> Result<Element> {
        let i = |key: u128| usize::try_from(key).unwrap_or(usize::MAX);
        let value = match generator {
            Generator::Primitive(k) => self.primitive.get(i(k)),
            Generator::Vertex(k) => self.groups.get(i(k)).map(|g| &self.primitive[g[0]]),
            Generator::Arrow(k) => self.arrows.get(i(k)).map(|a| &a.lift),
            Generator::Square([target, source, k]) => self
                .squares
                .get(i(target))
                .and_then(|row| row.get(i(source)))
                .and_then(|corner| corner.basis.get(i(k))),
        };
        Ok(fiber.element(
            Mode::Full,
            value
                .cloned()
                .ok_or_else(|| Error::plain("partition Gabriel generator index is out of range"))?,
        ))
    }
}

fn report_coordinates(entries: &Coordinates) -> String {
    let terms = entries
        .iter()
        .map(|(i, value)| format!("({i},{value})"))
        .collect::<Vec<_>>()
        .join(", ");
    format!("[{terms}]")
}

fn basis(index: usize) -> Coordinates {
    Coordinates::from([(index, ONE)])
}

fn sandwich(
    model: &mut Model,
    fiber: &Fiber,
    e: &Coordinates,
    value: &Coordinates,
    f: &Coordinates,
) -> Result<Coordinates> {
    let right = fiber.multiply(model, value, f)?;
    fiber.multiply(model, e, &right)
}

fn corner<'a>(
    model: &mut Model,
    fiber: &Fiber,
    e: &Coordinates,
    f: &Coordinates,
    values: impl Iterator<Item = &'a Coordinates>,
) -> Result<Span> {
    let mut span = Span::default();
    for value in values {
        let at = sandwich(model, fiber, e, value, f)?;
        span.insert(model, at)?;
    }
    Ok(span)
}

fn semisimple_corner(
    model: &mut Model,
    fiber: &Fiber,
    diagrams: &[Coordinates],
    e: &Coordinates,
    f: &Coordinates,
) -> Result<usize> {
    let mut span = Echelon::default();
    for diagram in diagrams {
        let value = sandwich(model, fiber, e, diagram, f)?;
        let quotient = fiber.project_coordinates(model, Mode::Semisimple, value)?;
        span.insert(model, integral_row(&quotient.coefficients)?, Vec::new())?;
    }
    Ok(span.rank())
}

fn minimal_polynomial(
    model: &mut Model,
    fiber: &Fiber,
    e: &Coordinates,
    x: &Coordinates,
) -> Result<Vec<Rational>> {
    let mut span = Echelon::default();
    let mut power = e.clone();
    for exponent in 0..=fiber.size() {
        let (row, denominator) = clear_denominators(&power)?;
        if let Some(relation) = span.insert(model, row, vec![(exponent, denominator)])? {
            let coefficients = normalized_relation(relation)?;
            if coefficients.last_key_value() != Some((&exponent, &ONE)) {
                return Err(Error::plain(
                    "partition minimal polynomial has an incorrect degree",
                ));
            }
            return Ok((0..=exponent)
                .map(|j| coefficients.get(&j).copied().unwrap_or(Rational::ZERO))
                .collect());
        }
        power = fiber.multiply(model, &power, x)?;
    }
    Err(Error::plain(
        "partition corner exceeded the minimal polynomial degree bound",
    ))
}

fn polynomial_on(
    model: &mut Model,
    fiber: &Fiber,
    e: &Coordinates,
    x: &Coordinates,
    polynomial: &[Rational],
) -> Result<Coordinates> {
    let mut value = Coordinates::new();
    for &coefficient in polynomial.iter().rev() {
        value = fiber.multiply(model, x, &value)?;
        add_multiple(&mut value, e, coefficient)?;
    }
    Ok(value)
}

fn split(
    model: &mut Model,
    fiber: &Fiber,
    diagrams: &[Coordinates],
    e: Coordinates,
    output: &mut Vec<Coordinates>,
) -> Result<()> {
    match semisimple_corner(model, fiber, diagrams, &e, &e)? {
        0 => {
            return Err(Error::plain(
                "nonzero partition idempotent has a zero semisimple corner",
            ))
        }
        1 => {
            output.push(e);
            return Ok(());
        }
        _ => (),
    }
    let candidates = diagrams
        .iter()
        .cloned()
        .chain((0..diagrams.len()).flat_map(|i| {
            (i + 1..diagrams.len()).map(move |j| Coordinates::from([(i, ONE), (j, ONE)]))
        }));
    for candidate in candidates {
        let x = sandwich(model, fiber, &e, &candidate, &e)?;
        let modulus = minimal_polynomial(model, fiber, &e, &x)?;
        let Some(projector) = polynomial::splitting_projector(model, &modulus, fiber.parameter.n)?
        else {
            continue;
        };
        let h = polynomial_on(model, fiber, &e, &x, &projector)?;
        if h.is_empty()
            || h == e
            || fiber.multiply(model, &h, &h)? != h
            || sandwich(model, fiber, &e, &h, &e)? != h
        {
            return Err(Error::plain(
                "partition polynomial CRT did not give a proper corner projector",
            ));
        }
        let mut complement = e;
        add_multiple(&mut complement, &h, NEGATIVE_ONE)?;
        split(model, fiber, diagrams, h, output)?;
        split(model, fiber, diagrams, complement, output)?;
        return Ok(());
    }
    Err(Error::plain(
        "partition primitive splitting exhausted its diagram and rational root search bound",
    ))
}

fn certify_primitives(
    model: &mut Model,
    fiber: &Fiber,
    identity: &Coordinates,
    primitive: &[Coordinates],
) -> Result<()> {
    let mut sum = Coordinates::new();
    for (i, e) in primitive.iter().enumerate() {
        add_multiple(&mut sum, e, ONE)?;
        for (j, f) in primitive.iter().enumerate() {
            let product = fiber.multiply(model, e, f)?;
            if (i == j && product != *e) || (i != j && !product.is_empty()) {
                return Err(Error::plain(
                    "partition primitive idempotents failed orthogonality",
                ));
            }
        }
    }
    if sum != *identity {
        return Err(Error::plain(
            "partition primitive idempotents did not sum to the identity",
        ));
    }
    Ok(())
}

fn simple_groups(
    model: &mut Model,
    fiber: &Fiber,
    diagrams: &[Coordinates],
    primitive: &[Coordinates],
) -> Result<Vec<Vec<usize>>> {
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for (i, e) in primitive.iter().enumerate() {
        let mut matched = false;
        for group in &mut groups {
            let f = &primitive[group[0]];
            match semisimple_corner(model, fiber, diagrams, e, f)? {
                0 => (),
                1 if !matched => {
                    group.push(i);
                    matched = true;
                }
                _ => {
                    return Err(Error::plain(
                        "partition primitives gave incompatible simple corners",
                    ))
                }
            }
        }
        if !matched {
            groups.push(vec![i]);
        }
    }
    Ok(groups)
}

struct Paths {
    quadratic: Vec<Quadratic>,
    relations: Vec<Coordinates>,
    cubic: Vec<[usize; 3]>,
    minimal_cubic: Vec<[usize; 3]>,
}

fn path_presentation(
    model: &mut Model,
    fiber: &Fiber,
    arrows: &[Arrow],
    squares: &[Vec<Span>],
) -> Result<Paths> {
    let mut paths = Paths {
        quadratic: Vec::new(),
        relations: Vec::new(),
        cubic: Vec::new(),
        minimal_cubic: Vec::new(),
    };
    let mut pivots = BTreeMap::<(usize, usize), Echelon>::new();
    for (a, left) in arrows.iter().enumerate() {
        for (b, right) in arrows
            .iter()
            .enumerate()
            .filter(|(_, right)| right.source == left.target)
        {
            let product = fiber.multiply(model, &right.lift, &left.lift)?;
            let (source, target) = (left.source, right.target);
            let span = &squares[target][source];
            let coordinates = span.coordinates(model, &product)?;
            let (row, denominator) = clear_denominators(&coordinates)?;
            if let Some(relation) = pivots.entry((source, target)).or_default().insert(
                model,
                row,
                vec![(paths.quadratic.len(), denominator)],
            )? {
                paths.relations.push(normalized_relation(relation)?);
            }
            paths.quadratic.push(Quadratic {
                arrows: [a, b],
                coordinates,
            });
            for (c, last) in arrows
                .iter()
                .enumerate()
                .filter(|(_, last)| last.source == right.target)
            {
                if !fiber.multiply(model, &last.lift, &product)?.is_empty() {
                    return Err(Error::plain("partition basic cubic path did not vanish"));
                }
                paths.cubic.push([a, b, c]);
            }
        }
    }
    for (j, row) in squares.iter().enumerate() {
        for (i, square) in row.iter().enumerate() {
            if pivots.get(&(i, j)).map_or(0, Echelon::rank) != square.basis.len() {
                return Err(Error::plain(
                    "partition quadratic paths did not span their radical-square corner",
                ));
            }
        }
    }
    paths.minimal_cubic = minimal_cubic(model, arrows, &paths)?;
    Ok(paths)
}

fn minimal_cubic(model: &mut Model, arrows: &[Arrow], paths: &Paths) -> Result<Vec<[usize; 3]>> {
    let indices: BTreeMap<_, _> = paths.cubic.iter().copied().zip(0..).collect();
    let mut span = Echelon::default();
    for relation in &paths.relations {
        let (&first, _) = relation.first_key_value().unwrap();
        let [a, b] = paths.quadratic[first].arrows;
        let (source, target) = (arrows[a].source, arrows[b].target);
        for (c, arrow) in arrows.iter().enumerate() {
            if arrow.source == target {
                let row = relation
                    .iter()
                    .map(|(&i, &coefficient)| {
                        let [a, b] = paths.quadratic[i].arrows;
                        (indices[&[a, b, c]], coefficient)
                    })
                    .collect();
                span.insert(model, integral_row(&row)?, Vec::new())?;
            }
            if arrow.target == source {
                let row = relation
                    .iter()
                    .map(|(&i, &coefficient)| {
                        let [a, b] = paths.quadratic[i].arrows;
                        (indices[&[c, a, b]], coefficient)
                    })
                    .collect();
                span.insert(model, integral_row(&row)?, Vec::new())?;
            }
        }
    }
    let mut independent = Vec::new();
    for (i, &path) in paths.cubic.iter().enumerate() {
        if span.insert(model, vec![(i, 1)], Vec::new())?.is_none() {
            independent.push(path);
        }
    }
    Ok(independent)
}
