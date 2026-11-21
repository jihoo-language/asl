use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::ops::RangeInclusive;

use crate::comodule::Generator;
use crate::error::{Error, Result, Span};
use crate::linear::xor;
use crate::model::{Model, TypeId};
use crate::q_koszul::{self, Indices};
use crate::syntax::Label;

const MAX_CELLS: usize = 300_000;

pub(crate) struct Rule {
    pub(crate) primitive: usize,
    pub(crate) span: Span,
    pub(crate) source: Label,
    pub(crate) targets: Vec<Label>,
}

pub(crate) struct Presentation {
    pub(crate) generators: Vec<(Label, usize)>,
    pub(crate) rules: Vec<Rule>,
}

pub(crate) struct Module {
    selected: Indices,
    ranks: BTreeMap<usize, usize>,
    images: HashMap<(usize, Generator), Box<[usize]>>,
}

pub(crate) fn generator(
    model: &mut Model,
    base: TypeId,
    degree: usize,
    index: usize,
) -> Result<Generator> {
    let rank = if let Some(module) = model.q_module(base) {
        module.rank(degree)
    } else {
        let ty = model.cohomology_type(base, degree)?;
        model.cohomology(ty)?.rank()
    };
    if index >= rank {
        return Err(Error::plain(format!(
            "qterm generator index {index} is out of range for {} in degree {degree} of rank {rank}",
            model.type_name(base),
        )));
    }
    Ok(Generator { degree, index })
}

pub(crate) fn term_bidegree(
    model: &mut Model,
    source: usize,
    target: usize,
    selected: &[usize],
    powers: &[usize],
) -> Result<(usize, i64)> {
    if selected.len() != powers.len() {
        return Err(Error::plain(
            "qterm requires one exponent per selected primitive",
        ));
    }
    model.spend(1 + powers.len())?;
    let mut length = 0usize;
    let mut degree = 0u128;
    for (&primitive, &power) in selected.iter().zip(powers) {
        length = length
            .checked_add(power)
            .ok_or_else(|| Error::plain("qterm cohomological degree is out of range"))?;
        // The sum fits u128 once the total exponent fits usize.
        degree += power as u128 * q_koszul::weight(primitive) as u128;
    }
    let degree_error = || Error::plain("qterm signed internal degree is out of range");
    let word = i128::try_from(degree).map_err(|_| degree_error())?;
    let signed = word
        .checked_add(target as i128 - source as i128)
        .ok_or_else(degree_error)?;
    let signed = i64::try_from(signed).map_err(|_| degree_error())?;
    Ok((length, signed))
}

impl Module {
    pub(crate) fn build(
        model: &mut Model,
        selected: Vec<usize>,
        presentation: Presentation,
    ) -> Result<Self> {
        let Presentation { generators, rules } = presentation;
        let selected = q_koszul::indices(model, selected)?;
        let mut ranks = BTreeMap::new();
        let mut names = HashMap::new();
        let mut coordinates = Vec::with_capacity(generators.len());
        let mut cells = 1 + 3 * generators.len();
        check_storage(cells)?;
        for (index, (name, degree)) in generators.iter().enumerate() {
            model.spend(1 + name.text.len())?;
            if names.insert(name.text.as_str(), index).is_some() {
                return Err(Error::new(
                    name.span,
                    format!("duplicate module generator '{}'", name.text),
                ));
            }
            let rank = ranks.entry(*degree).or_insert(0);
            coordinates.push(Generator {
                degree: *degree,
                index: *rank,
            });
            *rank += 1;
        }

        let mut actions = HashMap::new();
        let mut active = BTreeSet::new();
        for rule in rules {
            model.spend(1 + rule.targets.len())?;
            if selected.binary_search(&rule.primitive).is_err() {
                return Err(Error::new(
                    rule.span,
                    format!("Q_{} is not declared for this module", rule.primitive),
                ));
            }
            let source = find(&names, &rule.source)?;
            let expected = coordinates[source]
                .degree
                .checked_sub(q_koszul::weight(rule.primitive));
            let mut support = BTreeSet::new();
            for target in rule.targets {
                let index = find(&names, &target)?;
                if Some(coordinates[index].degree) != expected {
                    return Err(Error::new(
                        target.span,
                        format!(
                            "Q_{} on '{}' must drop degree by {}",
                            rule.primitive,
                            rule.source.text,
                            q_koszul::weight(rule.primitive),
                        ),
                    ));
                }
                if !support.remove(&index) {
                    support.insert(index);
                }
            }
            cells += 3 + support.len();
            check_storage(cells)?;
            if !support.is_empty() {
                active.insert(rule.primitive);
            }
            if actions
                .insert((rule.primitive, source), support.into_iter().collect())
                .is_some()
            {
                return Err(Error::new(
                    rule.span,
                    format!(
                        "duplicate Q_{} rule for '{}'",
                        rule.primitive, rule.source.text
                    ),
                ));
            }
        }

        let active: Vec<_> = active.into_iter().collect();
        for (source, (name, _)) in generators.iter().enumerate() {
            for (position, &left) in active.iter().enumerate() {
                model.spend(1)?;
                let left_image = image(&actions, left, source);
                if !apply(model, &actions, left, left_image)?.is_empty() {
                    return Err(Error::new(
                        name.span,
                        format!("Q_{left} squared is nonzero on '{}'", name.text),
                    ));
                }
                for &right in &active[position + 1..] {
                    model.spend(1)?;
                    let right_image = image(&actions, right, source);
                    let left_right = apply(model, &actions, left, right_image)?;
                    let right_left = apply(model, &actions, right, left_image)?;
                    model.spend(left_right.len().min(right_left.len()))?;
                    if left_right != right_left {
                        return Err(Error::new(
                            name.span,
                            format!("Q_{left} and Q_{right} do not commute on '{}'", name.text),
                        ));
                    }
                }
            }
        }

        let mut images = HashMap::new();
        let mut retained = 1 + 2 * ranks.len();
        for ((primitive, source), targets) in actions {
            if targets.is_empty() {
                continue;
            }
            let indices: Vec<_> = targets
                .iter()
                .map(|&target| coordinates[target].index)
                .collect();
            // A rule has one target degree, so declaration order induces
            // increasing indices within that degree.
            debug_assert!(indices.windows(2).all(|pair| pair[0] < pair[1]));
            retained += 4 + indices.len();
            check_storage(retained)?;
            model.spend(indices.len())?;
            images.insert((primitive, coordinates[source]), indices.into_boxed_slice());
        }
        model.allocate(retained)?;
        Ok(Self {
            selected,
            ranks,
            images,
        })
    }

    pub(crate) fn bound(&self) -> Option<usize> {
        self.ranks.keys().next_back().copied()
    }

    pub(crate) fn rank(&self, degree: usize) -> usize {
        self.ranks.get(&degree).copied().unwrap_or(0)
    }

    pub(crate) fn graded_ranks(&self, range: RangeInclusive<usize>) -> Vec<(usize, usize)> {
        if range.is_empty() {
            return Vec::new();
        }
        self.ranks
            .range(range)
            .map(|(&degree, &rank)| (degree, rank))
            .collect()
    }

    pub(crate) fn image(&self, primitive: usize, source: Generator) -> &[usize] {
        self.images
            .get(&(primitive, source))
            .map_or(&[], Box::as_ref)
    }

    pub(crate) fn check_selection(&self, selected: &[usize]) -> Result<()> {
        if let Some(index) = selected
            .iter()
            .find(|index| self.selected.binary_search(index).is_err())
        {
            return Err(Error::plain(format!(
                "Q_{index} is not declared for this module"
            )));
        }
        Ok(())
    }
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(Error::plain("module declaration exceeds the storage limit"));
    }
    Ok(())
}

fn find(names: &HashMap<&str, usize>, label: &Label) -> Result<usize> {
    names.get(label.text.as_str()).copied().ok_or_else(|| {
        Error::new(
            label.span,
            format!("unknown module generator '{}'", label.text),
        )
    })
}

fn image(
    actions: &HashMap<(usize, usize), Vec<usize>>,
    primitive: usize,
    source: usize,
) -> &[usize] {
    actions.get(&(primitive, source)).map_or(&[], Vec::as_slice)
}

fn apply(
    model: &mut Model,
    actions: &HashMap<(usize, usize), Vec<usize>>,
    primitive: usize,
    input: &[usize],
) -> Result<Vec<usize>> {
    let mut result = Vec::new();
    for &source in input {
        model.spend(1)?;
        result = xor(model, &result, image(actions, primitive, source))?;
    }
    Ok(result)
}
