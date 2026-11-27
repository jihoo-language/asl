use std::collections::{BTreeMap, HashSet, VecDeque};
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::Basis;
use crate::set_family::Diagram;
use crate::sheaf::{self, Matrix, Sheaf};
use crate::sheaf_hyper::{Map, NamedSheaf};
use crate::syntax::Label;

pub(crate) type Permutation = Rc<[usize]>;
pub(crate) type Coefficients = Vec<(usize, i64)>;

pub(crate) enum Query {
    Size,
    Characters,
    Element(usize),
    Stabilizer(usize),
    Maps,
    Orbitals,
}

#[derive(Clone, Eq, Hash, PartialEq)]
struct Pair {
    source: Permutation,
    target: Permutation,
}

impl Pair {
    fn identity(source: usize, target: usize) -> Self {
        Self {
            source: (0..source).collect(),
            target: (0..target).collect(),
        }
    }

    fn size(&self) -> usize {
        self.source.len() + self.target.len()
    }

    fn after(&self, model: &mut Model, before: &Self) -> Result<Self> {
        model.spend(1 + self.size())?;
        Ok(Self {
            source: before
                .source
                .iter()
                .map(|&image| self.source[image])
                .collect(),
            target: before
                .target
                .iter()
                .map(|&image| self.target[image])
                .collect(),
        })
    }

    fn inverse(&self, model: &mut Model) -> Result<Self> {
        model.spend(1 + self.size())?;
        Ok(Self {
            source: inverse(&self.source).expect("certified permutation"),
            target: inverse(&self.target).expect("certified permutation"),
        })
    }

    fn fixed(&self) -> (usize, usize) {
        (fixed(&self.source), fixed(&self.target))
    }
}

struct Transport {
    forward: Pair,
    back: Pair,
}

struct Orbits {
    labels: Vec<usize>,
    representatives: Vec<usize>,
    sizes: Vec<usize>,
}

impl Orbits {
    fn build(
        model: &mut Model,
        size: usize,
        retained: usize,
        generators: &[Pair],
        action: impl Fn(&Pair, usize) -> usize,
    ) -> Result<Self> {
        check_cells(retained.saturating_add(4usize.saturating_mul(size)))?;
        let mut orbits = Self {
            labels: vec![usize::MAX; size],
            representatives: Vec::new(),
            sizes: Vec::new(),
        };
        let mut queue = VecDeque::new();
        for i in 0..size {
            model.spend(1)?;
            if orbits.labels[i] != usize::MAX {
                continue;
            }
            let orbit = orbits.representatives.len();
            orbits.representatives.push(i);
            orbits.labels[i] = orbit;
            queue.push_back(i);
            let mut count = 0;
            while let Some(point) = queue.pop_front() {
                count += 1;
                model.spend(1 + generators.len())?;
                for generator in generators {
                    let next = action(generator, point);
                    if orbits.labels[next] == usize::MAX {
                        orbits.labels[next] = orbit;
                        queue.push_back(next);
                    }
                }
            }
            orbits.sizes.push(count);
        }
        Ok(orbits)
    }
}

pub(crate) struct Monodromy {
    pub(crate) source: Rc<str>,
    pub(crate) target: Rc<str>,
    group: Vec<Pair>,
    transports: Vec<Option<Transport>>,
    source_orbits: Vec<usize>,
    orbitals: Orbits,
}

impl Monodromy {
    pub(crate) fn build(
        model: &mut Model,
        source: &Label,
        target: &Label,
        basepoint: &Label,
    ) -> Result<Self> {
        let left = model.set_family(source)?;
        let right = model.set_family(target)?;
        if left.base != right.base {
            return Err(Error::plain(
                "joint monodromy requires set families on the same named site",
            ));
        }
        let root = model.point(left.base, basepoint)?;
        let site = model.sheaf_site_data(left.base)?;
        let covers = site.covers();
        let count = covers.len();
        let mut reverse = vec![Vec::new(); count];
        let family = LoopFamilies {
            left: &left,
            right: &right,
            source,
            target,
        };
        for (x, row) in covers.iter().enumerate() {
            model.spend(1 + row.len())?;
            for &y in row {
                reverse[y].push(x);
            }
        }
        let n = left.size(root);
        let m = right.size(root);
        let identity = Pair::identity(n, m);
        let weight = 4 + n + m;
        let mut stored = 8 + 4 * count + weight;
        let mut transports: Vec<Option<Transport>> = (0..count).map(|_| None).collect();
        transports[root] = Some(Transport {
            forward: identity.clone(),
            back: identity.clone(),
        });
        let mut tree = vec![None; count];
        let mut queue = VecDeque::from([root]);
        while let Some(parent) = queue.pop_front() {
            model.spend(1 + covers[parent].len() + reverse[parent].len())?;
            let neighbors = covers[parent]
                .iter()
                .map(|&x| (x, true))
                .chain(reverse[parent].iter().map(|&x| (x, false)));
            for (point, forward) in neighbors {
                if transports[point].is_some() {
                    continue;
                }
                let (x, y) = if forward {
                    (parent, point)
                } else {
                    (point, parent)
                };
                let edge = family.edge(model, x, y)?;
                let edge = if forward { edge } else { edge.inverse(model)? };
                let forward_map =
                    edge.after(model, &transports[parent].as_ref().unwrap().forward)?;
                let back = forward_map.inverse(model)?;
                stored += 2 * weight;
                check_cells(stored)?;
                transports[point] = Some(Transport {
                    forward: forward_map,
                    back,
                });
                tree[point] = Some((parent, forward));
                queue.push_back(point);
            }
        }
        drop(reverse);
        let mut generators = Vec::new();
        let mut unique = HashSet::new();
        for (x, row) in covers.iter().enumerate() {
            let Some(start) = &transports[x] else {
                continue;
            };
            model.spend(1 + row.len())?;
            for &y in row {
                if tree[y] == Some((x, true)) || tree[x] == Some((y, false)) {
                    continue;
                }
                let end = transports[y].as_ref().expect("same Hasse component");
                let edge = family.edge(model, x, y)?;
                let along = edge.after(model, &start.forward)?;
                let loop_action = end.back.after(model, &along)?;
                if loop_action != identity && unique.insert(loop_action.clone()) {
                    check_cells(stored + weight * (generators.len() + 1))?;
                    generators.push(loop_action);
                }
            }
        }
        drop(unique);
        let group = joint_image(
            model,
            identity,
            &generators,
            stored + weight * generators.len(),
        )?;
        stored += group.len() * weight;
        let retained = stored + weight * generators.len();
        let source_orbits =
            Orbits::build(model, n, retained, &generators, |pair, x| pair.source[x])?
                .representatives;
        let orbitals = Orbits::build(
            model,
            n * m,
            retained + source_orbits.len(),
            &generators,
            |pair, x| pair.target[x / n] * n + pair.source[x % n],
        )?;
        let orbital_cells = orbitals.labels.len() + 2 * orbitals.sizes.len();
        stored += orbital_cells + source_orbits.len();
        check_cells(stored)?;
        model.allocate(stored)?;
        Ok(Self {
            source: source.text.as_str().into(),
            target: target.text.as_str().into(),
            group,
            transports,
            source_orbits,
            orbitals,
        })
    }

    pub(crate) fn size(&self) -> usize {
        self.group.len()
    }

    pub(crate) fn element(&self, index: usize) -> Result<(&[usize], &[usize])> {
        let pair = self.group.get(index).ok_or_else(|| {
            Error::plain(format!(
                "joint monodromy element {index} is outside the image"
            ))
        })?;
        Ok((&pair.source, &pair.target))
    }

    pub(crate) fn characters(&self, model: &mut Model) -> Result<BTreeMap<(usize, usize), usize>> {
        let mut counts = BTreeMap::new();
        for pair in &self.group {
            model.spend(1 + pair.size())?;
            let position = pair.fixed();
            if !counts.contains_key(&position) {
                check_cells(5 * (counts.len() + 1))?;
            }
            *counts.entry(position).or_insert(0) += 1;
        }
        Ok(counts)
    }

    pub(crate) fn stabilizer(
        &self,
        model: &mut Model,
        index: usize,
    ) -> Result<(usize, usize, usize)> {
        let basis = &self.group[0];
        if index >= basis.source.len() {
            return Err(Error::plain(format!(
                "point stabilizer index {index} is outside the source fiber"
            )));
        }
        model.spend(basis.size() + self.group.len())?;
        let mut source_fixed = vec![true; basis.source.len()];
        let mut target_fixed = vec![true; basis.target.len()];
        let mut order = 0;
        for pair in &self.group {
            if pair.source[index] == index {
                order += 1;
                model.spend(pair.size())?;
                common_fixed(&mut source_fixed, &pair.source);
                common_fixed(&mut target_fixed, &pair.target);
            }
        }
        Ok((
            order,
            source_fixed.into_iter().filter(|&x| x).count(),
            target_fixed.into_iter().filter(|&x| x).count(),
        ))
    }

    pub(crate) fn natural_count(&self, model: &mut Model) -> Result<u128> {
        let mut counts = Vec::with_capacity(self.source_orbits.len());
        for &point in &self.source_orbits {
            let (_, _, target) = self.stabilizer(model, point)?;
            if target == 0 {
                return Ok(0);
            }
            counts.push(target);
        }
        counts.into_iter().try_fold(1u128, |out, count| {
            out.checked_mul(count as u128)
                .ok_or_else(|| Error::plain("equivariant map count exceeds u128"))
        })
    }

    pub(crate) fn orbital_sizes(&self) -> &[usize] {
        &self.orbitals.sizes
    }

    pub(crate) fn orbital_map(
        &self,
        model: &mut Model,
        source: NamedSheaf,
        target: NamedSheaf,
        coefficients: Coefficients,
    ) -> Result<Map> {
        if source.sheaf.prime != target.sheaf.prime {
            return Err(Error::plain(
                "orbital sheaf maps require the same prime field",
            ));
        }
        let prime = source.sheaf.prime;
        let mut weights = vec![0; self.orbitals.sizes.len()];
        for (index, value) in coefficients {
            model.spend(1)?;
            let Some(coefficient) = weights.get_mut(index) else {
                return Err(Error::plain(format!(
                    "orbital index {index} is outside the product action"
                )));
            };
            *coefficient = prime.add(*coefficient, sheaf::residue(prime, value));
        }
        let n = self.group[0].source.len();
        let zero = Rc::new(Matrix::default());
        let mut matrices = Vec::with_capacity(self.transports.len());
        let mut stored = self.transports.len();
        for transport in &self.transports {
            let Some(at) = transport else {
                matrices.push(zero.clone());
                continue;
            };
            let matrix = Matrix::from_action(model, n, |model, column| {
                model.spend(1 + at.back.target.len())?;
                let point = at.back.source[column];
                let entries = at
                    .back
                    .target
                    .iter()
                    .enumerate()
                    .filter_map(|(row, &covector)| {
                        let orbit = self.orbitals.labels[covector * n + point];
                        let coefficient = weights[orbit];
                        (coefficient != 0).then_some((row, coefficient))
                    })
                    .collect();
                Ok(entries)
            })?;
            stored += matrix.cells();
            check_cells(stored)?;
            matrices.push(Rc::new(matrix));
        }
        Map::from_matrices(model, source, target, matrices)
    }
}

struct LoopFamilies<'a> {
    left: &'a Diagram,
    right: &'a Diagram,
    source: &'a Label,
    target: &'a Label,
}

impl LoopFamilies<'_> {
    fn edge(&self, model: &mut Model, x: usize, y: usize) -> Result<Pair> {
        Ok(Pair {
            source: bijection(model, self.left, self.source, x, y)?,
            target: bijection(model, self.right, self.target, x, y)?,
        })
    }
}

fn bijection(
    model: &mut Model,
    family: &Diagram,
    name: &Label,
    x: usize,
    y: usize,
) -> Result<Permutation> {
    let size = family.size(x);
    let invalid = || {
        Error::plain(format!(
            "cover '{} < {}' is not bijective in set family '{}'",
            family.site.names[x], family.site.names[y], name.text
        ))
    };
    if size != family.size(y) {
        return Err(invalid());
    }
    model.spend(1 + 2 * size)?;
    let images: Permutation = (0..size).map(|i| family.apply(x, y, i)).collect();
    let mut used = vec![false; size];
    for &point in images.iter() {
        if std::mem::replace(&mut used[point], true) {
            return Err(invalid());
        }
    }
    Ok(images)
}

fn inverse(images: &[usize]) -> Option<Permutation> {
    let mut rows = vec![0; images.len()];
    for (source, &target) in images.iter().enumerate() {
        *rows.get_mut(target)? = source;
    }
    Some(rows.into())
}

fn joint_image(
    model: &mut Model,
    identity: Pair,
    generators: &[Pair],
    retained: usize,
) -> Result<Vec<Pair>> {
    let weight = 10 + identity.size();
    let mut seen = HashSet::from([identity.clone()]);
    let mut group = vec![identity];
    let mut index = 0;
    while index < group.len() {
        model.spend(1)?;
        for generator in generators {
            let candidate = generator.after(model, &group[index])?;
            model.spend(1 + candidate.size())?;
            if seen.insert(candidate.clone()) {
                check_cells(retained + (group.len() + 1) * weight)?;
                group.push(candidate);
            }
        }
        index += 1;
    }
    Ok(group)
}

fn fixed(images: &[usize]) -> usize {
    images.iter().enumerate().filter(|&(i, &j)| i == j).count()
}

fn common_fixed(common: &mut [bool], images: &[usize]) {
    for (i, (present, &image)) in common.iter_mut().zip(images).enumerate() {
        *present &= i == image;
    }
}

pub(crate) fn linearize(model: &mut Model, family: &Diagram, prime: Prime) -> Result<Sheaf> {
    let dimensions = (0..family.site.names.len())
        .map(|x| family.size(x))
        .collect();
    Sheaf::from_cover_maps(model, family.base, prime, dimensions, |model, x, y| {
        let matrix = Matrix::from_action(model, family.size(x), |_, i| {
            Ok(vec![(family.apply(x, y, i), 1)])
        })?;
        Ok(Rc::new(matrix))
    })
}

pub(crate) fn map_rank(model: &mut Model, map: &Map, point: usize) -> Result<usize> {
    let mut basis = Basis::default();
    for (_, column) in map.matrix(point).columns() {
        model.spend(1 + column.len())?;
        basis.insert(model, map.source.sheaf.prime, column.to_vec())?;
    }
    Ok(basis.rank())
}

fn check_cells(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        Err(Error::plain("joint monodromy exceeds its storage limit"))
    } else {
        Ok(())
    }
}
