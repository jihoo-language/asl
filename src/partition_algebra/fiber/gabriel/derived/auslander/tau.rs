use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::syntax::ParsedExpr;

use super::{Auslander, Nakayama, Parameter, Path};

mod basis;

const MAX_SUMMANDS: usize = 15;
type Mask = u16;
type Vector = Vec<i128>;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) enum Mode {
    Pair,
    Torsion,
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Grade {
    pub(crate) parameter: Parameter,
    pub(crate) mode: Mode,
}

impl Grade {
    pub(crate) fn render(self) -> String {
        let Parameter { rank, n } = self.parameter;
        let name = match self.mode {
            Mode::Pair => "RVPartitionTauPair",
            Mode::Torsion => "RVPartitionTauTorsion",
        };
        format!("{name}({rank},{n})")
    }
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Record {
    pub(crate) grade: Grade,
    mask: Mask,
}

impl Record {
    pub(crate) fn render(self) -> String {
        format!(
            "{}{:?}",
            self.grade.render(),
            bits(self.mask).collect::<Vec<_>>()
        )
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Literal {
    Empty,
    Module([u128; 2]),
    Shift(u128),
    Face([u128; 2]),
    Chamber(u128),
    Torsion(u128),
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Drop(u128),
    Mutate(u128),
    Fac,
    Recover,
}

impl Unary {
    pub(crate) fn input_mode(self) -> Mode {
        match self {
            Self::Recover => Mode::Torsion,
            _ => Mode::Pair,
        }
    }

    pub(crate) fn output_mode(self) -> Mode {
        match self {
            Self::Fac => Mode::Torsion,
            _ => Mode::Pair,
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Binary {
    Combine,
    Join,
    Meet,
}

impl Binary {
    pub(crate) fn mode(self) -> Mode {
        match self {
            Self::Combine => Mode::Pair,
            Self::Join | Self::Meet => Mode::Torsion,
        }
    }
}

pub(crate) enum Expression {
    Literal(Parameter, Literal),
    Unary(Unary, Box<ParsedExpr>),
    Binary(Binary, Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Copy)]
pub(crate) enum ObjectQuery {
    Summands,
    Vectors,
    Members,
    Contains([u128; 2]),
}

impl ObjectQuery {
    pub(crate) fn mode(self) -> Mode {
        match self {
            Self::Summands | Self::Vectors => Mode::Pair,
            Self::Members | Self::Contains(_) => Mode::Torsion,
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum PairQuery {
    Compatible,
    Inclusion,
    Wall,
}

impl PairQuery {
    pub(crate) fn mode(self) -> Mode {
        match self {
            Self::Compatible | Self::Wall => Mode::Pair,
            Self::Inclusion => Mode::Torsion,
        }
    }
}

pub(crate) enum Report {
    Summary(Parameter),
    Green(Parameter),
    Walls(Parameter),
    Object(ObjectQuery, Box<ParsedExpr>),
    Pair(PairQuery, Box<[ParsedExpr; 2]>),
}

pub(crate) struct Fan {
    parameter: Parameter,
    modules: Vec<[u128; 2]>,
    faces: Vec<Vec<Mask>>,
    vectors: Vec<Vector>,
    generated: Vec<Mask>,
    chambers: Vec<Chamber>,
    torsions: BTreeMap<Mask, usize>,
    mutations: BTreeMap<(usize, usize), usize>,
    walls: Vec<Wall>,
}

struct Chamber {
    face: Mask,
    torsion: Mask,
    c_vectors: Vec<Vector>,
}

struct Wall {
    small: usize,
    large: usize,
    brick: usize,
    normal: Vector,
    small_sign: i8,
}

impl Fan {
    pub(super) fn new(model: &mut Model, algebra: &Nakayama, category: &Auslander) -> Result<Self> {
        let rank = algebra.paths.len();
        let module_count = category.modules.len();
        let candidate_count = module_count + rank;
        if candidate_count > MAX_SUMMANDS {
            return Err(Error::plain(
                "partition tau-tilting candidate count exceeds 15",
            ));
        }
        let modules = (0..module_count).map(|i| category.key(i)).collect();
        let translates: Vec<_> = (0..module_count)
            .map(|i| category.meshes.get(&i).map(super::mesh::Mesh::translate))
            .collect();
        let mut projectives = vec![None; rank];
        for &i in &category.projectives {
            projectives[category.modules[i].head] = Some(i);
        }
        let projectives: Vec<_> = projectives
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| Error::plain("partition tau-tilting is missing a vertex projective"))?;
        let mut compatible = vec![0; candidate_count];
        for (a, row) in compatible.iter_mut().enumerate() {
            for b in 0..candidate_count {
                model.spend(1)?;
                let ok = match (a < module_count, b < module_count) {
                    (true, true) => {
                        translates[b].is_none_or(|t| category.homs[t][a].is_empty())
                            && translates[a].is_none_or(|t| category.homs[t][b].is_empty())
                    }
                    (true, false) => category.homs[a][projectives[b - module_count]].is_empty(),
                    (false, true) => category.homs[b][projectives[a - module_count]].is_empty(),
                    (false, false) => true,
                };
                if ok {
                    *row |= bit(b);
                }
            }
        }
        let faces = enumerate_faces(model, rank, &compatible)?;
        let generated = generated_masks(model, algebra, category)?;
        let mut vectors = Vec::new();
        for state in &category.modules {
            let mut vector = vec![0; rank];
            vector[state.head] = 1;
            let resolution = algebra.resolution(state.head as u128, state.length as u128)?;
            if let Some(syzygy) = resolution.at(1) {
                vector[syzygy.head] -= 1;
            }
            vectors.push(vector);
        }
        for i in 0..rank {
            let mut negative = vec![0; rank];
            negative[i] = -1;
            vectors.push(negative);
        }
        let mut chambers = Vec::new();
        let mut torsions = BTreeMap::new();
        for &face in &faces[rank] {
            let torsion = generated_by(face, &generated);
            let g: Vec<_> = bits(face).map(|k| vectors[k].as_slice()).collect();
            let c_vectors = basis::dual(model, &g)?;
            if torsions.insert(torsion, chambers.len()).is_some() {
                return Err(Error::plain(
                    "distinct partition chambers have the same torsion class",
                ));
            }
            chambers.push(Chamber {
                face,
                torsion,
                c_vectors,
            });
        }
        if faces.iter().flatten().any(|&f| {
            !chambers.iter().any(|c| subset(f, c.face))
                || !torsions.contains_key(&generated_by(f, &generated))
        }) {
            return Err(Error::plain(
                "partition tau-rigid face has no compatible chamber or torsion completion",
            ));
        }
        let mut fan = Self {
            parameter: category.parameter,
            modules,
            faces,
            vectors,
            generated,
            chambers,
            torsions,
            mutations: BTreeMap::new(),
            walls: Vec::new(),
        };
        fan.check_lattice(model)?;
        fan.check_walls(model, algebra, category)?;
        model.allocate(
            8 * fan.faces.iter().map(Vec::len).sum::<usize>() + 12 * fan.walls.len() * (rank + 1),
        )?;
        Ok(fan)
    }

    fn check_lattice(&self, model: &mut Model) -> Result<()> {
        let full = bit(self.modules.len()) - 1;
        if !self.torsions.contains_key(&0) || !self.torsions.contains_key(&full) {
            return Err(Error::plain(
                "partition torsion lattice lacks its top or bottom",
            ));
        }
        for &a in self.torsions.keys() {
            for &b in self.torsions.keys() {
                self.lattice(model, a, b)?;
            }
        }
        Ok(())
    }

    fn lattice(&self, model: &mut Model, a: Mask, b: Mask) -> Result<(Mask, Mask)> {
        let meet = a & b;
        let upper = a | b;
        model.spend(self.torsions.len())?;
        let join = self
            .torsions
            .keys()
            .filter(|&&t| subset(upper, t))
            .fold(bit(self.modules.len()) - 1, |m, &t| m & t);
        if !self.torsions.contains_key(&meet) || !self.torsions.contains_key(&join) {
            return Err(Error::plain(
                "partition torsion classes failed the computed join or meet",
            ));
        }
        Ok((meet, join))
    }

    fn check_walls(
        &mut self,
        model: &mut Model,
        algebra: &Nakayama,
        category: &Auslander,
    ) -> Result<()> {
        let rank = algebra.paths.len();
        let mut mutation_edges = BTreeSet::new();
        for &almost in &self.faces[rank - 1] {
            model.spend(self.chambers.len())?;
            let completions: Vec<_> = self
                .chambers
                .iter()
                .enumerate()
                .filter_map(|(i, c)| subset(almost, c.face).then_some(i))
                .collect();
            let [a, b] = completions.as_slice() else {
                return Err(Error::plain("almost complete partition tau-tilting pair does not have exactly two completions"));
            };
            for (a, b) in [(*a, *b), (*b, *a)] {
                let removed = removed_position(self.chambers[a].face, almost)?;
                if self.mutations.insert((a, removed), b).is_some() {
                    return Err(Error::plain("partition tau-tilting mutation is ambiguous"));
                }
            }
            mutation_edges.insert((*a, *b));
        }
        let mut hasse = BTreeSet::new();
        for (small, lower) in self.chambers.iter().enumerate() {
            for (large, upper) in self.chambers.iter().enumerate() {
                if lower.torsion == upper.torsion || !subset(lower.torsion, upper.torsion) {
                    continue;
                }
                model.spend(self.chambers.len())?;
                if self.chambers.iter().any(|middle| {
                    middle.torsion != lower.torsion
                        && middle.torsion != upper.torsion
                        && subset(lower.torsion, middle.torsion)
                        && subset(middle.torsion, upper.torsion)
                }) {
                    continue;
                }
                hasse.insert((small.min(large), small.max(large)));
                let perpendicular: Vec<_> = bits(upper.torsion)
                    .filter(|&i| {
                        category.homs[i][i].len() == 1
                            && bits(lower.torsion).all(|j| category.homs[i][j].is_empty())
                    })
                    .collect();
                let [brick] = perpendicular.as_slice() else {
                    return Err(Error::plain(
                        "partition torsion cover does not have a unique perpendicular brick",
                    ));
                };
                let module = category.modules[*brick];
                let mut normal = vec![0; rank];
                for length in 0..module.length {
                    normal[algebra.end(Path {
                        start: module.head,
                        length,
                    })] += 1;
                }
                let common = lower.face & upper.face;
                let left = &lower.c_vectors[removed_position(lower.face, common)?];
                let right = &upper.c_vectors[removed_position(upper.face, common)?];
                let negative: Vec<_> = normal.iter().map(|&x| -x).collect();
                let small_sign = match (left, right) {
                    (a, b) if *a == negative && *b == normal => -1,
                    (a, b) if *a == normal && *b == negative => 1,
                    _ => {
                        return Err(Error::plain(
                            "partition wall brick disagrees with the dual g-bases",
                        ))
                    }
                };
                self.walls.push(Wall {
                    small,
                    large,
                    brick: *brick,
                    normal,
                    small_sign,
                });
            }
        }
        if hasse != mutation_edges || self.mutations.len() != rank * self.chambers.len() {
            return Err(Error::plain(
                "partition mutation adjacency does not match the torsion Hasse covers",
            ));
        }
        Ok(())
    }

    fn record(&self, mode: Mode, mask: Mask) -> Record {
        Record {
            grade: Grade {
                parameter: self.parameter,
                mode,
            },
            mask,
        }
    }

    pub(crate) fn count(&self, mode: Mode) -> usize {
        match mode {
            Mode::Pair => self.faces.iter().map(Vec::len).sum(),
            Mode::Torsion => self.torsions.len(),
        }
    }

    pub(crate) fn values(&self, model: &mut Model, mode: Mode) -> Result<Vec<Value>> {
        model.allocate(self.count(mode))?;
        let records: Vec<_> = match mode {
            Mode::Pair => self.faces.iter().flatten().copied().collect(),
            Mode::Torsion => self.torsions.keys().copied().collect(),
        };
        Ok(records
            .into_iter()
            .map(|mask| Value::RVPartitionTau(self.record(mode, mask)))
            .collect())
    }

    fn find_module(&self, key: [u128; 2]) -> Result<usize> {
        self.modules
            .binary_search(&key)
            .map_err(|_| Error::plain("partition tau-tilting module top or length is out of range"))
    }

    fn is_face(&self, mask: Mask) -> bool {
        self.faces
            .get(mask.count_ones() as usize)
            .is_some_and(|row| row.binary_search(&mask).is_ok())
    }

    fn chamber(&self, mask: Mask) -> Result<usize> {
        self.chambers
            .binary_search_by_key(&mask, |c| c.face)
            .map_err(|_| Error::plain("partition operation requires a maximal tau-tilting chamber"))
    }

    pub(crate) fn literal(&self, kind: Literal) -> Result<Record> {
        let at = |i: u128| usize::try_from(i).unwrap_or(usize::MAX);
        let mask = match kind {
            Literal::Empty => 0,
            Literal::Module(key) => {
                let mask = bit(self.find_module(key)?);
                if !self.is_face(mask) {
                    return Err(Error::plain("partition module is not tau-rigid"));
                }
                mask
            }
            Literal::Shift(i) if at(i) < self.vectors[0].len() => bit(self.modules.len() + at(i)),
            Literal::Shift(_) => {
                return Err(Error::plain(
                    "partition shifted projective vertex is out of range",
                ))
            }
            Literal::Face([size, index]) => *self
                .faces
                .get(at(size))
                .and_then(|row| row.get(at(index)))
                .ok_or_else(|| Error::plain("partition tau-tilting face index is out of range"))?,
            Literal::Chamber(i) => {
                self.chambers
                    .get(at(i))
                    .ok_or_else(|| {
                        Error::plain("partition tau-tilting chamber index is out of range")
                    })?
                    .face
            }
            Literal::Torsion(i) => {
                let mask = self
                    .chambers
                    .get(at(i))
                    .ok_or_else(|| Error::plain("partition torsion class index is out of range"))?
                    .torsion;
                return Ok(self.record(Mode::Torsion, mask));
            }
        };
        Ok(self.record(Mode::Pair, mask))
    }

    pub(crate) fn unary(&self, op: Unary, value: Record) -> Result<Record> {
        let position = |i: u128| usize::try_from(i).unwrap_or(usize::MAX);
        match op {
            Unary::Drop(i) => {
                let k = bits(value.mask).nth(position(i)).ok_or_else(|| {
                    Error::plain("partition tau-tilting summand index is out of range")
                })?;
                Ok(self.record(Mode::Pair, value.mask & !bit(k)))
            }
            Unary::Mutate(i) => {
                let chamber = self.chamber(value.mask)?;
                let next = self.mutations.get(&(chamber, position(i))).ok_or_else(|| {
                    Error::plain("partition tau-tilting summand index is out of range")
                })?;
                Ok(self.record(Mode::Pair, self.chambers[*next].face))
            }
            Unary::Fac => Ok(self.record(Mode::Torsion, generated_by(value.mask, &self.generated))),
            Unary::Recover => {
                Ok(self.record(Mode::Pair, self.chambers[self.torsions[&value.mask]].face))
            }
        }
    }

    pub(crate) fn binary(
        &self,
        model: &mut Model,
        op: Binary,
        a: Record,
        b: Record,
    ) -> Result<Record> {
        let mask = match op {
            Binary::Combine => {
                let union = a.mask | b.mask;
                if !self.is_face(union) {
                    return Err(Error::plain(
                        "partition tau-rigid pairs have incompatible summands",
                    ));
                }
                union
            }
            Binary::Join => self.lattice(model, a.mask, b.mask)?.1,
            Binary::Meet => self.lattice(model, a.mask, b.mask)?.0,
        };
        Ok(self.record(op.mode(), mask))
    }

    fn summand(&self, k: usize) -> String {
        if let Some(&[i, l]) = self.modules.get(k) {
            format!("M({i},{l})")
        } else {
            format!("P({})[1]", k - self.modules.len())
        }
    }

    pub(crate) fn object_report(&self, query: ObjectQuery, value: Record) -> Result<String> {
        let mask = value.mask;
        match query {
            ObjectQuery::Summands => Ok(format!(
                "rvPartitionTauSummands = [{}]",
                bits(mask)
                    .map(|k| self.summand(k))
                    .collect::<Vec<_>>()
                    .join(", ")
            )),
            ObjectQuery::Vectors => {
                let g: Vec<_> = bits(mask).map(|k| &self.vectors[k]).collect();
                let mut output = format!("rvPartitionTauVectors = g={g:?}");
                if let Ok(chamber) = self.chamber(mask) {
                    output.push_str(&format!(", c={:?}", self.chambers[chamber].c_vectors));
                }
                Ok(output)
            }
            ObjectQuery::Members => Ok(format!(
                "rvPartitionTauMembers = {:?}",
                bits(mask).map(|k| self.modules[k]).collect::<Vec<_>>()
            )),
            ObjectQuery::Contains(key) => {
                let index = self.find_module(key)?;
                Ok(format!(
                    "rvPartitionTauContains[{},{}] = {}",
                    key[0],
                    key[1],
                    mask & bit(index) != 0
                ))
            }
        }
    }

    pub(crate) fn pair_report(&self, query: PairQuery, a: Record, b: Record) -> Result<String> {
        match query {
            PairQuery::Compatible => Ok(format!(
                "rvPartitionTauCompatible = {}",
                self.is_face(a.mask | b.mask)
            )),
            PairQuery::Inclusion => Ok(format!(
                "rvPartitionTauInclusion = {}",
                subset(a.mask, b.mask)
            )),
            PairQuery::Wall => {
                let a = self.chamber(a.mask)?;
                let b = self.chamber(b.mask)?;
                let wall = self
                    .walls
                    .iter()
                    .find(|w| (w.small, w.large) == (a, b) || (w.large, w.small) == (a, b))
                    .ok_or_else(|| {
                        Error::plain("partition tau-tilting chambers do not share a wall")
                    })?;
                Ok(format!("rvPartitionTauWall = {}", self.wall_string(wall)))
            }
        }
    }

    fn wall_string(&self, w: &Wall) -> String {
        let [i, l] = self.modules[w.brick];
        format!(
            "cover=({}, {}), brick=M({i},{l}), normal={:?}, signs=[{}, {}]",
            w.small, w.large, w.normal, w.small_sign, -w.small_sign
        )
    }

    pub(crate) fn stability_chamber(&self, index: u128) -> Result<(Vec<i128>, Vec<[u128; 2]>)> {
        let chamber = usize::try_from(index)
            .ok()
            .and_then(|i| self.chambers.get(i))
            .ok_or_else(|| Error::plain("partition tau-tilting chamber index is out of range"))?;
        Ok((
            self.interior(chamber.face),
            bits(chamber.torsion).map(|i| self.modules[i]).collect(),
        ))
    }

    pub(crate) fn stability_wall(&self, [a, b]: [u128; 2]) -> Result<(Vec<i128>, [u128; 2])> {
        let edge = self
            .walls
            .iter()
            .find(|w| {
                [w.small as u128, w.large as u128] == [a, b]
                    || [w.large as u128, w.small as u128] == [a, b]
            })
            .ok_or_else(|| Error::plain("partition tau-tilting chambers do not share a wall"))?;
        let common = self.chambers[edge.small].face & self.chambers[edge.large].face;
        Ok((self.interior(common), self.modules[edge.brick]))
    }

    fn interior(&self, face: Mask) -> Vec<i128> {
        let mut sum = vec![0; self.vectors[0].len()];
        for i in bits(face) {
            for (x, y) in sum.iter_mut().zip(&self.vectors[i]) {
                *x += y;
            }
        }
        sum
    }

    pub(crate) fn summary(&self) -> String {
        let Parameter { rank, n } = self.parameter;
        let counts: Vec<_> = self.faces.iter().map(Vec::len).collect();
        format!("rvPartitionTauTilting[{rank},{n}] = faces={counts:?}, chambers={}, torsionClasses={}, covers={}",self.chambers.len(),self.torsions.len(),self.walls.len())
    }

    pub(crate) fn walls(&self) -> String {
        let Parameter { rank, n } = self.parameter;
        format!(
            "rvPartitionTauWalls[{rank},{n}] = [{}]",
            self.walls
                .iter()
                .map(|w| self.wall_string(w))
                .collect::<Vec<_>>()
                .join("; ")
        )
    }

    pub(crate) fn green(&self, model: &mut Model) -> Result<String> {
        let start = self.torsions[&0];
        let top = self.torsions[&(bit(self.modules.len()) - 1)];
        let mut queue = VecDeque::from([start]);
        let mut previous = BTreeMap::from([(start, None)]);
        while let Some(source) = queue.pop_front() {
            if source == top {
                break;
            }
            model.spend(self.walls.len())?;
            for (edge, wall) in self
                .walls
                .iter()
                .enumerate()
                .filter(|(_, w)| w.small == source)
            {
                if let std::collections::btree_map::Entry::Vacant(slot) = previous.entry(wall.large)
                {
                    slot.insert(Some(edge));
                    queue.push_back(wall.large);
                }
            }
        }
        let mut chambers = vec![top];
        let mut bricks = Vec::new();
        let mut current = top;
        while current != start {
            let edge = previous.get(&current).and_then(|p| *p).ok_or_else(|| {
                Error::plain("partition torsion lattice has no saturated bottom-to-top path")
            })?;
            let wall = &self.walls[edge];
            current = wall.small;
            chambers.push(current);
            bricks.push(self.modules[wall.brick]);
        }
        chambers.reverse();
        bricks.reverse();
        let Parameter { rank, n } = self.parameter;
        Ok(format!(
            "rvPartitionTauGreen[{rank},{n}] = length={}, chambers={chambers:?}, bricks={bricks:?}",
            bricks.len()
        ))
    }
}

fn bit(i: usize) -> Mask {
    1 << i
}

fn bits(mask: Mask) -> impl Iterator<Item = usize> {
    (0..Mask::BITS as usize).filter(move |&i| mask & bit(i) != 0)
}

fn subset(a: Mask, b: Mask) -> bool {
    a & !b == 0
}

fn removed_position(face: Mask, common: Mask) -> Result<usize> {
    if (face & !common).count_ones() != 1 {
        return Err(Error::plain(
            "partition chamber wall does not remove exactly one summand",
        ));
    }
    Ok(bits(face).position(|i| common & bit(i) == 0).unwrap())
}

fn enumerate_faces(
    model: &mut Model,
    rank: usize,
    compatibility: &[Mask],
) -> Result<Vec<Vec<Mask>>> {
    let end = bit(compatibility.len());
    let mut faces = vec![Vec::new(); rank + 1];
    faces[0].push(0);
    let mut valid = vec![false; usize::from(end)];
    valid[0] = true;
    for mask in 1..end {
        model.spend(1)?;
        let k = mask.trailing_zeros() as usize;
        let prev = mask & !bit(k);
        if valid[usize::from(prev)] && subset(mask, compatibility[k]) {
            let size = mask.count_ones() as usize;
            if size > rank {
                return Err(Error::plain(
                    "partition tau-rigid face exceeded the simple rank",
                ));
            }
            valid[usize::from(mask)] = true;
            faces[size].push(mask);
        }
    }
    Ok(faces)
}

fn generated_masks(
    model: &mut Model,
    algebra: &Nakayama,
    category: &Auslander,
) -> Result<Vec<Mask>> {
    let mut result = Vec::new();
    for (source, a) in category.modules.iter().enumerate() {
        let mut mask = 0;
        for (target, b) in category.modules.iter().enumerate() {
            let mut images = BTreeSet::new();
            for f in category.basis(source, target) {
                for k in 0..a.length {
                    model.spend(1)?;
                    images.extend(category.image_layer(algebra, f, k));
                }
            }
            let surjective = images.len() == b.length;
            if images.contains(&0) != surjective {
                return Err(Error::plain(
                    "partition uniserial evaluation does not satisfy the head criterion",
                ));
            }
            if surjective {
                mask |= bit(target);
            }
        }
        result.push(mask);
    }
    Ok(result)
}

fn generated_by(face: Mask, masks: &[Mask]) -> Mask {
    bits(face)
        .filter_map(|k| masks.get(k))
        .fold(0, |x, y| x | y)
}
