use std::collections::{BTreeMap, BTreeSet};

use crate::error::{Error, Result};
use crate::lee_linear::Echelon;
use crate::model::Model;

use super::{Auslander, Basis, Nakayama, Path};

pub(super) struct Mesh {
    translate: usize,
    parts: Vec<Part>,
    incoming: usize,
    outgoing: usize,
}

struct Part {
    middle: usize,
    left: Basis,
    right: Basis,
    arrows: [usize; 2],
    sign: i128,
}

impl Mesh {
    pub(super) fn translate(&self) -> usize {
        self.translate
    }

    pub(super) fn render(&self, data: &Auslander) -> String {
        let middle: Vec<_> = self.parts.iter().map(|p| data.key(p.middle)).collect();
        let left: Vec<_> = self.parts.iter().map(|p| (p.arrows[0], p.sign)).collect();
        let right: Vec<_> = self.parts.iter().map(|p| (p.arrows[1], 1)).collect();
        let relation: Vec<_> = self.relation().collect();
        format!("translate={:?}, middle={middle:?}, left={left:?}, right={right:?}, mesh={relation:?}, incomingFactors={}, outgoingFactors={}",data.key(self.translate),self.incoming,self.outgoing)
    }

    fn relation(&self) -> impl Iterator<Item = ([usize; 2], i128)> + '_ {
        self.parts.iter().map(|p| (p.arrows, p.sign))
    }

    pub(super) fn map(&self, side: u128, term: u128) -> Result<(Basis, i128)> {
        let part = usize::try_from(term)
            .ok()
            .and_then(|i| self.parts.get(i))
            .ok_or_else(|| Error::plain("partition almost split middle index is out of range"))?;
        match side {
            0 => Ok((part.left, part.sign)),
            1 => Ok((part.right, 1)),
            _ => Err(Error::plain(
                "partition almost split map side must be 0 or 1",
            )),
        }
    }
}

pub(super) fn all_meshes(
    model: &mut Model,
    algebra: &Nakayama,
    data: &Auslander,
) -> Result<BTreeMap<usize, Mesh>> {
    let mut meshes = BTreeMap::new();
    for (i, top) in data.modules.iter().enumerate() {
        if data.projectives.contains(&i) {
            continue;
        }
        let next = algebra.end(Path {
            start: top.head,
            length: 1,
        });
        let translate = data.module([next as u128, top.length as u128])?;
        let middle = data.module([top.head as u128, top.length as u128 + 1])?;
        let mut parts = vec![part(data, (translate, middle, i), (1, 0), 1)?];
        if top.length > 1 {
            let middle = data.module([next as u128, top.length as u128 - 1])?;
            parts.push(part(data, (translate, middle, i), (0, 1), -1)?);
        }
        check_exactness(model, algebra, data, translate, i, &parts)?;
        let (incoming, outgoing) = check_factors(model, data, translate, i, &parts)?;
        meshes.insert(
            i,
            Mesh {
                translate,
                parts,
                incoming,
                outgoing,
            },
        );
    }
    Ok(meshes)
}

fn part(
    data: &Auslander,
    modules: (usize, usize, usize),
    shifts: (usize, usize),
    sign: i128,
) -> Result<Part> {
    let (source, middle, target) = modules;
    let (left_shift, right_shift) = shifts;
    let left = Basis {
        source,
        target: middle,
        shift: left_shift,
    };
    let right = Basis {
        source: middle,
        target,
        shift: right_shift,
    };
    // Irreducibility and distinct middle summands give both minimality
    // checks modulo the radical of the middle endomorphism ring.
    let arrow = |x| {
        data.arrows
            .binary_search(&x)
            .map_err(|_| Error::plain("partition almost split component is not an irreducible map"))
    };
    Ok(Part {
        middle,
        left,
        right,
        arrows: [arrow(left)?, arrow(right)?],
        sign,
    })
}

fn check_exactness(
    model: &mut Model,
    algebra: &Nakayama,
    data: &Auslander,
    translate: usize,
    target: usize,
    parts: &[Part],
) -> Result<()> {
    let mut composite = BTreeMap::<Basis, i128>::new();
    let mut middle_dimension = 0;
    let mut range = Echelon::default();
    for part in parts {
        if let Some(basis) = data.then(part.left, part.right) {
            *composite.entry(basis).or_default() += part.sign;
        }
        let length = data.modules[part.middle].length;
        middle_dimension += length;
        for t in 0..length {
            if let Some(j) = data.image_layer(algebra, part.right, t) {
                range.insert(model, vec![(j, 1)], Vec::new())?;
            }
        }
    }
    let nontrivial_middle = parts.len() == 2 && parts[0].middle == parts[1].middle;
    if nontrivial_middle || composite.values().any(|&c| c != 0) {
        return Err(Error::plain(
            "partition almost split maps do not form a minimal complex",
        ));
    }
    let mut image = Echelon::default();
    for t in 0..data.modules[translate].length {
        let mut column = Vec::new();
        let mut offset = 0;
        for part in parts {
            if let Some(j) = data.image_layer(algebra, part.left, t) {
                column.push((offset + j, part.sign));
            }
            offset += data.modules[part.middle].length;
        }
        image.insert(model, column, Vec::new())?;
    }
    if image.rank() != data.modules[translate].length
        || range.rank() != data.modules[target].length
        || image.rank() + range.rank() != middle_dimension
    {
        return Err(Error::plain(
            "partition almost split sequence failed the actual image and kernel ranks",
        ));
    }
    Ok(())
}

fn check_factors(
    model: &mut Model,
    data: &Auslander,
    translate: usize,
    target: usize,
    parts: &[Part],
) -> Result<(usize, usize)> {
    let radical = &data.powers[1];
    let mut incoming = 0;
    let mut outgoing = 0;
    for i in 0..data.modules.len() {
        let expected_in: BTreeSet<_> = data
            .basis(i, target)
            .filter(|b| radical.contains(b))
            .collect();
        let expected_out: BTreeSet<_> = data
            .basis(translate, i)
            .filter(|b| radical.contains(b))
            .collect();
        let mut actual_in = BTreeSet::new();
        let mut actual_out = BTreeSet::new();
        for part in parts {
            for basis in data.basis(i, part.middle) {
                model.spend(1)?;
                actual_in.extend(data.then(basis, part.right));
            }
            for basis in data.basis(part.middle, i) {
                model.spend(1)?;
                actual_out.extend(data.then(part.left, basis));
            }
        }
        if expected_in != actual_in || expected_out != actual_out {
            return Err(Error::plain(
                "partition almost split factors do not span the nonsplit map spaces",
            ));
        }
        incoming += expected_in.len();
        outgoing += expected_out.len();
    }
    Ok((incoming, outgoing))
}

struct Word {
    start: usize,
    end: usize,
    letters: Vec<usize>,
    image: Option<Basis>,
}

pub(super) fn certify_presentation(model: &mut Model, data: &Auslander) -> Result<Vec<usize>> {
    let relations: Vec<_> = data
        .meshes
        .values()
        .map(|mesh| mesh.relation().collect::<Vec<_>>())
        .collect();
    let mut words: Vec<_> = (0..data.modules.len())
        .map(|i| Word {
            start: i,
            end: i,
            letters: Vec::new(),
            image: Some(Basis {
                source: i,
                target: i,
                shift: 0,
            }),
        })
        .collect();
    let mut layers = Vec::new();
    for (degree, power) in data.powers.iter().enumerate() {
        let next = data.powers.get(degree + 1);
        let expected: BTreeSet<_> = power
            .iter()
            .copied()
            .filter(|b| next.is_none_or(|n| !n.contains(b)))
            .collect();
        let images: BTreeSet<_> = words
            .iter()
            .filter_map(|p| p.image)
            .filter(|b| next.is_none_or(|n| !n.contains(b)))
            .collect();
        let rank = ideal_rank(model, &words, &relations)?;
        let quotient = words.len() - rank;
        if images != expected || quotient != expected.len() {
            return Err(Error::plain("partition Auslander mesh path quotient does not present its actual endomorphism algebra"));
        }
        layers.push(quotient);
        if quotient == 0 {
            return Ok(layers);
        }
        let mut next_words = Vec::new();
        for word in words {
            for (j, &arrow) in data
                .arrows
                .iter()
                .enumerate()
                .filter(|(_, a)| a.source == word.end)
            {
                model.spend(2 + word.letters.len())?;
                if next_words.len() >= 4096 {
                    return Err(Error::plain(
                        "partition Auslander path layer exceeds the limit of 4096 words",
                    ));
                }
                let mut letters = word.letters.clone();
                letters.push(j);
                next_words.push(Word {
                    start: word.start,
                    end: arrow.target,
                    letters,
                    image: word.image.and_then(|f| data.then(f, arrow)),
                });
            }
        }
        words = next_words;
    }
    Err(Error::plain(
        "partition Auslander mesh presentation did not terminate at the radical bound",
    ))
}

fn ideal_rank(
    model: &mut Model,
    words: &[Word],
    relations: &[Vec<([usize; 2], i128)>],
) -> Result<usize> {
    let indices: BTreeMap<_, _> = words
        .iter()
        .enumerate()
        .map(|(i, w)| ((w.start, w.letters.as_slice()), i))
        .collect();
    let mut span = Echelon::default();
    for word in words {
        for (position, chunk) in word.letters.windows(2).enumerate() {
            model.spend(1 + relations.len())?;
            for relation in relations.iter().filter(|r| r[0].0 == chunk) {
                let mut row = Vec::new();
                let mut replacement = word.letters.clone();
                for &(path, coefficient) in relation {
                    replacement[position..position + 2].copy_from_slice(&path);
                    let index = *indices
                        .get(&(word.start, replacement.as_slice()))
                        .ok_or_else(|| {
                            Error::plain("partition mesh relation substituted a noncomposable path")
                        })?;
                    row.push((index, coefficient));
                }
                row.sort_unstable();
                span.insert(model, row, Vec::new())?;
            }
        }
    }
    Ok(span.rank())
}
