use std::rc::Rc;

use crate::complex::xor_sorted;
use crate::error::{Error, Result, Span};
use crate::model::{Model, Value};

const WORD_BITS: usize = u64::BITS as usize;

pub(crate) fn support(
    model: &mut Model,
    modulus: usize,
    residues: Vec<(usize, Span)>,
) -> Result<Value> {
    let mut entries = Vec::with_capacity(residues.len());
    for (residue, at) in residues {
        if residue == 0 || residue >= modulus {
            return Err(Error::new(
                at,
                format!("bar residue {residue} must be in 1..{modulus}"),
            ));
        }
        entries.push(residue);
    }
    if !entries.is_empty() {
        model.spend(entries.len() * (1 + entries.len().ilog2() as usize))?;
    }
    entries.sort_unstable();
    let mut odd = Vec::new();
    for residue in entries {
        if odd.last() == Some(&residue) {
            odd.pop();
        } else {
            odd.push(residue);
        }
    }
    model.allocate(odd.len())?;
    Ok(Value::Bar(odd.into()))
}

pub(crate) fn binomial(model: &mut Model, modulus: usize, index: usize, at: Span) -> Result<Value> {
    if index == 0 {
        return Err(Error::new(at, "binom index must be at least 1"));
    }
    let mut support = Vec::new();
    let mut residue = index;
    while residue < modulus {
        model.spend(1)?;
        model.allocate(1)?;
        support.push(residue);
        // Lucas's criterion requires every bit of index in the residue.
        // Incrementing and restoring those bits gives the next such residue.
        residue = (residue + 1) | index;
    }
    Ok(Value::Bar(support.into()))
}

pub(crate) fn add(model: &mut Model, left: &[usize], right: &[usize]) -> Result<Value> {
    model.spend(left.len() + right.len())?;
    let support = xor_sorted(left, right);
    model.allocate(support.len())?;
    Ok(Value::Bar(support.into()))
}

#[derive(Clone, Copy)]
pub(crate) enum Layout {
    Explicit,
    Repeated,
}

pub(crate) struct Cochain {
    pub(crate) support: Rc<[usize]>,
    pub(crate) span: Span,
}

pub(crate) fn evaluate(
    model: &mut Model,
    modulus: usize,
    arity: usize,
    layout: Layout,
    cochains: Vec<Cochain>,
) -> Result<Value> {
    if cochains.iter().all(|cochain| cochain.support.is_empty()) {
        return Ok(Value::CyclicClass(false));
    }

    let words = modulus.div_ceil(WORD_BITS);
    let cells = cochains.len().saturating_mul(words);
    model.allocate(cells)?;
    model.spend(cells)?;
    let mut bits = vec![0; cells];
    for (row, cochain) in bits.chunks_exact_mut(words).zip(&cochains) {
        model.spend(cochain.support.len())?;
        for &residue in cochain.support.iter() {
            row[residue / WORD_BITS] |= 1 << (residue % WORD_BITS);
        }
    }

    let mut offsets = Vec::new();
    if matches!(layout, Layout::Explicit) {
        let mut offset = 0;
        for gap in 1..arity {
            offsets.push(offset);
            offset += arity + 1 - gap;
        }
    }
    let system = DefiningSystem {
        modulus,
        layout,
        words,
        bits,
        offsets,
    };
    match layout {
        Layout::Explicit => {
            for gap in 1..arity {
                for start in 0..=arity - gap {
                    let index = system.row_index(start, gap);
                    system.check_explicit(model, start, gap, cochains[index].span)?;
                }
            }
        }
        Layout::Repeated => system.check_repeated(model, &cochains)?,
    }

    if modulus & 1 != 0 {
        // The periodic resolution has zero H^2 in odd order.
        return Ok(Value::CyclicClass(false));
    }

    // Pair the top obstruction with sum_g [g|1] by summing the first
    // factor over g and evaluating the second factor at 1 in each term.
    let mut class = false;
    model.spend(arity - 1)?;
    for split in 1..arity {
        let left = &cochains[system.row_index(0, split)].support;
        let right = &cochains[system.row_index(split, arity - split)].support;
        class ^= left.len() & 1 != 0 && right.first() == Some(&1);
    }
    Ok(Value::CyclicClass(class))
}

struct DefiningSystem {
    modulus: usize,
    layout: Layout,
    words: usize,
    bits: Vec<u64>,
    offsets: Vec<usize>,
}

impl DefiningSystem {
    fn row_index(&self, start: usize, gap: usize) -> usize {
        match self.layout {
            Layout::Explicit => self.offsets[gap - 1] + start,
            Layout::Repeated => gap - 1,
        }
    }

    fn cochain(&self, start: usize, gap: usize) -> &[u64] {
        let offset = self.row_index(start, gap) * self.words;
        &self.bits[offset..offset + self.words]
    }

    fn check_explicit(&self, model: &mut Model, start: usize, gap: usize, at: Span) -> Result<()> {
        let interval = self.cochain(start, gap);
        for g in 1..self.modulus {
            let constant = has_bit(interval, g);
            for block in 0..self.words {
                model.spend(gap)?;
                let delta = self.differential_block(interval, constant, g, block);
                let mut products = 0;
                for split in 1..gap {
                    if has_bit(self.cochain(start, split), g) {
                        products ^= self.cochain(start + split, gap - split)[block];
                    }
                }
                require_equation(start, gap, g, block, delta, products, at)?;
            }
        }
        Ok(())
    }

    fn check_repeated(&self, model: &mut Model, cochains: &[Cochain]) -> Result<()> {
        // The equations are P_(g+h) = P_g P_h in a truncated polynomial
        // ring. P_0 = 1 and the generator row give P_h = P_1^h and P_1^m = 1.
        let mut active_splits = Vec::new();
        for (index, member) in cochains.iter().enumerate() {
            let gap = index + 1;
            let interval = self.cochain(0, gap);
            let constant = has_bit(interval, 1);
            for block in 0..self.words {
                model.spend(1 + active_splits.len())?;
                let delta = self.differential_block(interval, constant, 1, block);
                let mut products = 0;
                for &split in &active_splits {
                    products ^= self.cochain(0, gap - split)[block];
                }
                require_equation(0, gap, 1, block, delta, products, member.span)?;
            }
            if constant {
                model.allocate(1)?;
                active_splits.push(gap);
            }
        }
        Ok(())
    }

    fn differential_block(&self, interval: &[u64], at_g: bool, g: usize, block: usize) -> u64 {
        let h_start = block * WORD_BITS;
        let width = (self.modulus - h_start).min(WORD_BITS);
        // Use a cyclic slice for f(g+h), including a partial last block.
        let sum_start = if g >= self.modulus - h_start {
            g - (self.modulus - h_start)
        } else {
            g + h_start
        };
        let constant = if at_g { low_mask(width) } else { 0 };
        constant ^ interval[block] ^ cyclic_block(interval, sum_start, width, self.modulus)
    }
}

fn require_equation(
    start: usize,
    gap: usize,
    g: usize,
    block: usize,
    delta: u64,
    products: u64,
    at: Span,
) -> Result<()> {
    let difference = delta ^ products;
    if difference == 0 {
        return Ok(());
    }
    let bit = difference.trailing_zeros() as usize;
    let h = block * WORD_BITS + bit;
    let end = start + gap;
    let left = (delta >> bit) & 1;
    let right = (products >> bit) & 1;
    let reason = if gap == 1 {
        "adjacent input is not a cocycle"
    } else {
        "defining-system equation fails"
    };
    Err(Error::new(
        at,
        format!("{reason} at interval [{start},{end}] for ({g},{h}): delta = {left}, products = {right}"),
    ))
}

fn has_bit(words: &[u64], residue: usize) -> bool {
    words[residue / WORD_BITS] & (1 << (residue % WORD_BITS)) != 0
}

fn low_mask(width: usize) -> u64 {
    u64::MAX >> (WORD_BITS - width)
}

fn contiguous_block(words: &[u64], start: usize, width: usize) -> u64 {
    let block = start / WORD_BITS;
    let offset = start % WORD_BITS;
    let mut bits = words[block] >> offset;
    let remaining = WORD_BITS - offset;
    if width > remaining {
        bits |= words[block + 1] << remaining;
    }
    bits & low_mask(width)
}

fn cyclic_block(words: &[u64], start: usize, width: usize, modulus: usize) -> u64 {
    let end_width = width.min(modulus - start);
    let mut bits = contiguous_block(words, start, end_width);
    if end_width < width {
        bits |= contiguous_block(words, 0, width - end_width) << end_width;
    }
    bits
}
