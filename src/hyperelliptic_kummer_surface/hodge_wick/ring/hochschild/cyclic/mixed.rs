use super::{grading_error, Algebra, Flavor, Key, Model, Rational, Result};

pub(super) fn connes_terms(
    model: &mut Model,
    key: &Key,
    power: i128,
    mut emit: impl FnMut(&mut Model, Key, Rational) -> Result<()>,
) -> Result<()> {
    if key.word[0] == 0 {
        return Ok(());
    }
    let arity = key.word.len() - 1;
    for cycle in 0..key.word.len() {
        model.spend(key.word.len())?;
        let mut word = Vec::with_capacity(key.word.len() + 1);
        word.push(0);
        word.extend_from_slice(&key.word[cycle..]);
        word.extend_from_slice(&key.word[..cycle]);
        let negative = !arity.is_multiple_of(2) && !cycle.is_multiple_of(2);
        emit(
            model,
            Key { power, word },
            Rational::integer(if negative { -1 } else { 1 }),
        )?;
    }
    Ok(())
}

pub(super) fn differential_terms(
    model: &mut Model,
    algebra: &Algebra,
    flavor: Flavor,
    key: &Key,
    mut emit: impl FnMut(&mut Model, Key, Rational) -> Result<()>,
) -> Result<()> {
    let arity = key.word.len() - 1;
    if arity > 0 {
        for (i, pair) in key.word.windows(2).enumerate() {
            for &(product, coefficient) in algebra.product(pair[0], pair[1]) {
                model.spend(1 + arity)?;
                let mut word = Vec::with_capacity(arity);
                word.extend_from_slice(&key.word[..i]);
                word.push(product);
                word.extend_from_slice(&key.word[i + 2..]);
                let coefficient = if !i.is_multiple_of(2) {
                    coefficient.negate()
                } else {
                    coefficient
                };
                emit(
                    model,
                    Key {
                        power: key.power,
                        word,
                    },
                    coefficient,
                )?;
            }
        }
        for &(product, coefficient) in algebra.product(key.word[arity], key.word[0]) {
            model.spend(1 + arity)?;
            let mut word = Vec::with_capacity(arity);
            word.push(product);
            word.extend_from_slice(&key.word[1..arity]);
            let coefficient = if !arity.is_multiple_of(2) {
                coefficient.negate()
            } else {
                coefficient
            };
            emit(
                model,
                Key {
                    power: key.power,
                    word,
                },
                coefficient,
            )?;
        }
    }
    if flavor != Flavor::Hochschild && !(flavor == Flavor::Cyclic && key.power == 0) {
        let power = key.power.checked_add(1).ok_or_else(grading_error)?;
        connes_terms(model, key, power, emit)?;
    }
    Ok(())
}
