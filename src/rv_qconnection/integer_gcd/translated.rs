use super::{
    arithmetic, certified_gcd, divides, gcd, gcd_images, height_bits, image, primes, reduce,
    remainder, residue, Model, Polynomial, Prime, Reconstruction, Result,
};
use crate::error::Error;

fn translate_image(
    model: &mut Model,
    polynomial: Vec<u32>,
    shift: i128,
    p: Prime,
) -> Result<Vec<u32>> {
    let shift = residue(shift, p);
    if shift == 0 || polynomial.len() <= 1 {
        return Ok(polynomial);
    }
    let (&leading, lower) = polynomial.split_last().unwrap();
    let mut result = Vec::with_capacity(polynomial.len());
    result.push(leading);
    for &coefficient in lower.iter().rev() {
        let size = result.len();
        model.spend(size)?;
        result.push(result[size - 1]);
        for i in (1..size).rev() {
            result[i] = p.add(p.mul(result[i], shift), result[i - 1]);
        }
        result[0] = p.add(p.mul(result[0], shift), coefficient);
    }
    Ok(result)
}

fn mapped(
    model: &mut Model,
    blocks: &[(usize, Polynomial)],
    first: usize,
    origin: usize,
    p: Prime,
) -> Result<Option<Vec<u32>>> {
    model.spend(blocks.len())?;
    if blocks
        .iter()
        .any(|(_, b)| residue(*b.last().unwrap(), p) == 0)
    {
        return Ok(None);
    }
    let (power, block) = &blocks[first];
    let anchor = image(model, block, p)?;
    let mut g = translate_image(model, anchor, origin as i128 - *power as i128, p)?;
    for (index, (power, block)) in blocks.iter().enumerate() {
        if g.len() == 1 {
            break;
        }
        if index == first {
            continue;
        }
        let shift = *power as i128 - origin as i128;
        // P is in the original block variable. Translate only the
        // smaller gcd to divide P there, then move the remainder back.
        model.spend(g.len())?;
        let divisor = translate_image(model, g.clone(), shift, p)?;
        let mut r = image(model, block, p)?;
        remainder(model, &mut r, &divisor, p)?;
        if !r.is_empty() {
            let r = translate_image(model, r, -shift, p)?;
            g = gcd_images(model, g, r, p)?;
        }
    }
    // If no Euclidean step touched the shortest block, it still needs
    // to be normalized before CRT.
    gcd_images(model, g, Vec::new(), p).map(Some)
}

fn ceil_log2(length: usize) -> usize {
    usize::BITS as usize - (length - 1).leading_zeros() as usize
}

fn translated_bits(polynomial: &[i128], shift: i128) -> usize {
    if shift == 0 || polynomial.len() == 1 {
        return height_bits(polynomial);
    }
    // The L1 height of (x+a)^j is exactly (1+|a|)^j.
    let offset = u128::BITS as usize - shift.unsigned_abs().leading_zeros() as usize;
    height_bits(polynomial) + ceil_log2(polynomial.len()) + (polynomial.len() - 1) * offset
}

fn certify(
    model: &mut Model,
    bound: usize,
    mut agrees: impl FnMut(&mut Model, Prime) -> Result<Option<bool>>,
) -> Result<bool> {
    let mut next = primes();
    let mut bits = 0;
    while bits < bound {
        let p = next(model)?;
        match agrees(model, p)? {
            Some(false) => return Ok(false),
            Some(true) => bits += p.get().ilog2() as usize,
            None => {}
        }
    }
    Ok(true)
}

fn divides_translated(
    model: &mut Model,
    block: &[i128],
    factor: &[i128],
    shift: i128,
) -> Result<bool> {
    if shift == 0 {
        return divides(model, block, factor);
    }
    let steps = block.len() + 1 - factor.len();
    // Each pseudo-remainder step scales the magnitude bound by at
    // most twice the shifted divisor height. A zero modular remainder
    // certifies zero pseudo-remainder when the leading term survives.
    let bound = steps
        .checked_mul(1 + translated_bits(factor, shift))
        .and_then(|height| height.checked_add(height_bits(block)))
        .ok_or_else(arithmetic)?;
    certify(model, bound, |model, p| {
        if residue(*factor.last().unwrap(), p) == 0 {
            return Ok(None);
        }
        let divisor = image(model, factor, p)?;
        let divisor = translate_image(model, divisor, shift, p)?;
        let mut block = image(model, block, p)?;
        remainder(model, &mut block, &divisor, p)?;
        Ok(Some(block.is_empty()))
    })
}

pub(crate) fn translated_gcd(
    model: &mut Model,
    blocks: &[(usize, Polynomial)],
    origin: usize,
) -> Result<Polynomial> {
    let (first, (_, block)) = blocks
        .iter()
        .enumerate()
        .min_by_key(|(_, (_, block))| block.len())
        .unwrap();
    let scale = blocks.iter().fold(0, |scale, (_, block)| {
        gcd(scale, block.last().unwrap().unsigned_abs())
    });
    certified_gcd(
        model,
        block.len(),
        scale,
        blocks.len() == 1,
        primes(),
        |model, p| mapped(model, blocks, first, origin, p),
        |model, candidate| {
            for (power, block) in blocks {
                if !divides_translated(model, block, candidate, *power as i128 - origin as i128)? {
                    return Ok(false);
                }
            }
            Ok(true)
        },
    )
}

fn verify_quotient(
    model: &mut Model,
    block: &[i128],
    factor: &[i128],
    shift: i128,
    quotient: &[i128],
) -> Result<bool> {
    // Bound every coefficient of F(x+shift) Q(x) - P(x) directly.
    let products = quotient.len().min(factor.len());
    let height = translated_bits(factor, shift) + height_bits(quotient) + ceil_log2(products);
    let bound = 1 + height.max(height_bits(block));
    certify(model, bound, |model, p| {
        let divisor = image(model, factor, p)?;
        let divisor = translate_image(model, divisor, shift, p)?;
        let quotient = image(model, quotient, p)?;
        let mut r = image(model, block, p)?;
        r.resize(block.len(), 0);
        for (i, &coefficient) in divisor.iter().enumerate() {
            if coefficient == 0 {
                continue;
            }
            model.spend(quotient.len())?;
            for (result, &other) in r[i..].iter_mut().zip(&quotient) {
                *result = p.add(*result, p.negate(p.mul(coefficient, other)));
            }
        }
        Ok(Some(r.iter().all(|&x| x == 0)))
    })
}

pub(crate) fn translated_quotient(
    model: &mut Model,
    block: &[i128],
    factor: &[i128],
    shift: i128,
) -> Result<Polynomial> {
    let size = block.len() + 1 - factor.len();
    let mut lift = Reconstruction::new(size)?;
    let mut next = primes();
    while lift.modulus.bits() <= 128 {
        let p = next(model)?;
        if residue(*factor.last().unwrap(), p) == 0 {
            continue;
        }
        let divisor = image(model, factor, p)?;
        let divisor = translate_image(model, divisor, shift, p)?;
        let mut r = image(model, block, p)?;
        let mut quotient = vec![0; size];
        reduce(model, &mut r, &divisor, p, |i, x| quotient[i] = x)?;
        if !r.is_empty() {
            return Err(Error::plain(
                "quantum differential fraction-free division was not exact",
            ));
        }
        lift.update(model, p, &quotient, 1)?;
    }
    if let Some(quotient) = lift.coefficients(model, 1)? {
        if verify_quotient(model, block, factor, shift, &quotient)? {
            return Ok(quotient);
        }
    }
    // Once M>2^128, a signed 128-bit quotient has unique centered
    // coordinates. Failure of conversion or certification is real.
    Err(arithmetic())
}
