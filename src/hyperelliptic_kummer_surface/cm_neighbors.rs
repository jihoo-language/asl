use std::cmp::Ordering;

use super::cm::{exact_divide, Form};
use super::cm_ring::{Class, Order};
use super::rational_sparse::{decimal, Integer, Magnitude, WORK_LIMBS};
use crate::error::{Error, Result};
use crate::model::{Model, MAX_SHOWN_BYTES};
use crate::odd_cobar::Prime;
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Statement {
    Address(usize, usize, Box<ParsedExpr>),
    Neighbors(usize, Box<ParsedExpr>),
    Counts(usize, Box<ParsedExpr>),
}

pub(crate) struct Expression(
    pub(crate) Label,
    pub(crate) usize,
    pub(crate) usize,
    pub(crate) Box<ParsedExpr>,
);

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
enum Direction {
    Ascending,
    Horizontal,
    Descending,
}

impl Direction {
    fn name(self) -> &'static str {
        match self {
            Self::Ascending => "ascending",
            Self::Horizontal => "horizontal",
            Self::Descending => "descending",
        }
    }
}

#[derive(Default)]
struct Counts {
    ascending: usize,
    horizontal: usize,
    descending: usize,
}

impl Counts {
    fn add(&mut self, direction: Direction) {
        match direction {
            Direction::Ascending => self.ascending += 1,
            Direction::Horizontal => self.horizontal += 1,
            Direction::Descending => self.descending += 1,
        }
    }

    fn report(&self) -> String {
        format!(
            "ascending={}, horizontal={}, descending={}",
            self.ascending, self.horizontal, self.descending
        )
    }
}

struct Neighbor {
    form: Form,
    discriminant: Magnitude,
    direction: Direction,
    content: u128,
}

pub(crate) fn prime(model: &mut Model, ell: usize) -> Result<Prime> {
    if ell == 2 {
        return Ok(Prime::validated(2));
    }
    if ell < 2 || ell % 2 == 0 {
        return Err(Error::plain("CM isogeny degree must be prime"));
    }
    Prime::new(model, ell)
}

pub(crate) fn check_address(ell: Prime, address: usize) -> Result<u32> {
    if address > ell.get() as usize {
        return Err(Error::plain(
            "CM neighbor address must be between 0 and ell inclusive",
        ));
    }
    Ok(address as u32)
}

pub(crate) fn target(
    model: &mut Model,
    source: Order,
    class: Class,
    ell: Prime,
    address: u32,
    target: Order,
) -> Result<Class> {
    let neighbor = Neighbor::new(model, source, class, ell, address)?;
    if target.magnitude() != neighbor.discriminant {
        return Err(Error::plain(format!(
            "CM neighbor discriminant -{} differs from the named target order",
            decimal(neighbor.discriminant),
        )));
    }
    Class::from_reduced(neighbor.form)
}

pub(crate) fn address_report(
    model: &mut Model,
    source: Order,
    class: Class,
    ell: Prime,
    address: u32,
) -> Result<String> {
    let neighbor = Neighbor::new(model, source, class, ell, address)?;
    let ell = ell.get();
    let label = if address == ell {
        "infinity".into()
    } else {
        address.to_string()
    };
    Ok(format!(
        "rvCMNeighborData[{ell},{address}] = address={label}, {}, content={}, multiplicity=1",
        neighbor.report(),
        neighbor.content,
    ))
}

pub(crate) fn neighbors_report(
    model: &mut Model,
    source: Order,
    class: Class,
    ell: Prime,
    with_forms: bool,
) -> Result<String> {
    let degree = ell.get();
    let count = degree as usize + 1;
    let sort_work = if with_forms {
        8 * WORK_LIMBS * (1 + count.ilog2() as usize)
    } else {
        0
    };
    model
        .spend(count.saturating_mul(8 * WORK_LIMBS + sort_work))
        .map_err(|_| Error::plain("CM neighbor enumeration exceeds the computation limit"))?;
    let mut counts = Counts::default();
    if !with_forms {
        for address in 0..=degree {
            let (_, _, direction) = restricted_form(model, class, ell, address)?;
            counts.add(direction);
        }
        return Ok(format!(
            "rvCMNeighborCounts[{degree}] = {}",
            counts.report()
        ));
    }
    let mut neighbors = Vec::with_capacity(count);
    for address in 0..=degree {
        let neighbor = Neighbor::new(model, source, class, ell, address)?;
        counts.add(neighbor.direction);
        neighbors.push(neighbor);
    }
    neighbors.sort_unstable_by(Neighbor::compare);
    let mut output = format!("rvCMNeighbors[{degree}] = {}, targets=[", counts.report());
    let mut rest = neighbors.as_slice();
    let mut separator = "";
    while let Some((first, _)) = rest.split_first() {
        let multiplicity = rest
            .iter()
            .take_while(|neighbor| first.compare(neighbor).is_eq())
            .count();
        let entry = format!(
            "{separator}{{{}, multiplicity={multiplicity}}}",
            first.report()
        );
        if output.len() + entry.len() + 1 > MAX_SHOWN_BYTES {
            return Err(Error::plain("CM neighbor report exceeds the display limit"));
        }
        output.push_str(&entry);
        separator = ", ";
        rest = &rest[multiplicity..];
    }
    output.push(']');
    Ok(output)
}

impl Neighbor {
    fn new(
        model: &mut Model,
        source: Order,
        class: Class,
        ell: Prime,
        address: u32,
    ) -> Result<Self> {
        let (mut form, content, direction) = restricted_form(model, class, ell, address)?;
        let l = u128::from(ell.get());
        let d = source.magnitude();
        let discriminant = match direction {
            Direction::Descending => d.multiply_u128(l * l).ok_or_else(arithmetic)?,
            Direction::Horizontal => d,
            Direction::Ascending => d
                .divide_exact(model, Magnitude::from_small(l * l))?
                .ok_or_else(arithmetic)?,
        };
        let divisor = Magnitude::from_small(content);
        form.a = form
            .a
            .divide_exact(model, divisor)?
            .ok_or_else(arithmetic)?;
        form.b = exact_divide(model, form.b, divisor)?;
        form.c = form
            .c
            .divide_exact(model, divisor)?
            .ok_or_else(arithmetic)?;
        form.reduce(model, |_, _| Ok(()))?;
        if form.discriminant()? != Some(discriminant) {
            return Err(Error::plain("CM neighbor has an inconsistent discriminant"));
        }
        Ok(Self {
            form,
            discriminant,
            direction,
            content,
        })
    }

    fn compare(&self, other: &Self) -> Ordering {
        self.discriminant
            .cmp(other.discriminant)
            .then_with(|| self.form.a.cmp(other.form.a))
            .then_with(|| match (self.form.b.negative, other.form.b.negative) {
                (true, false) => Ordering::Less,
                (false, true) => Ordering::Greater,
                (true, true) => other.form.b.magnitude.cmp(self.form.b.magnitude),
                (false, false) => self.form.b.magnitude.cmp(other.form.b.magnitude),
            })
            .then_with(|| self.form.c.cmp(other.form.c))
    }

    fn report(&self) -> String {
        let Form { a, b, c } = &self.form;
        let sign = if b.negative { "-" } else { "" };
        format!(
            "discriminant=-{}, form=[{}, {sign}{}, {}], direction={}",
            decimal(self.discriminant),
            decimal(*a),
            decimal(b.magnitude),
            decimal(*c),
            self.direction.name(),
        )
    }
}

fn restricted_form(
    model: &mut Model,
    Class([a, b, c]): Class,
    ell: Prime,
    address: u32,
) -> Result<(Form, u128, Direction)> {
    model.spend(8 * WORK_LIMBS)?;
    let l = i128::from(ell.get());
    let form = if address == ell.get() {
        Form {
            a: Magnitude::from_small(a as u128),
            b: Integer::product(l, b),
            c: Integer::product(l * l, c).magnitude,
        }
    } else {
        let s = i128::from(address);
        Form {
            a: Integer::product(l * l, a).magnitude,
            b: Integer::product(2 * a, s)
                .add(Integer::from_i128(b))?
                .multiply(l)?,
            c: Integer::product(a, s * s)
                .add(Integer::product(b, s))?
                .add(Integer::from_i128(c))?
                .magnitude,
        }
    };
    let content = form.content(model)?.small().ok_or_else(arithmetic)?;
    let l = l as u128;
    let direction = if content == 1 {
        Direction::Descending
    } else if content == l {
        Direction::Horizontal
    } else if content == l * l {
        Direction::Ascending
    } else {
        return Err(Error::plain("primitive CM neighbor has unexpected content"));
    };
    Ok((form, content, direction))
}

fn arithmetic() -> Error {
    Error::plain("CM neighbor exceeds the supported exact arithmetic range")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kronecker_counts_and_dual_edges_agree_for_all_sample_classes(
    ) -> std::result::Result<(), String> {
        verify_counts_and_duals().map_err(|error| error.to_string())
    }

    fn verify_counts_and_duals() -> Result<()> {
        let mut model = Model::new();
        for (d, ell, surface) in [
            (-23, 2, true),
            (-23, 3, true),
            (-23, 5, true),
            (-23, 23, true),
            (-3, 2, true),
            (-3, 3, true),
            (-3, 7, true),
            (-4, 2, true),
            (-4, 3, true),
            (-4, 5, true),
            (-12, 2, false),
            (-16, 2, false),
            (-27, 3, false),
            (-36, 3, false),
            (-92, 2, false),
            (-368, 2, false),
        ] {
            model.reset_work();
            let prime = prime(&mut model, ell)?;
            let order = Order::new(d)?;
            for class in order.elements(&mut model, 100)? {
                model.reset_work();
                let mut counts = Counts::default();
                for address in 0..=prime.get() {
                    let neighbor = Neighbor::new(&mut model, order, class, prime, address)?;
                    counts.add(neighbor.direction);
                    let target_d = neighbor
                        .discriminant
                        .signed(true)
                        .expect("small discriminant");
                    let target_order = Order::new(target_d)?;
                    let target = Class::from_reduced(neighbor.form)?;
                    assert!(target_order.elements(&mut model, 100)?.contains(&target));
                    let reverse_direction = match neighbor.direction {
                        Direction::Ascending => Direction::Descending,
                        Direction::Horizontal => Direction::Horizontal,
                        Direction::Descending => Direction::Ascending,
                    };
                    let mut found_dual = false;
                    for reverse in 0..=prime.get() {
                        let dual = Neighbor::new(&mut model, target_order, target, prime, reverse)?;
                        if dual.discriminant == order.magnitude()
                            && dual.direction == reverse_direction
                            && Class::from_reduced(dual.form)? == class
                        {
                            found_dual = true;
                            break;
                        }
                    }
                    assert!(
                        found_dual,
                        "discriminant {d}, class {:?}, degree {ell}, address {address}",
                        class.0
                    );
                }
                let expected_horizontal = if surface {
                    (1 + kronecker(d, ell as i128)) as usize
                } else {
                    0
                };
                let expected_ascending = usize::from(!surface);
                assert_eq!(counts.horizontal, expected_horizontal);
                assert_eq!(counts.ascending, expected_ascending);
                assert_eq!(
                    counts.descending,
                    ell + 1 - expected_ascending - expected_horizontal
                );
            }
        }
        Ok(())
    }

    fn kronecker(d: i128, ell: i128) -> i8 {
        if d % ell == 0 {
            return 0;
        }
        if ell == 2 {
            return if matches!(d.rem_euclid(8), 1 | 7) {
                1
            } else {
                -1
            };
        }
        if (1..ell).any(|i| i * i % ell == d.rem_euclid(ell)) {
            1
        } else {
            -1
        }
    }
}
