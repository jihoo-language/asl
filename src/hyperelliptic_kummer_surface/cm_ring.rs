use super::cm::{self, bezout_projection, exact_divide, positive, residue_mod, Form};
use super::rational_sparse::{Integer, Magnitude, WORK_LIMBS};
use crate::error::{Error, Result};
use crate::model::Model;
use crate::syntax::{Label, ParsedExpr};

pub(crate) const GROUP_LIMIT: usize = 4096;

pub(crate) enum Statement {
    Define(Label, i128),
    Order(Box<ParsedExpr>),
    Group(Label),
    SameGenus(Box<[ParsedExpr; 2]>),
}

pub(crate) enum Expression {
    Literal(Label, [i128; 3]),
    Principal(Label),
    Period(Label, Label),
    Compose(Box<[ParsedExpr; 2]>),
    Inverse(Box<ParsedExpr>),
    Power(i128, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) struct Order {
    discriminant: i128,
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Class(pub(crate) [i128; 3]);

impl Order {
    pub(crate) fn new(discriminant: i128) -> Result<Self> {
        if discriminant >= 0 || !matches!(discriminant.rem_euclid(4), 0 | 1) {
            return Err(Error::plain(
                "CM order discriminant must be negative and congruent to 0 or 1 modulo 4",
            ));
        }
        Ok(Self { discriminant })
    }

    pub(super) fn magnitude(self) -> Magnitude {
        Magnitude::from_small(self.discriminant.unsigned_abs())
    }

    fn parity(self) -> i128 {
        self.discriminant.rem_euclid(2)
    }

    pub(crate) fn cardinality_bound(self) -> usize {
        let max_a = (self.discriminant.unsigned_abs() / 3).isqrt();
        usize::try_from(max_a * (max_a + 1) / 2).unwrap_or(usize::MAX)
    }

    pub(crate) fn principal(self) -> Class {
        let b0 = self.parity();
        Class([1, b0, -((self.discriminant - b0) / 4)])
    }

    pub(crate) fn class(self, model: &mut Model, [a, b, c]: [i128; 3]) -> Result<Class> {
        if a <= 0 || c <= 0 {
            return Err(Error::plain("CM class requires a positive definite form"));
        }
        let form = Form {
            a: Magnitude::from_small(a as u128),
            b: Integer::from_i128(b),
            c: Magnitude::from_small(c as u128),
        };
        self.canonical(model, form)
    }

    pub(crate) fn period_class(self, model: &mut Model, period: &cm::Reduced) -> Result<Class> {
        let form = period.primitive_form(model, self.magnitude())?;
        Class::from_reduced(form)
    }

    pub(crate) fn inverse(self, model: &mut Model, Class([a, b, c]): Class) -> Result<Class> {
        self.class(model, [a, -b, c])
    }

    pub(crate) fn compose(self, model: &mut Model, left: Class, right: Class) -> Result<Class> {
        let Class([a, b, _]) = left;
        let Class([e, f, _]) = right;
        if a == 1 {
            return Ok(right);
        }
        if e == 1 {
            return Ok(left);
        }
        model.spend(20 * WORK_LIMBS)?;
        let b0 = self.parity();
        let k = (b + b0) / 2;
        let l = (f + b0) / 2;
        let omega_constant = (self.discriminant - b0) / 4;
        // Products of the generators a, omega-k and e, omega-l.
        let mut generators = [
            (Integer::product(a, e), Integer::from_i128(0)),
            (Integer::product(a, -l), Integer::from_i128(a)),
            (Integer::product(e, -k), Integer::from_i128(e)),
            (
                Integer::product(k, l).add(Integer::from_i128(omega_constant))?,
                Integer::from_i128(b0 - k - l),
            ),
        ];
        let mut content = Magnitude::default();
        for (x, y) in generators {
            content = content.gcd(model, x.magnitude)?.gcd(model, y.magnitude)?;
        }
        for (x, y) in &mut generators {
            *x = exact_divide(model, *x, content)?;
            *y = exact_divide(model, *y, content)?;
        }
        // The first generator is a constant, so working modulo it
        // gives an equivalent choice of R in the generator R+omega.
        let modulus = generators[0].0.magnitude;
        let r = bezout_projection(model, generators.map(|(x, y)| (y, x)), modulus)?.ok_or_else(
            || Error::plain("primitive CM ideal must project onto the omega coordinate"),
        )?;
        let mut a = Magnitude::default();
        for (x, y) in generators {
            let constant = x.add(positive(r).multiply_wide(y)?.negate())?;
            a = a.gcd(model, constant.magnitude)?;
        }
        let two_a = a.multiply_u128(2).ok_or_else(arithmetic)?;
        let residue = residue_mod(
            model,
            positive(r).multiply(-2)?.add(Integer::from_i128(-b0))?,
            two_a,
        )?;
        let b = if residue.cmp(a).is_gt() {
            positive(two_a.difference(residue)).negate()
        } else {
            positive(residue)
        };
        let c = self.last_coefficient(model, a, b)?;
        self.canonical(model, Form { a, b, c })
    }

    fn last_coefficient(self, model: &mut Model, a: Magnitude, b: Integer) -> Result<Magnitude> {
        let numerator = b.multiply_wide(b)?.add(positive(self.magnitude()))?;
        let denominator = a.multiply_u128(4).ok_or_else(arithmetic)?;
        Ok(exact_divide(model, numerator, denominator)?.magnitude)
    }

    fn canonical(self, model: &mut Model, mut form: Form) -> Result<Class> {
        if form.discriminant()? != Some(self.magnitude()) {
            return Err(Error::plain(
                "CM class form discriminant differs from the named order",
            ));
        }
        if form.content(model)?.small() != Some(1) {
            return Err(Error::plain("CM class requires a primitive form"));
        }
        form.reduce(model, |_, _| Ok(()))?;
        Class::from_reduced(form)
    }

    pub(crate) fn power(
        self,
        model: &mut Model,
        mut exponent: u128,
        mut base: Class,
    ) -> Result<Class> {
        let mut value = self.principal();
        while exponent != 0 {
            if exponent & 1 != 0 {
                value = self.compose(model, value, base)?;
            }
            exponent >>= 1;
            if exponent != 0 {
                base = self.compose(model, base, base)?;
            }
        }
        Ok(value)
    }

    pub(crate) fn signed_power(
        self,
        model: &mut Model,
        exponent: i128,
        base: Class,
    ) -> Result<Class> {
        let base = if exponent < 0 {
            self.inverse(model, base)?
        } else {
            base
        };
        self.power(model, exponent.unsigned_abs(), base)
    }

    pub(crate) fn class_order(self, model: &mut Model, element: Class) -> Result<usize> {
        if element == self.principal() {
            return Ok(1);
        }
        let mut size = 0;
        cm::each_primitive_reduced(model, self.discriminant.unsigned_abs(), |_| {
            size += 1;
            Ok(())
        })?;
        let primes = factor_order(model, size)?;
        self.order_from_size(model, element, size, &primes)
    }

    fn order_from_size(
        self,
        model: &mut Model,
        element: Class,
        mut size: usize,
        primes: &[(usize, u32)],
    ) -> Result<usize> {
        for &(p, _) in primes {
            while size % p == 0
                && self.power(model, (size / p) as u128, element)? == self.principal()
            {
                size /= p;
            }
        }
        Ok(size)
    }

    pub(crate) fn elements(self, model: &mut Model, limit: usize) -> Result<Vec<Class>> {
        let mut elements = Vec::new();
        cm::each_primitive_reduced(model, self.discriminant.unsigned_abs(), |form| {
            if elements.len() >= limit {
                return Err(Error::plain(
                    "CM class group exceeds the finite-value limit",
                ));
            }
            elements.push(Class(form));
            Ok(())
        })?;
        Ok(elements)
    }

    pub(crate) fn same_genus(
        self,
        model: &mut Model,
        x: Class,
        y: Class,
        limit: usize,
    ) -> Result<bool> {
        if x == y {
            return Ok(true);
        }
        let inverse = self.inverse(model, y)?;
        let relative = self.compose(model, x, inverse)?;
        for z in self.elements(model, limit)? {
            if self.compose(model, z, z)? == relative {
                return Ok(true);
            }
        }
        Ok(false)
    }

    pub(crate) fn group_report(
        self,
        model: &mut Model,
        name: &str,
        limit: usize,
    ) -> Result<String> {
        let elements = self.elements(model, limit)?;
        let size = elements.len();
        let primes = factor_order(model, size)?;
        let orders = elements
            .into_iter()
            .map(|element| self.order_from_size(model, element, size, &primes))
            .collect::<Result<Vec<_>>>()?;
        let mut factors: Vec<usize> = Vec::new();
        for (p, exponent) in primes {
            let primary = primary_factors(model, &orders, p, exponent)?;
            if primary.len() > factors.len() {
                factors.splice(0..0, std::iter::repeat_n(1, primary.len() - factors.len()));
            }
            for (factor, part) in factors.iter_mut().rev().zip(primary.into_iter().rev()) {
                *factor *= part;
            }
        }
        let genera = 1usize << factors.iter().filter(|&&factor| factor % 2 == 0).count();
        Ok(format!("rvCMClassGroup[{name}] = discriminant={}, order={size}, invariant_factors={factors:?}, genera={genera}", self.discriminant))
    }
}

impl Class {
    pub(super) fn from_reduced(form: Form) -> Result<Self> {
        Ok(Self([
            form.a.signed(false).ok_or_else(arithmetic)?,
            form.b.signed().ok_or_else(arithmetic)?,
            form.c.signed(false).ok_or_else(arithmetic)?,
        ]))
    }
}

fn factor_order(model: &mut Model, mut n: usize) -> Result<Vec<(usize, u32)>> {
    let mut factors = Vec::new();
    let mut p = 2;
    while p <= n / p {
        model.spend(1)?;
        if n % p == 0 {
            let mut e = 0;
            while n % p == 0 {
                n /= p;
                e += 1;
            }
            factors.push((p, e));
        }
        p += 1;
    }
    if n > 1 {
        factors.push((n, 1));
    }
    Ok(factors)
}

fn primary_factors(
    model: &mut Model,
    orders: &[usize],
    p: usize,
    exponent: u32,
) -> Result<Vec<usize>> {
    model.spend(orders.len() * exponent as usize)?;
    let mut primary = Vec::new();
    let mut power = 1;
    let mut previous = 0;
    for _ in 0..exponent {
        power *= p;
        let mut count = orders.iter().filter(|&&order| power % order == 0).count();
        let mut log = 0;
        while count > 1 && count % p == 0 {
            count /= p;
            log += 1;
        }
        if count != 1 || log < previous {
            return Err(arithmetic());
        }
        let rank = log - previous;
        if power == p {
            primary.resize(rank, 1);
        }
        for factor in primary.iter_mut().rev().take(rank) {
            *factor *= p;
        }
        previous = log;
    }
    Ok(primary)
}

fn arithmetic() -> Error {
    Error::plain("CM class group calculation exceeds the supported exact arithmetic range")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_small_ideal_product_is_closed_and_associative() -> std::result::Result<(), String> {
        verify_group_laws().map_err(|error| error.to_string())
    }

    fn verify_group_laws() -> Result<()> {
        let mut model = Model::new();
        for d in [-3, -4, -20, -23, -63, -84, -207, -368, -576] {
            let order = Order::new(d)?;
            let principal = order.principal();
            let classes = order.elements(&mut model, 100)?;
            for &x in &classes {
                model.reset_work();
                let inverse = order.inverse(&mut model, x)?;
                assert_eq!(order.compose(&mut model, x, principal)?.0, x.0);
                assert_eq!(order.compose(&mut model, x, inverse)?.0, principal.0);
                for &y in &classes {
                    let xy = order.compose(&mut model, x, y)?;
                    assert!(classes.contains(&xy), "discriminant {d}");
                    assert_eq!(xy.0, order.compose(&mut model, y, x)?.0);
                    for &z in &classes {
                        let yz = order.compose(&mut model, y, z)?;
                        assert_eq!(
                            order.compose(&mut model, xy, z)?.0,
                            order.compose(&mut model, x, yz)?.0,
                            "discriminant {d}, {:?}, {:?}, {:?}",
                            x.0,
                            y.0,
                            z.0
                        );
                    }
                }
            }
        }
        Ok(())
    }
}
