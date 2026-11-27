use crate::rv_cyclotomic as cyclo;

use super::{check_report, Error, Label, Model, Object, ParsedExpr, Rc, Result, Value};

mod field;
mod fusion;
mod jordan;

pub(crate) use fusion::verify_fusion;
pub(crate) use jordan::Push;

pub(crate) enum Command {
    Declare(Label, Label, Vec<usize>),
    Object(Label),
    Push(Label, Label, Label),
    PushReport(Label),
    PushResult(Label, Label),
    Fusion(Label, Label),
}

pub(crate) enum Expression {
    Coefficient(Label, usize, usize),
    Bulk(Label, usize, bool),
    Matrix(Label, usize, usize, usize, usize),
    Map(Label, Label, usize, usize, usize, Box<ParsedExpr>),
}

type Polynomial = Vec<Vec<i128>>;

pub(crate) struct Defect {
    pub(crate) mark: Rc<str>,
    pub(crate) d: usize,
    indices: Vec<usize>,
    p: Polynomial,
    q: Polynomial,
}

fn factor(
    model: &mut Model,
    d: usize,
    exponents: impl IntoIterator<Item = usize>,
) -> Result<Polynomial> {
    let mut polynomial = vec![vec![1]];
    for r in exponents {
        let mut product = vec![Vec::new(); polynomial.len() + 1];
        for (a, c) in polynomial.iter().enumerate() {
            cyclo::accumulate(model, &mut product[a + 1], c, false)?;
            let term = cyclo::shifted(model, d, c, r as i128)?;
            cyclo::accumulate(model, &mut product[a], &term, true)?;
        }
        polynomial = product;
    }
    Ok(polynomial)
}

impl Defect {
    pub(crate) fn marked(model: &mut Model, mark: Rc<str>, indices: Vec<usize>) -> Result<Self> {
        let definition = model.verbitsky_tate(&mark)?;
        let [d] = definition.powers.as_slice() else {
            return Err(Error::plain("LG defect requires a rank-one Tate marking"));
        };
        Self::new(model, mark, *d, indices)
    }

    pub(crate) fn new(
        model: &mut Model,
        mark: Rc<str>,
        d: usize,
        mut indices: Vec<usize>,
    ) -> Result<Self> {
        let field = model.cyclotomic(d)?;
        indices.sort_unstable();
        if indices.is_empty() || indices.len() >= d || indices.last().is_some_and(|&r| r >= d) {
            return Err(Error::plain(
                "LG defect requires a nonempty proper subset of root indices 0 through d−1",
            ));
        }
        if indices.windows(2).any(|r| r[0] == r[1]) {
            return Err(Error::plain("LG defect root index occurs twice"));
        }
        super::storage((d + 2).saturating_mul(3 + 2 * field.len()))?;
        let p = factor(model, d, indices.iter().copied())?;
        let q = factor(
            model,
            d,
            (0..d).filter(|r| indices.binary_search(r).is_err()),
        )?;
        let mut product = vec![Vec::new(); d + 1];
        for (a, p) in p.iter().enumerate() {
            for (b, q) in q.iter().enumerate() {
                model.spend(1)?;
                let term = cyclo::product_coefficients(model, d, p, q)?;
                cyclo::accumulate(model, &mut product[a + b], &term, false)?;
            }
        }
        if product[0] != [-1] || product[d] != [1] || product[1..d].iter().any(|c| !c.is_empty()) {
            return Err(Error::plain("LG defect square is not x^d−y^d"));
        }
        Ok(Self {
            mark,
            d,
            indices,
            p,
            q,
        })
    }

    pub(crate) fn size(&self) -> usize {
        self.indices.len()
            + self
                .p
                .iter()
                .chain(&self.q)
                .map(|c| 3 + 2 * c.len())
                .sum::<usize>()
            + 12
    }

    pub(crate) fn report(&self) -> Result<String> {
        let report = format!(
            "W=x^{}−y^{} support={:?} p={:?} q={:?} square=verified",
            self.d, self.d, self.indices, self.p, self.q
        );
        check_report(&report)?;
        Ok(report)
    }

    pub(crate) fn coefficient(
        &self,
        model: &mut Model,
        factor: usize,
        power: usize,
    ) -> Result<Value> {
        let polynomial = match factor {
            0 => &self.p,
            1 => &self.q,
            _ => return Err(Error::plain("LG defect factor index must be zero or one")),
        };
        let coefficient = polynomial
            .get(power)
            .ok_or_else(|| Error::plain("LG defect coefficient index is out of range"))?;
        cyclo::persist(model, coefficient.clone())
    }

    pub(crate) fn bulk(&self, model: &mut Model, degree: usize, reverse: bool) -> Result<Value> {
        let coefficient = self.transmission(model, degree, reverse)?;
        cyclo::persist(model, coefficient)
    }

    fn transmission(&self, model: &mut Model, degree: usize, reverse: bool) -> Result<Vec<i128>> {
        if degree >= self.d - 1 {
            return Err(Error::plain(
                "LG defect bulk degree must lie between zero and d−2",
            ));
        }
        // Reversing orientation exchanges which powers survive the bulk residue.
        let k = if reverse { self.d - 2 - degree } else { degree };
        let r = self.p.len() - 1;
        let s = self.q.len() - 1;
        let mut residue = Vec::new();
        for (a, p) in self.p.iter().enumerate() {
            let Some(b) = (k + 1).checked_sub(a).filter(|&b| b <= s) else {
                continue;
            };
            let coefficient = (a * (s - b)) as i128 - (b * (r - a)) as i128;
            if coefficient == 0 {
                continue;
            }
            let mut term = cyclo::product_coefficients(model, self.d, p, &self.q[b])?;
            for c in &mut term {
                *c = c
                    .checked_mul(coefficient)
                    .ok_or_else(|| Error::plain("LG defect residue coefficient is out of range"))?;
            }
            cyclo::accumulate(model, &mut residue, &term, false)?;
        }

        let mut monomials = vec![0; self.d];
        for &r in &self.indices {
            let exponent = (self.d - (r * (k + 1)) % self.d) % self.d;
            monomials[exponent] += 1;
        }
        let character = cyclo::reduced(model, self.d, monomials)?;
        let d = self.d as i128;
        for coefficient in &mut residue {
            if *coefficient % d != 0 {
                return Err(Error::plain(
                    "LG defect residue is not integral in the cyclotomic basis",
                ));
            }
            *coefficient /= d;
        }
        if residue != character {
            return Err(Error::plain(
                "LG defect derivative residue failed the cyclotomic character identity",
            ));
        }
        Ok(residue)
    }

    fn check_source(&self, source: &Object) -> Result<()> {
        if self.mark != source.mark || self.d != source.d {
            return Err(Error::plain(
                "LG defect and boundary brane require the same Tate marking",
            ));
        }
        Ok(())
    }
}
