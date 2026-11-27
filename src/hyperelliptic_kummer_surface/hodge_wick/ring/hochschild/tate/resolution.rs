use super::{check_report, storage, BTreeMap, Definition, Error, Mode, Model, Result, Write};

type Bimodule = BTreeMap<(usize, usize), i128>;

#[derive(Clone, Copy)]
enum Map {
    Difference,
    Norm,
}

impl Map {
    fn from_degree(degree: i128) -> Self {
        if degree.rem_euclid(2) == 0 {
            Self::Norm
        } else {
            Self::Difference
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Difference => "X-Y",
            Self::Norm => "norm",
        }
    }

    fn rank(self, d: usize) -> u128 {
        let d = d as u128;
        match self {
            // The quotient for the first map is A. The nonzero images
            // of the norm occur once in each total degree d-1,...,2d-2.
            Self::Difference => d * (d - 1),
            Self::Norm => d,
        }
    }
}

fn apply(
    model: &mut Model,
    d: usize,
    map: Map,
    a: usize,
    b: usize,
    mut emit: impl FnMut((usize, usize), i128) -> Result<()>,
) -> Result<()> {
    match map {
        Map::Difference => {
            model.spend(1)?;
            if a < d - 1 {
                emit((a + 1, b), 1)?;
            }
            if b < d - 1 {
                emit((a, b + 1), -1)?;
            }
        }
        Map::Norm => {
            let max = d - 1 - a;
            if b <= max {
                let count = max - b + 1;
                storage(count.saturating_mul(8))?;
                model.spend(count)?;
                for j in b..=max {
                    emit((a + j, b + (d - 1 - j)), 1)?;
                }
            }
        }
    }
    Ok(())
}

fn image(model: &mut Model, d: usize, map: Map, a: usize, b: usize) -> Result<Bimodule> {
    let mut result = Bimodule::new();
    apply(model, d, map, a, b, |term, coefficient| {
        result.insert(term, coefficient);
        Ok(())
    })?;
    Ok(result)
}

fn composition(model: &mut Model, d: usize, first: Map, second: Map) -> Result<Bimodule> {
    let mut result = Bimodule::new();
    for ((a, b), entry) in image(model, d, first, 0, 0)? {
        apply(model, d, second, a, b, |term, coefficient| {
            let value = result.entry(term).or_default();
            *value += entry * coefficient;
            if *value == 0 {
                result.remove(&term);
            }
            storage(result.len().saturating_mul(8))
        })?;
    }
    Ok(result)
}

fn validate_degree(mode: Mode, degree: i128) -> Result<()> {
    if mode == Mode::Ordinary && degree < 0 {
        Err(Error::plain(
            "ordinary Tate resolution maps have nonnegative degree",
        ))
    } else {
        Ok(())
    }
}

impl Definition {
    pub(crate) fn period_map(
        &self,
        model: &mut Model,
        index: usize,
        mode: Mode,
        degree: i128,
        a: usize,
        b: usize,
    ) -> Result<String> {
        validate_degree(mode, degree)?;
        let d = self.factor(index)?;
        if a >= d || b >= d {
            return Err(Error::plain(
                "Tate bimodule exponent is outside the truncated ring",
            ));
        }
        if mode == Mode::Ordinary && degree == 0 {
            return Ok(if a < d - b {
                format!("augmentation=[({}, 1)]", a + b)
            } else {
                String::from("augmentation=[]")
            });
        }
        let map = Map::from_degree(degree);
        let mut text = format!("{}=[", map.name());
        let mut separator = "";
        for ((x, y), value) in image(model, d, map, a, b)? {
            write!(text, "{separator}({x}, {y}, {value})").unwrap();
            check_report(&text)?;
            separator = ", ";
        }
        text.push(']');
        Ok(text)
    }

    pub(crate) fn tensor_period(
        &self,
        model: &mut Model,
        indices: &[usize],
        left: &[usize],
        right: &[usize],
    ) -> Result<String> {
        let count = self.powers.len();
        if indices.len() != count || left.len() != count || right.len() != count {
            return Err(Error::plain(
                "Tate tensor resolution lists need one entry per factor",
            ));
        }
        let mut report = "[".to_owned();
        let mut separator = "";
        let mut prefix_odd = false;
        for i in 0..count {
            let d = self.powers[i];
            if left[i] >= d || right[i] >= d {
                return Err(Error::plain(
                    "Tate bimodule exponent is outside the truncated ring",
                ));
            }
            let degree = indices[i];
            if degree == 0 {
                continue;
            }
            let map = Map::from_degree(degree as i128);
            let sign = if prefix_odd { -1 } else { 1 };
            apply(model, d, map, left[i], right[i], |(a, b), coefficient| {
                let mut target = indices.to_vec();
                let mut x = left.to_vec();
                let mut y = right.to_vec();
                target[i] -= 1;
                x[i] = a;
                y[i] = b;
                write!(
                    report,
                    "{separator}({target:?}, {x:?}, {y:?}, {})",
                    sign * coefficient
                )
                .unwrap();
                separator = ", ";
                check_report(&report)
            })?;
            prefix_odd ^= !degree.is_multiple_of(2);
        }
        report.push(']');
        Ok(report)
    }

    pub(crate) fn resolution_check(
        &self,
        model: &mut Model,
        index: usize,
        mode: Mode,
        degree: i128,
    ) -> Result<String> {
        validate_degree(mode, degree)?;
        let d = self.factor(index)?;
        for (a, b) in [(Map::Difference, Map::Norm), (Map::Norm, Map::Difference)] {
            if !composition(model, d, a, b)?.is_empty() {
                return Err(Error::plain(
                    "Verbitsky Tate period maps do not compose to zero",
                ));
            }
        }
        let through_degree_zero = mode == Mode::Ordinary && degree == 0;
        let dimension = if through_degree_zero { d } else { d - 1 };
        let map = Map::from_degree(degree);
        let (map, rank) = if through_degree_zero {
            ("augmentation", d as u128)
        } else {
            (map.name(), map.rank(d))
        };
        Ok(format!("map={map} rank={rank} period_compositions=0 HH={dimension} Hochschild_cohomology={dimension}"))
    }

    pub(crate) fn complex_map(
        &self,
        index: usize,
        mode: Mode,
        degree: i128,
        cochain: bool,
        exponent: usize,
    ) -> Result<String> {
        validate_degree(mode, degree)?;
        let d = self.factor(index)?;
        if exponent >= d {
            return Err(Error::plain(
                "Tate chain exponent is outside the truncated ring",
            ));
        }
        let active = if cochain {
            degree.rem_euclid(2) == 1
        } else {
            (mode == Mode::Stable || degree > 0) && degree.rem_euclid(2) == 0
        };
        Ok(if active && exponent == 0 {
            format!("[({}, {d})]", d - 1)
        } else {
            String::from("[]")
        })
    }
}
