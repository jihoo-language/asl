use super::operations::insertions;
use super::{
    grading_range, storage, work_limit, Algebra, BTreeMap, Cochain, Complex, Entry, Error, Form,
    Grade, Key, Label, Model, ParsedExpr, Rational, Rc, Result, Value, Write, MAX_SHOWN_BYTES,
};

pub(crate) enum Command {
    Declare(Label, usize, Vec<ParsedExpr>),
    Report(Label),
    Product(Label, usize, usize),
    Direction(Label, Vec<Entry>),
}

#[derive(Clone, Copy)]
pub(crate) enum Lookup {
    Coefficient,
    Residual,
    Obstruction,
}

pub(crate) struct Expression(pub(crate) Lookup, pub(crate) Label, pub(crate) Vec<usize>);

type Powers = Vec<usize>;
type Layer = BTreeMap<Powers, Rc<Cochain>>;

pub(crate) struct Jet {
    grade: Grade,
    algebra: Rc<Algebra>,
    shifts: Vec<i128>,
    requested: usize,
    coefficients: Vec<Layer>,
    residuals: Vec<Layer>,
    obstructions: Layer,
}

impl Jet {
    pub(crate) fn new(
        model: &mut Model,
        form: &Form,
        grade: &Grade,
        requested: usize,
        tangents: Vec<(i128, Rc<Cochain>)>,
    ) -> Result<Self> {
        let count = tangents.len();
        if count == 0 || requested == 0 {
            return Err(Error::plain(
                "Verbitsky associative jet needs a positive order and at least one parameter",
            ));
        }
        let terms = monomial_count(requested, count)?;
        storage(terms.saturating_mul(count.saturating_add(8)))?;
        model.spend(terms)?;
        let algebra = form.hochschild_algebra(model, grade)?;
        let mut jet = Self {
            grade: grade.clone(),
            algebra,
            shifts: tangents.iter().map(|(shift, _)| *shift).collect(),
            requested,
            coefficients: vec![Layer::new()],
            residuals: vec![Layer::new(), Layer::new()],
            obstructions: Layer::new(),
        };
        let mut first = Layer::new();
        let mut cells = 20 + count * 2;
        for (i, (shift, value)) in tangents.into_iter().enumerate() {
            if value.empty() {
                continue;
            }
            let initial_grade = Grade {
                arity: 2,
                shift,
                ..grade.clone()
            };
            if !form
                .hochschild_differential(model, &initial_grade, &value)?
                .empty()
            {
                return Err(Error::plain(format!(
                    "Verbitsky associative tangent {i} is not a Hochschild 2-cocycle"
                )));
            }
            let mut power = vec![0; count];
            power[i] = 1;
            cells = record_cells(cells, &power, &value)?;
            first.insert(power, value);
        }
        jet.coefficients.push(first);
        for total in 2..=requested {
            model.spend(total)?;
            let residuals = jet.associator(model, total)?;
            let mut corrections = Layer::new();
            let mut recorded = Layer::new();
            for (power, raw) in residuals {
                if raw.empty() {
                    continue;
                }
                let grade = jet.component(&power, 3)?;
                let complex = Complex::new(model, &jet.algebra, &grade)?;
                let obstruction = complex.project(model, &raw)?;
                if obstruction.empty() {
                    let correction = complex.primitive(model, &jet.algebra, &grade, &raw)?;
                    let correction_grade = grade.with_arity_offset(-1)?;
                    let check =
                        form.hochschild_differential(model, &correction_grade, &correction)?;
                    if !check.compare(model, &raw)?.is_eq() {
                        return Err(Error::plain(
                            "Verbitsky associative primitive failed its associativity equation",
                        ));
                    }
                    cells = record_cells(cells, &power, &correction)?;
                    corrections.insert(power.clone(), Rc::new(correction));
                } else {
                    cells = record_cells(cells, &power, &obstruction)?;
                    jet.obstructions.insert(power.clone(), Rc::new(obstruction));
                }
                cells = record_cells(cells, &power, &raw)?;
                recorded.insert(power, Rc::new(raw));
            }
            jet.residuals.push(recorded);
            if !jet.obstructions.is_empty() {
                // Every monomial of this total order was tested. Only
                // the earlier complete layers give an associative jet.
                break;
            }
            jet.coefficients.push(corrections);
        }
        model.allocate(cells + 8 * (jet.coefficients.len() + jet.residuals.len()))?;
        Ok(jet)
    }

    fn associator(&self, model: &mut Model, total: usize) -> Result<BTreeMap<Powers, Cochain>> {
        let mut residuals = BTreeMap::<Powers, Cochain>::new();
        let mut cells = 0;
        for i in 1..total {
            for (a, f) in &self.coefficients[i] {
                for (b, g) in &self.coefficients[total - i] {
                    model.spend(1 + self.shifts.len())?;
                    let power = a.iter().zip(b).map(|(&x, &y)| x + y).collect::<Vec<_>>();
                    let entry = residuals.entry(power).or_insert_with(|| {
                        cells += self.shifts.len() + 12;
                        Cochain::default()
                    });
                    let prior = entry.cells;
                    insertions(model, entry, f, g, false)?;
                    cells = cells - prior + entry.cells;
                    storage(cells)?;
                }
            }
        }
        Ok(residuals)
    }

    fn through(&self) -> usize {
        self.coefficients.len() - 1
    }

    fn visited(&self) -> usize {
        self.residuals.len() - 1
    }

    fn total(&self, power: &[usize]) -> Result<usize> {
        if power.len() != self.shifts.len() {
            return Err(Error::plain(
                "Verbitsky associative monomial has the wrong number of exponents",
            ));
        }
        power.iter().try_fold(0usize, |sum, &value| {
            sum.checked_add(value).ok_or_else(grading_range)
        })
    }

    fn component(&self, power: &[usize], arity: i128) -> Result<Grade> {
        let shift =
            power
                .iter()
                .zip(&self.shifts)
                .try_fold(0i128, |sum, (&exponent, &weight)| {
                    weight
                        .checked_mul(exponent as i128)
                        .and_then(|term| sum.checked_add(term))
                        .ok_or_else(grading_range)
                })?;
        Ok(Grade {
            arity,
            shift,
            ..self.grade.clone()
        })
    }

    pub(crate) fn lookup_grade(&self, operation: Lookup, power: &[usize]) -> Result<(Grade, bool)> {
        self.layer(operation, power)?;
        let arity = if matches!(operation, Lookup::Coefficient) {
            2
        } else {
            3
        };
        Ok((
            self.component(power, arity)?,
            matches!(operation, Lookup::Obstruction),
        ))
    }

    fn layer(&self, operation: Lookup, power: &[usize]) -> Result<Option<&Layer>> {
        let total = self.total(power)?;
        match operation {
            Lookup::Coefficient if (1..=self.through()).contains(&total) => {
                Ok(self.coefficients.get(total))
            }
            Lookup::Coefficient => Err(Error::plain(
                "Verbitsky associative coefficient is outside the certified positive orders",
            )),
            Lookup::Residual if (2..=self.visited()).contains(&total) => {
                Ok(self.residuals.get(total))
            }
            Lookup::Obstruction if (2..=self.visited()).contains(&total) => {
                Ok(Some(&self.obstructions))
            }
            _ => Err(Error::plain(
                "Verbitsky associative obstruction order has yet to be tested",
            )),
        }
    }

    pub(crate) fn lookup(
        &self,
        model: &mut Model,
        operation: Lookup,
        power: &[usize],
    ) -> Result<Value> {
        match self
            .layer(operation, power)?
            .and_then(|layer| layer.get(power))
        {
            Some(value) => Ok(Value::VerbitskyHochschild(Rc::clone(value))),
            None => Cochain::default().persist(model),
        }
    }

    pub(crate) fn report(&self, model: &mut Model) -> Result<String> {
        let mut text = format!(
            "weights={:?} requested={} associative_through={} corrections=[",
            self.shifts,
            self.requested,
            self.through()
        );
        let mut separator = "";
        for layer in self.coefficients.iter().skip(2) {
            for (power, value) in layer {
                write!(text, "{separator}({power:?}, {})", value.text(model)?).unwrap();
                check_report(&text)?;
                separator = ", ";
            }
        }
        text.push_str("] obstructions=[");
        separator = "";
        for (power, class) in &self.obstructions {
            write!(
                text,
                "{separator}({power:?}, w={}, {})",
                self.component(power, 3)?.shift,
                class.text(model)?
            )
            .unwrap();
            check_report(&text)?;
            separator = ", ";
        }
        text.push(']');
        check_report(&text)?;
        Ok(text)
    }

    pub(crate) fn product(&self, model: &mut Model, left: usize, right: usize) -> Result<String> {
        if left >= self.algebra.len() || right >= self.algebra.len() {
            return Err(Error::plain(
                "Verbitsky associative product index is outside the ring basis",
            ));
        }
        let mut report = "[".to_owned();
        let mut separator = "";
        let original = self.algebra.product(left, right);
        if !original.is_empty() {
            let zero = vec![0; self.shifts.len()];
            report_product(model, &mut report, &zero, original.iter().copied())?;
            separator = ", ";
        }
        let start = Key {
            word: vec![left, right],
            output: 0,
        };
        for layer in self.coefficients.iter().skip(1) {
            for (power, value) in layer {
                model.spend(1)?;
                let mut entries = value
                    .entries
                    .range(start.clone()..)
                    .take_while(|(key, _)| key.word == start.word)
                    .peekable();
                if entries.peek().is_some() {
                    report.push_str(separator);
                    report_product(
                        model,
                        &mut report,
                        power,
                        entries.map(|(key, &value)| (key.output, value)),
                    )?;
                    separator = ", ";
                }
            }
        }
        report.push(']');
        check_report(&report)?;
        Ok(report)
    }

    pub(crate) fn direction(&self, model: &mut Model, values: &[Entry]) -> Result<String> {
        if values.len() != self.shifts.len() {
            return Err(Error::plain(
                "Verbitsky associative direction has the wrong number of parameters",
            ));
        }
        let values = values
            .iter()
            .map(|&value| Rational::input(model, value))
            .collect::<Result<Vec<_>>>()?;
        if self.obstructions.is_empty() {
            return Ok(format!("associative_through={}", self.through()));
        }
        let mut groups = BTreeMap::<i128, Cochain>::new();
        for (power, value) in &self.obstructions {
            if power
                .iter()
                .zip(&values)
                .any(|(&exponent, &coordinate)| exponent > 0 && coordinate.is_zero())
            {
                continue;
            }
            let mut factor = Rational::integer(1);
            for (&exponent, &coordinate) in power.iter().zip(&values) {
                if exponent > 0 {
                    let value = rational_power(model, coordinate, exponent)?;
                    factor = factor.multiply(model, value)?;
                }
            }
            let group = groups.entry(self.component(power, 3)?.shift).or_default();
            for (key, &entry) in &value.entries {
                let coefficient = entry.multiply(model, factor)?;
                group.insert(model, key.clone(), coefficient)?;
            }
        }
        let obstructed = groups.values().any(|group| !group.empty());
        let mut report = format!("order={} obstructed={obstructed} classes=[", self.visited());
        let mut separator = "";
        for (shift, class) in groups.iter().filter(|(_, class)| !class.empty()) {
            write!(report, "{separator}({shift}, {})", class.text(model)?).unwrap();
            check_report(&report)?;
            separator = ", ";
        }
        report.push(']');
        check_report(&report)?;
        Ok(report)
    }
}

fn monomial_count(order: usize, count: usize) -> Result<usize> {
    let size = order as u128 + count as u128;
    let k = order.min(count) as u128;
    let mut count = 1u128;
    for i in 1..=k {
        count = count.checked_mul(size - k + i).ok_or_else(work_limit)? / i;
        storage(usize::try_from(count).map_err(|_| work_limit())?)?;
    }
    usize::try_from(count).map_err(|_| work_limit())
}

fn record_cells(cells: usize, powers: &[usize], value: &Cochain) -> Result<usize> {
    let next = cells + powers.len() + value.cells + 12;
    storage(next)?;
    Ok(next)
}

fn rational_power(model: &mut Model, mut base: Rational, mut exponent: usize) -> Result<Rational> {
    let mut value = Rational::integer(1);
    while exponent > 0 {
        if !exponent.is_multiple_of(2) {
            value = value.multiply(model, base)?;
        }
        exponent /= 2;
        if exponent > 0 {
            base = base.multiply(model, base)?;
        }
    }
    Ok(value)
}

fn report_product(
    model: &mut Model,
    report: &mut String,
    powers: &[usize],
    coordinates: impl IntoIterator<Item = (usize, Rational)>,
) -> Result<()> {
    write!(report, "({powers:?}, [").unwrap();
    let mut separator = "";
    for (index, coefficient) in coordinates {
        write!(report, "{separator}({index}, {})", coefficient.text(model)?).unwrap();
        check_report(report)?;
        separator = ", ";
    }
    report.push_str("])");
    check_report(report)
}

fn check_report(report: &str) -> Result<()> {
    if report.len() >= MAX_SHOWN_BYTES {
        Err(Error::plain(
            "Verbitsky associative jet exceeds the report limit",
        ))
    } else {
        Ok(())
    }
}
