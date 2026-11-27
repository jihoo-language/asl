use super::{
    integral, persist, persist_class, Braid, Cell, Chain, ChainSum, Class, Context, Error, Form,
    Grade, Model, Result, Saddle, Slice,
};
use crate::lee_linear::{Echelon, Vector};
use crate::rv_lens::gcd;

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Action {
    E,
    F,
    H,
    Minus,
    Zero,
    PlusTwice,
    D0,
    LeeZero,
    LeePlus,
}

impl Action {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::E => "khE",
            Self::F => "khF",
            Self::H => "khH",
            Self::Minus => "khVMinus",
            Self::Zero => "khVZero",
            Self::PlusTwice => "khWPlus",
            Self::D0 => "khD0",
            Self::LeeZero => "khLee0",
            Self::LeePlus => "khLeePlus",
        }
    }

    fn shifts(self) -> (i64, i64, i64) {
        match self {
            Self::E => (0, 2, 2),
            Self::F => (0, -2, -2),
            Self::H => (0, 0, 0),
            Self::Minus => (1, 0, -2),
            Self::Zero => (1, 2, 0),
            Self::PlusTwice | Self::LeePlus => (1, 4, 2),
            Self::LeeZero => (1, 4, 0),
            Self::D0 => (1, 0, 0),
        }
    }

    pub(crate) fn target(self, grade: Grade, form: Form) -> Result<(Grade, Form)> {
        let lee = matches!(self, Self::LeeZero | Self::LeePlus);
        let input_k = match form {
            Form::AnnularChain(k) => k,
            Form::AnnularClass(k) if !lee => k,
            _ => {
                return Err(Error::plain(format!(
                    "{} requires an annular Khovanov {}",
                    self.name(),
                    if lee { "chain" } else { "chain or class" }
                )));
            }
        };
        let (di, dj, dk) = self.shifts();
        let grade = Grade {
            j: grade
                .j
                .checked_add(dj)
                .ok_or_else(|| Error::plain("Khovanov quantum degree is out of range"))?,
            ..grade.shift(di)?
        };
        let k = input_k
            .checked_add(dk)
            .ok_or_else(|| Error::plain("Khovanov annular degree is out of range"))?;
        Ok((
            grade,
            if form.is_class() {
                Form::AnnularClass(k)
            } else {
                Form::AnnularChain(k)
            },
        ))
    }

    fn column(
        self,
        context: &mut Context<'_>,
        model: &mut Model,
        cell: &Cell,
        coefficient: i128,
        output: &mut ChainSum,
    ) -> Result<()> {
        let resolution = context.resolution(model, cell.state)?;
        match self {
            Self::E | Self::F => {
                // Circle labels are assigned by smallest vertex. Seam vertices
                // come first, in radial order, so the essential labels in
                // increasing order give the first seam encounters.
                model.spend(resolution.essential.len())?;
                let essential = resolution
                    .essential
                    .iter()
                    .enumerate()
                    .filter(|(_, essential)| **essential);
                for (depth, (circle, _)) in essential.enumerate() {
                    model.spend(1 + cell.xs.len().max(1).ilog2() as usize)?;
                    let place = cell.xs.binary_search(&circle);
                    if place.is_ok() == (self == Self::E) {
                        model.spend(cell.xs.len())?;
                        let mut xs = cell.xs.to_vec();
                        if let Ok(index) = place {
                            xs.remove(index);
                        } else {
                            xs.insert(place.unwrap_err(), circle);
                        }
                        let sign = if depth.is_multiple_of(2) { 1 } else { -1 };
                        output.add(
                            model,
                            Cell {
                                state: cell.state,
                                xs: xs.into(),
                            },
                            multiply(coefficient, sign)?,
                        )?;
                    }
                }
            }
            Self::H => {
                let weight = i128::from(resolution.k(model, &cell.xs)?);
                output.add(model, cell.clone(), multiply(coefficient, weight)?)?;
            }
            Self::Minus | Self::D0 | Self::LeeZero | Self::LeePlus => {
                let (_, _, shift) = self.shifts();
                let k = resolution.k(model, &cell.xs)? + shift;
                let saddle = if matches!(self, Self::LeeZero | Self::LeePlus) {
                    Saddle::Perturbation
                } else {
                    Saddle::Khovanov
                };
                context.edges_for(
                    model,
                    cell,
                    false,
                    saddle,
                    Some(Slice::Exact(k)),
                    |model, target, sign| output.add(model, target, multiply(coefficient, sign)?),
                )?;
            }
            Self::Zero | Self::PlusTwice => {
                let inner = if self == Self::Zero {
                    Self::Minus
                } else {
                    Self::Zero
                };
                let single = [(cell.clone(), coefficient)];
                let left = inner.apply(context, model, &single)?;
                let right = Self::E.apply(context, model, &single)?;
                Self::E.extend(context, model, &left, 1, output)?;
                inner.extend(context, model, &right, -1, output)?;
            }
        }
        Ok(())
    }

    fn extend(
        self,
        context: &mut Context<'_>,
        model: &mut Model,
        input: &[(Cell, i128)],
        scale: i128,
        output: &mut ChainSum,
    ) -> Result<()> {
        for (cell, coefficient) in input {
            self.column(context, model, cell, multiply(*coefficient, scale)?, output)?;
        }
        Ok(())
    }

    fn apply(
        self,
        context: &mut Context<'_>,
        model: &mut Model,
        input: &[(Cell, i128)],
    ) -> Result<Vec<(Cell, i128)>> {
        let mut output = ChainSum::default();
        self.extend(context, model, input, 1, &mut output)?;
        Ok(output.terms.into_iter().collect())
    }
}

fn multiply(a: i128, b: i128) -> Result<i128> {
    a.checked_mul(b).ok_or_else(integral::arithmetic)
}

fn relation(
    context: &mut Context<'_>,
    model: &mut Model,
    input: &[(Cell, i128)],
    name: &str,
    word_sum: &[(i128, &[Action])],
) -> Result<()> {
    let mut output = ChainSum::default();
    for &(scalar, word) in word_sum {
        let (first, rest) = word.split_first().expect("nonempty operator word");
        if let [second] = rest {
            let middle = second.apply(context, model, input)?;
            first.extend(context, model, &middle, scalar, &mut output)?;
        } else {
            debug_assert!(rest.is_empty());
            first.extend(context, model, input, scalar, &mut output)?;
        }
    }
    if output.terms.is_empty() {
        Ok(())
    } else {
        Err(Error::plain(format!(
            "annular Khovanov action failed {name}"
        )))
    }
}

impl Braid {
    pub(crate) fn action_chain(
        &self,
        model: &mut Model,
        action: Action,
        chain: &Chain,
    ) -> Result<Chain> {
        let output = action.apply(&mut Context::with_annular(self, true), model, chain)?;
        persist(model, output)
    }

    pub(crate) fn action_class(
        &self,
        model: &mut Model,
        action: Action,
        grade: Grade,
        k: i64,
        class: &Class,
    ) -> Result<Class> {
        if class.is_empty() {
            return Ok(Default::default());
        }
        let (target, target_form) = action.target(grade, Form::AnnularClass(k))?;
        let chain = self
            .group(model, grade, Some(k))?
            .representative(model, class)?;
        let image = self.action_chain(model, action, &chain)?;
        self.group(model, target, target_form.annular())?
            .project(model, &image)
    }

    pub(crate) fn action_check(&self, model: &mut Model) -> Result<usize> {
        use Action::{LeePlus, LeeZero, Minus, PlusTwice, Zero, D0, E, F, H};
        let mut context = Context::with_annular(self, true);
        context.enhanced_cells(model, |context, model, cell| {
            let input = [(cell, 1)];
            let mut verify =
                |name, terms: &[(i128, &[Action])]| relation(context, model, &input, name, terms);
            verify("[E,F]=H", &[(1, &[E, F]), (-1, &[F, E]), (-1, &[H])])?;
            verify("[H,E]=2E", &[(1, &[H, E]), (-1, &[E, H]), (-2, &[E])])?;
            verify("[H,F]=-2F", &[(1, &[H, F]), (-1, &[F, H]), (2, &[F])])?;
            for action in [E, F, H] {
                verify("even chain map", &[(1, &[D0, action]), (-1, &[action, D0])])?;
            }
            for action in [Minus, Zero, PlusTwice] {
                verify("odd chain map", &[(1, &[D0, action]), (1, &[action, D0])])?;
                let weight = i128::from(action.shifts().2);
                verify(
                    "odd weight",
                    &[(1, &[H, action]), (-1, &[action, H]), (-weight, &[action])],
                )?;
            }
            verify("[F,v-]=0", &[(1, &[F, Minus]), (-1, &[Minus, F])])?;
            verify(
                "[F,v0]=2v-",
                &[(1, &[F, Zero]), (-1, &[Zero, F]), (-2, &[Minus])],
            )?;
            verify(
                "[F,w+]=2v0",
                &[(1, &[F, PlusTwice]), (-1, &[PlusTwice, F]), (-2, &[Zero])],
            )?;
            verify("[E,w+]=0", &[(1, &[E, PlusTwice]), (-1, &[PlusTwice, E])])?;
            verify("w+=-2 Lee+", &[(1, &[PlusTwice]), (2, &[LeePlus])])?;
            for (a, b) in [
                (Minus, Minus),
                (Minus, Zero),
                (Zero, PlusTwice),
                (PlusTwice, PlusTwice),
            ] {
                verify("strict odd anticommutator", &[(1, &[a, b]), (1, &[b, a])])?;
            }
            verify(
                "{v-,w+}=2{d0,Lee0}",
                &[
                    (1, &[Minus, PlusTwice]),
                    (1, &[PlusTwice, Minus]),
                    (-2, &[D0, LeeZero]),
                    (-2, &[LeeZero, D0]),
                ],
            )?;
            verify(
                "v0^2=-{d0,Lee0}",
                &[(1, &[Zero, Zero]), (1, &[D0, LeeZero]), (1, &[LeeZero, D0])],
            )
        })
    }

    fn highest_relations(&self, model: &mut Model, grade: Grade, k: i64) -> Result<Vec<Vector>> {
        if !self.possible(grade, Form::AnnularClass(k)) {
            return Ok(Vec::new());
        }
        let source = self.group(model, grade, Some(k))?;
        if source.invariants.free == 0 {
            return Ok(Vec::new());
        }
        let (target_grade, target_form) = Action::E.target(grade, Form::AnnularClass(k))?;
        let target = self.group(model, target_grade, target_form.annular())?;
        let mut context = Context::with_annular(self, true);
        let mut echelon = Echelon::default();
        let mut kernel = Vec::new();
        let mut cells = 0;
        // Smith orders the free generators before the torsion generators.
        // Project the actual E-images first, then retain their free coordinates.
        for i in 0..source.invariants.free {
            let generator = source.generator(model, i)?;
            let cycle = source.representative(model, &generator)?;
            let image = Action::E.apply(&mut context, model, &cycle)?.into();
            let class = target.project(model, &image)?;
            let column = class
                .iter()
                .take_while(|&&(i, _)| i < target.invariants.free)
                .copied()
                .collect();
            if let Some(relation) = echelon.insert(model, column, vec![(i, 1)])? {
                cells += 3 * relation.len() + 1;
                integral::storage(cells)?;
                kernel.push(relation);
            }
        }
        Ok(kernel)
    }

    pub(crate) fn irreps(
        &self,
        model: &mut Model,
        i: i64,
        difference: i64,
    ) -> Result<Vec<(i64, usize)>> {
        let n = self.ends.len() as i64;
        let mut result = Vec::new();
        for k in (n.rem_euclid(2)..=n).step_by(2) {
            model.spend(1)?;
            let Some(j) = difference.checked_add(k) else {
                continue;
            };
            let dimension = self.highest_relations(model, Grade { i, j }, k)?.len();
            if dimension != 0 {
                result.push((k, dimension));
            }
        }
        Ok(result)
    }

    pub(crate) fn highest(
        &self,
        model: &mut Model,
        grade: Grade,
        k: i64,
        index: usize,
    ) -> Result<Chain> {
        if k < 0 {
            return Err(Error::plain("annular highest weight must be nonnegative"));
        }
        let relation = self
            .highest_relations(model, grade, k)?
            .into_iter()
            .nth(index)
            .ok_or_else(|| {
                Error::plain("annular highest-weight generator index is out of range")
            })?;
        let class = persist_class(model, relation.into_iter().collect())?;
        let cycle = self
            .group(model, grade, Some(k))?
            .representative(model, &class)?;
        let mut context = Context::with_annular(self, true);
        let mut powers = vec![cycle.to_vec()];
        let mut cells = cycle.iter().map(|(cell, _)| cell.cells()).sum::<usize>();
        loop {
            let next = Action::E.apply(&mut context, model, powers.last().unwrap())?;
            if next.is_empty() {
                break;
            }
            cells += 1 + next.iter().map(|(cell, _)| cell.cells()).sum::<usize>();
            integral::storage(cells)?;
            powers.push(next);
        }
        // The finite sl2 highest-weight projector has coefficients c_0=1,
        // c_r=-c_(r-1)/(r(k+r+1)) on F^r E^r. Clear a common denominator
        // and evaluate from the highest nonzero E-power by Horner's rule.
        let mut numerator = Vec::new();
        let mut absolute = 1;
        for (r, power) in powers.into_iter().enumerate().rev() {
            let mut sum = ChainSum::default();
            Action::F.extend(&mut context, model, &numerator, 1, &mut sum)?;
            let scalar = if r.is_multiple_of(2) {
                absolute
            } else {
                -absolute
            };
            for (cell, value) in power {
                sum.add(model, cell, multiply(value, scalar)?)?;
            }
            numerator = sum.terms.into_iter().collect();
            if r != 0 {
                absolute = multiply(absolute, r as i128 * (i128::from(k) + r as i128 + 1))?;
            }
        }
        let content = numerator
            .iter()
            .map(|(_, value)| value.unsigned_abs())
            .fold(0, gcd);
        if content == 0 {
            return Err(Error::plain(
                "annular highest-weight projector killed its class",
            ));
        }
        let content = i128::try_from(content).map_err(|_| integral::arithmetic())?;
        for (_, value) in &mut numerator {
            *value /= content;
        }
        persist(model, numerator)
    }
}
