use crate::cobar::{self, Monomial, Word};
use crate::error::Result;
use crate::model::Model;

#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub(crate) enum Profile {
    #[default]
    Full,
    An(usize),
}

impl Profile {
    pub(crate) fn can_restrict_to(self, target: Self) -> bool {
        match (self, target) {
            (Self::Full, _) => true,
            (Self::An(source), Self::An(target)) => target <= source,
            (Self::An(_), Self::Full) => false,
        }
    }

    pub(crate) fn basis(
        self,
        model: &mut Model,
        s: usize,
        t: usize,
        limit: usize,
    ) -> Result<Vec<Word>> {
        cobar::basis_profile(model, self, s, t, limit)
    }

    pub(crate) fn differential_word(self, model: &mut Model, word: &Word) -> Result<Vec<Word>> {
        self.differential_support(model, std::slice::from_ref(word))
    }

    pub(crate) fn differential_support(
        self,
        model: &mut Model,
        words: &[Word],
    ) -> Result<Vec<Word>> {
        cobar::differential_profile(model, self, words)
    }

    pub(crate) fn restrict(self, model: &mut Model, words: &[Word]) -> Result<Vec<Word>> {
        let mut projected = Vec::new();
        for word in words {
            model.spend(cobar::word_cells(word))?;
            if word.iter().all(|factor| {
                factor
                    .iter()
                    .enumerate()
                    .all(|(index, &power)| self.allows_power(index, power))
            }) {
                projected.push(word.clone());
            }
        }
        Ok(projected)
    }

    pub(crate) fn allows_power(self, index: usize, power: usize) -> bool {
        match self {
            Self::Full => true,
            Self::An(n) => power <= exponent_bound(n, index),
        }
    }

    pub(crate) fn bounds(self, max_factor: usize) -> Option<Bounds> {
        let Self::An(n) = self else {
            return None;
        };
        let width = n.saturating_add(1);
        // The first forbidden power is xi_1^(2^(n+1)). Below that
        // degree the profile and the full algebra have identical bases.
        if width >= usize::BITS as usize || max_factor < 1usize << width {
            return None;
        }
        let mut variables = Vec::with_capacity(width);
        let mut max_degree = 0usize;
        // The xi_1 exponent range spans every generator weight, so the
        // allowed monomials realize each degree through the top degree.
        for index in 0..width {
            let weight = (1usize << (index + 1)) - 1;
            let exponent = exponent_bound(n, index);
            variables.push(Variable {
                weight,
                exponent,
                lower_capacity: max_degree,
            });
            max_degree = max_degree.saturating_add(weight.saturating_mul(exponent));
        }
        Some(Bounds {
            variables,
            max_degree,
        })
    }
}

fn exponent_bound(n: usize, index: usize) -> usize {
    let Some(bits) = n
        .checked_sub(index)
        .map(|difference| difference.saturating_add(1))
    else {
        return 0;
    };
    if bits >= usize::BITS as usize {
        usize::MAX
    } else {
        (1usize << bits) - 1
    }
}

struct Variable {
    weight: usize,
    exponent: usize,
    lower_capacity: usize,
}

pub(crate) struct Bounds {
    variables: Vec<Variable>,
    pub(crate) max_degree: usize,
}

impl Bounds {
    // Fill degrees in increasing lexicographic order, from the rightmost
    // factor, with each positive factor bounded by the profile's top degree.
    pub(crate) fn distribute(&self, degrees: &mut [usize], mut excess: usize) {
        for degree in degrees.iter_mut().rev() {
            let placed = excess.min(self.max_degree - 1);
            *degree = 1 + placed;
            excess -= placed;
        }
        debug_assert_eq!(excess, 0);
    }

    pub(crate) fn next_degrees(&self, degrees: &mut [usize]) -> bool {
        let mut excess = 0;
        for index in (0..degrees.len()).rev() {
            if degrees[index] < self.max_degree && excess > 0 {
                degrees[index] += 1;
                self.distribute(&mut degrees[index + 1..], excess - 1);
                return true;
            }
            excess += degrees[index] - 1;
        }
        false
    }

    pub(crate) fn monomials(
        &self,
        model: &mut Model,
        degree: usize,
        limit: usize,
        cells: &mut usize,
    ) -> Result<Vec<Monomial>> {
        let variables = self
            .variables
            .partition_point(|variable| variable.weight <= degree);
        let mut catalogue = Monomials {
            variables: &self.variables[..variables],
            exponents: vec![0; variables],
            terms: Vec::new(),
            limit,
            cells,
        };
        catalogue.enumerate(model, variables, degree)?;
        Ok(catalogue.terms)
    }
}

struct Monomials<'a> {
    variables: &'a [Variable],
    exponents: Vec<usize>,
    terms: Vec<Monomial>,
    limit: usize,
    cells: &'a mut usize,
}

impl Monomials<'_> {
    fn enumerate(
        &mut self,
        model: &mut Model,
        remaining_variables: usize,
        degree: usize,
    ) -> Result<()> {
        model.spend(1)?;
        let Some(index) = remaining_variables.checked_sub(1) else {
            debug_assert_eq!(degree, 0);
            if self.terms.len() == self.limit {
                return Err(cobar::basis_limit(self.limit));
            }
            model.spend(self.exponents.len())?;
            let last = self
                .exponents
                .iter()
                .rposition(|&power| power != 0)
                .unwrap();
            *self.cells += last + 2;
            cobar::check_storage(*self.cells)?;
            self.terms.push(Monomial::from(&self.exponents[..=last]));
            return Ok(());
        };

        let Variable {
            weight,
            exponent,
            lower_capacity,
        } = self.variables[index];
        // Every exponent outside this interval is impossible, regardless of
        // the choices for lower generators. The final generator is xi_1.
        let first = degree.saturating_sub(lower_capacity).div_ceil(weight);
        let last = exponent.min(degree / weight);
        for power in first..=last {
            self.exponents[index] = power;
            self.enumerate(model, index, degree - power * weight)?;
        }
        Ok(())
    }
}
