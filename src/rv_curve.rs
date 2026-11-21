use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BTreeSet, BinaryHeap};

use crate::binary_field::BinaryField;
use crate::error::{Error, Result};
use crate::linear::Basis;
use crate::model::Model;
use crate::q_koszul::{self, Indices};
use crate::rv_ideal::{self, Ideal};
use crate::rv_poly::{self, Monomial, Polynomial};
use crate::{rv_rees, rv_scheme};

type Row = Vec<Polynomial>;

pub(crate) struct Curve {
    pub(crate) selected: Indices,
    pub(crate) orders: Vec<usize>,
    degree: usize,
    columns: usize,
    independent: Vec<usize>,
    rows: Vec<Row>,
}

impl Curve {
    pub(crate) fn build(
        model: &mut Model,
        selected: Indices,
        forms: Vec<(usize, Polynomial)>,
    ) -> Result<Self> {
        if selected.len() != 2 {
            return Err(Error::plain(
                "rvcurve requires exactly two parameter variables",
            ));
        }
        let Some(&(degree, _)) = forms.first() else {
            return Err(Error::plain(
                "rvcurve requires at least one coordinate form",
            ));
        };
        if degree == 0 || forms.iter().any(|(other, _)| *other != degree) {
            return Err(Error::plain(
                "rvcurve coordinate forms must have a common positive degree",
            ));
        }
        let columns = forms.len();
        model.spend(columns)?;
        rv_ideal::check_storage(columns)?;
        let mut cells = 0;
        let mut affine = Vec::with_capacity(columns);
        for (_, form) in forms {
            model.spend(rv_poly::cells(&form))?;
            cells += 1 + 2 * form.len();
            rv_ideal::check_storage(cells)?;
            // Exponent vectors are (s,t). On s=1 the exponent of t reverses
            // the order of the homogeneous support.
            affine.push(
                form.iter()
                    .rev()
                    .map(|term| [term[1]].into())
                    .collect::<Polynomial>(),
            );
        }
        if !affine
            .iter()
            .any(|p| p.last().is_some_and(|term| term[0] == degree))
        {
            return Err(Error::plain(
                "rvcurve has a common projective basepoint at [0:1]",
            ));
        }
        if !is_one(&affine_gcd(model, &affine)?) {
            return Err(Error::plain(
                "rvcurve has common projective basepoints on the chart s != 0",
            ));
        }

        let independent = coordinate_basis(model, &affine)?;
        let (orders, rows) = order_sequence(model, &affine, independent.len(), None)?;
        model.allocate(
            1 + orders.len()
                + independent.len()
                + rows.iter().map(|row| row_cells(row)).sum::<usize>(),
        )?;
        Ok(Self {
            selected,
            orders,
            degree,
            columns,
            independent,
            rows,
        })
    }

    pub(crate) fn frobenius_orders(
        &self,
        model: &mut Model,
        exponent: usize,
    ) -> Result<Vec<usize>> {
        let (_, orders, _) = self.frobenius_rows(model, exponent)?;
        Ok(orders)
    }

    pub(crate) fn frobenius_wronskian(
        &self,
        model: &mut Model,
        exponent: usize,
    ) -> Result<(usize, Polynomial)> {
        let (q, orders, rows) = self.frobenius_rows(model, exponent)?;
        let sum = orders.iter().map(|&order| order as u128).sum::<u128>();
        let total = (q as u128 + orders.len() as u128) * self.degree as u128;
        let degree = total
            .checked_sub(2 * sum)
            .ok_or_else(|| Error::plain("Frobenius Wronskian degree is negative"))?;
        let degree = usize::try_from(degree)
            .map_err(|_| Error::plain("Frobenius Wronskian degree is out of range"))?;
        let determinant = rv_poly::determinant_polynomials(model, rows)?;
        let polynomial = homogenize(
            model,
            &determinant,
            degree,
            "Frobenius Wronskian exceeds its homogeneous degree",
        )?;
        Ok((degree, polynomial))
    }

    fn frobenius_rows(
        &self,
        model: &mut Model,
        exponent: usize,
    ) -> Result<(usize, Vec<usize>, Vec<Row>)> {
        if !(1..usize::BITS as usize).contains(&exponent) {
            return Err(Error::plain("Frobenius extension exponent must be positive and within the machine exponent range"));
        }
        let q = 1usize << exponent;
        self.degree
            .checked_mul(q)
            .ok_or_else(|| Error::plain("Frobenius row degree is out of range"))?;
        let series: Row = self
            .independent
            .iter()
            .map(|&column| self.rows[0][column].clone())
            .collect();
        let mut cells = 1 + series.len();
        let mut top = Vec::with_capacity(series.len());
        for polynomial in &series {
            let cost = rv_poly::cells(polynomial);
            cells += cost;
            rv_ideal::check_storage(cells)?;
            model.spend(cost)?;
            top.push(polynomial.iter().map(|term| [term[0] * q].into()).collect());
        }
        let (orders, rows) = order_sequence(model, &series, series.len(), Some(&top))?;
        let mut matrix = Vec::with_capacity(series.len());
        matrix.push(top);
        matrix.extend(rows);
        Ok((q, orders, matrix))
    }

    fn plucker_degree(&self, k: usize) -> Result<usize> {
        let orders = self.orders.get(..=k).ok_or_else(|| {
            Error::plain(format!(
                "osculating order {k} is out of range for coordinate rank {}",
                self.orders.len()
            ))
        })?;
        let sum = orders.iter().map(|&order| order as u128).sum::<u128>();
        let total = (orders.len() as u128) * (self.degree as u128);
        let result = total
            .checked_sub(2 * sum)
            .ok_or_else(|| Error::plain("osculating Plucker degree is negative"))?;
        usize::try_from(result)
            .map_err(|_| Error::plain("osculating Plucker degree is out of range"))
    }

    fn affine_minor(&self, model: &mut Model, columns: &[usize]) -> Result<Polynomial> {
        if columns.len() == 1 {
            return Ok(self.rows[0][columns[0]].clone());
        }
        let size = columns.len();
        rv_ideal::check_storage(size.saturating_mul(size))?;
        model.spend(size * size)?;
        let matrix = self.rows[..size]
            .iter()
            .map(|row| columns.iter().map(|&column| row[column].clone()).collect())
            .collect();
        rv_poly::determinant_polynomials(model, matrix)
    }

    pub(crate) fn minor(
        &self,
        model: &mut Model,
        k: usize,
        mut columns: Vec<usize>,
    ) -> Result<(usize, Polynomial)> {
        let degree = self.plucker_degree(k)?;
        if columns.len() != k + 1 {
            return Err(Error::plain(
                "rvOscMinor requires exactly k+1 coordinate positions",
            ));
        }
        model.spend(columns.len() * (1 + columns.len().ilog2() as usize))?;
        columns.sort_unstable();
        if columns.last().is_some_and(|&index| index >= self.columns) {
            return Err(Error::plain(
                "rvOscMinor coordinate position is out of range",
            ));
        }
        if columns.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(Error::plain(
                "rvOscMinor requires distinct coordinate positions",
            ));
        }
        let affine = self.affine_minor(model, &columns)?;
        let polynomial = homogenize(
            model,
            &affine,
            degree,
            "affine minor exceeds its Plucker degree",
        )?;
        Ok((degree, polynomial))
    }

    pub(crate) fn image(&self, model: &mut Model, k: usize, selected: Indices) -> Result<Ideal> {
        self.plucker_degree(k)?;
        if !choose_matches(self.columns, k + 1, selected.len()) {
            return Err(Error::plain(
                "rvosculating requires one new variable per lexicographic Plucker coordinate",
            ));
        }
        if selected
            .iter()
            .any(|index| self.selected.binary_search(index).is_ok())
        {
            return Err(Error::plain(
                "rvosculating requires fresh coordinate labels",
            ));
        }
        let mut columns: Vec<_> = (0..=k).collect();
        let mut polynomials = Vec::with_capacity(selected.len());
        let mut cells = 0;
        loop {
            model.spend(columns.len())?;
            let polynomial = self.affine_minor(model, &columns)?;
            cells += rv_poly::cells(&polynomial);
            rv_ideal::check_storage(cells)?;
            polynomials.push(polynomial);
            if !next_columns(&mut columns, self.columns) {
                break;
            }
        }

        normalize_map(model, &mut polynomials)?;
        // The common homogeneous degree makes affine dehomogenization
        // injective on each homogeneous relation after substitution.
        // Track the relation degree with a separate scaling parameter.
        let kernel = rv_rees::graph(model, std::iter::empty(), 1, &polynomials)?;
        let kernel = rv_ideal::eliminate_prefix(model, kernel, 1)?;
        let width = selected.len();
        let kernel = rv_scheme::saturate_coordinates(model, &kernel, width, 0..width)?;
        Ideal::derived(model, selected, kernel)
    }
}

pub(crate) fn order_at(
    model: &mut Model,
    polynomial: &Polynomial,
    modulus: usize,
    point: &[usize],
) -> Result<Option<usize>> {
    let field = BinaryField::new(model, modulus)?;
    let &[s, t] = point else {
        return Err(Error::plain("rvOrdAt requires two field coordinates"));
    };
    field.check_element(s)?;
    field.check_element(t)?;
    if s == 0 && t == 0 {
        return Err(Error::plain(
            "the zero vector does not define a projective point",
        ));
    }
    let Some(first) = polynomial.first() else {
        return Ok(None);
    };
    if s == 0 {
        return Ok(Some(first[0]));
    }
    if t == 0 {
        return Ok(Some(polynomial.last().unwrap()[1]));
    }

    let parameter = field.multiply(t, field.inverse(model, s)?);
    model.spend(rv_poly::cells(polynomial))?;
    let mut terms = polynomial
        .iter()
        .rev()
        .map(|term| {
            let power = term[1];
            Ok((power, field.power(model, parameter, power as u128)?))
        })
        .collect::<Result<Vec<_>>>()?;

    // Rescale the local coordinate so the point is X=1. If the largest
    // power of two present is b, write g(X)=A(X)+X^b B(X). At X=1+Z,
    // g is (A+B)(1+Z)+Z^b B(1+Z). The first summand has degree below b.
    let mut order = 0;
    while terms.len() > 1 {
        model.spend(1 + terms.len())?;
        let bit = 1usize << terms.last().unwrap().0.ilog2();
        let split = terms.partition_point(|&(power, _)| power < bit);
        let mut upper = terms.split_off(split);
        for (power, _) in &mut upper {
            *power -= bit;
        }
        let low = add_bands(&terms, &upper);
        terms = if low.is_empty() {
            order += bit;
            upper
        } else {
            low
        };
    }
    Ok(Some(order))
}

fn add_bands(mut left: &[(usize, usize)], mut right: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut sum = Vec::with_capacity(left.len() + right.len());
    while let (Some((&a, tail_a)), Some((&b, tail_b))) = (left.split_first(), right.split_first()) {
        match a.0.cmp(&b.0) {
            Ordering::Less => {
                sum.push(a);
                left = tail_a;
            }
            Ordering::Greater => {
                sum.push(b);
                right = tail_b;
            }
            Ordering::Equal => {
                let coefficient = a.1 ^ b.1;
                if coefficient != 0 {
                    sum.push((a.0, coefficient));
                }
                left = tail_a;
                right = tail_b;
            }
        }
    }
    sum.extend_from_slice(left);
    sum.extend_from_slice(right);
    sum
}

fn homogenize(
    model: &mut Model,
    affine: &Polynomial,
    degree: usize,
    degree_error: &str,
) -> Result<Polynomial> {
    model.spend(rv_poly::cells(affine))?;
    affine
        .iter()
        .rev()
        .map(|term| {
            let t = term[0];
            let s = degree
                .checked_sub(t)
                .ok_or_else(|| Error::plain(degree_error))?;
            Ok([s, t].into())
        })
        .collect()
}

fn row_cells(row: &[Polynomial]) -> usize {
    1 + row.iter().map(|p| rv_poly::cells(p)).sum::<usize>()
}

fn is_one(polynomial: &[Monomial]) -> bool {
    polynomial.len() == 1 && polynomial[0][0] == 0
}

fn affine_gcd(model: &mut Model, polynomials: &[Polynomial]) -> Result<Polynomial> {
    let mut gcd: Polynomial = [].into();
    for polynomial in polynomials.iter().filter(|p| !p.is_empty()) {
        let mut other = polynomial.clone();
        while !other.is_empty() {
            let remainder = binary_remainder(model, &gcd, &other)?;
            gcd = other;
            other = remainder;
        }
        if is_one(&gcd) {
            break;
        }
    }
    Ok(gcd)
}

fn binary_remainder(
    model: &mut Model,
    numerator: &[Monomial],
    denominator: &[Monomial],
) -> Result<Polynomial> {
    let leading = denominator.last().expect("nonzero divisor")[0];
    if leading == 0 {
        return Ok([].into());
    }
    model.spend(1 + numerator.len() * (1 + numerator.len().max(1).ilog2() as usize))?;
    let mut powers: BTreeSet<_> = numerator.iter().map(|term| term[0]).collect();
    while let Some(&highest) = powers.last() {
        let Some(shift) = highest.checked_sub(leading) else {
            break;
        };
        powers.pop_last();
        model.spend(1 + denominator.len() * (1 + powers.len().max(1).ilog2() as usize))?;
        for term in &denominator[..denominator.len() - 1] {
            let power = term[0] + shift;
            if !powers.remove(&power) {
                powers.insert(power);
                rv_ideal::check_storage(1 + 2 * powers.len())?;
            }
        }
    }
    Ok(powers.into_iter().map(|power| [power].into()).collect())
}

fn coordinate_basis(model: &mut Model, polynomials: &[Polynomial]) -> Result<Vec<usize>> {
    // Translation f(T+X) is invertible over F2(T), so the full Hasse matrix
    // has rank equal to the F2 span of the coordinate polynomials.
    let mut span = Basis::default();
    let mut independent = Vec::new();
    for (column, polynomial) in polynomials.iter().enumerate() {
        let support: Vec<_> = polynomial.iter().map(|term| term[0]).collect();
        if span.insert(model, support)? {
            independent.push(column);
        }
    }
    Ok(independent)
}

fn order_sequence(
    model: &mut Model,
    polynomials: &[Polynomial],
    rank: usize,
    leading: Option<&[Polynomial]>,
) -> Result<(Vec<usize>, Vec<Row>)> {
    let needed = rank - usize::from(leading.is_some());
    if needed == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    let columns = polynomials.len();
    let mut exponents = BTreeMap::<usize, Vec<usize>>::new();
    let mut cells = 0;
    for (column, polynomial) in polynomials.iter().enumerate() {
        for term in polynomial.iter() {
            model.spend(1 + exponents.len().max(1).ilog2() as usize)?;
            let coordinates = exponents.entry(term[0]).or_default();
            cells += if coordinates.is_empty() { 4 } else { 1 };
            rv_ideal::check_storage(cells)?;
            coordinates.push(column);
        }
    }

    // Merge the increasing submask lists of the exponents actually present.
    // Lucas parity excludes every other derivative without visiting it.
    let mut queue: BinaryHeap<_> = exponents
        .keys()
        .map(|&power| Reverse((0usize, power)))
        .collect();
    let index_cells = cells + 2 * queue.len();
    rv_ideal::check_storage(index_cells)?;
    model.spend(1 + queue.len())?;
    let mut echelon = RationalBasis::default();
    let mut rows = Vec::with_capacity(needed);
    let mut orders = Vec::with_capacity(needed);
    let mut retained = 0;
    if let Some(row) = leading {
        retained = row_cells(row);
        if !echelon.insert(model, row, index_cells + retained)? {
            return Err(Error::plain("Frobenius row vanishes"));
        }
    }
    while let Some(&Reverse((order, _))) = queue.peek() {
        model.spend(1 + columns)?;
        let mut row = vec![Vec::new(); columns];
        let mut cells = 1 + columns;
        while let Some(&Reverse((candidate, power))) = queue.peek() {
            if candidate != order {
                break;
            }
            queue.pop();
            let coordinates = &exponents[&power];
            model.spend(1 + queue.len().max(1).ilog2() as usize + coordinates.len())?;
            cells += 2 * coordinates.len();
            rv_ideal::check_storage(index_cells + cells + retained + echelon.cells)?;
            let term: Monomial = [power - order].into();
            for &column in coordinates {
                row[column].push(term.clone());
            }
            let next = order.wrapping_sub(power) & power;
            if next > order {
                queue.push(Reverse((next, power)));
            }
        }
        let row: Row = row.into_iter().map(Into::into).collect();
        if echelon.insert(model, &row, index_cells + retained + cells)? {
            orders.push(order);
            retained += row_cells(&row);
            rv_ideal::check_storage(index_cells + retained + echelon.cells)?;
            rows.push(row);
            if rows.len() == needed {
                return Ok((orders, rows));
            }
        }
    }
    Err(Error::plain(
        "Hasse order search did not reach the coordinate rank",
    ))
}

#[derive(Default)]
struct RationalBasis {
    rows: Vec<(usize, Row)>,
    cells: usize,
}

impl RationalBasis {
    fn insert(&mut self, model: &mut Model, row: &[Polynomial], retained: usize) -> Result<bool> {
        model.spend(1 + row.len())?;
        let mut row = row.to_vec();
        let mut cells = row_cells(&row);
        let mut previous: Option<&[Monomial]> = None;
        // Bareiss updates keep entries as minors rather than introducing
        // nested rational denominators in the elimination over F2(T).
        for (column, pivot_row) in &self.rows {
            let pivot = &pivot_row[*column];
            let factor = row[*column].clone();
            for (index, current) in row.iter_mut().enumerate() {
                model.spend(1)?;
                let value = if index == *column {
                    [].into()
                } else {
                    rv_poly::bareiss_update(
                        model,
                        pivot,
                        current,
                        &factor,
                        &pivot_row[index],
                        previous,
                    )?
                };
                cells = cells - rv_poly::cells(current) + rv_poly::cells(&value);
                rv_ideal::check_storage(cells + self.cells + retained)?;
                *current = value;
            }
            previous = Some(pivot);
        }
        let Some(column) = row.iter().position(|p| !p.is_empty()) else {
            return Ok(false);
        };
        self.cells += cells + 1;
        rv_ideal::check_storage(self.cells + retained)?;
        self.rows.push((column, row));
        Ok(true)
    }
}

fn choose_matches(n: usize, k: usize, expected: usize) -> bool {
    let space = k.min(n - k);
    let mut count = 1u128;
    for i in 1..=space {
        count = count * (n - i + 1) as u128 / i as u128;
        if count > expected as u128 {
            return false;
        }
    }
    count == expected as u128
}

fn next_columns(columns: &mut [usize], bound: usize) -> bool {
    let length = columns.len();
    for i in (0..length).rev() {
        if columns[i] < bound - (length - i) {
            let first = columns[i] + 1;
            for (offset, column) in columns[i..].iter_mut().enumerate() {
                *column = first + offset;
            }
            return true;
        }
    }
    false
}

fn normalize_map(model: &mut Model, polynomials: &mut [Polynomial]) -> Result<()> {
    let gcd = affine_gcd(model, polynomials)?;
    if gcd.is_empty() {
        return Err(Error::plain("the selected osculating minors all vanish"));
    }
    let mut cells = polynomials.iter().map(|p| rv_poly::cells(p)).sum::<usize>();
    if !is_one(&gcd) {
        for polynomial in polynomials.iter_mut() {
            let reduced = rv_poly::exact_divide(model, polynomial, &gcd)?;
            cells = cells - rv_poly::cells(polynomial) + rv_poly::cells(&reduced);
            rv_ideal::check_storage(cells)?;
            *polynomial = reduced;
        }
    }
    model.spend(cells)?;
    let step = polynomials
        .iter()
        .flat_map(|p| p.iter())
        .fold(0, |step, term| q_koszul::gcd(step, term[0]));
    // Substitution T -> T^step is injective. Removing this common exponent
    // scaling preserves the ideal even for an inseparable parametrization.
    if step > 1 {
        for polynomial in polynomials {
            *polynomial = polynomial
                .iter()
                .map(|term| [term[0] / step].into())
                .collect();
        }
    }
    Ok(())
}
