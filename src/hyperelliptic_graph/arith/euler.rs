use crate::error::{Error, Result};
use crate::hyperelliptic::weil_polynomial;
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_unram::Ring;

use super::{Arithmetic, Chart};

impl Chart {
    fn infinity_correction(&self, degree: usize) -> i128 {
        if self.polynomial.len().is_multiple_of(2) {
            0
        } else if self.nonsquare && !degree.is_multiple_of(2) {
            1
        } else {
            -1
        }
    }
}

fn counted_traces(
    model: &mut Model,
    prime: Prime,
    degree: usize,
    selected: &[&Chart],
) -> Result<Vec<i128>> {
    if selected.is_empty() {
        return Ok(Vec::new());
    }
    let field = counted_field(model, prime, degree)?;
    let polynomials: Vec<_> = selected
        .iter()
        .map(|chart| chart.polynomial.as_slice())
        .collect();
    let traces = field.residue_polynomial_traces(model, &polynomials)?;
    // Add the projective points at infinity to the affine character count.
    Ok(traces
        .into_iter()
        .zip(selected)
        .map(|(trace, chart)| trace + chart.infinity_correction(degree))
        .collect())
}

pub(super) struct Local {
    toric: usize,
    root_number: i32,
    graph: Vec<i128>,
    abelian: Vec<i128>,
    full: Vec<i128>,
    components: Vec<Component>,
}

struct Component {
    traces: Vec<i128>,
    factor: Vec<i128>,
}

impl Local {
    pub(super) fn compute(model: &mut Model, curve: &Arithmetic) -> Result<Self> {
        let prime = curve.cover.grade.prime;
        let toric = curve.cover.cycles.len();
        let (vertices, edges) = fixed_counts(model, curve, 1)?;
        // Trace additivity on the graph incidence complex gives the sum
        // of the eigenvalues on H_1 after subtracting H_0.
        let trace = 1 + edges as i128 - vertices as i128;
        let split_twice = toric as i128 + trace;
        let split = usize::try_from(split_twice / 2)
            .ok()
            .filter(|&split| split <= toric && split_twice == 2 * split as i128)
            .ok_or_else(|| Error::plain("invalid Frobenius trace on hyperelliptic cycles"))?;
        let mut graph = vec![1];
        for j in 0..toric {
            let sign = if j < split { -1 } else { 1 };
            graph = multiply(model, &graph, &[1, sign])?;
        }

        let charts = &curve.charts;
        let mut traces = vec![Vec::new(); charts.len()];
        let max_genus = charts.iter().filter_map(|c| c.genus).max().unwrap_or(0);
        for degree in 1..=max_genus {
            let selected: Vec<_> = charts
                .iter()
                .enumerate()
                .filter(|(_, chart)| chart.genus.is_some_and(|g| degree <= g))
                .collect();
            let chosen: Vec<_> = selected.iter().map(|(_, chart)| *chart).collect();
            for ((i, _), trace) in selected
                .into_iter()
                .zip(counted_traces(model, prime, degree, &chosen)?)
            {
                traces[i].push(trace);
            }
        }
        let mut components = Vec::with_capacity(charts.len());
        let mut abelian = vec![1];
        for sums in traces {
            let factor = weil_polynomial(model, u64::from(prime.get()), &sums)?;
            abelian = multiply(model, &abelian, &factor)?;
            components.push(Component {
                traces: sums,
                factor,
            });
        }
        let full = multiply(model, &graph, &abelian)?;
        let expected = 2 * curve.cover.grade.genus() - toric;
        if full.len() != 1 + expected || full.last() == Some(&0) {
            return Err(Error::plain(
                "semistable hyperelliptic Euler factor has an incorrect degree",
            ));
        }
        let root_number = if split.is_multiple_of(2) { 1 } else { -1 };
        let cells = 8
            + 2 * (graph.len() + abelian.len() + full.len())
            + components
                .iter()
                .map(|c| 4 + 2 * (c.traces.len() + c.factor.len()))
                .sum::<usize>();
        model.allocate(cells)?;
        Ok(Self {
            toric,
            root_number,
            graph,
            abelian,
            full,
            components,
        })
    }

    pub(super) fn report(&self, model: &mut Model) -> Result<String> {
        model.spend(self.graph.len() + self.abelian.len() + self.full.len())?;
        Ok(format!(
            "rvHyperLocal = graph={:?}, abelian={:?}, polynomial={:?}, conductor={}, root_number={}",
            self.graph, self.abelian, self.full, self.toric, self.root_number,
        ))
    }

    pub(super) fn factors(&self, model: &mut Model, charts: &[Chart]) -> Result<String> {
        model.spend(
            charts.iter().map(Chart::cells).sum::<usize>()
                + self
                    .components
                    .iter()
                    .map(|c| 1 + c.traces.len() + c.factor.len())
                    .sum::<usize>(),
        )?;
        let data: Vec<_> = charts
            .iter()
            .zip(&self.components)
            .enumerate()
            .map(|(i, (chart, component))| {
                (
                    i,
                    chart.genus,
                    &chart.polynomial,
                    &component.traces,
                    &component.factor,
                )
            })
            .collect();
        Ok(format!("rvHyperLocalFactors = {data:?}"))
    }
}

pub(super) fn nodal(model: &mut Model, curve: &Arithmetic, degree: usize) -> Result<String> {
    let prime = curve.cover.grade.prime;
    let q = field_size(prime, degree)?;
    let base_points = q.checked_add(1).ok_or_else(range)?;
    let (vertices, edges) = fixed_counts(model, curve, degree)?;
    let positive: Vec<_> = curve
        .charts
        .iter()
        .enumerate()
        .filter(|(_, chart)| chart.genus.is_some_and(|g| g > 0))
        .collect();
    let chosen: Vec<_> = positive.iter().map(|(_, chart)| *chart).collect();
    let mut traces = vec![0; curve.charts.len()];
    for ((i, _), trace) in positive
        .into_iter()
        .zip(counted_traces(model, prime, degree, &chosen)?)
    {
        traces[i] = trace;
    }
    let mut points = 0i128;
    for (i, vertex) in curve.cover.vertices.iter().enumerate() {
        model.spend(1)?;
        if degree.is_multiple_of(2) || curve.vertices[i] == i {
            // A fixed geometrically rational component has q+1 points.
            let component_points = base_points
                .checked_sub(traces[vertex.base])
                .ok_or_else(range)?;
            points = points.checked_add(component_points).ok_or_else(range)?;
        }
    }
    // The two branches of a fixed node remain individually fixed because
    // the graph action preserves its two distinct base endpoints.
    points = points.checked_sub(edges as i128).ok_or_else(range)?;
    let trace = q
        .checked_mul(vertices as i128)
        .and_then(|x| x.checked_add(1))
        .and_then(|x| x.checked_sub(points))
        .ok_or_else(range)?;
    Ok(format!(
        "rvHyperNodalPoints[{degree}] = {points}, fixed_vertices={vertices}, fixed_edges={edges}, trace={trace}"
    ))
}

fn field_size(prime: Prime, degree: usize) -> Result<i128> {
    if degree == 0 {
        return Err(Error::plain(
            "hyperelliptic point-count extension degree must be positive",
        ));
    }
    let degree = u32::try_from(degree).map_err(|_| range())?;
    i128::from(prime.get())
        .checked_pow(degree)
        .ok_or_else(range)
}

fn counted_field(model: &mut Model, prime: Prime, degree: usize) -> Result<Ring> {
    let count = usize::try_from(field_size(prime, degree)?).map_err(|_| range())?;
    // Charge an initial field scan before constructing a polynomial modulus.
    model.spend(count)?;
    Ring::field_of_degree(model, prime, degree)
}

fn fixed_counts(model: &mut Model, curve: &Arithmetic, degree: usize) -> Result<(usize, usize)> {
    let permanent = degree.is_multiple_of(2);
    let mut count_fixed = |permutation: &[usize]| -> Result<usize> {
        model.spend(permutation.len())?;
        Ok(if permanent {
            permutation.len()
        } else {
            permutation
                .iter()
                .enumerate()
                .filter(|&(i, &j)| i == j)
                .count()
        })
    };
    Ok((count_fixed(&curve.vertices)?, count_fixed(&curve.edges)?))
}

fn multiply(model: &mut Model, lhs: &[i128], rhs: &[i128]) -> Result<Vec<i128>> {
    model.spend(lhs.len() * rhs.len())?;
    let mut result = vec![0i128; lhs.len() + rhs.len() - 1];
    for (i, &a) in lhs.iter().enumerate() {
        for (j, &b) in rhs.iter().enumerate() {
            let value = a.checked_mul(b).ok_or_else(range)?;
            result[i + j] = result[i + j].checked_add(value).ok_or_else(range)?;
        }
    }
    Ok(result)
}

fn range() -> Error {
    Error::plain("semistable hyperelliptic count or coefficient exceeds signed 128-bit range")
}
