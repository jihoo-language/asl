use std::rc::Rc;

use crate::complex::{for_each_boundary, OrderComplex, Simplex};
use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::{self, Basis, Coordinate, QuotientBasis};
use crate::set_family::Diagram;
use crate::syntax::{Label, ParsedExpr};

use super::postnikov::{Cochain, Postnikov};

pub(crate) enum Command {
    Declare(Label, Label, Label, Vec<(Label, Label, usize)>),
    Data(Label),
    Edge(Label, Label, Label),
    Obstruction(Label, Box<ParsedExpr>),
    Check(Label, Box<[ParsedExpr; 2]>),
    Coefficient(Vec<Label>, Box<ParsedExpr>),
    CyclePair(Box<ParsedExpr>, Vec<usize>),
}

pub(crate) enum Expression {
    Literal(Label, usize, Vec<usize>),
    Pullback(Label, Box<ParsedExpr>),
    Lift(Label, Box<ParsedExpr>),
    H2Basis(Label, usize),
    Coboundary(Box<ParsedExpr>),
    BoundaryWitness(Box<ParsedExpr>),
    LiftDifference(Label, Box<[ParsedExpr; 3]>),
}

pub(super) type Arrows = Vec<Vec<Option<usize>>>;
type Rows = Vec<Vec<Coordinate>>;

pub(crate) struct Flat {
    pub(crate) base: TypeId,
    pub(crate) post_name: Rc<str>,
    pub(crate) post: Rc<Postnikov>,
    arrows: Arrows,
    chains: [Vec<Simplex>; 5],
    rows: [Rows; 4],
    d0: Basis,
    h2: QuotientBasis,
    d2: Basis,
    betti: [usize; 4],
}

impl Flat {
    pub(crate) fn new(
        model: &mut Model,
        base: TypeId,
        post_name: Rc<str>,
        post: Rc<Postnikov>,
        transitions: Vec<(Label, Label, usize)>,
    ) -> Result<Self> {
        let site = model.poset(base).expect("named finite space");
        let (order, covers) = Diagram::site_orders(model, &site)?;
        let n = site.names.len();
        model.allocate(2 * n * n + 3 * n)?;
        let mut arrows = vec![vec![None; n]; n];
        for (source, target, g) in transitions {
            let x = model.point(base, &source)?;
            let y = model.point(base, &target)?;
            if !site.strict_less(x, y) {
                return Err(Error::new(
                    source.span,
                    "flat edge requires x < y in its space",
                ));
            }
            post.crossed
                .object(g)
                .map_err(|error| error.at(source.span))?;
            if arrows[x][y].replace(g).is_some() {
                return Err(Error::new(source.span, "duplicate flat edge label"));
            }
        }
        for &x in &order {
            for &y in &covers[x] {
                if arrows[x][y].is_none() {
                    return Err(Error::plain(format!(
                        "missing flat label on cover '{} < {}'",
                        site.names[x], site.names[y]
                    )));
                }
            }
        }
        // Every saturated path ends with a cover. Visiting intermediate
        // vertices in topological order lets us check all such products.
        for &x in &order {
            for &y in &order {
                if !site.strict_less(x, y) {
                    continue;
                }
                let first = arrows[x][y].expect("reachable through a cover");
                for &z in &covers[y] {
                    model.spend(1)?;
                    let next = arrows[y][z].expect("given cover");
                    let product = post.crossed.quotient_mul(first, next);
                    let value = arrows[x][z].get_or_insert(product);
                    if *value != product {
                        return Err(Error::plain(format!(
                            "flat labels disagree along '{} < {} < {}'",
                            site.names[x], site.names[y], site.names[z]
                        )));
                    }
                }
            }
        }
        let points: Vec<_> = (0..n).collect();
        let mut complex = OrderComplex::induced(model, &site, &points)?;
        let mut chains: [Vec<Simplex>; 5] = Default::default();
        chains[0] = complex.vertices();
        let mut rows: [Rows; 4] = Default::default();
        let mut all_columns = Vec::with_capacity(4);
        let p = post.prime;
        for (degree, d) in rows.iter_mut().enumerate() {
            chains[degree + 1] = complex.extend(model, &chains[degree])?;
            let (faces, columns) = twisted_differential(
                model,
                post.as_ref(),
                &arrows,
                &chains[degree],
                &chains[degree + 1],
            )?;
            *d = faces;
            all_columns.push(columns);
        }
        let [cols0, cols1, cols2, cols3]: [Rows; 4] =
            all_columns.try_into().expect("four cochain differentials");
        let mut d0 = Basis::default();
        let mut d1 = Basis::default();
        let mut d2 = Basis::default();
        let mut d3 = Basis::default();
        for (columns, basis) in [(cols0, &mut d0), (cols1, &mut d1), (cols2.clone(), &mut d2)] {
            for (i, column) in columns.into_iter().enumerate() {
                basis.insert_with_witness(model, p, column, vec![(i, 1)])?;
            }
        }
        for column in cols3 {
            d3.insert(model, p, column)?;
        }
        let ranks = [0, d0.rank(), d1.rank(), d2.rank(), d3.rank()];
        let mut betti = [0; 4];
        for (degree, value) in betti.iter_mut().enumerate() {
            *value = chains[degree]
                .len()
                .checked_sub(ranks[degree] + ranks[degree + 1])
                .ok_or_else(|| Error::plain("twisted cochain ranks are inconsistent"))?;
        }
        let mut cols2 = cols2.into_iter().enumerate().peekable();
        let h2 = QuotientBasis::build(model, p, d1, chains[2].len(), |_, index| {
            while cols2.peek().is_some_and(|(i, _)| *i < index) {
                cols2.next();
            }
            Ok(cols2.next().expect("degree-two cochain column").1)
        })?;
        debug_assert_eq!(h2.rank(), betti[2]);
        let chain_cells = chains
            .iter()
            .enumerate()
            .map(|(degree, layer)| (degree + 4) * layer.len())
            .sum::<usize>();
        let row_cells = rows
            .iter()
            .enumerate()
            .map(|(degree, d)| (2 * (degree + 2) + 3) * d.len())
            .sum::<usize>();
        model.allocate(chain_cells + row_cells + d0.cells() + d2.cells() + 10)?;
        Ok(Self {
            base,
            post_name,
            post,
            arrows,
            chains,
            rows,
            d0,
            h2,
            d2,
            betti,
        })
    }

    pub(crate) fn prime(&self) -> Prime {
        self.post.prime
    }

    pub(crate) fn size(&self, degree: usize) -> usize {
        self.chains[degree].len()
    }

    pub(crate) fn edge(&self, x: usize, y: usize) -> Result<usize> {
        if x == y {
            return Ok(0);
        }
        self.arrows[x][y].ok_or_else(|| Error::plain("flat edge requires x <= y in its space"))
    }

    pub(crate) fn literal(
        &self,
        model: &mut Model,
        degree: usize,
        dense: Vec<usize>,
    ) -> Result<Cochain> {
        self.coordinates(model, degree, dense).map(Cochain)
    }

    fn coordinates(
        &self,
        model: &mut Model,
        degree: usize,
        dense: Vec<usize>,
    ) -> Result<Vec<Coordinate>> {
        let size = self.size(degree);
        let p = self.prime().get() as usize;
        if dense.len() != size || dense.iter().any(|&x| x >= p) {
            return Err(Error::plain(format!(
                "degree-{degree} requires {size} residues from 0 through {}",
                p - 1
            )));
        }
        model.spend(size)?;
        Ok(dense
            .into_iter()
            .enumerate()
            .filter(|&(_, x)| x != 0)
            .map(|(i, x)| (i, x as u32))
            .collect())
    }

    pub(crate) fn coefficient(
        &self,
        model: &mut Model,
        input: &Cochain,
        points: &[usize],
    ) -> Result<u32> {
        let chains = &self.chains[points.len() - 1];
        model.spend(points.len() * (1 + chains.len().max(1).ilog2() as usize))?;
        let index = chains
            .binary_search_by(|chain| chain.as_slice().cmp(points))
            .map_err(|_| Error::plain("flat cochain coefficient needs a strict chain"))?;
        Ok(input.at(index))
    }

    pub(crate) fn pullback(
        &self,
        model: &mut Model,
        degree: usize,
        input: &Cochain,
    ) -> Result<Cochain> {
        pullback_layer(
            model,
            &self.post,
            &self.arrows,
            &self.chains[degree],
            degree,
            input,
        )
    }

    pub(crate) fn coboundary(
        &self,
        model: &mut Model,
        degree: usize,
        input: &Cochain,
    ) -> Result<Cochain> {
        let mut result = Vec::new();
        let p = self.prime();
        if input.0.is_empty() {
            return Ok(Cochain(result));
        }
        let rows = &self.rows[degree];
        model.spend(rows.len() * (degree + 2) * (1 + input.0.len().ilog2() as usize))?;
        for (i, row) in rows.iter().enumerate() {
            let value = row
                .iter()
                .fold(0, |sum, &(j, value)| p.add(sum, p.mul(value, input.at(j))));
            if value != 0 {
                result.push((i, value));
            }
        }
        Ok(Cochain(result))
    }

    pub(crate) fn boundary_witness(
        &self,
        model: &mut Model,
        degree: usize,
        input: &Cochain,
    ) -> Result<Cochain> {
        let p = self.prime();
        let result = match degree {
            2 => self.h2.boundary_witness(model, p, input.0.clone())?,
            _ => {
                let basis = if degree == 1 { &self.d0 } else { &self.d2 };
                let (rest, witness) = basis.reduce_with_witness(model, p, input.0.clone())?;
                rest.is_empty().then_some(witness)
            }
        };
        result
            .map(Cochain)
            .ok_or_else(|| Error::plain("local cochain is not a coboundary"))
    }

    pub(crate) fn cycle_pair(
        &self,
        model: &mut Model,
        degree: usize,
        input: &Cochain,
        dense: Vec<usize>,
    ) -> Result<u32> {
        let cycle = self.coordinates(model, degree, dense)?;
        let p = self.prime();
        if degree > 0 {
            model.spend(self.size(degree - 1) + cycle.len() * (degree + 1))?;
            let mut boundary = vec![0; self.size(degree - 1)];
            // Dual boundary is the transpose of the local coboundary.
            for &(i, value) in &cycle {
                for &(j, coefficient) in &self.rows[degree - 1][i] {
                    boundary[j] = p.add(boundary[j], p.mul(value, coefficient));
                }
            }
            if boundary.into_iter().any(|x| x != 0) {
                return Err(Error::plain("supplied dual local chain is not a cycle"));
            }
        }
        model.spend(cycle.len() * (1 + input.0.len().max(1).ilog2() as usize))?;
        Ok(cycle
            .into_iter()
            .fold(0, |sum, (i, value)| p.add(sum, p.mul(value, input.at(i)))))
    }

    pub(crate) fn h2_basis(&self, model: &mut Model, index: usize) -> Result<Cochain> {
        if index >= self.h2.rank() {
            return Err(Error::plain("local H^2 basis index is out of range"));
        }
        Ok(Cochain(self.h2.representative(
            model,
            self.prime(),
            &[(index, 1)],
        )?))
    }

    fn obstruction(&self, model: &mut Model, k: &Cochain) -> Result<Cochain> {
        self.post.require_cocycle(model, k)?;
        self.pullback(model, 3, k)
    }

    pub(crate) fn lift(&self, model: &mut Model, k: &Cochain) -> Result<Cochain> {
        let obstruction = self.obstruction(model, k)?;
        let (rest, witness) = self
            .d2
            .reduce_with_witness(model, self.prime(), obstruction.0)?;
        if !rest.is_empty() {
            return Err(Error::plain(
                "flat lift is obstructed by a nonzero pulled-back Postnikov class",
            ));
        }
        Ok(Cochain(witness))
    }

    pub(crate) fn check_lift(&self, model: &mut Model, k: &Cochain, b: &Cochain) -> Result<()> {
        let obstruction = self.obstruction(model, k)?;
        self.match_lift(model, &obstruction, b)
    }

    fn match_lift(&self, model: &mut Model, obstruction: &Cochain, b: &Cochain) -> Result<()> {
        if self.coboundary(model, 2, b)? != *obstruction {
            return Err(Error::plain(
                "supplied lift does not have the pulled-back associator as coboundary",
            ));
        }
        Ok(())
    }

    pub(crate) fn lift_difference(
        &self,
        model: &mut Model,
        k: &Cochain,
        a: &Cochain,
        b: &Cochain,
    ) -> Result<Cochain> {
        let obstruction = self.obstruction(model, k)?;
        self.match_lift(model, &obstruction, a)?;
        self.match_lift(model, &obstruction, b)?;
        let p = self.prime();
        let diff = prime_linear::add_scaled(model, p, &b.0, &a.0, p.negate(1))?;
        let class = self
            .h2
            .project(model, p, diff)?
            .ok_or_else(|| Error::plain("difference of checked lifts must be closed"))?;
        Ok(Cochain(self.h2.representative(model, p, &class)?))
    }

    pub(crate) fn data(&self, name: &str) -> String {
        let chain_counts = self.chains.each_ref().map(Vec::len);
        let p = self.prime().get();
        format!(
            "rvFlatData[{name}] = prime={p}, chains={chain_counts:?}, betti={:?}",
            self.betti
        )
    }

    pub(crate) fn obstruction_report(
        &self,
        model: &mut Model,
        name: &str,
        input: &Cochain,
    ) -> Result<String> {
        let obstruction = self.obstruction(model, input)?;
        let prime = self.prime();
        let pulled_support = obstruction.0.len();
        let rest = self.d2.reduce(model, prime, obstruction.0)?;
        let obstructed = !rest.is_empty();
        let (classes, card) = if obstructed {
            ("0".to_string(), "0".to_string())
        } else {
            let p = prime.get();
            let [h0, h1, h2, _] = self.betti;
            (
                power_of_prime(p, h2 as i64),
                power_of_prime(p, h0 as i64 - h1 as i64 + h2 as i64),
            )
        };
        Ok(format!(
            "rvFlatObstruction[{name}] = obstructed={obstructed}, pulled_support={pulled_support}, reduced={rest:?}, lift_classes={classes}, fixed_cardinality={card}"
        ))
    }
}

pub(super) fn twisted_differential(
    model: &mut Model,
    post: &Postnikov,
    arrows: &Arrows,
    source: &[Simplex],
    target: &[Simplex],
) -> Result<(Rows, Rows)> {
    let mut rows = Vec::with_capacity(target.len());
    let mut columns = vec![Vec::new(); source.len()];
    let minus = post.prime.negate(1);
    for_each_boundary(model, source, target, |_, index, faces| {
        let simplex = &target[index];
        let g = arrows[simplex[0]][simplex[1]].expect("strict edge");
        let mut row = Vec::with_capacity(faces.len());
        for (i, coordinate) in faces.into_iter().enumerate() {
            let coefficient = match i {
                0 => post.action[g],
                i if i.is_multiple_of(2) => 1,
                _ => minus,
            };
            row.push((coordinate, coefficient));
            columns[coordinate].push((index, coefficient));
        }
        rows.push(row);
        Ok(())
    })?;
    Ok((rows, columns))
}

pub(super) fn pullback_layer(
    model: &mut Model,
    post: &Postnikov,
    arrows: &Arrows,
    chains: &[Simplex],
    degree: usize,
    input: &Cochain,
) -> Result<Cochain> {
    model.spend(chains.len() * degree)?;
    let mut result = Vec::new();
    let mut labels = [0; 3];
    for (i, chain) in chains.iter().enumerate() {
        for (j, edge) in chain.windows(2).enumerate() {
            labels[j] = arrows[edge[0]][edge[1]].expect("strict edge");
        }
        let value = post.coefficient(model, input, &labels[..degree])?;
        if value != 0 {
            result.push((i, value));
        }
    }
    Ok(Cochain(result))
}

fn power_of_prime(p: u32, exponent: i64) -> String {
    let n = exponent.unsigned_abs() as u32;
    let value = (p as u128)
        .checked_pow(n)
        .map_or_else(|| format!("{p}^{n}"), |x| x.to_string());
    if exponent < 0 {
        format!("1/{value}")
    } else {
        value
    }
}

pub(crate) fn persist(model: &mut Model, c: Cochain) -> Result<Value> {
    model.allocate(2 * c.0.len() + 1)?;
    Ok(Value::FlatCochain(Rc::new(c)))
}
