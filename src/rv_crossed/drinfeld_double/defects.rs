use std::rc::Rc;

use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::qz::Scalar;
use crate::model::{Model, Value};
use crate::rv_cyclotomic;
use crate::rv_lens;
use crate::syntax::{Label, ParsedExpr};

use super::boundaries::Boundary;
use super::representations::{Image, Irrep};
use super::{characters, negate_phase, representations, Double, Quotient, Simple};

mod balanced;

pub(crate) enum Command {
    Declare(Label, [Label; 2]),
    Data(Label),
    Coset(Label, usize),
    Simple(Label, usize),
    Phase(usize, Box<ParsedExpr>),
    Image([usize; 2], Box<ParsedExpr>),
    Fuse(Box<[ParsedExpr; 2]>),
    Projector([usize; 3], Box<[ParsedExpr; 2]>),
}

pub(crate) enum Expression {
    Basis(Label, usize),
    Theta(Label, [usize; 3]),
    Entry([usize; 3], Box<ParsedExpr>),
    Character(usize, Box<ParsedExpr>),
    Opposite(Label, Box<ParsedExpr>),
}

pub(crate) struct Category {
    pub(crate) conductor: usize,
    pub(crate) boundaries: [Rc<Boundary>; 2],
    boundary_names: [Rc<str>; 2],
    intersection: Vec<usize>,
    group: Quotient,
    cosets: Vec<Coset>,
    coset_index: Vec<usize>,
    simples: Vec<Defect>,
}

struct Coset {
    support: Vec<usize>,
    theta: Vec<Vec<usize>>,
    radical: Vec<usize>,
    first: usize,
    count: usize,
    dimension: usize,
}

struct Defect {
    coset: usize,
    simple: Simple,
    representation: Irrep,
}

impl Category {
    pub(crate) fn new(
        model: &mut Model,
        double: &Double,
        boundary_names: [Rc<str>; 2],
        boundaries: [Rc<Boundary>; 2],
    ) -> Result<Self> {
        let h = boundaries[0].members()?;
        let k = boundaries[1].members()?;
        let annulus = boundaries[0].annulus(model, &boundaries[1])?;
        let intersection = h
            .iter()
            .filter(|x| k.binary_search(x).is_ok())
            .copied()
            .collect::<Vec<_>>();
        let group = double.group.restrict(model, &intersection)?;
        let (supports, coset_index) = double_cosets(model, &double.group, h, k)?;
        let exponent = characters::exponent(model, &group)?;
        let mut common = 1;
        let mut multipliers = Vec::with_capacity(supports.len());
        for support in &supports {
            let mut rows = Vec::new();
            for &x in &intersection {
                let mut row = Vec::new();
                for &y in &intersection {
                    let q = multiplier(model, double, &boundaries, support[0], x, y)?;
                    let denominator = q.root_order().ok_or_else(conductor_limit)?;
                    common = (common / rv_lens::gcd(common, denominator))
                        .checked_mul(denominator)
                        .ok_or_else(conductor_limit)?;
                    row.push(q);
                }
                rows.push(row);
            }
            multipliers.push(rows);
        }
        let conductor = common
            .checked_mul(exponent as u128)
            .and_then(|x| usize::try_from(x).ok())
            .ok_or_else(conductor_limit)?;
        rv_cyclotomic::check_order(conductor).map_err(|_| conductor_limit())?;
        let mut cosets = Vec::with_capacity(supports.len());
        let mut simples = Vec::new();
        for (support, multiplier) in supports.into_iter().zip(multipliers) {
            let g = support[0];
            let theta = multiplier
                .into_iter()
                .map(|row| {
                    row.into_iter()
                        .map(|q| {
                            q.root_index(conductor)
                                .expect("conductor divisible by denominator")
                        })
                        .collect()
                })
                .collect::<Vec<Vec<_>>>();
            characters::verify_multiplier(model, &group, &theta, conductor)?;
            let radical = characters::radical(model, &theta)?;
            let sectors = characters::extend(model, &group, g, &theta, &radical, conductor)?;
            let irreps =
                representations::induce_all(model, &group, &theta, &radical, &sectors, conductor)?;
            orthogonality(model, &sectors, &radical, group.order(), conductor)?;
            let dimension = sectors[0].dimension;
            let first = simples.len();
            let count = sectors.len();
            let index = cosets.len();
            model.allocate(
                support.len()
                    + theta.len() * (theta.len() + 3)
                    + radical.len()
                    + 9
                    + count * (2 * group.order() + 5),
            )?;
            cosets.push(Coset {
                support,
                theta,
                radical,
                first,
                count,
                dimension,
            });
            for (simple, representation) in sectors.into_iter().zip(irreps) {
                simples.push(Defect {
                    coset: index,
                    simple,
                    representation,
                });
            }
        }
        if simples.len() as i128 != annulus {
            return Err(Error::plain(format!(
                "defect-sector count {} disagrees with bulk annulus overlap {annulus}",
                simples.len()
            )));
        }
        model.allocate(intersection.len() * (intersection.len() + 5) + coset_index.len() + 20)?;
        Ok(Self {
            conductor,
            boundaries,
            boundary_names,
            intersection,
            group,
            cosets,
            coset_index,
            simples,
        })
    }

    pub(crate) fn count(&self) -> usize {
        self.simples.len()
    }

    pub(crate) fn check_index(&self, index: usize) -> Result<()> {
        if index >= self.count() {
            return Err(Error::plain(
                "double boundary defect simple index is out of range",
            ));
        }
        Ok(())
    }

    fn simple(&self, index: usize) -> Result<&Defect> {
        self.check_index(index)?;
        Ok(&self.simples[index])
    }

    fn coset(&self, index: usize) -> Result<&Coset> {
        self.cosets
            .get(index)
            .ok_or_else(|| Error::plain("double boundary defect coset index is out of range"))
    }

    fn local(&self, x: usize) -> Result<usize> {
        self.intersection
            .binary_search(&x)
            .map_err(|_| Error::plain("defect action is defined on the boundary intersection"))
    }

    pub(crate) fn data(&self, name: &str) -> String {
        let [a, b] = self.boundary_names.each_ref();
        let dimension_sum: usize = self.simples.iter().map(|d| d.simple.dimension.pow(2)).sum();
        format!(
            "rvDoubleDefectData[{name}] = boundaries=[{a},{b}], intersection={:?}, cosets={}, conductor={}, simples={}, sum_squared_fiber_dimensions={dimension_sum}, annulus_verified=true",
            self.intersection,self.cosets.len(),self.conductor,self.count()
        )
    }

    pub(crate) fn coset_report(&self, name: &str, index: usize) -> Result<String> {
        let coset = self.coset(index)?;
        let radical = coset
            .radical
            .iter()
            .map(|&x| self.intersection[x])
            .collect::<Vec<_>>();
        Ok(format!(
            "rvDoubleDefectCoset[{name},{index}] = representative={}, support={:?}, radical={radical:?}, first_simple={}, simples={}, fiber_dimension={}",
            coset.support[0],coset.support,coset.first,coset.count,coset.dimension
        ))
    }

    pub(crate) fn simple_report(&self, name: &str, index: usize) -> Result<String> {
        let defect = self.simple(index)?;
        let simple = &defect.simple;
        let phases = self
            .intersection
            .iter()
            .copied()
            .zip(&simple.phases)
            .filter_map(|(x, p)| p.map(|p| (x, p)))
            .collect::<Vec<_>>();
        let isotropic = defect
            .representation
            .isotropic()
            .iter()
            .map(|&x| self.intersection[x])
            .collect::<Vec<_>>();
        let reps = defect
            .representation
            .representatives()
            .iter()
            .map(|&x| self.intersection[x])
            .collect::<Vec<_>>();
        let dimension = simple.dimension;
        let support = self.cosets[defect.coset].support.len();
        let fp = if Rc::ptr_eq(&self.boundaries[0], &self.boundaries[1]) {
            format!(", fp_dimension={dimension}")
        } else {
            String::new()
        };
        Ok(format!(
            "rvDoubleDefectSimple[{name},{index}] = coset={}, representative={}, fiber_dimension={dimension}, support_size={support}, vector_space_dimension={}, radical_phases={phases:?}, isotropic={isotropic:?}, coset_representatives={reps:?}{fp}",
            defect.coset,simple.flux,support*dimension,
        ))
    }

    fn scalar(&self, model: &mut Model, phase: usize) -> Result<Scalar> {
        Scalar::literal(model, phase as i128, self.conductor as u128)
    }

    pub(crate) fn theta(&self, model: &mut Model, [coset, x, y]: [usize; 3]) -> Result<Scalar> {
        self.scalar(
            model,
            self.coset(coset)?.theta[self.local(x)?][self.local(y)?],
        )
    }

    pub(crate) fn phase(&self, index: usize, x: usize) -> Result<Option<usize>> {
        Ok(self.simple(index)?.simple.phases[self.local(x)?])
    }

    pub(crate) fn image(&self, index: usize, [x, column]: [usize; 2]) -> Result<Image> {
        let defect = self.simple(index)?;
        defect.representation.action[self.local(x)?]
            .get(column)
            .copied()
            .ok_or_else(|| Error::plain("defect representation basis index is out of range"))
    }

    pub(crate) fn entry(
        &self,
        model: &mut Model,
        index: usize,
        [x, row, column]: [usize; 3],
    ) -> Result<Value> {
        self.image(index, [x, row])?;
        let image = self.image(index, [x, column])?;
        if image.target == row {
            self.monomial(model, image.phase, 1)
        } else {
            Ok(rv_lens::zero())
        }
    }

    pub(crate) fn character(&self, model: &mut Model, index: usize, x: usize) -> Result<Value> {
        let Some(phase) = self.phase(index, x)? else {
            return Ok(rv_lens::zero());
        };
        self.monomial(model, phase, self.simple(index)?.simple.dimension as i128)
    }

    fn monomial(&self, model: &mut Model, phase: usize, weight: i128) -> Result<Value> {
        let mut polynomial = vec![0; phase + 1];
        polynomial[phase] = weight;
        let numerator = rv_cyclotomic::reduced(model, self.conductor, polynomial)?;
        let result = rv_lens::fraction(model, &numerator, 1)?;
        rv_lens::rational_value(model, result)
    }

    pub(crate) fn check_opposite(&self, other: &Self) -> Result<()> {
        if !Rc::ptr_eq(&self.boundaries[0], &other.boundaries[1])
            || !Rc::ptr_eq(&self.boundaries[1], &other.boundaries[0])
        {
            return Err(Error::plain(
                "opposite defect category must reverse the same two nominal boundaries",
            ));
        }
        Ok(())
    }

    pub(crate) fn opposite(
        &self,
        model: &mut Model,
        double: &Double,
        other: &Self,
        index: usize,
    ) -> Result<usize> {
        self.check_opposite(other)?;
        let defect = self.simple(index)?;
        let g = defect.simple.flux;
        let inverse = double.group.inv(g);
        let target = &other.cosets[other.coset_index[inverse]];
        let representative = target.support[0];
        let translate = double.group.mul(g, representative);
        let [h, k] = correction_pair(double, &self.boundaries, translate)?;
        let ih = double.group.mul(inverse, h);
        let shifts = self
            .intersection
            .iter()
            .map(|&x| {
                let first = tensor_scalar(model, double, x, g, inverse)?.negate();
                let left = commutator(model, &self.boundaries[0], x, h)?;
                let right = commutator(model, &self.boundaries[1], x, k)?;
                let step_h = tensor_scalar(model, double, x, inverse, h)?;
                // On the reversed bimodule, H translates by right action
                // and K by left action, so their tensor phases are ordered differently.
                let step_k = tensor_scalar(model, double, x, k, ih)?;
                add_terms(model, [first, left, right, step_h, step_k])
            })
            .collect::<Result<Vec<_>>>()?;
        let source = &self.cosets[defect.coset];
        for (x, &a) in shifts.iter().enumerate() {
            for (y, &b) in shifts.iter().enumerate() {
                let original = self.scalar(model, source.theta[x][y])?.negate();
                let xy = shifts[self.group.mul(x, y)].negate();
                let expected = add_terms(model, [original, a, b, xy])?;
                if expected.root_index(other.conductor) != Some(target.theta[x][y]) {
                    return Err(Error::plain(
                        "opposite boundary multiplier rephasing failed",
                    ));
                }
            }
        }
        let phases = defect
            .simple
            .phases
            .iter()
            .zip(shifts)
            .map(|(p, shift)| {
                p.map(|p| {
                    let value = self.scalar(model, p)?.negate().add(model, shift)?;
                    value
                        .root_index(other.conductor)
                        .ok_or_else(conductor_limit)
                })
                .transpose()
            })
            .collect::<Result<Vec<_>>>()?;
        (target.first..target.first + target.count)
            .find(|&i| {
                other.simples[i].simple.phases == phases
                    && other.simples[i].simple.dimension == defect.simple.dimension
            })
            .ok_or_else(|| Error::plain("opposite projective character did not match a simple"))
    }
}

fn correction_pair(
    double: &Double,
    boundaries: &[Rc<Boundary>; 2],
    target: usize,
) -> Result<[usize; 2]> {
    let left = boundaries[0].members()?;
    let right = boundaries[1].members()?;
    for &h in left {
        for &k in right {
            if double.group.mul(h, k) == target {
                return Ok([h, k]);
            }
        }
    }
    Err(Error::plain(
        "opposite double coset has no subgroup translation",
    ))
}

fn commutator(model: &mut Model, boundary: &Boundary, x: usize, h: usize) -> Result<Scalar> {
    boundary
        .cochain([x, h])?
        .add(model, boundary.cochain([h, x])?.negate())
}

fn tensor_scalar(
    model: &mut Model,
    double: &Double,
    g: usize,
    x: usize,
    y: usize,
) -> Result<Scalar> {
    Scalar::literal(
        model,
        double.tensor_phase(g, x, y) as i128,
        double.conductor as u128,
    )
}

fn multiplier(
    model: &mut Model,
    double: &Double,
    boundaries: &[Rc<Boundary>; 2],
    g: usize,
    x: usize,
    y: usize,
) -> Result<Scalar> {
    let bulk = tensor_scalar(model, double, g, x, y)?;
    add_terms(
        model,
        [
            boundaries[0].cochain([x, y])?,
            boundaries[1].cochain([x, y])?.negate(),
            bulk,
        ],
    )
}

fn add_terms(model: &mut Model, terms: impl IntoIterator<Item = Scalar>) -> Result<Scalar> {
    let mut value = Scalar::zero();
    for term in terms {
        if term != Scalar::zero() {
            value = value.add(model, term)?;
        }
    }
    Ok(value)
}

fn double_cosets(
    model: &mut Model,
    group: &Quotient,
    h: &[usize],
    k: &[usize],
) -> Result<(Vec<Vec<usize>>, Vec<usize>)> {
    model.spend(group.order() + h.len() * k.len())?;
    let mut hk = vec![false; group.order()];
    for &x in h {
        for &y in k {
            hk[group.mul(x, y)] = true;
        }
    }
    let hk = hk
        .into_iter()
        .enumerate()
        .filter_map(|(i, yes)| yes.then_some(i))
        .collect::<Vec<_>>();
    let mut attached = vec![usize::MAX; group.order()];
    let mut cosets = Vec::new();
    for g in 0..group.order() {
        if attached[g] != usize::MAX {
            continue;
        }
        model.spend(hk.len())?;
        let mut support = hk.iter().map(|&x| group.mul(g, x)).collect::<Vec<_>>();
        support.sort_unstable();
        for &x in &support {
            attached[x] = cosets.len();
        }
        cosets.push(support);
    }
    Ok((cosets, attached))
}

fn orthogonality(
    model: &mut Model,
    simples: &[Simple],
    radical: &[usize],
    n: usize,
    conductor: usize,
) -> Result<()> {
    if simples.iter().map(|s| s.dimension.pow(2)).sum::<usize>() != n {
        return Err(Error::plain(
            "defect dimensions do not square-sum to the intersection order",
        ));
    }
    for (i, a) in simples.iter().enumerate() {
        for (j, b) in simples.iter().enumerate().skip(i) {
            let mut polynomial = vec![0; conductor];
            model.spend(radical.len() + conductor)?;
            for &x in radical {
                let phase = (a.phases[x].unwrap() + negate_phase(conductor, b.phases[x].unwrap()))
                    % conductor;
                polynomial[phase] += (a.dimension * b.dimension) as i128;
            }
            let reduced = rv_cyclotomic::reduced(model, conductor, polynomial)?;
            if (i == j && reduced != [n as i128]) || (i != j && !reduced.is_empty()) {
                return Err(Error::plain(
                    "defect projective characters are not orthonormal",
                ));
            }
        }
    }
    Ok(())
}

fn conductor_limit() -> Error {
    Error::plain("boundary defect multiplier exceeds the cyclotomic conductor limit")
}
