use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_cyclotomic;
use crate::rv_lens;

use super::super::representations::{self, IntertwinerOrbit};
use super::{add_terms, commutator, conductor_limit, tensor_scalar, Category, Double, Image};

struct Balanced {
    coset: usize,
    conductor: usize,
    scale: usize,
    tensor: Vec<Vec<Image>>,
}

impl Category {
    fn balanced(&self, model: &mut Model, double: &Double, [a, b]: [usize; 2]) -> Result<Balanced> {
        if !std::rc::Rc::ptr_eq(&self.boundaries[0], &self.boundaries[1]) {
            return Err(Error::plain(
                "balanced defect self-fusion requires one nominal boundary",
            ));
        }
        let (a, b) = (self.simple(a)?, self.simple(b)?);
        let (g, h) = (a.simple.flux, b.simple.flux);
        let product = double.group.mul(g, h);
        let coset = self.coset_index[product];
        let block = &self.cosets[coset];
        let t = double
            .group
            .mul(double.group.inv(product), block.support[0]);
        self.local(t)?;
        // Balancing identifies all decompositions of the grade gh with
        // the chosen M_g ⊗ N_h. Transport on this fiber is the tensor
        // operator with α_x(g,h). Right translation by t then adds the
        // character of t and α_x(gh,t) to reach the catalog's representative.
        let mut conductor = self.conductor as u128;
        let mut shifts = Vec::with_capacity(self.intersection.len());
        for &x in &self.intersection {
            let tensor = tensor_scalar(model, double, x, g, h)?;
            let translation = commutator(model, &self.boundaries[0], x, t)?;
            let mixed = tensor_scalar(model, double, x, product, t)?;
            let phase = add_terms(model, [tensor, translation, mixed])?;
            let denominator = phase.root_order().ok_or_else(conductor_limit)?;
            conductor = (conductor / rv_lens::gcd(conductor, denominator))
                .checked_mul(denominator)
                .ok_or_else(conductor_limit)?;
            shifts.push(phase);
        }
        let conductor = usize::try_from(conductor).map_err(|_| conductor_limit())?;
        rv_cyclotomic::check_order(conductor).map_err(|_| conductor_limit())?;
        let scale = conductor / self.conductor;
        let phases = shifts
            .iter()
            .map(|p| p.root_index(conductor).expect("common conductor"))
            .collect::<Vec<_>>();
        let left = rescale(&a.representation.action, scale);
        let right = rescale(&b.representation.action, scale);
        let tensor = representations::tensor_action(model, &left, &right, &phases, conductor)?;
        let theta = block
            .theta
            .iter()
            .map(|row| row.iter().map(|x| scale * x).collect())
            .collect::<Vec<_>>();
        representations::verify_action(model, &self.group, &theta, &tensor, conductor)?;
        Ok(Balanced {
            coset,
            conductor,
            scale,
            tensor,
        })
    }

    pub(crate) fn self_fusion(
        &self,
        model: &mut Model,
        double: &Double,
        indices: [usize; 2],
    ) -> Result<String> {
        let product = self.balanced(model, double, indices)?;
        let block = &self.cosets[product.coset];
        let mut channels = Vec::new();
        let mut dimension = 0;
        for c in block.first..block.first + block.count {
            let result = self.simple(c)?;
            let multiplicity = product.orbits(model, &result.representation.action)?.len();
            if multiplicity != 0 {
                dimension += multiplicity * result.simple.dimension;
                channels.push((c, multiplicity));
            }
        }
        if dimension != product.tensor[0].len() {
            return Err(Error::plain(
                "balanced defect fusion fails fiber-dimension conservation",
            ));
        }
        let [a, b] = indices;
        Ok(format!(
            "rvDoubleDefectFuse[{a},{b}] = coset={}, fiber_dimension={dimension}, sectors={channels:?}",product.coset
        ))
    }

    pub(crate) fn projector(
        &self,
        model: &mut Model,
        double: &Double,
        indices: [usize; 2],
        [c, row, col]: [usize; 3],
    ) -> Result<String> {
        let product = self.balanced(model, double, indices)?;
        let defect = self.simple(c)?;
        let size = product.tensor[0].len() * defect.simple.dimension;
        if row >= size || col >= size {
            return Err(Error::plain(
                "defect fusion projector coordinate is out of range",
            ));
        }
        let mut coefficient = rv_lens::Rational::zero();
        if defect.coset == product.coset {
            for orbit in product.orbits(model, &defect.representation.action)? {
                let find = |j| orbit.iter().find_map(|&(i, p)| (i == j).then_some(p));
                if let (Some(a), Some(b)) = (find(row), find(col)) {
                    let phase = (a + product.conductor - b) % product.conductor;
                    let mut monomial = vec![0; phase + 1];
                    monomial[phase] = 1;
                    let reduced = rv_cyclotomic::reduced(model, product.conductor, monomial)?;
                    coefficient = rv_lens::fraction(model, &reduced, orbit.len() as u128)?;
                    break;
                }
            }
        }
        Ok(format!(
            "rvDoubleDefectFusionProjector[{c},{row},{col}] = RVCycloQ({}){{{:?} / {}}}",
            product.conductor, coefficient.numerator, coefficient.denominator,
        ))
    }
}

impl Balanced {
    fn orbits(&self, model: &mut Model, action: &[Vec<Image>]) -> Result<Vec<IntertwinerOrbit>> {
        let target = rescale(action, self.scale);
        representations::intertwiner_orbits(model, &self.tensor, &target, self.conductor)
    }
}

fn rescale(action: &[Vec<Image>], scale: usize) -> Vec<Vec<Image>> {
    action
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| Image {
                    target: x.target,
                    phase: x.phase * scale,
                })
                .collect()
        })
        .collect()
}
