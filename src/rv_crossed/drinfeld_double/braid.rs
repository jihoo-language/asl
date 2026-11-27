use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_lens;
use crate::syntax::{Label, ParsedExpr};

use super::representations::Image;
use super::{negate_phase, Double, Monomial};

const MAX_WORD: usize = 4096;
pub(crate) type Colors = [usize; 3];

pub(crate) enum Command {
    Representation(Label, usize),
    MatrixColumn(Label, [usize; 3]),
    Image(usize, Box<ParsedExpr>),
}

pub(crate) enum Expression {
    Word(Label, Colors, Vec<i128>),
    Inverse(Box<ParsedExpr>),
    Compose(Box<[ParsedExpr; 2]>),
    Trace(Box<ParsedExpr>),
    Entry([usize; 2], Box<ParsedExpr>),
    Associator(Label, Colors),
    Hopf(Label, [usize; 2]),
    MatrixEntry(Label, [usize; 4]),
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Operator(pub(crate) Vec<Image>);

impl Operator {
    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(2 + 2 * self.0.len())?;
        Ok(Value::DoubleBraid(Rc::new(self)))
    }

    pub(crate) fn image(&self, column: usize) -> Result<Image> {
        self.0
            .get(column)
            .copied()
            .ok_or_else(|| Error::plain("double braid basis index is out of range"))
    }

    pub(crate) fn compose(
        &self,
        model: &mut Model,
        other: &Self,
        conductor: usize,
    ) -> Result<Self> {
        model.spend(self.0.len())?;
        let data = self
            .0
            .iter()
            .map(|a| {
                let b = other.0[a.target];
                Image {
                    target: b.target,
                    phase: (a.phase + b.phase) % conductor,
                }
            })
            .collect();
        Ok(Self(data))
    }

    pub(crate) fn inverse(&self, model: &mut Model, conductor: usize) -> Result<Self> {
        model.spend(self.0.len())?;
        let mut inverse = vec![
            Image {
                target: 0,
                phase: 0
            };
            self.0.len()
        ];
        for (i, a) in self.0.iter().enumerate() {
            inverse[a.target] = Image {
                target: i,
                phase: negate_phase(conductor, a.phase),
            };
        }
        Ok(Self(inverse))
    }

    pub(crate) fn trace(&self, model: &mut Model, double: &Double) -> Result<Value> {
        model.spend(self.0.len())?;
        let mut counts = vec![0; double.conductor];
        for (i, image) in self.0.iter().enumerate() {
            if i == image.target {
                counts[image.phase] += 1;
            }
        }
        double.rational(model, counts, 1)
    }

    pub(crate) fn entry(
        &self,
        model: &mut Model,
        double: &Double,
        row: usize,
        column: usize,
    ) -> Result<Value> {
        self.image(row)?;
        let image = self.image(column)?;
        if image.target == row {
            double.unit_phase(model, image.phase)
        } else {
            Ok(rv_lens::zero())
        }
    }
}

impl Double {
    pub(crate) fn representation(&self, name: &str, sector: usize) -> Result<String> {
        self.check_index(sector)?;
        Ok(self.representations[sector].data(name, sector))
    }

    pub(crate) fn matrix_column(&self, [sector, x, column]: [usize; 3]) -> Result<Image> {
        let simple = self.simple(sector)?;
        if x >= self.group.order() || column >= simple.dimension {
            return Err(Error::plain(
                "projective matrix element or basis index is out of range",
            ));
        }
        Ok(self.representations[sector].action[x][column])
    }

    pub(crate) fn matrix_entry(
        &self,
        model: &mut Model,
        [sector, x, row, column]: [usize; 4],
    ) -> Result<Value> {
        self.matrix_column([sector, x, row])?;
        let image = self.matrix_column([sector, x, column])?;
        if image.target == row {
            self.unit_phase(model, image.phase)
        } else {
            Ok(rv_lens::zero())
        }
    }

    fn unit_phase(&self, model: &mut Model, phase: usize) -> Result<Value> {
        self.value(model, Monomial { phase, weight: 1 }, 1)
    }

    fn associator_phase(&self, [a, b, c]: Colors) -> usize {
        let [g, h, j] = [a, b, c].map(|i| self.simples[i].flux);
        let n = self.group.order();
        negate_phase(self.conductor, self.cocycle[(g * n + h) * n + j])
    }

    pub(crate) fn associator(&self, model: &mut Model, colors: Colors) -> Result<Value> {
        self.check_colors(colors)?;
        self.unit_phase(model, self.associator_phase(colors))
    }

    pub(crate) fn hopf(&self, model: &mut Model, [a, b]: [usize; 2]) -> Result<Value> {
        let (x, y) = (self.simple(a)?, self.simple(b)?);
        let (first, last) = (
            &self.representations[a].action[y.flux],
            &self.representations[b].action[x.flux],
        );
        model.spend(first.len() * last.len())?;
        let mut counts = vec![0; self.conductor];
        for (i, u) in first.iter().enumerate() {
            for (j, v) in last.iter().enumerate() {
                if u.target == i && v.target == j {
                    counts[(u.phase + v.phase) % self.conductor] += 1;
                }
            }
        }
        self.rational(model, counts, self.group.order() as u128)
    }

    pub(crate) fn check_colors(&self, colors: Colors) -> Result<()> {
        for i in colors {
            self.check_index(i)?;
        }
        Ok(())
    }

    pub(crate) fn braid(
        &self,
        model: &mut Model,
        mut colors: Colors,
        word: &[i128],
    ) -> Result<(Colors, Operator)> {
        self.check_colors(colors)?;
        if word.len() > MAX_WORD {
            return Err(Error::plain(
                "double braid word exceeds its length limit of 4096",
            ));
        }
        let mut dimensions = colors.map(|i| self.simples[i].dimension);
        let size = dimensions.iter().product();
        let mut images: Vec<_> = (0..size).map(|target| Image { target, phase: 0 }).collect();
        for &crossing in word {
            let i = match crossing {
                1 | -1 => 0,
                2 | -2 => 1,
                _ => {
                    return Err(Error::plain(
                        "double braid generators must be 1, -1, 2 or -2",
                    ));
                }
            };
            let [a, b] = [colors[i], colors[i + 1]];
            let [g, h] = [a, b].map(|j| self.simples[j].flux);
            let mut swapped = colors;
            swapped.swap(i, i + 1);
            let mut next_dimensions = dimensions;
            next_dimensions.swap(i, i + 1);
            let associator = if i == 0 {
                0
            } else {
                (self.associator_phase(colors)
                    + negate_phase(self.conductor, self.associator_phase(swapped)))
                    % self.conductor
            };
            model.spend(1 + size * (1 + dimensions[i]))?;
            for image in &mut images {
                let mut coords = decode(dimensions, image.target);
                let action = if crossing > 0 {
                    let action = self.representations[b].action[g][coords[i + 1]];
                    coords[i + 1] = action.target;
                    action
                } else {
                    let action =
                        self.representations[a].inverse_image(h, coords[i], self.conductor);
                    coords[i] = action.target;
                    action
                };
                coords.swap(i, i + 1);
                image.target = encode(next_dimensions, coords);
                image.phase = (image.phase + action.phase + associator) % self.conductor;
            }
            colors = swapped;
            dimensions = next_dimensions;
        }
        Ok((colors, Operator(images)))
    }
}

fn encode([_, d1, d2]: [usize; 3], [a, b, c]: [usize; 3]) -> usize {
    (a * d1 + b) * d2 + c
}

fn decode([_, d1, d2]: [usize; 3], col: usize) -> [usize; 3] {
    [col / (d1 * d2), (col / d2) % d1, col % d2]
}
