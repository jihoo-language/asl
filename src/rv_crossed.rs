use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod postnikov;

const MAX_ORDER: usize = 128;

pub(crate) enum Command {
    Cyclic(Label, usize),
    Group(Label, Vec<Vec<usize>>),
    Crossed {
        name: Label,
        h: Label,
        e: Label,
        boundary: Vec<usize>,
        action: Vec<Vec<usize>>,
    },
    Data(Label),
    QuotientMul(Label, usize, usize),
}

pub(crate) enum Expression {
    Element(Label, usize),
    Multiply(Box<[ParsedExpr; 2]>),
    Inverse(Box<ParsedExpr>),
    Boundary(Label, Box<ParsedExpr>),
    Act(Label, Box<[ParsedExpr; 2]>),
    KernelAct(Label, usize, Box<ParsedExpr>),
    Associator(Label, [usize; 3]),
}

pub(crate) struct Group {
    products: Vec<Vec<usize>>,
    inverses: Vec<usize>,
}

impl Group {
    pub(crate) fn cyclic(model: &mut Model, order: usize) -> Result<Self> {
        check_group_order(order)?;
        model.spend(order * order)?;
        let products = (0..order)
            .map(|a| (0..order).map(|b| (a + b) % order).collect())
            .collect();
        let inverses = (0..order).map(|a| (order - a) % order).collect();
        Ok(Self { products, inverses })
    }

    pub(crate) fn new(model: &mut Model, products: Vec<Vec<usize>>) -> Result<Self> {
        let n = products.len();
        check_group_order(n)?;
        model.spend(n * n)?;
        if products
            .iter()
            .any(|row| row.len() != n || row.iter().any(|&i| i >= n))
        {
            return Err(Error::plain(
                "finite group requires a square table of in-range elements",
            ));
        }
        if (0..n).any(|i| products[0][i] != i || products[i][0] != i) {
            return Err(Error::plain("finite group table must have identity 0"));
        }
        let inverses = products
            .iter()
            .enumerate()
            .map(|(i, row)| {
                row.iter()
                    .enumerate()
                    .position(|(j, &product)| product == 0 && products[j][i] == 0)
                    .ok_or_else(|| {
                        Error::plain(
                            "finite group table must give a two-sided inverse for every element",
                        )
                    })
            })
            .collect::<Result<Vec<_>>>()?;
        model.spend(n * n * n)?;
        for (a, row) in products.iter().enumerate() {
            for (b, &ab) in row.iter().enumerate() {
                for (c, &bc) in products[b].iter().enumerate() {
                    if products[ab][c] != products[a][bc] {
                        return Err(Error::plain(format!(
                            "finite group table is not associative at ({a},{b},{c})"
                        )));
                    }
                }
            }
        }
        Ok(Self { products, inverses })
    }

    pub(crate) fn order(&self) -> usize {
        self.inverses.len()
    }

    pub(crate) fn cells(&self) -> usize {
        self.order() * (self.order() + 4)
    }

    pub(crate) fn element(&self, i: usize) -> Result<usize> {
        if i >= self.order() {
            return Err(Error::plain("finite group element index is out of range"));
        }
        Ok(i)
    }

    pub(crate) fn mul(&self, a: usize, b: usize) -> usize {
        self.products[a][b]
    }

    pub(crate) fn inv(&self, a: usize) -> usize {
        self.inverses[a]
    }
}

fn check_group_order(n: usize) -> Result<()> {
    if !(1..=MAX_ORDER).contains(&n) {
        return Err(Error::plain("finite group order must be between 1 and 128"));
    }
    Ok(())
}

pub(crate) struct Crossed {
    pub(crate) h_name: Rc<str>,
    pub(crate) e_name: Rc<str>,
    h: Rc<Group>,
    e: Rc<Group>,
    boundary: Vec<usize>,
    action: Vec<Vec<usize>>,
    kernel: Vec<usize>,
    classes: Vec<usize>,
    standard: Choices,
}

pub(super) struct Choices {
    section: Vec<usize>,
    lifts: Vec<Vec<usize>>,
}

impl Crossed {
    pub(crate) fn new(
        model: &mut Model,
        groups: [(Rc<str>, Rc<Group>); 2],
        boundary: Vec<usize>,
        action: Vec<Vec<usize>>,
    ) -> Result<Self> {
        let [(h_name, h), (e_name, e)] = groups;
        verify(model, &h, &e, &boundary, &action)?;
        let kernel = (0..h.order()).filter(|&a| boundary[a] == 0).collect();
        let mut taken = vec![usize::MAX; e.order()];
        let mut image = vec![None; e.order()];
        for (a, &b) in boundary.iter().enumerate() {
            image[b].get_or_insert(a);
        }
        let representatives: Vec<_> = image
            .iter()
            .enumerate()
            .filter_map(|(i, entry)| entry.is_some().then_some(i))
            .collect();
        let mut section = Vec::new();
        model.spend(e.order() * e.order())?;
        for e0 in 0..e.order() {
            if taken[e0] != usize::MAX {
                continue;
            }
            let quotient = section.len();
            section.push(e0);
            for &x in &representatives {
                taken[e.mul(x, e0)] = quotient;
            }
        }
        let lifts = section
            .iter()
            .map(|&a| {
                section
                    .iter()
                    .map(|&b| {
                        let ab = e.mul(a, b);
                        let defect = e.mul(ab, e.inv(section[taken[ab]]));
                        image[defect].expect("section defect belongs to the image subgroup")
                    })
                    .collect()
            })
            .collect();
        Ok(Self {
            h_name,
            e_name,
            h,
            e,
            boundary,
            action,
            kernel,
            classes: taken,
            standard: Choices { section, lifts },
        })
    }

    pub(crate) fn cells(&self) -> usize {
        self.h.order() * (self.e.order() + 2)
            + self.e.order() * 4
            + self.order() * (self.order() + 4)
    }

    pub(crate) fn boundary(&self, a: usize) -> usize {
        self.boundary[a]
    }

    pub(crate) fn act(&self, g: usize, a: usize) -> usize {
        self.action[g][a]
    }

    pub(crate) fn kernel_act(&self, g: usize, a: usize) -> Result<usize> {
        self.object(g)?;
        if self.boundary(a) != 0 {
            return Err(Error::plain(
                "skeletal action requires an element of the crossed-module kernel",
            ));
        }
        Ok(self.act(self.standard.section[g], a))
    }

    pub(crate) fn data(&self, name: &str) -> String {
        format!(
            "rvCrossedData[{name}] = kernel={:?}, quotient_representatives={:?}",
            self.kernel, self.standard.section
        )
    }

    pub(crate) fn order(&self) -> usize {
        self.standard.section.len()
    }

    pub(crate) fn object(&self, g: usize) -> Result<()> {
        if g >= self.order() {
            return Err(Error::plain(
                "skeletal quotient object index is out of range",
            ));
        }
        Ok(())
    }

    pub(crate) fn quotient_mul(&self, g: usize, h: usize) -> usize {
        let standard = &self.standard.section;
        self.classes[self.e.mul(standard[g], standard[h])]
    }

    pub(crate) fn standard_associator(&self, [g, h, l]: [usize; 3]) -> Result<usize> {
        for x in [g, h, l] {
            self.object(x)?;
        }
        self.associator(&self.standard, [g, h, l])
    }

    fn check_choices(
        &self,
        model: &mut Model,
        section: Vec<usize>,
        lifts: Vec<Vec<usize>>,
    ) -> Result<Choices> {
        let q = self.order();
        if section.len() != q
            || section[0] != 0
            || section
                .iter()
                .enumerate()
                .any(|(g, &e)| self.classes.get(e) != Some(&g))
        {
            return Err(Error::plain(
                "Postnikov section must give a normalized representative of every quotient object",
            ));
        }
        if lifts.len() != q || lifts.iter().any(|row| row.len() != q) {
            return Err(Error::plain(
                "Postnikov lifts must form a square quotient table",
            ));
        }
        model.spend(q * q)?;
        for (g, row) in lifts.iter().enumerate() {
            for (h, &lift) in row.iter().enumerate() {
                let quotient = self.quotient_mul(g, h);
                let target = self.e.mul(section[g], section[h]);
                if ((g == 0 || h == 0) && lift != 0)
                    || self
                        .boundary
                        .get(lift)
                        .is_none_or(|&image| self.e.mul(image, section[quotient]) != target)
                {
                    return Err(Error::plain(format!(
                        "Postnikov lift at ({g},{h}) fails the normalized section equation"
                    )));
                }
            }
        }
        Ok(Choices { section, lifts })
    }

    fn associator(&self, choices: &Choices, [g, h, l]: [usize; 3]) -> Result<usize> {
        let gh = self.quotient_mul(g, h);
        let hl = self.quotient_mul(h, l);
        let lift = &choices.lifts;
        let mut value = self.act(choices.section[g], lift[h][l]);
        for factor in [lift[g][hl], self.h.inv(lift[gh][l]), self.h.inv(lift[g][h])] {
            value = self.h.mul(value, factor);
        }
        if self.boundary(value) != 0 {
            return Err(Error::plain("crossed-module associator left the kernel"));
        }
        Ok(value)
    }
}

fn verify(
    model: &mut Model,
    h: &Group,
    e: &Group,
    boundary: &[usize],
    action: &[Vec<usize>],
) -> Result<()> {
    let (n, m) = (h.order(), e.order());
    if boundary.len() != n || boundary[0] != 0 || boundary.iter().any(|&g| g >= m) {
        return Err(Error::plain(
            "crossed-module boundary must map H to E and preserve identity",
        ));
    }
    if action.len() != m
        || action
            .iter()
            .any(|row| row.len() != n || row.iter().any(|&x| x >= n))
    {
        return Err(Error::plain(
            "crossed-module action needs one in-range row for each E element",
        ));
    }
    model.spend(n * n + m * n * n + m * m * n + n * m + n * n)?;
    for a in 0..n {
        for b in 0..n {
            if boundary[h.mul(a, b)] != e.mul(boundary[a], boundary[b]) {
                return Err(Error::plain(
                    "crossed-module boundary must be a group homomorphism",
                ));
            }
        }
    }
    if action[0].iter().copied().enumerate().any(|(i, x)| i != x) {
        return Err(Error::plain(
            "crossed-module identity action must be trivial",
        ));
    }
    for (g, row) in action.iter().enumerate() {
        let mut touched = vec![false; n];
        for &x in row {
            if std::mem::replace(&mut touched[x], true) {
                return Err(Error::plain(
                    "crossed-module action rows must be automorphisms",
                ));
            }
        }
        for (a, &a_image) in row.iter().enumerate() {
            for (b, &b_image) in row.iter().enumerate() {
                if row[h.mul(a, b)] != h.mul(a_image, b_image) {
                    return Err(Error::plain(
                        "crossed-module action must preserve multiplication",
                    ));
                }
            }
            let conjugate = e.mul(e.mul(g, boundary[a]), e.inv(g));
            if boundary[a_image] != conjugate {
                return Err(Error::plain(
                    "first Peiffer identity (boundary equivariance) fails",
                ));
            }
        }
        for (g2, second) in action.iter().enumerate() {
            let product = &action[e.mul(g, g2)];
            if second.iter().zip(product).any(|(&x, &y)| row[x] != y) {
                return Err(Error::plain(
                    "crossed-module action must respect multiplication in E",
                ));
            }
        }
    }
    for (a, &a_image) in boundary.iter().enumerate() {
        for (b, &action_image) in action[a_image].iter().enumerate() {
            if action_image != h.mul(h.mul(a, b), h.inv(a)) {
                return Err(Error::plain(
                    "second Peiffer identity (boundary acts by conjugation) fails",
                ));
            }
        }
    }
    Ok(())
}
