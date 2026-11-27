use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::sheaf::{self, Cochain, Cochains, Generator, Matrix, Open, Place, Sheaf, Term};
use crate::sheaf_cone::{self, Apex};
use crate::sheaf_hyper::{Complex, NamedSheaf};
use crate::sheaf_image::Site;
use crate::sheaf_leray::Input;
use crate::syntax::{Label, MAX_IDENTIFIER_BYTES};

use super::Disk;

struct DiskApex {
    alpha0: Vec<Vec<Term>>,
    alpha1: Vec<Vec<Term>>,
    internal: Rc<Matrix>,
}

impl Apex for DiskApex {
    fn height(&self) -> Result<usize> {
        Ok(2)
    }

    fn size(&self, degree: usize) -> usize {
        match degree {
            0 => self.alpha0.len(),
            1 => self.alpha1.len(),
            _ => 0,
        }
    }

    fn column(&self, degree: usize, index: usize) -> &[Term] {
        match degree {
            0 => &self.alpha0[index],
            1 => &self.alpha1[index],
            _ => unreachable!("two-term apex"),
        }
    }

    fn differential(&self, degree: usize) -> Option<Rc<Matrix>> {
        (degree == 0).then(|| self.internal.clone())
    }
}

impl DiskApex {
    fn build(model: &mut Model, disk: &Disk, link: &Sheaf, edge: [usize; 2]) -> Result<Self> {
        let mut alpha0 = Vec::with_capacity(disk.nearby);
        let mut alpha1 = Vec::<Vec<Term>>::with_capacity(disk.vanishing);
        let mut cells = 0;
        let chains: Vec<Rc<[usize]>> = (0..4).map(|point| [point].into()).collect();
        let edge: Rc<[usize]> = edge.into();
        for index in 0..disk.nearby {
            model.spend(1 + chains.len())?;
            let terms: Vec<_> = chains
                .iter()
                .map(|chain| (term(0, chain, index), 1))
                .collect();
            cells += 1 + 4 * terms.len();
            sheaf::check_cells(cells)?;
            alpha0.push(terms);
        }
        for index in 0..disk.vanishing {
            let column = disk.v.column(index);
            model.spend(1 + column.len())?;
            cells += 1 + 4 * column.len();
            sheaf::check_cells(cells)?;
            alpha1.push(
                column
                    .iter()
                    .map(|&(i, x)| (term(0, &edge, i), x))
                    .collect(),
            );
        }
        for (i, column) in alpha0.iter().enumerate() {
            let differential = link.differential_support(model, column)?;
            let mut image = Cochain::default();
            for &(j, x) in disk.u.column(i) {
                for (term, value) in &alpha1[j] {
                    image.add(model, disk.prime, term.clone(), disk.prime.mul(x, *value))?;
                }
            }
            if differential != image.into_terms() {
                return Err(Error::plain(
                    "disk attachment must commute with the actual Roos differential",
                ));
            }
        }
        Ok(Self {
            alpha0,
            alpha1,
            internal: disk.u.clone(),
        })
    }
}

pub(crate) struct Realization {
    pub(crate) name: Rc<str>,
    pub(crate) nearby: NamedSheaf,
    pub(crate) ambient: Rc<Complex>,
    pub(crate) costalk_name: Rc<str>,
    pub(crate) costalk: Rc<Complex>,
    disk_prime: crate::odd_cobar::Prime,
    apex_chain: Rc<[usize]>,
    edge_chain: Rc<[usize]>,
    edge_indices: Vec<usize>,
}

impl Realization {
    pub(crate) fn build(
        model: &mut Model,
        name: &Label,
        disk: &Disk,
        site: Site,
        apex: usize,
        roles: [usize; 4],
    ) -> Result<Self> {
        let [a, b, c, d] = roles;
        let link = model.poset(site.source).expect("named open site");
        if link.names.len() != 4
            || roles.iter().any(|&i| i >= 4)
            || roles
                .iter()
                .enumerate()
                .any(|(i, x)| roles[..i].contains(x))
        {
            return Err(Error::plain(
                "disk realization requires four distinct points a,b,c,d in its link",
            ));
        }
        model.spend(16)?;
        for x in 0..4 {
            for y in 0..4 {
                if link.strict_less(x, y) != ((x == a || x == b) && (y == c || y == d)) {
                    return Err(Error::plain(
                        "disk realization requires exactly the four links from a,b to c,d",
                    ));
                }
            }
        }
        let name_b = derived_name(name, "nearby")?;
        let name_c = derived_name(name, "costalk")?;
        model.check_new_sheaf(&name_b)?;
        model.check_new_sheaf_complex(name)?;
        model.check_new_sheaf_complex(&name_c)?;
        let id = Rc::new(Matrix::identity(model, disk.nearby)?);
        let diagram = Sheaf::from_cover_maps(
            model,
            site.source,
            disk.prime,
            vec![disk.nearby; 4],
            |_, x, y| {
                Ok(if (x, y) == (b, d) {
                    disk.t.clone()
                } else {
                    id.clone()
                })
            },
        )?;
        let nearby = NamedSheaf {
            name: name_b.text.into(),
            sheaf: Rc::new(diagram),
        };
        let input = Input::Sheaf(nearby.clone());
        let attachment = DiskApex::build(model, disk, &nearby.sheaf, [b, d])?;
        let ambient = Rc::new(sheaf_cone::strict_ambient(
            model,
            name,
            &input,
            &site,
            apex,
            &attachment,
        )?);
        let costalk = Rc::new(disk.complex(model, &name_c, true)?);
        let b_image = site.image(b);
        let mask: Vec<_> = (0..4)
            .map(|u| site.contains_in_future(b_image, u))
            .collect();
        let open = Open::new(model, &nearby.sheaf, &mask)?;
        let layer = open.layer(model, 1, MAX_PIVOT_CELLS)?;
        let edge: Rc<[usize]> = [b, d].into();
        let edge_indices = (0..disk.nearby)
            .map(|j| layer.index(model, &term(0, &edge, j)))
            .collect::<Result<Vec<_>>>()?;
        model.allocate(9 + edge_indices.len())?;
        Ok(Self {
            name: name.text.as_str().into(),
            nearby,
            ambient,
            costalk_name: name_c.text.into(),
            costalk,
            disk_prime: disk.prime,
            apex_chain: [apex].into(),
            edge_chain: [apex, b_image].into(),
            edge_indices,
        })
    }

    pub(crate) fn lift(&self, model: &mut Model, degree: usize, raw: &[Term]) -> Result<Value> {
        let mut result = Cochain::default();
        for (generator, value) in raw {
            match degree {
                1 => result.add(
                    model,
                    self.disk_prime,
                    term(1, &self.apex_chain, generator.index),
                    *value,
                )?,
                2 => result.add(
                    model,
                    self.disk_prime,
                    term(1, &self.edge_chain, self.edge_indices[generator.index]),
                    self.disk_prime.negate(*value),
                )?,
                _ => unreachable!("nonzero disk costalk cochain is in degree one or two"),
            }
        }
        sheaf::persist(model, result.into_terms())
    }
}

fn term(vertical: usize, chain: &Rc<[usize]>, index: usize) -> Generator {
    Generator {
        vertical,
        place: Place::Chain(chain.clone()),
        index,
    }
}

fn derived_name(name: &Label, suffix: &str) -> Result<Label> {
    let text = format!("{}__{suffix}", name.text);
    if text.len() > MAX_IDENTIFIER_BYTES {
        return Err(Error::new(
            name.span,
            "derived disk name exceeds the identifier length limit",
        ));
    }
    Ok(Label {
        text,
        span: name.span,
    })
}
