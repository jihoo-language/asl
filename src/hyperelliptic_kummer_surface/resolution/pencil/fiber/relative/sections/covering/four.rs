use std::rc::Rc;

use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::nikulin::binary_rref;
use crate::model::Model;
use crate::syntax::{Label, ParsedExpr};

use super::{cassels::Pairing, Descent, MordellImage};

pub(crate) enum Expression {
    Pair(Label, Box<[ParsedExpr; 2]>),
    Contains(Label, Box<ParsedExpr>),
    Basis(Label, usize, bool),
}

pub(crate) struct Certificate {
    image: Rc<MordellImage>,
    coordinates: Coordinates,
    form: Vec<u32>,
    source: Vec<u32>,
    symplectic: Vec<u32>,
    rank: usize,
    radical: Vec<u32>,
    order: Option<u128>,
    torsion_places: Option<[String; 2]>,
}

impl Certificate {
    pub(crate) fn new(
        model: &mut Model,
        image: Rc<MordellImage>,
        witnesses: &[Rc<Pairing>],
    ) -> Result<Self> {
        let pencil = image.pencil();
        let selmer = pencil.selmer_space(model)?;
        let size = selmer.sha_dimension(image.rank())?;
        if witnesses.len() != size {
            return Err(Error::plain(format!("four-Selmer certificate requires {size} Cassels witnesses for a basis modulo the global image")));
        }
        let mut coordinates = Coordinates::default();
        let mut source = Vec::new();
        for (index, witness) in witnesses.iter().enumerate() {
            if !Rc::ptr_eq(pencil, &witness.cover.pencil) {
                return Err(Error::plain(
                    "four-Selmer witness must use the certificate's named pencil",
                ));
            }
            let class = selmer
                .encode_local(model, &witness.cover.b)?
                .ok_or_else(inconsistent)?;
            let vector = image.reduce(model, class)?;
            coordinates.insert(vector, 1 << index)?;
            source.push(class);
        }
        let classes = source
            .iter()
            .chain(image.rows())
            .map(|&vector| selmer.decode(model, vector))
            .collect::<Result<Vec<_>>>()?;
        let mut form = Vec::new();
        for witness in witnesses {
            let row = witness.bits(model, &classes)?;
            if row[size..].iter().any(|&entry| entry) {
                return Err(Error::plain(
                    "Cassels evaluation did not annihilate the certified global Mordell-Weil image",
                ));
            }
            let mask = row[..size]
                .iter()
                .enumerate()
                .fold(0, |mask, (i, &entry)| mask | (u32::from(entry) << i));
            form.push(mask);
        }
        let (normal, rank) = normal_form(model, &form)?;
        let symplectic: Vec<_> = normal
            .iter()
            .map(|&vector| combination(&source, vector))
            .collect();
        for (&lift, &vector) in symplectic.iter().zip(&normal) {
            if coordinates.reduce(image.reduce(model, lift)?, 0) != (0, vector) {
                return Err(inconsistent());
            }
        }
        let mut radical = image.rows().to_vec();
        radical.extend_from_slice(&symplectic[rank..]);
        if binary_rref(&mut radical, selmer.width()).len() != radical.len() {
            return Err(inconsistent());
        }
        for &lift in &radical {
            let (remainder, vector) = coordinates.reduce(image.reduce(model, lift)?, 0);
            if remainder != 0 || combination(&form, vector) != 0 {
                return Err(inconsistent());
            }
        }
        let torsion_places = pencil.split_torsion_certificate(model)?;
        let order =
            if torsion_places.is_some() {
                let exponent = 2 * image.rank() + 2 + size + (size - rank);
                Some(1u128.checked_shl(exponent as u32).ok_or_else(|| {
                    Error::plain("four-Selmer order exceeds the exact integer range")
                })?)
            } else {
                None
            };
        let place_cells = torsion_places
            .as_ref()
            .map_or(0, |places| places.iter().map(String::len).sum());
        model.allocate(
            106 + form.len() + source.len() + symplectic.len() + radical.len() + place_cells,
        )?;
        Ok(Self {
            image,
            coordinates,
            form,
            source,
            symplectic,
            rank,
            radical,
            order,
            torsion_places,
        })
    }

    pub(crate) fn pencil_name(&self) -> &Label {
        &self.image.pencil_name
    }

    pub(crate) fn pairing(
        &self,
        model: &mut Model,
        left: &Descent,
        right: &Descent,
    ) -> Result<u32> {
        let left = self.coordinate(model, left)?.ok_or_else(not_selmer)?;
        let right = self.coordinate(model, right)?.ok_or_else(not_selmer)?;
        model.spend(self.form.len())?;
        Ok(pair(&self.form, left, right))
    }

    pub(crate) fn contains(&self, model: &mut Model, class: &Descent) -> Result<bool> {
        let Some(vector) = self.coordinate(model, class)? else {
            return Ok(false);
        };
        model.spend(self.form.len())?;
        Ok(combination(&self.form, vector) == 0)
    }

    pub(crate) fn basis(&self, model: &mut Model, index: usize, quotient: bool) -> Result<Descent> {
        let basis = if quotient {
            &self.symplectic
        } else {
            &self.radical
        };
        let vector = *basis.get(index).ok_or_else(|| {
            Error::plain(format!(
                "four-Selmer basis index {index} is outside dimension {}",
                basis.len()
            ))
        })?;
        self.image
            .pencil()
            .selmer_space(model)?
            .decode(model, vector)
    }

    fn coordinate(&self, model: &mut Model, class: &Descent) -> Result<Option<u32>> {
        let selmer = self.image.pencil().selmer_space(model)?;
        let Some(vector) = selmer.encode_local(model, class)? else {
            return Ok(None);
        };
        model.spend(selmer.width())?;
        let quotient = self.image.reduce(model, vector)?;
        let (remainder, coordinate) = self.coordinates.reduce(quotient, 0);
        if remainder != 0 {
            return Err(inconsistent());
        }
        Ok(Some(coordinate))
    }

    pub(crate) fn report(&self, name: &str) -> String {
        let size = self.form.len();
        let sha_divisible = size - self.rank;
        let mut factors = Vec::new();
        if sha_divisible > 0 {
            factors.push(format!("(Z/4)^{sha_divisible}"));
        }
        if self.rank > 0 {
            factors.push(format!("(Z/2)^{}", self.rank));
        }
        let sha_four = if factors.is_empty() {
            "0".into()
        } else {
            factors.join(" + ")
        };
        let order = match (self.order, &self.torsion_places) {
            (Some(order), Some([additive, multiplicative])) => format!(
                "Sel4 order={order}, torsion=(Z/2)^2 certified by {additive} and {multiplicative}"
            ),
            _ => "Sel4 order unavailable without exact torsion certificate".into(),
        };
        let matrix: Vec<Vec<_>> = self
            .form
            .iter()
            .map(|row| (0..size).map(|j| (row >> j) & 1).collect())
            .collect();
        let pairs: Vec<_> = self.symplectic[..self.rank]
            .chunks_exact(2)
            .map(|pair| (pair[0], pair[1]))
            .collect();
        format!(
            "rvK3Selmer4[{name}] = dim Sel2={}, dim E(K)/2={}, quotient source masks={:?}, Cassels matrix={matrix:?}, rank={}, symplectic pairs={pairs:?}, quotient radical lifts={:?}, dim image Sel4->Sel2={}, image basis masks={:?}, dim Sha[2] intersect 2Sha={sha_divisible}, Sha[4]={sha_four}, {order}",
            self.image.rows().len()+size, self.image.rows().len(), self.source, self.rank,
            &self.symplectic[self.rank..], self.radical.len(), self.radical,
        )
    }
}

#[derive(Default)]
struct Coordinates([Option<(u32, u32)>; 32]);

impl Coordinates {
    fn insert(&mut self, vector: u32, tag: u32) -> Result<()> {
        let (row, row_tag) = self.reduce(vector, tag);
        if row == 0 {
            return Err(Error::plain(
                "four-Selmer witnesses are dependent modulo the global Mordell-Weil image",
            ));
        }
        self.0[row.trailing_zeros() as usize] = Some((row, row_tag));
        Ok(())
    }

    fn reduce(&self, mut vector: u32, mut tag: u32) -> (u32, u32) {
        for (i, entry) in self.0.iter().enumerate() {
            if let Some((row, row_tag)) = entry {
                if vector & (1 << i) != 0 {
                    vector ^= row;
                    tag ^= row_tag;
                }
            }
        }
        (vector, tag)
    }
}

fn combination(basis: &[u32], vector: u32) -> u32 {
    basis
        .iter()
        .enumerate()
        .filter(|(i, _)| vector & (1 << i) != 0)
        .fold(0, |sum, (_, row)| sum ^ row)
}

fn pair(matrix: &[u32], left: u32, right: u32) -> u32 {
    (combination(matrix, left) & right).count_ones() % 2
}

fn normal_form(model: &mut Model, form: &[u32]) -> Result<(Vec<u32>, usize)> {
    let size = form.len();
    model.spend(size.pow(4))?;
    for (i, &row) in form.iter().enumerate() {
        if row & (1 << i) != 0 || (0..i).any(|j| ((row >> j) & 1) != ((form[j] >> i) & 1)) {
            return Err(Error::plain(
                "Cassels evaluations failed the alternating matrix check",
            ));
        }
    }
    let mut remaining: Vec<_> = (0..size).map(|i| 1u32 << i).collect();
    let mut basis = Vec::new();
    while let Some((i, j)) = remaining.iter().enumerate().find_map(|(i, &u)| {
        (i + 1..remaining.len())
            .find(|&j| pair(form, u, remaining[j]) == 1)
            .map(|j| (i, j))
    }) {
        let v = remaining.remove(j);
        let u = remaining.remove(i);
        for w in &mut remaining {
            *w ^= (pair(form, *w, v) * u) ^ (pair(form, *w, u) * v);
        }
        basis.extend([u, v]);
    }
    let rank = basis.len();
    basis.extend(remaining);
    if binary_rref(&mut basis.clone(), size).len() != size
        || binary_rref(&mut form.to_vec(), size).len() != rank
    {
        return Err(inconsistent());
    }
    for (i, &u) in basis.iter().enumerate() {
        for (j, &v) in basis.iter().enumerate() {
            let expected = u32::from(i < rank && (i ^ 1) == j);
            if pair(form, u, v) != expected {
                return Err(inconsistent());
            }
        }
    }
    Ok((basis, rank))
}

fn not_selmer() -> Error {
    Error::plain("four-Selmer pairing requires two everywhere locally soluble classes")
}

fn inconsistent() -> Error {
    Error::plain("four-Selmer certificate failed its exact binary linear algebra checks")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn degenerate_normal_form_retains_noncoordinate_radical_lifts() {
        let mut model = Model::new();
        let form = [0b0110, 0b0101, 0b0011, 0];
        let (basis, rank) =
            normal_form(&mut model, &form).unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(rank, 2);
        assert_eq!(basis, [1, 2, 7, 8]);
        assert_eq!(
            basis[rank..]
                .iter()
                .map(|&v| combination(&form, v))
                .collect::<Vec<_>>(),
            [0, 0]
        );
        assert!(normal_form(&mut model, &[1]).is_err());
        assert!(normal_form(&mut model, &[2, 0]).is_err());
    }

    #[test]
    fn tracked_elimination_recovers_the_input_basis_coordinates() {
        let mut coordinates = Coordinates::default();
        let rows = [0b1100, 0b0110, 0b0011];
        for (i, row) in rows.into_iter().enumerate() {
            coordinates
                .insert(row, 1 << i)
                .unwrap_or_else(|error| panic!("{error}"));
        }
        for mask in 0..8 {
            assert_eq!(coordinates.reduce(combination(&rows, mask), 0), (0, mask));
        }
        assert!(coordinates.insert(rows[0] ^ rows[2], 8).is_err());
        assert_ne!(coordinates.reduce(0b10000, 0).0, 0);
    }
}
