use crate::error::{Error, Result};
use crate::model::Model;

use super::super::{
    tree::{self, Node, Tree},
    Ambient, Balanced, Category, Module,
};
use super::branch::Bulk;
use super::movement::{Central, Layout, Map};

pub(super) fn unit(
    model: &mut Model,
    ambient: &Ambient<'_>,
    category: &Category,
    simple: usize,
) -> Result<String> {
    let object = Bulk::new(model, ambient, category, simple)?;
    let regular = Bulk::new(model, ambient, category, 0)?;
    let map = Central::new(model, ambient, [&object, &object], &regular.module)?;
    let [source, target] = &map.products;
    let mut count = 0;
    for (g, basis) in source.bases.iter().enumerate() {
        model.spend(basis.len())?;
        for (i, pair) in basis.iter().enumerate() {
            let [u, h] = pair.grades;
            let left_unit = object
                .module
                .grade_action(ambient, 1, h, u, pair.columns[0]);
            let crossing = map.maps[0].0[g][i];
            let output = target.bases[g][crossing.target];
            let [k, v] = output.grades;
            let right_unit = object
                .module
                .grade_action(ambient, 0, k, v, output.columns[1]);
            if left_unit.target != right_unit.target
                || left_unit.phase != ambient.sum([crossing.phase, right_unit.phase])
            {
                return Err(Error::plain(
                    "central map does not equal the regular left and right unitors",
                ));
            }
            count += 1;
        }
    }
    Ok(format!(
        "rvDoubleCentralUnit[{simple}] = columns={count}, conductor={}, verified=true",
        ambient.conductor
    ))
}

pub(super) fn whisker(
    model: &mut Model,
    ambient: &Ambient<'_>,
    bulk: [&Bulk; 3],
    defects: [&Module; 2],
) -> Result<String> {
    let [bh, bk, bl] = bulk;
    let [m, n] = defects;
    let arrays = [[&bh.module, m, n], [m, &bk.module, n], [m, n, &bl.module]];
    // Left-associated and right-associated groupings for the three
    // positions of the bulk-induced strand.
    let pairs = arrays
        .iter()
        .map(|inputs| {
            let [l, r] = &tree::TRIPLE;
            Ok([
                Tree::build(model, ambient, inputs, l)?,
                Tree::build(model, ambient, inputs, r)?,
            ])
        })
        .collect::<Result<Vec<_>>>()?;
    let assoc = pairs
        .iter()
        .zip(&arrays)
        .map(|([l, r], inputs)| tree::rotation(model, ambient, inputs, l, r, [0, 1, 2, 3]).map(Map))
        .collect::<Result<Vec<_>>>()?;
    let c_m = Central::new(model, ambient, [bh, bk], m)?;
    let c_n = Central::new(model, ambient, [bk, bl], n)?;
    verify_inner_basis(&pairs[0][0], 0, &c_m.products[0])?;
    verify_inner_basis(&pairs[1][0], 0, &c_m.products[1])?;
    verify_inner_basis(&pairs[1][1], 1, &c_n.products[0])?;
    verify_inner_basis(&pairs[2][1], 1, &c_n.products[1])?;
    let first = c_m.whisker(
        model,
        ambient,
        layout(&pairs[0][0]),
        layout(&pairs[1][0]),
        0,
    )?;
    let second = c_n.whisker(
        model,
        ambient,
        layout(&pairs[1][1]),
        layout(&pairs[2][1]),
        1,
    )?;
    let direct = Central::new(model, ambient, [bh, bl], layout(&pairs[0][1]).factors[1])?;
    let [source, target] = &direct.products;
    if source.bases != layout(&pairs[0][1]).product.bases
        || target.bases != layout(&pairs[2][0]).product.bases
        || &source.module != pairs[0][1].module()
        || &target.module != pairs[2][0].module()
    {
        return Err(Error::plain(
            "central coherence paths differ in their outer balanced bases",
        ));
    }
    let composite = assoc[0]
        .inverse(ambient)
        .followed(ambient, &first)
        .followed(ambient, &assoc[1])
        .followed(ambient, &second)
        .followed(ambient, &assoc[2].inverse(ambient));
    if composite != direct.maps[0] {
        return Err(Error::plain(
            "central half-braiding fails the whiskering hexagon",
        ));
    }
    let (phases, moves) = direct.maps[0].statistics();
    Ok(format!("rvDoubleCentralWhisker = dimension={}, conductor={}, nonzero_phases={phases}, moved_columns={moves}, verified=true",
        source.module.dimension(),ambient.conductor))
}

fn verify_inner_basis(tree: &Tree<'_>, side: usize, expected: &Balanced) -> Result<()> {
    let Node::Join(branch) = &tree.node else {
        unreachable!("threefold product");
    };
    if layout(&branch.children[side]).product.bases != expected.bases {
        return Err(Error::plain(
            "central whiskering changed the inner balanced basis",
        ));
    }
    Ok(())
}

fn layout<'a>(tree: &'a Tree<'_>) -> Layout<'a> {
    let Node::Join(branch) = &tree.node else {
        unreachable!("central coherence tree has three inputs");
    };
    Layout {
        product: &branch.product,
        factors: branch.children.each_ref().map(Tree::module),
    }
}
