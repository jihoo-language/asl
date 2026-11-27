use crate::error::{Error, Result};
use crate::model::Model;

use super::graded::Ambient;
use super::tree::{self, Tree};
use super::{check_categories, decompose, sector_counts, Category, Double, Operand};

pub(crate) fn associate(
    model: &mut Model,
    double: &Double,
    target: &Category,
    inputs: [(&Category, Operand); 3],
    entry: Option<[usize; 2]>,
) -> Result<String> {
    let [(a, x), (b, y), (c, z)] = inputs;
    check_categories(&[a, b, c], target)?;
    let ambient = Ambient::new(
        model,
        double,
        &[a, b, c, target],
        [x.conductor(), y.conductor(), z.conductor()]
            .into_iter()
            .flatten(),
    )?;
    let modules = [
        x.module(model, &ambient, a)?,
        y.module(model, &ambient, b)?,
        z.module(model, &ambient, c)?,
    ];
    let [left, right] = tree::TRIPLE.map(|shape| Tree::build(model, &ambient, &modules, &shape));
    let [left, right] = [left?, right?];
    let images = tree::rotation(model, &ambient, &modules, &left, &right, [0, 1, 2, 3])?;
    let first = sector_counts(&decompose(model, &ambient, target, left.module())?);
    let second = sector_counts(&decompose(model, &ambient, target, right.module())?);
    if first != second {
        return Err(Error::plain(
            "relative associator changed the simple-sector decomposition",
        ));
    }
    if let Some([grade, column]) = entry {
        let image = images
            .get(grade)
            .and_then(|row| row.get(column))
            .ok_or_else(|| {
                Error::plain("relative associator grade or basis column is out of range")
            })?;
        Ok(format!(
            "rvDoubleRelativeAssocImage[{grade},{column}] = target={}, phase={}, conductor={}",
            image.target, image.phase, ambient.conductor
        ))
    } else {
        let nonzero = images.iter().flatten().filter(|x| x.phase != 0).count();
        Ok(format!(
            "rvDoubleRelativeAssociator = dimension={}, conductor={}, nonzero_phases={nonzero}, sectors={first:?}, verified=true",
            left.module().dimension(),ambient.conductor,
        ))
    }
}
