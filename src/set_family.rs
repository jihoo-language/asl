use std::ops::{ControlFlow, Range};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Poset, TypeId, Value};
use crate::syntax::{Label, MAX_IDENTIFIER_BYTES};

pub(crate) const MAX_OCCURRENCES: usize = 4096;
pub(crate) const MAX_SOLUTIONS: usize = 4096;

pub(crate) type Component = (Label, Vec<usize>);
pub(crate) type Transition = (Label, Label, Vec<usize>);
type Arrow = Rc<[usize]>;
type Arrows = Vec<Vec<Option<Arrow>>>;

pub(crate) struct Diagram {
    pub(crate) base: TypeId,
    pub(crate) site: Rc<Poset>,
    pub(crate) order: Vec<usize>,
    offsets: Vec<usize>,
    arrows: Arrows,
}

impl Diagram {
    pub(crate) fn declare(
        model: &mut Model,
        base: TypeId,
        fibers: Vec<(Label, usize)>,
        transitions: Vec<Transition>,
    ) -> Result<Self> {
        let site = model.poset(base).expect("named space");
        let count = site.names.len();
        let mut sizes = vec![None; count];
        for (point, size) in fibers {
            let x = model.point(base, &point)?;
            if sizes[x].replace(size).is_some() {
                return Err(Error::new(point.span, "duplicate set-family fiber"));
            }
        }
        let sizes = sizes
            .into_iter()
            .zip(&site.names)
            .map(|(size, name)| {
                size.ok_or_else(|| Error::plain(format!("missing set-family fiber at '{name}'")))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut arrows = vec![vec![None; count]; count];
        for (source, target, entries) in transitions {
            let x = model.point(base, &source)?;
            let y = model.point(base, &target)?;
            if !site.strict_less(x, y) {
                return Err(Error::new(
                    source.span,
                    "set-family transition requires x < y in its site",
                ));
            }
            Self::check_entries(model, sizes[x], sizes[y], &entries)
                .map_err(|error| error.at(source.span))?;
            if arrows[x][y].replace(Arrow::from(entries)).is_some() {
                return Err(Error::new(source.span, "duplicate set-family transition"));
            }
        }

        let (order, covers) = Self::site_orders(model, &site)?;
        let offsets = Self::build_offsets(&sizes)?;
        let mut cells = arrows
            .iter()
            .flatten()
            .flatten()
            .map(|arrow| arrow.len())
            .sum::<usize>();
        check_cells(cells)?;
        for &x in &order {
            for &y in &covers[x] {
                if arrows[x][y].is_none() {
                    return Err(Error::plain(format!(
                        "missing set-family transition on '{} < {}'",
                        site.names[x], site.names[y]
                    )));
                }
            }

            for &y in &order {
                if !site.strict_less(x, y) {
                    continue;
                }
                let first = arrows[x][y]
                    .as_ref()
                    .expect("reachable through a cover")
                    .clone();
                for &z in &covers[y] {
                    let next = arrows[y][z].as_ref().ok_or_else(|| {
                        Error::plain(format!(
                            "missing set-family transition on '{} < {}'",
                            site.names[y], site.names[z]
                        ))
                    })?;
                    model.spend(first.len() + 1)?;
                    if let Some(expected) = &arrows[x][z] {
                        if !first.iter().map(|&i| next[i]).eq(expected.iter().copied()) {
                            return Err(Error::plain(format!(
                                "set-family transitions disagree along '{} < {} < {}'",
                                site.names[x], site.names[y], site.names[z]
                            )));
                        }
                    } else {
                        cells += first.len();
                        check_cells(cells)?;
                        arrows[x][z] = Some(first.iter().map(|&i| next[i]).collect());
                    }
                }
            }
        }

        model.allocate(cells + count * count + 2 * count + 1)?;
        Ok(Self {
            base,
            site,
            order,
            offsets,
            arrows,
        })
    }

    pub(crate) fn generated(
        model: &mut Model,
        base: TypeId,
        sizes: Vec<usize>,
        mut transition: impl FnMut(&mut Model, usize, usize) -> Result<Vec<usize>>,
    ) -> Result<Self> {
        let site = model.poset(base).expect("named space");
        let count = site.names.len();
        let offsets = Self::build_offsets(&sizes)?;
        let (order, _) = Self::site_orders(model, &site)?;
        let mut cells = count * count + 2 * count + 1;
        let mut arrows = vec![vec![None; count]; count];
        model.spend(count * count)?;
        for (x, row) in arrows.iter_mut().enumerate() {
            for (y, arrow) in row.iter_mut().enumerate() {
                if site.strict_less(x, y) {
                    cells += sizes[x];
                    check_cells(cells)?;
                    let entries = transition(model, x, y)?;
                    Self::check_entries(model, sizes[x], sizes[y], &entries)?;
                    *arrow = Some(Arrow::from(entries));
                }
            }
        }
        model.allocate(cells)?;
        Ok(Self {
            base,
            site,
            order,
            offsets,
            arrows,
        })
    }

    fn site_orders(model: &mut Model, site: &Poset) -> Result<(Vec<usize>, Vec<Vec<usize>>)> {
        let n = site.names.len();
        let mut lower = vec![0; n];
        let mut uppers = vec![Vec::new(); n];
        model.spend(n * n)?;
        for (x, upper) in uppers.iter_mut().enumerate() {
            for (y, degree) in lower.iter_mut().enumerate() {
                if site.strict_less(x, y) {
                    upper.push(y);
                    *degree += 1;
                }
            }
        }
        let mut order: Vec<_> = (0..n).collect();
        order.sort_unstable_by_key(|&x| (lower[x], x));
        for upper in &mut uppers {
            upper.sort_unstable_by_key(|&y| (lower[y], y));
        }
        let mut covers = vec![Vec::new(); n];
        for (x, upper) in uppers.iter().enumerate() {
            let mut covered = vec![false; n];
            for &y in upper {
                model.spend(1)?;
                if covered[y] {
                    continue;
                }
                covers[x].push(y);
                model.spend(uppers[y].len())?;
                for &z in &uppers[y] {
                    covered[z] = true;
                }
            }
        }
        Ok((order, covers))
    }

    fn build_offsets(sizes: &[usize]) -> Result<Vec<usize>> {
        let mut offsets = vec![0];
        let mut count = 0usize;
        for &size in sizes {
            count = count
                .checked_add(size)
                .filter(|&n| n <= MAX_OCCURRENCES)
                .ok_or_else(|| Error::plain("set-family exceeds the limit of 4096 occurrences"))?;
            offsets.push(count);
        }
        Ok(offsets)
    }

    pub(crate) fn size(&self, point: usize) -> usize {
        self.offsets[point + 1] - self.offsets[point]
    }

    pub(crate) fn occurrence(&self, point: usize, index: usize) -> usize {
        self.offsets[point] + index
    }

    fn fiber(&self, point: usize) -> Range<usize> {
        self.offsets[point]..self.offsets[point + 1]
    }

    pub(crate) fn count(&self) -> usize {
        *self.offsets.last().unwrap()
    }

    pub(crate) fn apply(&self, x: usize, y: usize, index: usize) -> usize {
        if x == y {
            index
        } else {
            self.arrows[x][y].as_ref().expect("comparable points")[index]
        }
    }

    pub(crate) fn components<'a>(&self, point: usize, entries: &'a [usize]) -> &'a [usize] {
        &entries[self.fiber(point)]
    }

    fn check_entries(
        model: &mut Model,
        source: usize,
        target: usize,
        entries: &[usize],
    ) -> Result<()> {
        if entries.len() != source {
            return Err(Error::plain(format!(
                "set map expects {source} images, found {}",
                entries.len()
            )));
        }
        model.spend(entries.len() + 1)?;
        if let Some(index) = entries.iter().find(|&&i| i >= target) {
            return Err(Error::plain(format!(
                "set map image {index} is outside its target fiber of size {target}"
            )));
        }
        Ok(())
    }

    pub(crate) fn identity(&self, model: &mut Model) -> Result<Value> {
        model.spend(self.count())?;
        model.allocate(self.count())?;
        Ok(Value::SetNat(
            (0..self.site.names.len())
                .flat_map(|x| 0..self.size(x))
                .collect(),
        ))
    }

    fn visit_sections(
        &self,
        model: &mut Model,
        points: &[usize],
        emit: impl FnMut(&[usize]) -> Result<ControlFlow<()>>,
    ) -> Result<bool> {
        let mut positions = vec![None; self.site.names.len()];
        model.spend(positions.len() + points.len())?;
        for (index, &point) in points.iter().enumerate() {
            positions[point] = Some(index);
        }
        let domains: Vec<_> = points.iter().map(|&x| self.size(x)).collect();
        let order: Vec<_> = self.order.iter().filter_map(|&x| positions[x]).collect();
        solve(
            model,
            &domains,
            &order,
            |model, x, i, push| {
                model.spend(self.site.names.len())?;
                push.extend(
                    self.arrows[points[x]]
                        .iter()
                        .zip(&positions)
                        .filter_map(|(map, &position)| Some((position?, map.as_ref()?[i]))),
                );
                Ok(())
            },
            emit,
        )
    }

    pub(crate) fn section_count(&self, model: &mut Model) -> Result<usize> {
        let points: Vec<_> = (0..self.site.names.len()).collect();
        let mut count = 0;
        self.visit_sections(model, &points, |_| {
            count += 1;
            Ok(ControlFlow::Continue(()))
        })?;
        Ok(count)
    }

    pub(crate) fn sections_on(
        &self,
        model: &mut Model,
        points: &[usize],
    ) -> Result<Vec<Vec<usize>>> {
        let mut rows = Vec::new();
        self.visit_sections(model, points, |row| {
            if rows.len() == MAX_SOLUTIONS {
                return Err(Error::plain(
                    "set-family exceeds the limit of 4096 sections",
                ));
            }
            check_cells((rows.len() + 1) * (row.len() + 1))?;
            rows.push(row.to_vec());
            Ok(ControlFlow::Continue(()))
        })?;
        let sort_steps = rows.len().max(1).ilog2() as usize + 1;
        model.spend(rows.len() * (points.len() + 1) * sort_steps)?;
        rows.sort_unstable();
        Ok(rows)
    }

    fn visit_natural(
        &self,
        model: &mut Model,
        target: &Self,
        emit: impl FnMut(&[usize]) -> Result<ControlFlow<()>>,
    ) -> Result<bool> {
        let sources: Vec<_> = (0..self.site.names.len())
            .flat_map(|x| (0..self.size(x)).map(move |i| (x, i)))
            .collect();
        let domains: Vec<_> = sources.iter().map(|&(x, _)| target.size(x)).collect();
        let order: Vec<_> = self.order.iter().flat_map(|&x| self.fiber(x)).collect();
        solve(
            model,
            &domains,
            &order,
            |model, point, value, push| {
                let (x, i) = sources[point];
                model.spend(self.site.names.len())?;
                push.extend(self.arrows[x].iter().enumerate().filter_map(|(y, map)| {
                    map.as_ref()
                        .map(|images| (self.occurrence(y, images[i]), target.apply(x, y, value)))
                }));
                Ok(())
            },
            emit,
        )
    }

    pub(crate) fn has_natural(&self, model: &mut Model, target: &Self) -> Result<bool> {
        self.visit_natural(model, target, |_| Ok(ControlFlow::Break(())))
    }

    pub(crate) fn natural_elements(&self, model: &mut Model, target: &Self) -> Result<Vec<Value>> {
        let mut rows = Vec::new();
        self.visit_natural(model, target, |row| {
            if rows.len() == MAX_SOLUTIONS {
                return Err(Error::plain(
                    "type exceeds the enumeration limit of 4096 elements",
                ));
            }
            check_cells((rows.len() + 1) * (row.len() + 1))?;
            rows.push(Arrow::from(row));
            Ok(ControlFlow::Continue(()))
        })?;
        let cells = rows.len() * (self.count() + 1);
        let sort_steps = rows.len().max(1).ilog2() as usize + 1;
        model.spend(cells * sort_steps)?;
        model.allocate(cells)?;
        rows.sort_unstable();
        Ok(rows.into_iter().map(Value::SetNat).collect())
    }

    pub(crate) fn natural(
        &self,
        model: &mut Model,
        target: &Self,
        components: Vec<Component>,
    ) -> Result<Value> {
        let mut images = vec![None; self.site.names.len()];
        for (point, entries) in components {
            let x = model.point(self.base, &point)?;
            Self::check_entries(model, self.size(x), target.size(x), &entries)
                .map_err(|error| error.at(point.span))?;
            if images[x].replace(entries).is_some() {
                return Err(Error::new(
                    point.span,
                    "duplicate natural transformation component",
                ));
            }
        }
        let mut rows = Vec::with_capacity(self.count());
        for (image, name) in images.into_iter().zip(&self.site.names) {
            rows.extend(image.ok_or_else(|| {
                Error::plain(format!(
                    "missing natural transformation component at '{name}'"
                ))
            })?);
        }
        self.certify(model, target, rows)
    }

    pub(crate) fn certify(
        &self,
        model: &mut Model,
        target: &Self,
        rows: Vec<usize>,
    ) -> Result<Value> {
        debug_assert_eq!(rows.len(), self.count());
        for x in 0..self.site.names.len() {
            let image = self.components(x, &rows);
            Self::check_entries(model, self.size(x), target.size(x), image)?;
            model.spend(self.site.names.len())?;
            for (y, arrow) in self.arrows[x].iter().enumerate() {
                if let Some(arrow) = arrow {
                    model.spend(arrow.len())?;
                    let upper = self.components(y, &rows);
                    if !arrow
                        .iter()
                        .zip(image)
                        .all(|(&i, &j)| upper[i] == target.apply(x, y, j))
                    {
                        return Err(Error::plain(format!(
                            "natural transformation does not commute along '{} < {}'",
                            self.site.names[x], self.site.names[y]
                        )));
                    }
                }
            }
        }
        model.allocate(rows.len())?;
        Ok(Value::SetNat(rows.into()))
    }

    pub(crate) fn compose(
        &self,
        model: &mut Model,
        middle: &Self,
        outer: &[usize],
        inner: &[usize],
    ) -> Result<Value> {
        model.spend(self.count())?;
        model.allocate(self.count())?;
        let rows = (0..self.site.names.len())
            .flat_map(|x| {
                let outer = middle.components(x, outer);
                self.components(x, inner).iter().map(|&i| outer[i])
            })
            .collect();
        Ok(Value::SetNat(rows))
    }
}

pub(crate) struct Elements {
    pub(crate) family: Rc<str>,
    pub(crate) points: Vec<(usize, usize)>,
    projection: Rc<[Value]>,
}

impl Elements {
    pub(crate) fn build(
        model: &mut Model,
        name: &Label,
        family: Rc<str>,
        diagram: &Diagram,
    ) -> Result<(TypeId, Self)> {
        model.check_named_space(name, diagram.count())?;
        let points: Vec<_> = (0..diagram.site.names.len())
            .flat_map(|x| (0..diagram.size(x)).map(move |i| (x, i)))
            .collect();
        let count = points.len();
        model.spend(count * count)?;
        model.allocate(count * count + 3 * count + 1)?;
        let mut names: Vec<_> = points
            .iter()
            .map(|&(x, i)| format!("{}_{i}", diagram.site.names[x]))
            .collect();
        if names.iter().any(|name| name.len() > MAX_IDENTIFIER_BYTES) {
            names = points
                .iter()
                .map(|&(x, i)| format!("point_{x}_{i}"))
                .collect();
        }
        let order = points
            .iter()
            .map(|&(x, i)| {
                points
                    .iter()
                    .map(|&(y, j)| {
                        (x == y || diagram.site.strict_less(x, y)) && diagram.apply(x, y, i) == j
                    })
                    .collect()
            })
            .collect();
        let projection = points.iter().map(|&(x, _)| Value::Point(x)).collect();
        let ty = model.insert_derived_space(name, names, order);
        Ok((
            ty,
            Self {
                family,
                points,
                projection,
            },
        ))
    }

    pub(crate) fn projection(&self) -> Value {
        Value::Map(self.projection.clone())
    }

    pub(crate) fn section(
        &self,
        model: &mut Model,
        diagram: &Diagram,
        index: usize,
    ) -> Result<Value> {
        if index >= MAX_SOLUTIONS {
            return Err(Error::plain(
                "set section index exceeds the limit of 4096 sections",
            ));
        }
        let points: Vec<_> = (0..diagram.site.names.len()).collect();
        let sections = diagram.sections_on(model, &points)?;
        let row = sections
            .get(index)
            .ok_or_else(|| Error::plain("set section index is outside its range"))?;
        model.allocate(row.len())?;
        Ok(Value::Map(
            row.iter()
                .enumerate()
                .map(|(x, &i)| Value::Point(diagram.occurrence(x, i)))
                .collect(),
        ))
    }

    pub(crate) fn map(&self, model: &mut Model, rows: &[usize], target: &Diagram) -> Result<Value> {
        model.spend(rows.len())?;
        model.allocate(rows.len())?;
        Ok(Value::Map(
            self.points
                .iter()
                .zip(rows)
                .map(|(&(x, _), &image)| Value::Point(target.occurrence(x, image)))
                .collect(),
        ))
    }
}

pub(crate) fn check_cells(cells: usize) -> Result<()> {
    if cells > 1_000_000 {
        return Err(Error::plain(
            "set-family exceeds the finite-value storage limit",
        ));
    }
    Ok(())
}

struct Branch {
    point: usize,
    next: usize,
    start: usize,
}

// A choice fixes every successor; inconsistent branches are stopped immediately.
// Follow a fixed topological order so skipping a forced vertex loses no condition.
pub(crate) fn solve(
    model: &mut Model,
    domains: &[usize],
    order: &[usize],
    mut force: impl FnMut(&mut Model, usize, usize, &mut Vec<(usize, usize)>) -> Result<()>,
    mut emit: impl FnMut(&[usize]) -> Result<ControlFlow<()>>,
) -> Result<bool> {
    const UNSET: usize = usize::MAX;
    if domains.contains(&0) {
        return Ok(false);
    }
    let mut row = vec![UNSET; domains.len()];
    let mut trail = Vec::new();
    let mut branches = Vec::<Branch>::new();
    let mut pending = Vec::new();
    let mut place = 0;
    loop {
        while place < order.len() {
            model.spend(1)?;
            if row[order[place]] == UNSET {
                branches.push(Branch {
                    point: place,
                    next: 0,
                    start: trail.len(),
                });
                break;
            }
            place += 1;
        }
        if place == order.len() {
            model.spend(row.len() + 1)?;
            if emit(&row)?.is_break() {
                return Ok(true);
            }
        }

        'search: loop {
            let Some(branch) = branches.last_mut() else {
                return Ok(false);
            };
            for point in trail.drain(branch.start..) {
                row[point] = UNSET;
            }
            let x = order[branch.point];
            if branch.next == domains[x] {
                branches.pop();
                continue;
            }
            model.spend(1)?;
            let value = branch.next;
            branch.next += 1;
            pending.clear();
            pending.push((x, value));
            force(model, x, value, &mut pending)?;
            for &(point, choice) in &pending {
                if row[point] == UNSET {
                    row[point] = choice;
                    trail.push(point);
                } else if row[point] != choice {
                    continue 'search;
                }
            }
            place = branch.point + 1;
            break;
        }
    }
}
