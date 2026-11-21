// Degree-zero maps of disconnected spaces act by their homology matrices.
space Two { a, b };
let swap = tab Two -> Two { a => Two.b, b => Two.a };
let first = \x : Two => Two.a;
let second = \x : Two => Two.b;
let u = emap(swap);
let a = emap(first);
let b = emap(second);
assert ecomp(u, u) == eid[Two];
assert ecomp(u, a) == b;
assert ecomp(a, u) == a;
assert ecomp(u, a) == emap(\x : Two => swap(first(x)));
let sum = add(a, b);
assert ecomp(sum, sum) == zero[MapExt(Two, Two, 0, 0)];
mapextdim[Two, Two, 0, 0];

// Degree −1 selects the circle's H_1, and degree 1 inserts it.
space Circle { a < ab, b < ab, b < bc, c < bc, c < ca, a < ca };
let select = mapextgen[Circle, One, 0, -1, 0];
let insert = mapextgen[One, Circle, 0, 1, 0];
assert ecomp(select, insert) == eid[One];
assert add(ecomp(insert, select), eid[Circle]) == mapextgen[Circle, Circle, 0, 0, 0];
mapextdim[Circle, One, 0, -1];
show ecomp(insert, select);
