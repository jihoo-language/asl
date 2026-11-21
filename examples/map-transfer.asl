space Two { a, b };
space Three { a, b, c };

// The matrix paths use h0, h1, h2, h0 in source-to-target order.
let f = add(mapextgen[One, Two, 1, 1, 0], mapextgen[One, Two, 1, 1, 1]);
let g = add(mapextgen[Two, Three, 1, 2, 0], mapextgen[Two, Three, 1, 2, 4]);
let h = mapextgen[Three, S, 1, 4, 0];
let k = mapextgen[S, One, 1, 1, 0];
let h2 = mapextgen[One, One, 1, 4, 0];
let fourth = eHigher(k, h, g, f);
assert fourth == ecomp(h2, h2);
show fourth;
assert eHigher(g, eid[Two], f) == zero[MapExt(One, Three, 1, 3)];
