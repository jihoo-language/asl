// Two components contribute two copies of the sphere Adams E2 page.
space Two { a, b };
let first = \u : One => Two.a;
let second = \u : One => Two.b;
let collapse = \u : Two => One.pt;
let unit = extgen[0, 0, 0];
let h0 = extgen[1, 1, 0];
let point = xextgen[One, 0, 0, 0];
let a = xmap(first, point);
let b = xmap(second, point);

xextdim[Two, 0, 0];
xextdim[Two, 1, 1];
show xaction(h0, add(a, b));
assert xaction(unit, a) == a;
assert xmap(collapse, a) == point;
assert xmap(collapse, add(a, b)) == zero[XExt(One, 0, 0)];
assert xmap(collapse, xaction(h0, b)) == xaction(h0, point);
