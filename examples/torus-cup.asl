// A characteristic cochain of one edge of this hexagon generates H^1.
space Circle { a < ab, b < ab, b < bc, c < bc, c < ca, a < ca };
let u = cocycle[Circle, 1] { a < ab };

// Moving the mark across the vertex ab adds a coboundary.
assert u == cocycle[Circle, 1] { b < ab };

let first = \p : Circle * Circle => fst(p);
let second = \p : Circle * Circle => snd(p);
let a = pull(first, u);
let b = pull(second, u);
let area = cup(a, b);
show area;

assert area == gen[Circle * Circle, 2, 0];
assert cup(a, a) == zero[Circle * Circle, 2];
assert cup(b, b) == zero[Circle * Circle, 2];

// Swapping the coordinates exchanges the loops and preserves area over F2.
let swap = \p : Circle * Circle => (snd(p), fst(p));
assert pull(swap, a) == b;
assert pull(swap, area) == cup(b, a);
assert pull(swap, area) == area;
