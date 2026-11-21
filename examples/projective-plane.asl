// The face poset of a six-vertex triangulation of the real projective plane.
// Its thirty-one points are six vertices, fifteen edges and ten faces.
space RP2 {
  a < ab, b < ab, a < ac, c < ac, a < ad, d < ad,
  a < ae, e < ae, a < af, f < af, b < bc, c < bc,
  b < bd, d < bd, b < be, e < be, b < bf, f < bf,
  c < cd, d < cd, c < ce, e < ce, c < cf, f < cf,
  d < de, e < de, d < df, f < df, e < ef, f < ef,

  ab < abc, ac < abc, bc < abc,
  ab < abd, ad < abd, bd < abd,
  ac < ace, ae < ace, ce < ace,
  ad < adf, af < adf, df < adf,
  ae < aef, af < aef, ef < aef,
  bc < bcf, bf < bcf, cf < bcf,
  bd < bde, be < bde, de < bde,
  be < bef, bf < bef, ef < bef,
  cd < cde, ce < cde, de < cde,
  cd < cdf, cf < cdf, df < cdf
};

// The square and the boundary of an octahedron meet at just the vertex o.
// Their union is an order complex for the wedge of a circle and a 2-sphere.
space Wedge {
  o < east, v < east, o < west, v < west,
  o < m, o < n, t < m, t < n, m < p, n < p, m < q, n < q
};

betti RP2;
betti Wedge;

// Both spaces have a unique nonzero class in degree one. Square it.
let squareRP = \u : H(RP2, 1) => cup(u, u);
let squareWedge = \u : H(Wedge, 1) => cup(u, u);
let x = squareRP(gen[RP2, 1, 0]);
let y = squareWedge(gen[Wedge, 1, 0]);
show x;
show y;
assert x == gen[RP2, 2, 0];
assert y == zero[Wedge, 2];

// The suspensions have the same Betti numbers and all positive-degree cup
// products vanish. Sq^1 still distinguishes their degree-two generators.
let alpha = gen[RP2, 1, 0];
let beta = gen[Wedge, 1, 0];
let test = \theta : A(1) => \v : H(Susp(RP2), 2) => act(theta, v);
betti Susp(RP2);
betti Susp(Wedge);
show test(sq[1])(susp(alpha));
show act(sq[1], susp(beta));
assert cup(susp(alpha), susp(alpha)) == zero[Susp(RP2), 4];
assert test(sq[1])(susp(alpha)) == susp(x);
assert act(sq[1], susp(beta)) == susp(y);

// Stability transports the action to the million-fold suspension of RP2
// with a disjoint basepoint. The result is a nonzero degree-1000002 class.
let distant = \v : Stable(RP2, 1000000, 1000001) => act(sq[1], v);
let result = distant(shift[1000000](alpha));
show result;
assert result == shift[1000000](x);
