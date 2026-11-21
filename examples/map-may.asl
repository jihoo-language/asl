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

// The target coaction detects Sq^1 taking the degree-one class to the top.
let top = mayClass[0, 0](mword[One, 0, 0, RP2, 2, 0](cword{}));
let topD = mayClass[0, 0](mword[One, 0, 0, RP2, 1, 0](cword{[1]}));
assert mayD(top) == topD;
assert topD == mapmaygen[One, RP2, 0, 1, 2, 0, 0];
let wedgeTop = mayClass[0, 0](mword[One, 0, 0, Wedge, 2, 0](cword{}));
assert mayD(wedgeTop) == zero[MapMay(One, Wedge, 0, 1, 2, 0)];

// Transposing the source coaction detects the same link in signed degree -1.
let down = mayClass[0, 0](mword[RP2, 1, 0, One, 0, 0](cword{}));
let downD = mayClass[0, 0](mword[RP2, 2, 0, One, 0, 0](cword{[1]}));
assert mayD(down) == downD;
assert downD == mapmaygen[RP2, One, 0, 1, -1, 0, 0];
let wedgeDown = mayClass[0, 0](mword[Wedge, 1, 0, One, 0, 0](cword{}));
assert mayD(wedgeDown) == zero[MapMay(Wedge, One, 0, 1, -1, 0)];

// The two nonzero summands in the Leibniz identity cancel.
let h0 = mapmaygen[One, One, 0, 1, 1, 0, 0];
assert ecomp(downD, top) == h0;
assert ecomp(down, topD) == h0;
assert mayD(ecomp(down, top)) == add(ecomp(downD, top), ecomp(down, topD));

let idRP = mayClass[0, 0](eRep(eid[RP2]));
let idOne = mayClass[0, 0](eRep(eid[One]));
assert ecomp(idRP, top) == top;
assert ecomp(top, idOne) == top;
assert mayClass[0, 0](mayRep(downD)) == downD;

// Actual maps give filtration-zero cycles, and their classes compose on later pages.
let intoRP = \point : One => RP2.a;
let collapseRP = \point : RP2 => One.pt;
let intoMay = mayClass[2, 0](eRep(emap(intoRP)));
let collapseMay = mayClass[2, 0](eRep(emap(collapseRP)));
let composed = mayClass[2, 0](eRep(emap(\point : One => collapseRP(intoRP(point)))));
assert ecomp(collapseMay, intoMay) == composed;
assert composed == mayClass[2, 0](eRep(eid[One]));

mapmaydim[One, RP2, 0, 0, 2, 0];
mapmaydim[One, RP2, 1, 0, 2, 0];
mapmaydim[One, Wedge, 1, 0, 2, 0];
mapmaydim[RP2, One, 1, 0, -1, 0];
mapmaydim[Wedge, One, 1, 0, -1, 0];
mapmaydim[Wedge, One, 9, 0, -1, 0];
mapextdim[Wedge, One, 0, -1];

// In signed degree -1, the stabilized slices of this length-one group
// have one class each in filtrations zero and one.
mapmaydim[Susp(Wedge), One, 9, 1, -1, 0];
mapmaydim[Susp(Wedge), One, 9, 1, -1, 1];
mapextdim[Susp(Wedge), One, 1, -1];
