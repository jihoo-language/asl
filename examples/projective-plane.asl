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

betti RP2;
