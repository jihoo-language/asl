// Over One, differential and composition use the ordinary cobar algebra.
let a = mword[One, 0, 0, One, 0, 0](cword{[0, 1]});
let b = mword[One, 0, 0, One, 0, 0](cword{[0, 0, 1]});
assert eD(eD(a)) == zero[MapCobar(One, One, 3, 3)];
assert eD(ecomp(b, a)) == add(ecomp(eD(b), a), ecomp(b, eD(a)));
assert eClass(eD(a)) == zero[MapExt(One, One, 2, 3)];
show eD(a);

// Changing only the coefficient degree can give a signed internal shift.
space Circle { a < ab, b < ab, b < bc, c < bc, c < ca, a < ca };
let insert = mword[One, 0, 0, Circle, 1, 0](cword{});
let select = mword[Circle, 1, 0, One, 0, 0](cword{});
assert ecomp(select, insert) == eRep(eid[One]);
assert eClass(eRep(eid[Circle])) == eid[Circle];
show select;

// A triangulated real projective plane has linked H_1 and H_2 generators.
space RP2 {
  a < ab, b < ab, a < ac, c < ac, a < ad, d < ad,
  a < ae, e < ae, a < af, f < af, b < bc, c < bc,
  b < bd, d < bd, b < be, e < be, b < bf, f < bf,
  c < cd, d < cd, c < ce, e < ce, c < cf, f < cf,
  d < de, e < de, d < df, f < df, e < ef, f < ef,
  ab < abc, ac < abc, bc < abc, ab < abd, ad < abd, bd < abd,
  ac < ace, ae < ace, ce < ace, ad < adf, af < adf, df < adf,
  ae < aef, af < aef, ef < aef, bc < bcf, bf < bcf, cf < bcf,
  bd < bde, be < bde, de < bde, be < bef, bf < bef, ef < bef,
  cd < cde, ce < cde, de < cde, cd < cdf, cf < cdf, df < cdf
};
let degreeZero = mword[RP2, 0, 0, RP2, 0, 0](cword{});
let degreeOne = mword[RP2, 1, 0, RP2, 1, 0](cword{});
let degreeTwo = mword[RP2, 2, 0, RP2, 2, 0](cword{});
let linked = mword[RP2, 2, 0, RP2, 1, 0](cword{[1]});
assert eD(degreeOne) == linked;
assert eD(degreeTwo) == linked;
let identity = eRep(emap(\x : RP2 => x));
assert identity == add(degreeZero, add(degreeOne, degreeTwo));
assert eD(identity) == zero[MapCobar(RP2, RP2, 1, 0)];
assert eClass(identity) == eid[RP2];
