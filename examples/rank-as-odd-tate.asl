// The wild deck action on H^1_dR of y^3 - y = x^7.
rvasodd Curve[3] {(7,1)};
rvASOddTateJordan[Curve];
rvASOddTatedim[Curve,0];
rvASOddTatedim[Curve,-1];

let x = rvASOddTategen[Curve,0,1];
let dual = rvASOddTategen[Curve,-1,3];
let y = rvASOddTategen[Curve,0,4];
show rvASOddTateF(y);
show rvASOddTateV(y);
show rvASOddTatePair(x,dual);
show rvASOddTatePair(dual,x);

assert rvASOddTateV(x) == zero[ASOddTate(Curve,0)];
assert rvASOddTateV(dual) == zero[ASOddTate(Curve,-1)];
assert rvASOddTateF(rvASOddTateF(y)) == zero[ASOddTate(Curve,0)];
assert rvASOddTateClass[0](rvASOddTateRep(y)) == y;

// Nonzero pairing against dual, which V kills, certifies x lies outside im F.
assert rvASOddTatePair(x,dual) == rvPrimeElt[3,1];
