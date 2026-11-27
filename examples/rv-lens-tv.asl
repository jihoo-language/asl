// In the doubled SU(2) level-two theory, equal |H_1| can give different values.
let x = rvLensTV[2,2,8,1];
let y = rvLensTV[2,2,8,3];
assert x == zero[RVCycloQ(16)];
assert y == rvCycloRat[1](rvCyclo[16,0]);
assert x == rvWzwTV[2,2]{S,T[8],S};
assert y == rvLensTV[2,2,8,-3];
show x;
show y;
// L(2,1) has the nonrational value (2-sqrt(2))/4.
let expected = rvCycloRat[4](add(rvCycloScale[2](rvCyclo[16,0]),add(neg(rvCyclo[16,2]),rvCyclo[16,6])));
let exact = rvLensTV[2,2,2,1];
assert exact == expected;
assert exact == rvWzwTV[2,2]{S,T[-2],S};
show exact;
// The vacuum normalization and L(p,1) at the abelian levels.
show rvLensTV[2,2,0,-1];
show rvLensTV[2,1,5,1];
show rvLensTV[2,1,6,1];
show rvLensTV[2,1,8,1];
show rvLensTV[3,1,5,1];
show rvLensTV[3,1,6,1];
