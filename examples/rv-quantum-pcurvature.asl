// In characteristic two, the Gr(2,4) spectrum vanishes but the operator does not.
let one2 = rvQGrPBasis[2,2,4,[]];
let v2 = rvQGrCurvature(one2);
let w2 = rvQGrCurvature(v2);
assert w2 == add(rvQGrPBasis[2,2,4,[2]],rvQGrPBasis[2,2,4,[1,1]]);
assert rvQGrCurvature(w2) == zero[RVQGrP(2,2,4)];
assert v2 == rvQGrCurvatureSatake(one2);
show v2;
show w2;
rvQGrCurvatureChar[2,2,4];
// At odd primes the characteristic polynomial has nonconstant q coefficients.
let one3 = rvQGrPBasis[3,2,4,[]];
let v3 = rvQGrCurvatureCol[3,2,4,[]];
assert v3 == rvQGrCurvatureSatake(one3);
assert rvQGrCurvature(rvQGrShift[1](one3)) == rvQGrShift[1](v3);
assert rvQGrCurvature(rvQGrNabla(one3)) == rvQGrNabla(v3);
show v3;
rvQGrCurvatureCoeff[3,2,4,2];
rvQGrCurvatureHorizontal[3,2,4];
rvQGrCurvatureCoeff[5,3,6,14];
rvQGrCurvatureCoeff[5,3,6,8];
rvQGrCurvatureCoeff[5,3,6,2];
// The prime can divide the ambient dimension. Duality preserves the connection.
let comp = rvQGrPBasis[3,4,6,[2,1]];
assert rvQGrCurvature(comp) == rvQGrCurvatureSatake(comp);
assert rvQGrTranspose(rvQGrCurvature(comp)) == rvQGrCurvature(rvQGrTranspose(comp));
rvQGrCurvatureCoeff[5,1,5,0];
