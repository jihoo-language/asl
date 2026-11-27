// A genuinely positive-degree multiplicity in Gr(3,6).
let x = rvQGrBasis[3,6,[2,1]];
let y = rvQGrBasis[3,6,[3,2,1]];
let z = rvQGrMul(x,y);
show z;
rvQGrCoeff[1,[2,1]](z);
assert rvQGrClassical(z) == rvQGrBasis[3,6,[3,3,3]];
rvQGrGW[3,6,1,[2,1],[3,2,1],[3,2,1]];
// A determinant class contributes a quantum term in Gr(2,4).
let u = rvQGrOne[2,4];
assert rvQGrMul(rvQGrBasis[2,4,[2]],rvQGrBasis[2,4,[1,1]]) == rvQGrShift[1](u);
rvQGrGW[2,4,1,[2],[1,1],[2,2]];
// Projective space and its complementary Grassmannian both satisfy h^37=q.
let h = rvQGrBasis[1,37,[1]];
let h36 = rvQGrBasis[1,37,[36]];
let q = rvQGrShift[1](rvQGrOne[1,37]);
assert rvQGrMul(h,h36) == q;
assert rvQGrMul(rvQGrTranspose(h),rvQGrTranspose(h36)) == rvQGrTranspose(q);
show rvQGrTranspose(q);
