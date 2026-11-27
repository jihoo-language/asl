rvunram B[3,7]{0,1};
rvcyclotomiclocal K[B,1];
let l = rvEisensteinPi[K,14];
let one = rvEisensteinElem[K,14]([1]);
let z = add(one,l);
let u2 = add(one,rvEisensteinMul(l,l));

rvcubicalgebra D[K,8](z,u2);
rvCubicAlgData[D];
let x = rvCubicAlgElem[D]([1],[2],[],[1],[],[],[],[],[1]);
rvCubicReducedChar(x);
show rvCubicReducedNorm(x);
assert rvCubicAlgMul(x,rvCubicInverse(x)) == rvCubicAlgBasis[D,0,0];

let b = add(one,neg(z));
rvcubicalgebra A[K,8](z,b);
rvcubicsplit S[A]([1],[-1]);
rvCubicAlgData[A];
rvCubicSplitData[S];
let u = rvCubicAlgBasis[A,1,0];
let v = rvCubicAlgBasis[A,0,1];
let t = add(rvCubicAlgBasis[A,0,0],neg(add(u,v)));
let m = rvCubicSplitMat[S](t);
show m;
assert rvCubicSplitDet(m) == rvCubicReducedNorm(t);
assert rvCubicAlgMul(t,rvCubicAdj(t)) == zero[RVCubicAlg(A)];
