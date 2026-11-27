rvunram B[3,7]{0,1};
rvcyclotomiclocal K[B,1];
rvcubicdiv D[K,2,8];

let one = rvCubicDivElem[D]([1]);
let z = zero[RVCubicDiv(D)];
let theta = rvCubicDivTheta[D];
let k0 = rvCubicTreeMatrix(one,z,z,one);
let k1 = rvCubicTreeMatrix(one,one,one,theta);
rvcubicschottky G[D]{[1,k0],[2,k1]};

rvCubicSchottkyData[G];
let a = rvCubicSchottkyWord[G,[1]];
let b = rvCubicSchottkyWord[G,[2]];
let ab = rvCubicSchottkyMul(a,b);
let ba = rvCubicSchottkyMul(b,a);
let w = rvCubicSchottkyWord[G,[1,2,-1]];
let o = rvCubicTreeRoot[D];
assert rvCubicSchottkyConjugate(w,b) == S.hi;
assert rvCubicSchottkyCyclic(w) == b;
assert rvCubicTreeLca(rvCubicSchottkyApply(ab,o),rvCubicSchottkyApply(ba,o)) == o;
rvCubicSchottkyLengths(w);
rvCubicSchottkyPath(w);
rvCubicSchottkyAxisCheck(w);
rvCubicSchottkyPingPong[G,1](rvCubicSchottkyDirection[G,2]);

rvCubicSchottkyCore[G];
rvCubicSchottkyNB[G,0,1];
rvCubicSchottkyWeighted[G,1,2];
rvCubicIhara[G];
rvCubicIharaCoeff[G,2];
rvCubicSchottkyPrimitive[G,8];
rvCubicSchottkyEnumerate[G,8];
