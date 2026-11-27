// F3 Frobenius swaps the two roots of x²+1.
rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
rvjacbase L[C,W];
let zeroW = rvUnramElem[W]{};
let oneW = rvUnramElem[W]{1};
let twoW = neg(oneW);
let a = rvUnramElem[W]{0,1};
let r0 = rvJacPoint[L](zeroW,zeroW);
let r1 = rvJacPoint[L](oneW,zeroW);
let r2 = rvJacPoint[L](twoW,zeroW);
let r3 = rvJacPoint[L](a,zeroW);
let e1 = r0;
let e2 = add(add(r0,r1),r2);
let f1 = r1;
let f2 = add(add(r0,r1),r3);
rvjactorsionbasis B[L,2](twoW){e1,e2,f1,f2};
rvJacRichelotKernels[B,C];

// The first stable plane has a smooth target. The other two split geometrically.
let k0 = rvJacRichelotKernel[B,C,0,0];
let k1 = rvJacRichelotKernel[B,C,0,1];
let k2 = rvJacRichelotKernel[B,C,0,2];
assert k0 == e1;
assert k1 == e2;
assert k2 == add(e1,e2);
assert rvJacRichelotKernel[B,C,1,1] == f1;
rvrichelotfromlevel D[B,C,0];
assert rvRichelotPush[D](k0) == zero[RVJac(D)];
assert rvRichelotPush[D](k1) == zero[RVJac(D)];
assert rvRichelotPush[D](k2) == zero[RVJac(D)];
rvRichelotData[D];
