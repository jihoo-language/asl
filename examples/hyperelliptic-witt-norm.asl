// A virtual negative line changes parity under quadratic tensor induction.
let odd = neg(rvWittSpectrumLine[1]);
let induced = rvWittSpectrumNorm[2](odd);
assert induced == rvWittSpectrumLine[-1];
show induced;
rvWittSpectrumNormGhost[2,1](odd);
rvWittSpectrumNormGhost[2,2](odd);

// Restrict a resolved surface from F9 to F3. Its full determinant has degree
// 576, but the norm ghost can be evaluated directly.
rvunram F[3,1]{0,1};
rvunram V[3,1]{1,0,1};
let z=rvUnramElem[F]{};
let o=rvUnramElem[F]{1};
rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
rvjacbase L[C,V];
let zv=rvUnramElem[V]{};
let ov=rvUnramElem[V]{1};
let a=rvUnramElem[V]{0,1};
let e1=rvJacPoint[L](zv,zv);
let f1=rvJacPoint[L](ov,zv);
let e2=add(add(e1,f1),rvJacPoint[L](neg(ov),zv));
let f2=add(add(e1,f1),rvJacPoint[L](a,zv));
rvjactorsionbasis B[L,2](neg(ov)){e1,e2,f1,f2};
rvk3autom Idle[B]{};
let extensionSpectrum = rvK3FrobHilbSpectrum[Idle,C,1,2];
rvWittSpectrumNormGhost[2,1](extensionSpectrum);
rvKummerHilbWeilRestrictPoints[C,1,2,1];
rvWittSpectrumNormGhost[2,2](extensionSpectrum);
rvKummerHilbWeilRestrictPoints[C,1,2,2];
