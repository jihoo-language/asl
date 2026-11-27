// Determinant factors are written in ascending powers of T.
let reduced = rvWittSpectrum{[1,0,-1],[1,-1]};
assert reduced == rvWittSpectrumLine[-1];
show reduced;
rvWittSpectrumGhost[1](reduced);
let pair = rvWittSpectrum{[1,-5,6]};
rvWittSpectrumGhost[2](pair);
show rvWittSpectrumExterior[2](pair);
assert add(pair,neg(pair)) == zero[RVWittSpectrum];

// A semilinear Specht summand and the Tate shift on a Hilbert–Chow diagonal.
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
rvk3autom Moving[B]{translate(8)};
rvhilbstratum Diagonal{(2,1)};
rvspechtstratum Intrinsic[Diagonal]{};
let middle = rvK3SpechtSpectrum[Moving,C,Intrinsic,1,2];
let shifted = rvK3StratumSpectrum[Moving,C,Diagonal,1,4];
assert shifted == rvWittSpectrumTensor(rvWittSpectrumLine[3],middle);
assert rvWittSpectrumAdams[2](middle) == rvK3SpechtSpectrum[Moving,C,Intrinsic,2,2];
rvWittSpectrumData(middle);
rvWittSpectrumGhost[2](middle);
rvWittSpectrumCoeff[2](shifted);
