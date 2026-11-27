// The spectrum of a three-cycle has determinant 1-T^3.
let cycle = rvWittSpectrumVerschiebung[3](rvWittSpectrumLine[1]);
show cycle;
rvWittSpectrumGhost[6](cycle);
rvWittSpectrumNecklace[3](cycle);
rvWittSpectrumNecklace[6](cycle);

// A marked translation moved by Frobenius on the resolved Kummer surface.
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
rvk3autom Idle[B]{};
rvK3FrobHilbOrbits[Idle,C,1,2];
rvK3FrobHilbOrbits[Moving,C,1,2];
rvK3FrobHilbOrbits[Moving,C,2,2];
let factor = rvK3FrobHilbSpectrum[Moving,C,1,1];
rvWittSpectrumData(factor);
rvWittSpectrumNecklace[2](factor);
assert rvWittSpectrumAdams[2](factor) == rvK3FrobHilbSpectrum[Moving,C,1,2];
