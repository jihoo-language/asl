rvunram W[3,1]{1,0,1};
let z=rvUnramElem[W]{};
let o=rvUnramElem[W]{1};
let i=rvUnramElem[W]{0,1};
rvhyperelliptic C[W]{z,neg(o),z,z,z,o};
let e1=rvJacPoint[C](z,z);
let f1=rvJacPoint[C](o,z);
let e2=add(add(e1,f1),rvJacPoint[C](neg(o),z));
let f2=add(add(e1,f1),rvJacPoint[C](i,z));
rvjactorsionbasis B[C,2](neg(o)){e1,e2,f1,f2};

rvog6 M(rvAbMukai[B,1,-1]{0,0,0,0,0,0});
let e=rvOG6[M,0](rvAbMukai[B,0,0]{1,0,0,0,0,0});
let f=rvOG6[M,0](rvAbMukai[B,0,0]{0,1,0,0,0,0});
let eDual=rvOG6[M,0](rvAbMukai[B,0,0]{0,0,0,0,0,-1});
let fDual=rvOG6[M,0](rvAbMukai[B,0,0]{0,0,0,0,1,0});
let hE=rvOG6[M,0](rvAbMukai[B,0,0]{0,0,1,0,0,0});
let hF=rvOG6[M,0](rvAbMukai[B,0,0]{0,0,0,-1,0,0});
let h=add(hE,hF);
let u=add(f,fDual);
rvog6cusp P(h,e,f);
rvog6nilpotentii A[P,3];
rvog6nilpotentiii B3[P,0](u);

// Left regular monodromy and conjugation on the marked polarized lattice.
rvOG6SpinData[A];
rvOG6SpinData[B3];
let s=rvOG6SpinLog[B3];
assert rvOG6CliffordMul(s,s)==zero[RVOG6Clifford(P,0)];
let t=rvOG6SpinExp[B3,2];
let tInverse=rvOG6SpinExp[B3,-2];
assert rvOG6CliffordMul(t,tInverse)==rvOG6CliffordBasis[P,0];
let c=rvOG6CliffordVector[P](eDual);
let left=rvOG6CliffordMul(t,c);
assert rvOG6CliffordMul(left,tInverse)==rvOG6CliffordVector[P](rvOG6NilpotentExp[B3,2](eDual));
