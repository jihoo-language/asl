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
let fDual=rvOG6[M,0](rvAbMukai[B,0,0]{0,0,0,0,1,0});
let hE=rvOG6[M,0](rvAbMukai[B,0,0]{0,0,1,0,0,0});
let hF=rvOG6[M,0](rvAbMukai[B,0,0]{0,0,0,-1,0,0});
let h=add(hE,hF);
rvog6cusp P(h,e,f);
rvog6nilpotentii A[P,3];
rvog6nilpotentiii B3[P,0](add(f,fDual));
rvog6nilpotentiii Mixed[P,0](add(f,add(fDual,fDual)));

// Integral divided powers and Smith spectra of the logarithm.
rvKSExteriorSmith[A,2,2];
rvKSExteriorSmith[B3,32,4];
rvKSExteriorSmith[Mixed,2,2];
rvog6nilpotentii Huge[P,-170141183460469231731687303715884105728];
rvKSExteriorSmith[Huge,32,16];
let blade=rvKSExteriorBasis[P,3];
let divided=rvKSExteriorDivided[B3,2](blade);
assert rvKSExteriorPower[B3,2](blade)==rvKSExteriorScale[2](divided);
