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

// Integral mapping torus and actual finite differences.
let x=rvKSExteriorBasis[P,65];
let difference=rvKSExteriorDifference[A,1,1](x);
assert difference==add(rvKSExteriorExp[A,1](x),neg(x));
rvKSExteriorCoeff[130](difference);
rvKSExteriorCoeff[130](rvKSExteriorPower[A,1](x));
rvKSExteriorGroupSmith[A,2,-2,2];
rvKSMappingTorus[A,-170141183460469231731687303715884105728,2];
rvKSMappingTorus[B3,0,32];
