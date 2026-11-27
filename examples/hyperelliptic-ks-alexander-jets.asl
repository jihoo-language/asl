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

// Completed Alexander exponents and finite character jets.
rvKSAlexander[B3,1,2,Q];
rvKSAlexander[B3,1,2,2];
rvKSAlexander[B3,1,32,2];
rvKSAlexanderJets[B3,1,2,2,Q];
rvKSAlexanderJets[B3,1,2,2,2];

let u=rvKSExteriorBasis[P,2];
let w=rvKSExteriorBasis[P,1];
let v=rvKSExteriorBasis[P,18446744073709551614];
let pre=rvKSExteriorBasis[P,18446744073709551613];
rvKSAlexanderWitness[A,1,Q](rvKSExteriorScale[3](u),neg(w))(rvKSExteriorScale[3](v),pre);
rvKSAlexanderWitness[A,1,2](rvKSExteriorScale[3](u),neg(w))(rvKSExteriorScale[3](v),pre);
