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

// Trace projector in the full matrix algebra.
rvKSVolumeSquare[P];
let ie=rvKSEmbed[P](e);
let ieDual=rvKSEmbed[P](eDual);
rvKSEndTracePair(ie,ieDual);
assert rvKSProject(ie)==ie;
let matrix=rvKSEndUnit[P,0,50];
let projection=rvKSProject(matrix);
assert rvKSProject(projection)==projection;
rvKSEndCoeff[0,50](projection);
rvKSEndTracePair(add(matrix,neg(projection)),ie);

let S=rvKSEndLeft(rvOG6SpinLog[B3]);
let bracket=add(rvKSEndCompose(S,ieDual),neg(rvKSEndCompose(ieDual,S)));
assert bracket==rvKSEmbed[P](rvOG6NilpotentPower[B3,1](eDual));
