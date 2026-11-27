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
let u=add(f,fDual);
rvog6cusp P(h,e,f);
rvog6nilpotentii A[P,3];
rvog6nilpotentiii B3[P,0](u);

// Exterior monodromy on cohomology.
rvKSExteriorData[A,1];
rvKSExteriorData[B3,2];
rvKSExteriorData[A,32];
let x=rvKSExteriorEmbed(rvOG6CliffordBasis[P,0]);
let y=rvKSExteriorBasis[P,2];
let xdy=rvKSExteriorCup(x,rvKSExteriorPower[A,1](y));
let dxy=rvKSExteriorCup(rvKSExteriorPower[A,1](x),y);
let xy=rvKSExteriorCup(x,y);
assert rvKSExteriorPower[A,1](xy)==add(dxy,xdy);
assert rvKSExteriorToSpin(rvKSExteriorPower[A,1](x))==rvOG6SpinLog[A];
assert rvKSExteriorExp[B3,4](xy)==rvKSExteriorCup(rvKSExteriorExp[B3,4](x),rvKSExteriorExp[B3,4](y));
assert rvKSExteriorExp[B3,-4](rvKSExteriorExp[B3,4](xy))==xy;
let top=rvKSExteriorBasis[P,18446744073709551615];
assert rvKSExteriorExp[A,4](top)==top;
