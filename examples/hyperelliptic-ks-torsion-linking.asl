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

// Complementary Wang torsion classes in degrees 2 and 64.
let w=rvKSExteriorBasis[P,1];
let u=rvKSExteriorBasis[P,2];
let v=rvKSExteriorBasis[P,18446744073709551614];
let pre=rvKSExteriorBasis[P,18446744073709551613];

// Linking in the oriented mapping torus.
assert rvKSExteriorDifference[A,1,1](w)==rvKSExteriorScale[3](u);
assert rvKSExteriorDifference[A,1,1](pre)==rvKSExteriorScale[3](v);
let value=rvKSLink[A,1,3,3](u,w,v,pre);
assert value==rvQZ[1,3];
assert rvKSLink[A,1,3,3](v,pre,u,w)==value;
rvQZOrder(value);
rvKSLinkCertificate[A,1,3,3](u,w,v,pre);
rvKSLinkCertificate[A,-2,6,6](u,neg(w),v,neg(pre));
