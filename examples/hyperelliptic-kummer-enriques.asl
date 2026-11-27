// Over F3 the two branch points at ±i are conjugate. Mask 2 is even and fixed.
rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let z=rvUnramElem[F]{};
let o=rvUnramElem[F]{1};
rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
rvjacbase L[C,W];
let zw=rvUnramElem[W]{};
let ow=rvUnramElem[W]{1};
let i=rvUnramElem[W]{0,1};
let e1=rvJacPoint[L](zw,zw);
let f1=rvJacPoint[L](ow,zw);
let e2=add(add(e1,f1),rvJacPoint[L](neg(ow),zw));
let f2=add(add(e1,f1),rvJacPoint[L](i,zw));
rvjactorsionbasis B[L,2](neg(ow)){e1,e2,f1,f2};

rvKummerEnriquesPoints[B,C,2,1];
rvKummerEnriquesPoints[B,C,2,2];
rvKummerEnriquesPoints[B,L,2,1];
rvKummerEnriquesZeta[B,C,2];
rvKummerEnriquesZeta[B,L,2];
