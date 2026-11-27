// The marking lives over F9. The even mask 2 descends to the curve over F3.
rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
rvjacbase L[C,W];
let zw = rvUnramElem[W]{};
let ow = rvUnramElem[W]{1};
let i = rvUnramElem[W]{0,1};
let e1 = rvJacPoint[L](zw,zw);
let f1 = rvJacPoint[L](ow,zw);
let e2 = add(add(e1,f1),rvJacPoint[L](neg(ow),zw));
let f2 = add(add(e1,f1),rvJacPoint[L](i,zw));
rvjactorsionbasis B[L,2](neg(ow)){e1,e2,f1,f2};

rvKummerTranslateElliptic[B,0,1,0,0];
rvKummerTranslateElliptic[B,9,1,0,0];
rvKummerTranslateElliptic[B,9,1,1,0];
rvKummerTranslateElliptic[B,9,2,0,0];

rvKummerEnriquesElliptic[B,C,2,1,0,0];
rvKummerEnriquesElliptic[B,C,2,2,2,-2];
rvKummerEnriquesElliptic[B,C,2,4,0,0];
