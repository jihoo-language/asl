rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let z = rvUnramElem[F]{};
let one = rvUnramElem[F]{1};
let two = rvUnramElem[F]{2};
rvhyperelliptic C[F]{z,one,z,z,z,one};
rvjacbase L[C,W];
let p = rvJacPoint[C](two,one);

rvJacHilbPoints[C,0,1];
rvJacHilbPoints[C,3,1];
rvGenKummerPoints[C,2];
rvKummerK3Points[C,1];
rvGenKummerPoints[C,3];
rvGenKummerFiber[C,3](p);
rvGenKummerPoints[C,4];
rvGenKummerFiber[C,4](p);

// The geometric two-torsion becomes rational in W. The base Frobenius
// determines the torsion order for each extension degree.
let zw = rvUnramElem[W]{};
let mw = rvUnramElem[W]{2};
let t0 = rvJacPoint[L](zw,zw);
let t1 = rvJacPoint[L](rvUnramElem[W]{1,1},zw);
let t2 = rvJacPoint[L](rvUnramElem[W]{1,2},zw);
let t3 = rvJacPoint[L](rvUnramElem[W]{2,1},zw);
rvjactorsionbasis B[L,2](mw){t0,add(add(t0,t1),t2),t1,add(add(t0,t1),t3)};
rvGenKummerLevel[B,C,1];
rvGenKummerLevel[B,C,2];
rvGenKummerPoints[L,2];
rvKummerK3Points[C,2];
