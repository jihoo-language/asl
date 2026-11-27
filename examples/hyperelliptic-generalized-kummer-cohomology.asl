// The generalized Kummer fourfold has genuine odd cohomology.
rvGenKummerBetti[3,2];
rvGenKummerBetti[3,3];
rvGenKummerBetti[3,4];
rvJacHilbBetti[2,3];

rvunram F[5,1]{0,1};
let zF = rvUnramElem[F]{};
let oF = rvUnramElem[F]{1};
rvhyperelliptic C[F]{zF,oF,zF,zF,zF,oF};
rvunram W[5,1]{3,0,1};
rvjacbase S[C,W];
let o = rvUnramElem[W]{1};
let a = rvUnramElem[W]{0,1};
let e1 = rvJacScale[2](rvJacPoint[S](o,a));
let e2 = add(e1,rvJacScale[2](rvJacPoint[S](rvUnramElem[W]{1,2},a)));
let f1 = rvJacScale[4](rvJacPoint[S](rvUnramElem[W]{1,1},a));
let f2 = rvJacScale[4](rvJacPoint[S](rvUnramElem[W]{1,3},a));
rvjactorsionbasis B[S,3](rvUnramElem[W]{2,2}){e1,e2,f1,f2};

// Frobenius over F25 acts by -5 on H^1 of the Jacobian.
rvGenKummerTrace[B,C,2,2];
rvGenKummerTrace[B,C,2,3];
rvGenKummerTrace[B,C,2,4];
rvJacHilbTrace[C,2,2,3];
rvGenKummerLevel[B,C,2];
