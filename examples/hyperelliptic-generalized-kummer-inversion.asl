// The complex Kummer surface and fourfold have different inversion actions.
rvGenKummerInversionEuler[2,1];
rvGenKummerInversionEuler[2,2];
rvGenKummerInversionEuler[3,3];
rvGenKummerInversionHodge[3,3,2,1];

// Full geometric 3-torsion of y²=x⁵+x over F25.
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

// The shift e1 moves under F5 Frobenius. The square retains a nonzero shift.
assert rvJacFrob[C,1](e1) == neg(e1);
let fixed = add(f2,rvJacFrob[C,1](f2));
assert rvJacFrob[C,1](fixed) == fixed;
assert rvJacTorsionCoord[B,1](fixed) == rvPrimeElt[3,1];
rvGenKummerInversionPoints[B,C,1](e1);
rvGenKummerInversionPoints[B,C,1](fixed);
rvGenKummerInversionPoints[B,C,2](e1);
rvGenKummerInversionPoints[B,C,2](fixed);

// Over F25 the odd cohomology detects inversion even when torsion is rational.
rvGenKummerInversionFrob[B,S,1,3](e1);
rvGenKummerInversionPoints[B,S,1](e1);
