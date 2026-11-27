// Frobenius has both fixed nodes and node cycles of length three.
rvunram F[3,1]{0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
let t = rvUnramElem[F]{2};
rvhyperelliptic C[F]{z,t,z,t,o,o};
rvKummerHilbPoints[C,2,1];
rvKummerHilbPoints[C,3,1];
rvKummerHilbTrace[C,2,1,2];
rvKummerHilbTrace[C,2,1,4];
rvKummerHilbTrace[C,2,1,6];
rvKummerHilbBetti[2,4];
rvKummerHilbBetti[3,4];

// Frobenius on the superspecial surface over F25 is scalar on H².
rvunram W[5,1]{2,0,1};
let zw = rvUnramElem[W]{};
let ow = rvUnramElem[W]{1};
rvhyperelliptic S[W]{zw,ow,zw,zw,zw,ow};
rvKummerHilbTrace[S,2,1,4];
rvKummerHilbPoints[S,2,1];
