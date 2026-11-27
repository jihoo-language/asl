// Here the cross quartic has the scaled cyclotomic factor t²+3t+9.
rvunram F[3,1]{0,1};
let z = rvUnramElem[F]{};
let one = rvUnramElem[F]{1};
let two = rvUnramElem[F]{2};
rvhyperelliptic C[F]{z,two,z,two,one,one};
rvKummerK3Zeta[C];
rvKummerK3Points[C,1];
rvKummerK3Points[C,3];
rvKummerK3Picard[C,1];
rvKummerK3Picard[C,3];

// On this supersingular example the Jacobian Frobenius over F25 is −5.
rvunram W[5,1]{2,0,1};
let zw = rvUnramElem[W]{};
let ow = rvUnramElem[W]{1};
rvhyperelliptic S[W]{zw,ow,zw,zw,zw,ow};
rvKummerK3Zeta[S];
rvKummerK3Points[S,1];
rvKummerK3Picard[S,1];
