rvunram F[3,1]{0,1};
let z = rvUnramElem[F]{};
let one = rvUnramElem[F]{1};
let two = rvUnramElem[F]{2};

// Ordinary Jacobian with ranks 10 and 20 over F3 and F27.
rvhyperelliptic C[F]{z,two,z,two,one,one};
rvKummerK3ArtinTate[C,1];
rvKummerK3ArtinTate[C,3];

// A supersingular Jacobian with a nonzero rank-one Cartier-Manin matrix.
// Geometric rank is 22. Full rational rank first occurs over F81.
rvhyperelliptic N[F]{one,z,z,z,z,one};
rvKummerK3ArtinTate[N,1];
rvKummerK3ArtinTate[N,2];
rvKummerK3ArtinTate[N,4];
rvKummerK3ArtinTate[N,8];

// A superspecial Jacobian with scalar Frobenius over F25.
rvunram W[5,1]{2,0,1};
let zw = rvUnramElem[W]{};
let ow = rvUnramElem[W]{1};
rvhyperelliptic S[W]{zw,ow,zw,zw,zw,ow};
rvKummerK3ArtinTate[S,1];
