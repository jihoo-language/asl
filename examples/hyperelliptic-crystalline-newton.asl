rvunram F[3,1]{0,1};
let z = rvUnramElem[F]{};
let one = rvUnramElem[F]{1};
rvhyperelliptic O[F]{z,neg(one),z,neg(one),one,one};
rvhyperelliptic P[F]{z,one,one,one,one,one};
rvhyperelliptic S[F]{one,z,z,z,z,one};

rvKummerK3Height[O];
rvKummerK3Height[P];
rvKummerK3Height[S];
rvKummerHilbNewton[O,2,4];
rvKummerHilbNewton[P,2,4];
rvKummerHilbNewton[S,2,4];
rvJacHilbNewton[O,2,3];

// At length three use a field of characteristic prime to the length.
rvunram W[5,1]{0,1};
let zw = rvUnramElem[W]{};
let ow = rvUnramElem[W]{1};
rvhyperelliptic A[W]{zw,ow,zw,ow,zw,ow};
rvhyperelliptic B[W]{zw,ow,zw,zw,zw,ow};
rvGenKummerNewton[A,3,3];
rvGenKummerNewton[B,3,3];
