// The diagonal stratum of length two and the two-double-point stratum of length four.
rvhilbstratum Double{(2,1)};
rvhilbstratum TwoDoubles{(2,1),(2,1)};
rvHilbStratumData[Double];
rvHilbStratumData[TwoDoubles];
rvK3StratumHodge[Double,1,1];
rvJacStratumHodge[Double,2,1];
rvGenKummerStratumHodge[Double,1,1,1];
rvGenKummerStratumHodge[Double,2,1,1];
rvGenKummerStratumHodge[TwoDoubles,2,2,2];
rvGenKummerStratumHodge[TwoDoubles,4,2,2];

// Two integral Nakajima classes of the same grade with different mode lengths.
rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
rvjacbase L[C,W];
let zw = rvUnramElem[W]{};
let ow = rvUnramElem[W]{1};
let a = rvUnramElem[W]{0,1};
let r0 = rvJacPoint[L](zw,zw);
let r1 = rvJacPoint[L](ow,zw);
let r2 = rvJacPoint[L](neg(ow),zw);
let ra = rvJacPoint[L](a,zw);
rvjactorsionbasis B[L,2](neg(ow)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),ra)};
let diagonal = rvNakajimaMonomial[B]{(2,0,1)};
let open = rvNakajimaMonomial[B]{(1,0,1),(1,1,1)};
let part = rvNakajimaStratum[Double](add(diagonal,open));
assert part == diagonal;
rvNakajimaPair(part,rvNakajimaMonomial[B]{(2,23,1)});

// Frobenius has eight fixed node labels here. A moving shift can remove them all.
let moving = add(add(r0,r1),ra);
rvGenKummerStratumFrob[B,C,Double,1,2](zero[RVJac(L)]);
rvGenKummerStratumFrob[B,C,Double,1,2](moving);
