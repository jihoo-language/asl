rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let z = rvUnramElem[F]{};
let one = rvUnramElem[F]{1};
let two = rvUnramElem[F]{2};
rvhyperelliptic C[F]{z,one,z,one,two,one};
rvjacbase CL[C,W];
rvrichelot DL[CL]{[0,1],[1],[1]}{[0,2],[1],[1]}{[],[1]};
rvrichelotdescend D[C,DL];
rvRichelotDescent[D];

// The x=1 point lifts in F9 with a Frobenius sign change.
let k = rvKummer[C](z,one,one,one);
assert rvKummerLiftable(k) == S.lo;
let extension = rvKummerExtend[CL](k);
let divisor = rvKummerLift(extension);
assert rvJacFrob[C,1](divisor) == neg(divisor);
assert rvKummerFrob[C,1](extension) == extension;
assert rvKummerDescend[C](rvJacKummerPoint(divisor)) == k;

let image = rvKummerRichelotDesc[D](k);
assert rvKummerLiftable(image) == S.lo;
assert rvKummerExtend[DL](image) == rvJacKummerPoint(rvRichelotPush[DL](divisor));
show image;
