rvunram F[5,1]{0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
let t = rvUnramElem[F]{2};
let m = rvUnramElem[F]{4};
rvhyperelliptic C[F]{o,o,z,z,z,o};

let point = rvJacPoint[C](z,o);
let double = add(point,point);
let image = rvJacKummerPoint(double);
assert image == rvJacKummerPoint(neg(double));
assert image == rvKummer[C](o,z,z,m);
assert image == rvKummer[C](t,z,z,rvUnramElem[F]{3});
assert rvKummerCoord[3](image) == m;
assert rvKummerLiftable(image) == S.hi;
assert rvKummerLift(image) == double;

let twist = rvKummer[C](z,o,o,o);
assert rvKummerLiftable(twist) == S.lo;
show image;
show twist;
