let u = rvSphBasis[[0,1]];
let a = rvSatake(u);
let b = rvSatake(rvSphBasis[[0,2]]);
let det = rvSatake(rvSphBasis[[1,1]]);

// q=v^2. The mixed coefficient of S_(0,2) is v^2-1.
rvSatakeCoeff[[1,0]](a);
rvSatakeCoeff[[0,1]](a);
rvSatakeCoeff[[1,1]](b);
rvSatakeCoeff[[2,0]](b);
rvSatakeCoeff[[1,1]](det);
assert rvSatake(rvSphQ[2]) == rvSatV[2,2];
assert rvSatakeMul(a,a) == rvSatake(rvSphMul(u,u));

// A negative central translation multiplies by a Laurent determinant.
let negative = rvSphBasis[[-3,-3]];
let translated = rvSatake(rvSphMul(negative,u));
assert translated == rvSatake(rvSphBasis[[-3,-2]]);
assert translated == rvSatakeMul(rvSatake(negative),a);
rvSatakeAt[2,[2,3]](translated);
rvSatakeAt[2,[2,3]](rvSatakeMul(rvSatV[2,-2],a));

let x = rvSphBasis[[0,1,1]];
let y = rvSphBasis[[0,0,1]];
assert rvSatake(rvSphMul(x,y)) == rvSatakeMul(rvSatake(x),rvSatake(y));
rvSatakeAt[2,[2,3,5]](rvSatake(rvSphMul(x,y)));
rvSatakeCoeff[[1,1,0]](rvSatake(rvSphBasis[[0,0,2]]));
