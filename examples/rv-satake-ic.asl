let adjoint = rvSphIC[[0,1,2]];
let q = rvSphQ[3];
let i = rvSphId[3];
let standard = add(rvSphBasis[[0,1,2]],rvSphMul(add(q,i),rvSphBasis[[1,1,1]]));
assert adjoint == standard;
rvSphCoeff[[1,1,1]](adjoint);
rvSatakeCoeff[[1,1,1]](rvSatSchur[[0,1,2]]);
rvSatakeCoeff[[1,1,1]](rvSatake(adjoint));
assert rvSatake(adjoint) == rvSatakeMul(rvSatV[3,4],rvSatSchur[[0,1,2]]);

// The rank-three adjoint character has a repeated summand in its square.
let q2 = rvSphMul(q,q);
let q4 = rvSphMul(q2,q2);
let left = rvSphMul(adjoint,adjoint);
let highest = rvSphIC[[0,2,4]];
let middle = rvSphMul(q,add(rvSphIC[[1,1,4]],rvSphIC[[0,3,3]]));
let repeated = rvSphMul(q2,rvSphIC[[1,2,3]]);
let bottom = rvSphMul(q4,rvSphIC[[2,2,2]]);
assert left == add(add(highest,middle),add(add(repeated,repeated),bottom));
rvSatakeAt[2,[2,3,5]](rvSatake(left));

let u = rvSphBasis[[0,1]];
let icExpansion = add(rvSphIC[[0,2]],rvSphMul(rvSphQ[2],rvSphIC[[1,1]]));
assert rvSphMul(u,u) == icExpansion;
assert rvSatakeInverse(rvSatake(icExpansion)) == icExpansion;
