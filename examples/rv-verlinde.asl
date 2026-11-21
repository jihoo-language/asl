// The adjoint tensor square changes with the affine level.
let adj2 = rvFusionBasis[2,[0,1,2]];
let adj3 = rvFusionBasis[3,[0,1,2]];
let adj4 = rvFusionBasis[4,[0,1,2]];
let square2 = rvFusionMul(adj2,adj2);
let square3 = rvFusionMul(adj3,adj3);
let square4 = rvFusionMul(adj4,adj4);
assert square2 == add(rvFusionOne[3,2],adj2);
show square2;
show square3;
show square4;
rvFusionCoeff[[0,2,4]](square4);
// The fundamental class at level one has order three under multiplication.
let f = rvFusionBasis[1,[-3,-3,-2]];
let dual = rvFusionDual(f);
assert rvFusionMul(f,f) == dual;
assert rvFusionMul(f,dual) == rvFusionOne[3,1];
assert rvFusionDual(dual) == f;
