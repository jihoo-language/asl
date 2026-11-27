// The SU(3) level-two adjoint has golden-ratio quantum dimension.
let adj = rvFusionBasis[2,[0,1,2]];
let q = rvQDim(adj);
let one = rvCyclo[15,0];
assert rvCycloMul(q,q) == add(q,one);
assert rvCycloConj(q) == q;
show q;
// The same nontrivial torus label gives conjugate values for the two fundamentals.
let f = rvFusionBasis[2,[0,0,1]];
let a = rvVerlinde[[0,0,1]](f);
let b = rvVerlinde[[0,0,1]](rvFusionDual(f));
assert rvCycloConj(a) == b;
assert rvCycloMul(a,b) == rvVerlinde[[0,0,1]](rvFusionMul(f,rvFusionDual(f)));
show a;
show b;
// Sewing handles uses the ordinary product in the level-two fusion ring.
show rvFusionHandle[3,2];
rvConformalDim[3,2,0,[]];
rvConformalDim[3,2,1,[]];
rvConformalDim[3,2,2,[]];
rvConformalDim[3,2,3,[]];
rvConformalDim[3,2,4,[]];
rvConformalDim[3,2,5,[]];
rvConformalDim[3,2,0,[[0,0,1],[0,1,1]]];
rvConformalDim[3,2,0,[[0,0,1],[0,0,1]]];
