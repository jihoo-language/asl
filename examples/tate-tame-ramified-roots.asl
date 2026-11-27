rvunram Q[3,4]{0,1};
let quadratic = rvTatePoly[3,4](-3,0,1);
rvTateTame[Q,0](quadratic);
let p3 = rvTameElem[Q,2,8]([],[1]);
assert rvTameMul(p3,p3) == rvTameElem[Q,2,8]([3]);
assert rvTameEval(p3,quadratic) == zero[RVTame(Q,2,8)];
assert rvTameInertia(p3) == neg(p3);

rvunram U[5,4]{1,1,1};
let cubic = rvTatePoly[5,4](-30,0,0,1);
rvTateTame[U,0](cubic);
let p5 = rvTameElem[U,3,12]([],[1]);
let zeta = rvTameTeich[U,3,12](rvUnramElem[U]{0,6});
assert rvTameInertia(p5) == rvTameMul(zeta,p5);
assert rvTameFrob[1](rvTameInertia(p5)) == rvTameInertia(rvTameInertia(p5));
let root = rvTameMul(rvTameElem[U,3,12]([111]),p5);
rvTameVal(root);
rvTameVal(rvTameEval(root,cubic));
