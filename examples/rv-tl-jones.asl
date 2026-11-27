// Three fundamental strands of SU(2) at level 3, total height 1.
let a = rvTLPath[3,[0,1,0,1]];
let b = rvTLPath[3,[0,1,2,1]];
let ua = rvTLBraid[[1,-2]](a);
let ub = rvTLBraid[[1,-2]](b);
let aa = rvTLCoord[[0,1,0,1]](ua);
let ba = rvTLCoord[[0,1,2,1]](ua);
let ab = rvTLCoord[[0,1,0,1]](ub);
let bb = rvTLCoord[[0,1,2,1]](ub);
let one = rvCycloRat[1](rvCyclo[20,0]);
let trace = rvTLTrace[3,3,1,[1,-2]];
assert trace == add(aa,bb);
assert add(rvCycloQMul(aa,bb),neg(rvCycloQMul(ba,ab))) == one;
let three = rvCycloRat[1](rvCycloScale[3](rvCyclo[20,0]));
assert add(add(rvCycloQMul(trace,trace),neg(rvCycloQMul(three,trace))),one) == zero[RVCycloQ(20)];
show trace;
show rvCycloQGalois[3](trace);

// The trefoil and figure-eight in the variable t = A^(-4).
// The exact coefficients below are reduced modulo the cyclotomic polynomial.
show rvTLJones[3,2,[1,1,1]];
show rvTLJones[3,3,[1,-2,1,-2]];
assert rvTLJones[3,2,[1,1,1]] == rvTLJones[3,3,[1,1,1,2]];
assert rvTLJones[3,2,[1,1,1]] == rvTLJones[3,3,[1,1,1,-2]];
