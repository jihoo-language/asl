// The regularized unit period of Gr(2,5) is the Apéry sequence for zeta(2).
let l = rvQGrAnnihilator[2,5];
let f = rvOreFactorial[5](l);
let d = rvOreEulerFactor(f);
let b = rvOreEulerReduce(f);
assert rvOreMul(d,b) == f;
let x = rvOre[0,1];
let y = add(rvOre[0,1],rvOre[0,0]);
let xm = add(rvOre[0,1],rvOre[0,0,-1]);
let q = rvOre[1,0];
let a = add(add(rvOre[0,2,11],rvOre[0,1,11]),rvOre[0,0,3]);
assert b == add(rvOreMul(x,x),neg(add(rvOreMul(q,a),rvOreMul(rvOre[2,0],rvOreMul(y,y)))));
// The factor has roots at 0 and 1. The normalized vector period makes both
// corresponding coefficients of B G vanish.
assert d == rvOreMul(rvOre[0,5],rvOreMul(xm,rvOreMul(xm,xm)));
show b;
rvQGrRegularizedRhs[2,5];
rvQGrRegularized[2,5,5];

// Projective space has a nonzero forcing term after Euler cancellation.
let p = rvOreFactorial[37](rvQGrAnnihilator[1,37]);
let g = rvOreEulerReduce(p);
assert g == add(rvOre[0,0],neg(q));
assert p == rvOreMul(rvOreEulerFactor(p),g);
show g;
rvQGrRegularizedRhs[1,37];

// The cyclic order is 14 before transformation and the reduced order is 4.
let w = rvOreFactorial[6](rvQGrAnnihilator[3,6]);
let r = rvOreEulerReduce(w);
assert w == rvOreMul(rvOreEulerFactor(w),r);
rvOreCoeff[0,4](r);
rvQGrRegularized[3,6,4];
