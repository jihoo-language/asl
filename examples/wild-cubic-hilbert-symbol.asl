rvunram B[3,7]{0,1};
rvcyclotomiclocal K[B,1];
let one = rvEisensteinElem[K,14]([1]);
let l = rvEisensteinPi[K,14];
let z = add(one,l);
let l2 = rvEisensteinMul(l,l);
let u2 = add(one,l2);
let u3 = add(one,rvEisensteinMul(l2,l));

rvCubicHilbertCertificate[K];
rvCubicHilbertInvariant(rvCubicClass(z),rvCubicClass(u2));
assert rvCubicHilbert(rvCubicClass(u3),rvCubicClass(l)) == rvCyclo[3,1];
rvCubicNormSubgroup(rvCubicClass(u3));

let a = rvCubicLocal[K,8,-2](z);
let b = rvCubicOneMinus[4](a);
show rvCubicClass(a);
show rvCubicClass(b);
assert rvCubicAlgebraSplits(rvCubicClass(a),rvCubicClass(b)) == S.hi;
