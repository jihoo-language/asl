rvunram B[3,7]{0,1};
rvcyclotomiclocal K[B,1];
rvcubicdiv D[K,1,4];
rvcubicdiv E[K,2,7];

rvCubicDivData[D];
let p = rvCubicDivPi[D];
let theta = rvCubicDivTheta[D];
let one = rvCubicDivElem[D]([1]);
let x = add(one,p);
let y = add(one,rvCubicDivMul(theta,p));
rvCubicDivCommutatorCheck(x,y);
rvCubicDivVal(add(rvCubicDivMul(rvCubicDivMul(x,x),x),neg(one)));
assert rvCubicDivMul(x,rvCubicDivInverse(x)) == one;
assert rvCubicDivMul(p,theta) == rvCubicDivMul(rvCubicDivConjugatePi(theta),p);

rvCubicDivData[E];
let q = rvCubicDivPi[E];
let w = rvCubicDivTheta[E];
let e = rvCubicDivElem[E]([1]);
rvCubicDivCommutatorCheck(add(e,q),add(e,rvCubicDivMul(w,q)));
rvCubicDivVal(q);
show rvCubicDivNorm[3](q);

let lambda = rvEisensteinPi[K,14];
let lambda2 = rvEisensteinMul(lambda,lambda);
let u3 = add(rvEisensteinElem[K,14]([1]),rvEisensteinMul(lambda,lambda2));
rvcubicalgebra A[K,8](u3,lambda2);
let v = rvCubicDivFromSymbol[E](rvCubicAlgBasis[A,0,1]);
assert v == rvCubicDivMul(q,q);
