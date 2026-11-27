rvunram Q3[3,4]{0,1};
rveisenstein K[Q3]{-3,0,0,1};
let pi = rvEisensteinPi[K,12];
assert rvEisensteinEval(pi,rvTatePoly[3,4](-3,0,0,1)) == zero[RVEisenstein(K,12)];
rvRamificationPolygon[K];
rvEisensteinVal(rvEisensteinEval(pi,rvTatePoly[3,4](0,0,3)));
rvEisensteinTrace(pi);
rvEisensteinNorm(pi);

rvcyclotomiclocal C[Q3,2];
rvEisensteinPolynomial[C];
rvRamificationPolygon[C];
rvHerbrand[C,18];
let q = rvEisensteinPi[C,18];
assert rvEisensteinSigma[2](rvEisensteinSigma[5](q)) == q;
assert rvEisensteinEval(rvEisensteinSigma[-1](q),rvTatePoly[3,3](3,9,18,21,15,6,1))
    == zero[RVEisenstein(C,18)];

rvunram Q2[2,3]{0,1};
rvcyclotomiclocal D[Q2,3];
rvHerbrand[D,9];
rvcyclotomiclocal E[Q3,3];
rvHerbrand[E,19];
