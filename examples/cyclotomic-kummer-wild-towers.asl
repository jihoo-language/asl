rvunram Q[3,5]{0,1};
rvcyclotomiclocal K[Q,1];

rvkummerlocal First[K,1]([1]);
rvKummerLocalData[First,24];
rvKummerLocalPolynomial[First,24];

rvkummerlocal Second[K,2]([1]);
rvKummerLocalData[Second,24];
rvKummerLocalPolynomial[Second,24];
let pi = rvKummerLocalPi[Second,24];
let y = rvKummerLocalY[Second,24];
let one = rvKummerLocalElem[Second,24](([1]));
let l = rvKummerLocalBase[Second,24](rvEisensteinPi[K,8]);
let w = add(y,neg(one));
assert rvKummerLocalMul(pi,l) == rvKummerLocalMul(w,w);
assert rvKummerLocalEval(pi,rvTatePoly[3,4](3,18,66,30,18,3,1)) == zero[RVKummerLocal(Second,24)];
assert rvKummerLocalSigma[-1](rvKummerLocalSigma[1](y)) == y;
rvKummerLocalTrace(pi);
rvKummerLocalNorm(pi);
rvKummerLocalNorm(y);

rvkummerlocal Boundary[K,3]([1]);
rvKummerLocalBoundary[Boundary];
