rvPeriodDim[12,plus];
rvPeriodDim[12,minus];
let p = rvPeriodPoly[12,plus](0,0,1,0,-3,0,3,0,-1,0,0);
let q = rvPeriodPoly[12,minus](0,4,0,-25,0,42,0,-25,0,4,0);
let t = rvHeckeT[12,2];
assert rvPeriodHecke(t,p) == rvPeriodScale[-24,1](p);
assert rvPeriodHecke(t,q) == rvPeriodScale[-24,1](q);
rvPeriodCoeff[2](rvPeriodHecke(t,p));

let a = rvHeckeT[24,2];
rvPeriodHeckeMatrix[24,plus](a);
rvPeriodHeckeMatrix[24,minus](a);
rvHeckeChar(a);
