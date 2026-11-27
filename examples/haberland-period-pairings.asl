let p = rvPeriodPoly[12,plus](0,0,1,0,-3,0,3,0,-1,0,0);
let q = rvPeriodPoly[12,minus](0,4,0,-25,0,42,0,-25,0,4,0);
rvPeriodPair(p,q);
rvPeriodPair(q,p);

let x = rvPeriodIntegralBasis[24,plus,0];
let y = rvPeriodIntegralBasis[24,minus,0];
let t = rvHeckeT[24,2];
rvPeriodPair(rvPeriodHecke(t,x),y);
rvPeriodPair(x,rvPeriodHecke(t,y));
rvPeriodPairing[24];
