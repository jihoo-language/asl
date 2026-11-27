rvunram B[3,7]{0,1};
rvcyclotomiclocal K[B,1];
rvcubicdiv D[K,2,5];

let c = rvCubicBaseChamber[D];
let u = rvCubicChamberPanel[0,0](c);
let v = rvCubicChamberPanel[1,0](u);
rvCubicGallery(c,v);
assert rvCubicOmega[1](c) == c;

let t0 = rvCubicIwT[D,[0]];
let t1 = rvCubicIwT[D,[1]];
let t01 = rvCubicIwT[D,[0,1]];
let omega = rvCubicIwOmega[D,1];
let delta = rvCubicIwDelta(c);
assert rvCubicHeckeMul(omega,t0) == rvCubicHeckeMul(t1,omega);
assert rvCubicIwAct(t01,delta) == rvCubicIwAct(t1,rvCubicIwAct(t0,delta));
rvCubicIwKernel(t01,c,v);
rvCubicIwVectorCoeff(rvCubicIwAct(t01,delta),v);

let z = rvCubicIwZ[D];
assert rvCubicHeckeMul(z,t0) == rvCubicHeckeMul(t0,z);
let e = rvCubicIwE[D];
assert rvCubicHeckeMul(e,e) == e;
let a2 = rvCubicSph[D,2];
let a22 = rvCubicHeckeMul(a2,a2);
assert rvCubicIwCompress(a22) == rvCubicHeckeMul(rvCubicIwCompress(a2),rvCubicIwCompress(a2));
assert rvCubicSatake(a22) == rvCubicHeckeMul(rvCubicSatake(a2),rvCubicSatake(a2));
show rvCubicSatake(a2);
rvCubicSatakeCounts(rvCubicTreeAt[D,2,0]);
