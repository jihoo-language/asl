rvunram B[3,7]{0,1};
rvcyclotomiclocal K[B,1];
rvcubicdiv D[K,2,5];

let root = rvCubicTreeRoot[D];
let left = rvCubicTreeAt[D,1,0];
let right = rvCubicTreeAt[D,1,27];
let f = add(rvCubicBoundaryCylinder(left),neg(rvCubicBoundaryCylinder(right)));
let refined = rvCubicBoundaryRefine[3](f);
assert refined == f;
assert rvCubicBoundaryCoarsen[1](refined) == f;
assert rvCubicBoundaryDifference[1](f) == f;
show refined;
rvCubicBoundaryMass(right);
rvCubicBoundaryValue(f,left);
rvCubicPoisson[3,1](f,left);
rvCubicPoissonEigen[3,1](f,left);

rvCubicKnappSpectrum[3,1,1];
rvCubicKnappSpectrum[3,1,2];
rvCubicKnappSpectrum[3,1,3];
rvCubicKnappKernel[3,1](left,left);
let a = rvCubicKnappStein[3,1](f);
assert a == rvCubicBoundaryScale[10,3](f);
assert rvCubicKnappStein[9,1](a) == f;
let moved = rvCubicPrincipalOmega[3,1,1](f);
assert rvCubicKnappStein[3,1](moved) == rvCubicPrincipalOmega[9,1,1](a);
rvCubicPoissonCovariantOmega[3,1,1](f,root);
