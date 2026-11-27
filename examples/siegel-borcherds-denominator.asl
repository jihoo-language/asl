// Coordinates are doubled. The Weyl vector is (1/2, 1/2, 1/2).
let rho=rvSiegelWeight[1,1,1];
let e=rvSiegelWeight[2,0,2];
assert add(rho,neg(rho)) == zero[RVSiegelWeight];
assert rvSiegelReflect[1](rvSiegelReflect[1](rho)) == rho;
rvSiegelPair(rho,e);
rvSiegelRoot(e);
rvSiegelRoot(rvSiegelWeight[2,2,2]);
rvSiegelDelta5At(rho);
rvSiegelDelta5At(rvSiegelReflect[1](rho));
rvIgusaPhi10[1,1,1];
