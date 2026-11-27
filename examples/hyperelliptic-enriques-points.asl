rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
rvjacbase L[C,W];
let zw = rvUnramElem[W]{};
let ow = rvUnramElem[W]{1};
let i = rvUnramElem[W]{0,1};
let e1 = rvJacPoint[L](zw,zw);
let f1 = rvJacPoint[L](ow,zw);
let e2 = add(add(e1,f1),rvJacPoint[L](neg(ow),zw));
let f2 = add(add(e1,f1),rvJacPoint[L](i,zw));
rvjactorsionbasis B[L,2](neg(ow)){e1,e2,f1,f2};
rvenriques Y[C,L,B,2];

card (\y:RVEnriques(Y)=>S.hi);
card (\y:RVEnriques(Y)=>rvEnriquesSplit(y));
assert (\y:RVEnriques(Y)=>rvEnriquesOrbit[Y](rvEnriquesRepresentative(y))) == (\y:RVEnriques(Y)=>y);
assert (\y:RVEnriques(Y)=>rvEnriquesOrbit[Y](rvKummerK3Translate[B,2](rvKummerSwitch(rvEnriquesRepresentative(y))))) == (\y:RVEnriques(Y)=>y);

// This exceptional point retains its tangent direction in the round trip.
let x = rvKummerExceptional(zero[RVKummer(C)],z,z,o,z);
let y = rvEnriquesOrbit[Y](x);
assert rvEnriquesSplit(y) == S.hi;
assert rvEnriquesLift(y) == x;
assert rvEnriquesRepresentative(y) == rvKummerK3Extend[L](x);
show y;

let xTwisted = rvKummerSmooth(rvKummer[L](ow,ow,rvUnramElem[W]{2,1},rvUnramElem[W]{1,2}));
let yTwisted = rvEnriquesOrbit[Y](xTwisted);
assert rvEnriquesSplit(yTwisted) == S.lo;
assert rvEnriquesRepresentative(yTwisted) == xTwisted;
assert rvEnriquesOrbit[Y](rvKummerK3Frob[C,1](xTwisted)) == yTwisted;
show yTwisted;
