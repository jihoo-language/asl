// A non-scalar cocycle with a nonzero internal Bockstein modulo rho^2.
let raw = add(rmcword[0,3]{([1],[])},add(rmcword[1,2]{([1],[1])},rmcword[1,2]{([],[0,1])}));
let z = extClass(rhoReduce[2](raw));
assert rhoTruncLift(z) == raw;
assert rhoReduce[2](rhoTruncLift(extRep(z))) == extRep(z);
let d = rhoConnectRep(z);
assert cobarD(raw) == rhoShift[2](d);
assert cobarD(d) == zero[RMCobar(3,4,0)];
let b = rhoBeta(z);
assert b == rhoReduce[2](rhoConnect(z));
assert rhoBeta(b) == zero[RMRhoExt(2,3,6,2)];
assert rhoConnect(rhoReduce[1](z)) == rhoShift[1](rhoConnect(z));
show b;

let tc = mextgen[0,0,-1,0];
let t = rhoReduce[1](tc);
assert realComplex(t) == tc;
assert realComplex(extRep(t)) == extRep(tc);
assert realComplex(yoneda(t,t)) == yoneda(tc,tc);
assert realComplex(rhoBeta(t)) == mextgen[1,1,0,0];
let u = extClass(rhoReduce[2](rmcword[0,2]{}));
assert rhoBeta(yoneda(u,z)) == add(yoneda(rhoBeta(u),z),yoneda(u,rhoBeta(z)));
assert rhoShift[2](rhoConnect(z)) == zero[RMExt(2,2,-2)];

// The scalar differential is invisible modulo rho^32 but not rho^33.
let tau32 = rmcword[0,32]{};
let a = extClass(rhoReduce[32](tau32));
rmrhoextdim[32,0,0,-32];
rmrhoextdim[33,0,0,-32];
let target = rmcword[0,0]{([],[32])};
assert rhoConnectRep(a) == target;
assert rhoConnect(a) == extClass(target);
assert rhoShift[32](rhoConnectRep(a)) == cobarD(tau32);
rhoorder(rhoConnect(a));
assert rhoBeta(a) == rhoReduce[32](rhoConnect(a));
assert rhoBeta(rhoBeta(a)) == zero[RMRhoExt(32,2,64,32)];
assert rhoConnect(rhoReduce[31](a)) == rhoShift[1](rhoConnect(a));
assert rhoReduce[33](cobarD(tau32)) == cobarD(rhoReduce[33](tau32));
show rhoConnect(a);
show rhoBeta(a);
