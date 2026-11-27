rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
rvjacbase L[C,W];
let zero = rvUnramElem[W]{};
let one = rvUnramElem[W]{1};
let i = rvUnramElem[W]{0,1};
let r0 = rvJacPoint[L](zero,zero);
let r1 = rvJacPoint[L](one,zero);
let r2 = rvJacPoint[L](neg(one),zero);
let ri = rvJacPoint[L](i,zero);
rvjactorsionbasis B[L,2](neg(one)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),ri)};

let root = rvK3FromKummerDiv(rvKummerNodeDiv[B,0]);
let trope = rvK3FromKummerDiv(rvKummerTropeDiv[B,0]);
let g = rvCliffordVector(root);
let t = rvCliffordVector(trope);
let unit = rvCliffordOne[B];
let omega = rvCliffordVolume[B];
let x = rvK3LatticeBasis[B,0];

assert rvCliffordMul(g,g) == neg(unit);
assert add(rvCliffordMul(g,t),rvCliffordMul(t,g)) == unit;
assert rvCliffordAct(g,x) == rvK3Reflect(root,x);
assert rvCliffordAct(rvCliffordMul(g,t),x) == rvK3Reflect(root,rvK3Reflect(trope,x));
assert rvCliffordMul(omega,omega) == unit;
assert rvCliffordMul(omega,g) == neg(rvCliffordMul(g,omega));
assert rvCliffordAct(omega,x) == neg(x);

let plus = rvCliffordScale[1,2](add(unit,omega));
let minus = rvCliffordScale[1,2](add(unit,neg(omega)));
assert rvCliffordMul(plus,plus) == plus;
assert rvCliffordMul(minus,minus) == minus;
assert rvCliffordMul(plus,minus) == zero[RVK3Clifford(B,0)];
assert rvCliffordMul(plus,rvCliffordMul(g,t)) == rvCliffordMul(rvCliffordMul(g,t),plus);

rvCliffordCoeff[0](rvCliffordMul(g,g));
rvCliffordCoeff[4194303](omega);
rvCliffordCoeff[4194303](plus);
rvCliffordCoeff[4194303](minus);
