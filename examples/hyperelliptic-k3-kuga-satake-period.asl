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

let u = rvK3Lattice[B]{[0,2,0,0,2,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]};
let v = rvK3Lattice[B]{[3,0,0,0,0,-4],[0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1]};
rvk3period P(u,v);
rvk3period Shear(u,add(u,v));
rvk3period Opposite(v,u);
rvk3period Scaled(add(u,u),add(v,add(v,v)));

let omega = rvK3PeriodBlade[P];
let unit = rvCliffordOne[B];
assert rvCliffordMul(omega,omega) == rvCliffordScale[-40,1](unit);
assert omega == rvCliffordScale[2,1](rvCliffordMul(rvCliffordVector(u),rvCliffordVector(v)));
assert rvK3PeriodBlade[Shear] == omega;
assert rvK3PeriodBlade[Opposite] == neg(omega);
assert rvK3PeriodBlade[Scaled] == rvCliffordScale[6,1](omega);

let plus = rvK3SpinorBasis[B,0,0];
let minus = rvK3SpinorBasis[B,1,37];
assert rvK3PeriodCM[P](plus) == rvK3SpinAct(omega,plus);
assert rvK3PeriodCM[P](rvK3PeriodCM[P](plus)) == rvK3SpinorScale[-40,1](plus);
assert rvK3PeriodCM[P](rvK3PeriodCM[P](minus)) == rvK3SpinorScale[-40,1](minus);
assert rvK3PeriodCM[Opposite](minus) == neg(rvK3PeriodCM[P](minus));

rvK3PeriodReport[P];
rvK3PeriodReport[Scaled];
rvK3PeriodIs11[P](rvK3FromKummerDiv(rvKummerNodeDiv[B,0]));
rvK3PeriodIs11[P](rvK3FromKummerDiv(rvKummerNodeDiv[B,15]));
rvK3PeriodIs11[P](u);
