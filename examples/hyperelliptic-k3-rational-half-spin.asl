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

let plus = rvK3SpinorBasis[B,0,0];
let minus = rvK3SpinorBasis[B,1,1023];
let another = rvK3SpinorBasis[B,0,1023];
let unit = rvCliffordOne[B];
let vol = rvCliffordVolume[B];

assert rvK3SpinAct(vol,plus) == plus;
assert rvK3SpinAct(vol,minus) == neg(minus);
let pplus = rvCliffordScale[1,2](add(unit,vol));
let pminus = rvCliffordScale[1,2](add(unit,neg(vol)));
assert rvK3SpinAct(pplus,plus) == plus;
assert rvK3SpinAct(pplus,minus) == zero[RVK3Spinor(B,1)];
assert rvK3SpinAct(pminus,minus) == minus;
assert rvK3SpinAct(pminus,plus) == zero[RVK3Spinor(B,0)];

let root = rvK3FromKummerDiv(rvKummerNodeDiv[B,15]);
let x = rvK3Lattice[B]{[1,0,0,0,0,0],[0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1]};
let g = rvCliffordVector(root);
let h = rvCliffordVector(x);
let twisted = rvK3SpinAct(g,plus);
assert twisted == neg(rvK3SpinorBasis[B,1,0]);
assert rvK3SpinAct(g,twisted) == neg(plus);
assert rvK3SpinAct(rvCliffordMul(g,h),minus) == rvK3SpinAct(g,rvK3SpinAct(h,minus));
assert rvK3SpinAct(add(rvCliffordMul(g,h),rvCliffordMul(h,g)),plus) == neg(plus);

let combination = add(rvK3SpinorScale[2,3](plus),rvK3SpinorScale[-5,7](another));
assert rvK3SpinAct(g,combination) == add(rvK3SpinorScale[2,3](twisted),rvK3SpinorScale[-5,7](rvK3SpinAct(g,another)));
assert rvK3SpinAct(add(g,h),combination) == add(rvK3SpinAct(g,combination),rvK3SpinAct(h,combination));
assert rvK3SpinAct(rvCliffordScale[2,3](g),combination) == rvK3SpinorScale[2,3](rvK3SpinAct(g,combination));
rvK3SpinCoeff[0](rvK3SpinAct(pplus,combination));
rvK3SpinCoeff[1023](rvK3SpinAct(pplus,combination));
rvK3SpinCoeff[1023](rvK3SpinAct(pminus,minus));
