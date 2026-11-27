// All branch points of x^5−x are rational over F9.
rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
rvjacbase D[C,W];
let zw = rvUnramElem[W]{};
let ow = rvUnramElem[W]{1};
let tw = rvUnramElem[W]{2};
let i = rvUnramElem[W]{0,1};
let t0 = rvJacPoint[D](zw,zw);
let t1 = rvJacPoint[D](ow,zw);
let t2 = rvJacPoint[D](tw,zw);
let ti = rvJacPoint[D](i,zw);
rvjactorsionbasis B[D,2](tw){t0,add(add(t0,t1),t2),t1,add(add(t0,t1),ti)};

let h = rvKummerHyperplane[B];
let e = rvKummerNodeDiv[B,0];
let t = rvKummerTropeDiv[B,0];
let u = rvKummerTropeDiv[B,1];
let zero = zero[RVKummerDiv(B)];
assert rvKummerDivFrob[C,1](rvKummerNodeDiv[B,8]) == rvKummerNodeDiv[B,10];
assert rvKummerDivFrob[C,1](rvKummerTropeDiv[B,8]) == rvKummerTropeDiv[B,10];
rvKummerDivPair(h,h);
rvKummerDivPair(e,e);
rvKummerDivPair(t,t);
rvKummerDivPair(h,t);
rvKummerDivPair(e,t);
rvKummerDivPair(t,u);

// H^4, H²E², and delta^4 on Hilb², followed by a mixed Hilb³ intersection.
rvKummerFujiki[2,0,0,0](h,h);
rvKummerFujiki[2,2,0,0](h,e);
rvKummerFujiki[2,2,1,1](zero,zero);
rvKummerFujiki[3,3,0,1](h,t);
