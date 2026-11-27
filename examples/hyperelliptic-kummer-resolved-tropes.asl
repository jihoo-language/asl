rvunram W[3,1]{1,0,1};
let z=rvUnramElem[W]{};
let o=rvUnramElem[W]{1};
let i=rvUnramElem[W]{0,1};
let m=neg(o);
rvhyperelliptic C[W]{z,m,z,z,z,o};
let t0=rvJacPoint[C](z,z);
let t1=rvJacPoint[C](o,z);
let tm=rvJacPoint[C](m,z);
let ti=rvJacPoint[C](i,z);
let e2=add(add(t0,t1),tm);
let f2=add(add(t0,t1),ti);
rvjactorsionbasis B[C,2](m){t0,e2,t1,f2};

// The strict tropes with masks 0 and 13 meet different directions above ti.
let origin=rvKummerTropeParam[B,0](z,o);
let translated=rvKummerK3Translate[B,13](origin);
assert translated==rvKummerTropeParam[B,13](z,o);
assert rvKummerBlowdown(translated)==rvJacKummerPoint(ti);
show rvKummerOnTrope[B,13](translated);
show rvKummerOnTrope[B,0](translated);
