rvunram W[3,1]{1,0,1};
let zeroW = rvUnramElem[W]{};
let oneW = rvUnramElem[W]{1};
let i = rvUnramElem[W]{0,1};
rvhyperelliptic C[W]{zeroW,neg(oneW),zeroW,zeroW,zeroW,oneW};

// The five finite branch points are 0, ±1, ±i.
let t0 = rvJacPoint[C](zeroW,zeroW);
let t1 = rvJacPoint[C](oneW,zeroW);
let tm = rvJacPoint[C](neg(oneW),zeroW);
let ti = rvJacPoint[C](i,zeroW);
let e1 = t0;
let e2 = add(add(t0,t1),tm);
let f1 = t1;
let f2 = add(add(t0,t1),ti);
rvjactorsionbasis B[C,2](neg(oneW)){e1,e2,f1,f2};

// Mask 13 is e1 + f1 + f2, the class of the branch point at i.
let image = rvKummerTranslate[B,13](zero[RVKummer(C)]);
assert image == rvJacKummerPoint(ti);
assert rvKummerTranslate[B,13](image) == zero[RVKummer(C)];
assert rvKummerThetaEntry[B,13,2,3]
    == rvUnramMul(i,rvKummerThetaEntry[B,13,1,3]);
assert rvKummerThetaFactor[B,1,4] == neg(rvKummerThetaFactor[B,4,1]);
assert rvKummerTropeCoeff[B,0,0] == oneW;
show image;
show rvKummerThetaEntry[B,13,1,3];
show rvKummerThetaEntry[B,13,2,3];
