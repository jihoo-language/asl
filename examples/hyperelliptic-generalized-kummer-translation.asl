// A fourfold's middle character and fixed locus under an order-three translation.
rvGenKummerTranslateChar[3,3,4];
rvGenKummerTranslateFixed[3,3];
rvGenKummerTranslateFixed[4,2];

// Full 2-torsion over F9. Frobenius over F3 sends f2 to f2+e2
// and fixes the other basis elements.
rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
rvjacbase D[C,W];
let zero = rvUnramElem[W]{};
let one = rvUnramElem[W]{1};
let two = rvUnramElem[W]{2};
let i = rvUnramElem[W]{0,1};
let e1 = rvJacPoint[D](zero,zero);
let f1 = rvJacPoint[D](one,zero);
let e2 = add(add(e1,f1),rvJacPoint[D](two,zero));
let f2 = add(add(e1,f1),rvJacPoint[D](i,zero));
rvjactorsionbasis B[D,2](two){e1,e2,f1,f2};
assert rvJacFrob[C,1](f2) == add(f2,e2);

// The affine torsion equation is solvable for e2 and insoluble for e1 or f2.
rvGenKummerTranslatePoints[B,C,1](e2);
rvGenKummerTranslatePoints[B,C,1](e1);
rvGenKummerTranslateFrob[B,C,1,2](f2);
rvGenKummerTranslateFrob[B,C,1,2](zero[RVJac(D)]);
