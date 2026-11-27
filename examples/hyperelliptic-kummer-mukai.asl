rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let z = rvUnramElem[F]{};
let one = rvUnramElem[F]{1};
rvhyperelliptic C[F]{z,neg(one),z,z,z,one};
rvjacbase L[C,W];
let zw = rvUnramElem[W]{};
let ow = rvUnramElem[W]{1};
let tw = rvUnramElem[W]{2};
let i = rvUnramElem[W]{0,1};
let t0 = rvJacPoint[L](zw,zw);
let t1 = rvJacPoint[L](ow,zw);
let t2 = rvJacPoint[L](tw,zw);
let ti = rvJacPoint[L](i,zw);
let e2 = add(add(t0,t1),t2);
let f2 = add(add(t0,t1),ti);
rvjactorsionbasis B[L,2](tw){t0,e2,t1,f2};

// Mukai calculations below.
let h = rvKummerHyperplane[B];
let t = rvKummerTropeDiv[B,0];
let e = rvKummerNodeDiv[B,0];
let embeddedT = rvK3FromKummerDiv(t);
assert rvK3ToKummerDiv(embeddedT) == t;
let structure = rvMukaiIdeal[B,0];
let line = rvMukaiLine(h);
assert line == rvMukaiTensor(structure,h);
assert rvMukaiTensor(rvMukaiTensor(structure,h),t) == rvMukaiLine(add(h,t));
rvMukaiPair(structure,line);
let ideal = rvMukaiIdeal[B,3];
rvMukaiPair(ideal,ideal);
rvMukaiPair(ideal,line);
rvMukaiPair(rvMukaiTensor(ideal,t),rvMukaiTensor(line,t));

// The incident exceptional curve and trope give an A2 pair of roots.
let a = rvMukai[0,0](rvK3FromKummerDiv(e));
let b = rvMukai[0,0](embeddedT);
let x = rvMukaiLine(t);
assert rvMukaiTwist(a,rvMukaiTwist(b,rvMukaiTwist(a,x)))
    == rvMukaiTwist(b,rvMukaiTwist(a,rvMukaiTwist(b,x)));
assert rvMukaiTwist(a,rvMukaiTwist(a,x)) == x;

let halfExceptional = rvMukai[1,1](zero[RVK3Lattice(B)]);
rvMukaiHilbFujiki[2,2](halfExceptional,halfExceptional);
rvMukaiHilbFujiki[3,3](rvMukai[0,0](rvMukaiC(line)),rvMukai[1,2](embeddedT));
