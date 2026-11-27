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

let u = rvK3Lattice[B]{[0,2,0,0,4,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]};
let v = rvK3Lattice[B]{[2,0,0,0,2,-6],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]};
rvk3period P(u,v);
rvk3period Conjugate(u,neg(v));
rvk3period Shear(add(u,v),v);

let aa = add(add(u,u),add(v,add(v,v)));
let bb = add(u,add(add(v,v),add(v,v)));
rvk3period Nonprimitive(add(aa,aa),bb);

rvK3CMForm[P];
rvK3CMForm[Conjugate];
rvK3CMEquivalent[P,Conjugate];
rvK3CMEquivalent[P,Shear];
rvK3CMEquivalent[P,Nonprimitive];
rvK3CMPeriods[P];
rvK3CMClassNumber[P];

let first = rvK3CMBasis[Nonprimitive,0];
let second = rvK3CMBasis[Nonprimitive,1];
rvK3LatticePair(first,first);
rvK3LatticePair(first,second);
rvK3LatticePair(second,second);
