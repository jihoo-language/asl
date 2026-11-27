rvunram W[3,1]{1,0,1};
let z=rvUnramElem[W]{};
let o=rvUnramElem[W]{1};
let i=rvUnramElem[W]{0,1};
rvhyperelliptic C[W]{z,neg(o),z,z,z,o};
let e1=rvJacPoint[C](z,z);
let f1=rvJacPoint[C](o,z);
let e2=add(add(e1,f1),rvJacPoint[C](neg(o),z));
let f2=add(add(e1,f1),rvJacPoint[C](i,z));
rvjactorsionbasis B[C,2](neg(o)){e1,e2,f1,f2};

let k3Zero=zero[RVK3Lattice(B)];
let v=rvMukai[0,1](k3Zero);
rvmukaimoduli Partner(v);
let w=rvMukaiModuliDual[Partner];
assert w==rvMukai[-1,0](k3Zero);
rvMukaiModuliLattice[Partner];
rvMukaiPair(v,w);
rvMukaiPair(w,w);

let h=rvK3FromKummerDiv(rvKummerHyperplane[B]);
let x=rvMukai[0,0](h);
let a=rvMukaiQuotient[Partner](x);
let b=rvMukaiQuotient[Partner](add(x,v));
assert a==b;
assert rvMukaiQuotientLift(a)==x;
assert rvMukaiQuotient[Partner](v)==zero[RVMukaiQuotient(Partner)];
rvMukaiQuotientPair(a,b);
rvMukaiModuliDiv[Partner](a);
rvMukaiModuliFujiki[Partner,1](a,add(x,v));
rvMukaiModuliRR[Partner](b);

// This isotropic vector has a nonzero middle component and rank two.
let c=rvK3Lattice[B]{[0,2,0,0,4,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]};
let u=rvMukai[2,2](c);
rvmukaimoduli Other(u);
let dual=rvMukaiModuliDual[Other];
rvMukaiPair(dual,dual);
rvMukaiPair(u,dual);
let basis=rvMukaiModuliBasis[Other,0];
assert rvMukaiQuotientLift(rvMukaiQuotient[Other](add(basis,u)))==basis;
