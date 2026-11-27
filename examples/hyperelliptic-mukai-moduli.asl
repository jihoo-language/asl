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

// Hilbert cube and the half exceptional class.
let zeroK3=zero[RVK3Lattice(B)];
let v=rvMukaiIdeal[B,3];
rvmukaimoduli Hilb(v);
let delta=rvMukai[1,2](zeroK3);
rvMukaiModuliLattice[Hilb];
rvMukaiModuliDiv[Hilb](delta);
rvMukaiModuliFujiki[Hilb,0](delta,delta);
rvMukaiModuliRR[Hilb](zero[RVMukai(B)]);
rvMukaiModuliRR[Hilb](add(delta,delta));

// Both outside coordinates are even. The middle class forces primitivity.
let h=rvK3FromKummerDiv(rvKummerHyperplane[B]);
let threeh=add(add(h,h),h);
let v2=rvMukai[2,-4](threeh);
rvmukaimoduli Mod(v2);
let g=rvMukaiModuliGlue[Mod];
rvMukaiModuliLattice[Mod];
rvMukaiPair(v2,rvMukaiModuliBasis[Mod,0]);
rvMukaiPair(g,v2);
rvMukaiModuliDiv[Mod](g);
