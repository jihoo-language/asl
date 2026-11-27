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

let a=rvK3LatticeBasis[B,0];
let zline=zero[RVK3Lattice(B)];
let determinant=rvHilbCohSchur[2,2]{1,1}(a);
let c1=rvHilbCohChern[2,1](a);
assert determinant==rvHilbCohScale[1,2](rvHilbCohCup(c1,c1));
assert rvHilbCohSchur[3,0]{2,1}(a)==rvHilbCohScale[8,1](rvHilbCohUnit[B,3]);
rvHilbExteriorEuler[2,2](a,zline);
rvHilbExteriorEuler[2,1](a,a);
