rvunram F[13,1]{0,1};
let z=rvUnramElem[F]{};
let o=rvUnramElem[F]{1};
let two=rvUnramElem[F]{2};
let three=rvUnramElem[F]{3};
let nine=rvUnramElem[F]{9};
let eleven=rvUnramElem[F]{11};
rvhyperelliptic C[F]{z,eleven,two,nine,three,o};
let r0=rvJacPoint[C](z,z);
let r1=rvJacPoint[C](o,z);
let r2=rvJacPoint[C](two,z);
let r3=rvJacPoint[C](three,z);
rvjactorsionbasis B[C,2](neg(o)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
rvk3pencil P[B,0,1];

let p=rvK3PencilParam[P](o,o);
rvK3FiberQuartic[P](p);
rvK3FiberJacobian[P](p);
show rvK3FiberJ[P](p);
rvK3PencilFiber[P](p);
rvK3FiberPoints[P,1](p);
rvK3FiberPoints[P,2](p);
rvk3fiberjac E[P](p);
rvECPointCount[E];

let central=rvK3PencilParam[P](z,o);
rvK3FiberQuartic[P](central);
rvK3FiberJacobian[P](central);
