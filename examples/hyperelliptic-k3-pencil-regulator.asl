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


let p=rvK3PencilSectionPoint[P]{[0],[1],[0],[6]}{[1]}{[0],[0],[7],[6],[3]}{[1]};
let q=rvK3PencilSectionPoint[P]{[12],[0],[12],[6]}{[1]}{[8],[8],[10],[3],[6]}{[1]};
rvK3PencilSectionHeight[P](p);
rvK3PencilSectionHeight[P](q);
rvK3PencilSectionPair[P](p,q);
rvK3PencilRegulator[P]{p,q};
