rvunram F[23,1]{0,1};
let z=rvUnramElem[F]{};
let o=rvUnramElem[F]{1};
let two=rvUnramElem[F]{2};
let four=rvUnramElem[F]{4};
let nine=rvUnramElem[F]{9};
rvhyperelliptic C[F]{z,four,neg(four),o,neg(two),o};
let r0=rvJacPoint[C](z,z);
let r1=rvJacPoint[C](o,z);
let r2=rvJacPoint[C](two,z);
let r3=rvJacPoint[C](nine,z);
rvjactorsionbasis B[C,2](neg(o)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
rvk3pencil P[B,0,14];

let p=rvK3PencilSectionPoint[P]{[0],[18],[2],[5]}{[1]}{[0],[0],[18],[11],[10],[5],[6]}{[1]};
let q=rvK3PencilSectionPoint[P]{[0],[18],[4],[5]}{[1]}{[0],[0],[14],[7],[18],[5],[6]}{[1]};
rvk3pencilmw2 M[P]{p,q};
rvK3PencilMW2[M];

let delta=rvK3PencilDescent2(p);
rvk3pencilcover RationalCover[P](delta);
let witness=rvK3PencilCoverPoint[RationalCover]{[2],[1]}{[0],[6],[18]}{[15],[19]}{[6]};
assert rvK3PencilCoverMap(witness)==p;
assert rvK3PencilCoverLift[RationalCover](p)==witness;
show witness;
rvK3PencilCoverStatus[RationalCover,M];

let a=rvK3PencilKummer2[P]{[0],[-10],[1]}{[0],[-10],[1]}{[1]};
let b=rvK3PencilKummer2[P]{[150],[-85],[5]}{[0],[-2],[1]}{[0],[-75],[5]};
rvk3pencilcover A[P](a);
rvk3pencilcover Bcover[P](b);
rvk3pencilcover Sum[P](add(a,b));

rvk3cassels T[A]{
  0{[5]}{[0],[13],[4]}{[1]},
  1{[11],[0],[3]}{[0],[4],[8],[3]}{[1]},
  2{[0],[19],[20],[1],[17]}{[0],[16],[18],[1],[17]}{[0],[3],[2]}
};
rvK3CasselsTangents[T];
rvK3CasselsPair[T](b);
rvK3CasselsPair[T](a);
rvK3CasselsPair[T](delta);

rvk3cassels Second[Bcover]{
  0{[4],[12]}{[5],[3]}{[1]},
  1{}{[6],[3]}{[1]},
  2{[0],[8]}{[7],[21]}{[1]}
};
rvk3selmer4 R[M]{T,Second};
rvK3Selmer4[R];
show rvK3Selmer4Pair[R](a,b);
assert rvK3Selmer4Contains[R](delta)==S.hi;
assert rvK3Selmer4Contains[R](a)==S.lo;
show rvK3Selmer4ImageBasis[R,0];
