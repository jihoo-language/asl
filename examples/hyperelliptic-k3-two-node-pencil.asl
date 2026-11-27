rvunram F[5,1]{0,1};
let z=rvUnramElem[F]{};
let o=rvUnramElem[F]{1};
let two=rvUnramElem[F]{2};
let three=rvUnramElem[F]{3};
rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
let r0=rvJacPoint[C](z,z);
let r1=rvJacPoint[C](o,z);
let r2=rvJacPoint[C](two,z);
let r3=rvJacPoint[C](three,z);
rvjactorsionbasis B[C,2](neg(o)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};

rvk3pencil P[B,0,1];
let f=rvK3PencilDiv[P];
rvKummerDivPair(f,f);
let fourNodes=add(add(rvKummerNodeDiv[B,4],rvKummerNodeDiv[B,7]),add(rvKummerNodeDiv[B,13],rvKummerNodeDiv[B,15]));
assert f==add(rvKummerDivScale[2](rvKummerTropeDiv[B,0]),fourNodes);

let p=rvK3PencilParam[P](o,two);
let x=rvK3PencilSection[P,4](p);
assert rvK3PencilMap[P](x)==p;
assert rvKummerOnTrope[B,4](x)==S.hi;
show p;
show x;
rvK3PencilFiber[P](p);
rvK3PencilCensus[P];
rvKummerK3Points[C,1];
