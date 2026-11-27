rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let z=rvUnramElem[F]{};
let o=rvUnramElem[F]{1};
rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
rvjacbase L[C,W];
let zw=rvUnramElem[W]{};
let ow=rvUnramElem[W]{1};
let i=rvUnramElem[W]{0,1};
let e1=rvJacPoint[L](zw,zw);
let f1=rvJacPoint[L](ow,zw);
let e2=add(add(e1,f1),rvJacPoint[L](neg(ow),zw));
let f2=add(add(e1,f1),rvJacPoint[L](i,zw));
rvjactorsionbasis B[L,2](neg(ow)){e1,e2,f1,f2};

rvk3autom Idle[B]{};
rvKummerK3Points[C,1];
rvK3FrobFixed[Idle,C,1];

rvk3autom Enriques[B]{switch(2)};
rvK3FrobFixed[Enriques,C,1];
rvK3FrobFixed[Enriques,C,2];
rvK3FrobZeta[Enriques,C];

rvk3autom Noncommuting[B]{translate(8)};
rvK3FrobFixed[Noncommuting,C,2];
rvK3FrobFixed[Noncommuting,C,4];
rvK3FrobZeta[Noncommuting,C];

rvk3autom Hyperbolic[B]{project(0),project(1),project(2)};
rvK3FrobFixed[Hyperbolic,C,1];
rvK3FrobFixed[Hyperbolic,C,2];
rvK3FrobZeta[Hyperbolic,C];
