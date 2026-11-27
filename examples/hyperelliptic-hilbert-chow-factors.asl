rvunram F[3,1]{0,1};
rvunram V[3,1]{1,0,1};
let z=rvUnramElem[F]{};
let o=rvUnramElem[F]{1};
rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
rvjacbase L[C,V];
let zv=rvUnramElem[V]{};
let ov=rvUnramElem[V]{1};
let a=rvUnramElem[V]{0,1};
let e1=rvJacPoint[L](zv,zv);
let f1=rvJacPoint[L](ov,zv);
let e2=add(add(e1,f1),rvJacPoint[L](neg(ov),zv));
let f2=add(add(e1,f1),rvJacPoint[L](a,zv));
rvjactorsionbasis B[L,2](neg(ov)){e1,e2,f1,f2};
rvk3autom Moving[B]{translate(8)};
rvk3autom Idle[B]{};

rvhilbstratum Diagonal{(2,1)};
rvhilbstratum Distinct{(1,2)};
// The square of Moving composed after Frobenius has a different shift.
rvK3StratumFrobTrace[Moving,C,Diagonal,1,4];
rvK3StratumFrobTrace[Moving,C,Diagonal,2,4];
rvK3StratumFrobTrace[Moving,C,Distinct,1,4];
rvK3FrobHilbTrace[Moving,C,2,1,4];

// In degree two the diagonal factor is 1−3T.
rvK3StratumEulerFactor[Moving,C,Diagonal,1,2,0];
rvK3StratumEulerFactor[Moving,C,Diagonal,1,2,1];
rvK3StratumEulerFactor[Moving,C,Diagonal,1,2,2];

// The degree-four diagonal factor starts with 1+6T−81T².
rvK3StratumEulerFactor[Moving,C,Diagonal,1,4,1];
rvK3StratumEulerFactor[Moving,C,Diagonal,1,4,2];
rvK3StratumEulerFactor[Idle,C,Diagonal,1,4,2];
