// Full 2-torsion of the Kummer surface over F9, with a moving translation over F3.
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

rvhilbstratum Square{(1,2)};
rvspechtstratum Symmetric[Square]{};
rvspechtstratum Alternating[Square]{(1,[1,1])};
rvK3SpechtFrobTrace[Moving,C,Symmetric,1,4];
rvK3SpechtFrobTrace[Moving,C,Alternating,1,4];

rvhilbstratum Triple{(2,3)};
rvspechtstratum Standard[Triple]{(2,[2,1])};
rvK3SpechtFrobTrace[Moving,C,Standard,1,4];
rvK3SpechtEulerFactor[Moving,C,Standard,1,2,2];

// Intrinsic degree 2 of this row diagram appears in Hilbert degree 6.
rvhilbstratum Weighted{(2,2)};
rvspechtstratum Row[Weighted]{(2,[2])};
rvK3SpechtFrobTrace[Moving,C,Row,1,2];
rvK3StratumFrobTrace[Moving,C,Weighted,1,6];
