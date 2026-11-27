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

let scaledUnit=rvNakajimaMonomial[B]{(1,0,2)};
let delta=rvHilbSquareBasis[B,2,22];
let boundary=rvNakajimaBoundary(scaledUnit);
assert rvNakajimaHilbSquare(boundary)==rvHilbSquareScale[2,1](delta);
assert rvNakajimaHilbSquare(rvNakajimaBoundary(boundary))==rvHilbSquareScale[2,1](rvHilbSquareCup(delta,delta));

let u=rvNakajimaMonomial[B]{(1,0,3)};
let d1=rvNakajimaBoundary(u);
let d2=rvNakajimaBoundary(d1);
let d3=rvNakajimaBoundary(d2);
let d4=rvNakajimaBoundary(d3);
let d5=rvNakajimaBoundary(d4);
let d6=rvNakajimaBoundary(d5);
let point=rvNakajimaMonomial[B]{(1,23,3)};
assert d6==rvNakajimaScale[-5760](point);
let zeroDiv=zero[RVKummerDiv(B)];
rvKummerFujiki[3,6,0,1](zeroDiv,zeroDiv);
