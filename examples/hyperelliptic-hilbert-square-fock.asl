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

let unit=rvHilbSquareBasis[B,0,0];
let delta=rvHilbSquareBasis[B,2,22];
let delta2=rvHilbSquareCup(delta,delta);
let pointFock=rvNakajimaMonomial[B]{(1,0,1),(1,23,1)};
let supportPoint=rvNakajimaHilbSquare(pointFock);
assert supportPoint==add(rvHilbSquareScale[1,24](rvHilbSquareC2[B]),rvHilbSquareScale[1,8](delta2));
assert rvHilbSquareNakajima[1](supportPoint)==pointFock;
rvHilbSquareIntegral(rvHilbSquareCup(supportPoint,supportPoint));
rvHilbSquareIntegral(rvHilbSquareCup(supportPoint,delta2));

let exceptional=rvNakajimaMonomial[B]{(2,0,1)};
let punctual=rvNakajimaMonomial[B]{(2,23,1)};
assert rvNakajimaHilbSquare(exceptional)==rvHilbSquareScale[2,1](delta);
rvNakajimaPair(exceptional,punctual);
rvHilbSquareIntegral(rvHilbSquareCup(rvNakajimaHilbSquare(exceptional),rvNakajimaHilbSquare(punctual)));

rvHilbSquareNakajimaDen(unit);
let scaledUnit=rvHilbSquareNakajima[2](unit);
assert scaledUnit==rvNakajimaCreate[1,0](rvNakajimaCreate[1,0](rvNakajimaVac[B]));
assert rvNakajimaHilbSquare(scaledUnit)==rvHilbSquareScale[2,1](unit);
