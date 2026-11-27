rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
rvjacbase L[C,W];
let zeroW = rvUnramElem[W]{};
let oneW = rvUnramElem[W]{1};
let a = rvUnramElem[W]{0,1};
let r0 = rvJacPoint[L](zeroW,zeroW);
let r1 = rvJacPoint[L](oneW,zeroW);
let r2 = rvJacPoint[L](neg(oneW),zeroW);
let ra = rvJacPoint[L](a,zeroW);
rvjactorsionbasis B[L,2](neg(oneW)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),ra)};

let vacuum = rvNakajimaVac[B];
assert rvNakajimaAnn[2,1](vacuum) == zero[RVNakajima(B,-2,-4)];
let x = rvNakajimaCreate[2,2](vacuum);
let left = rvNakajimaAnn[2,1](rvNakajimaCreate[2,2](x));
let right = rvNakajimaCreate[2,2](rvNakajimaAnn[2,1](x));
assert add(left,neg(right)) == rvNakajimaScale[-2](x);

rvNakajimaPair(rvNakajimaCreate[1,0](vacuum),rvNakajimaCreate[1,23](vacuum));
rvNakajimaPair(rvNakajimaCreate[2,0](vacuum),rvNakajimaCreate[2,23](vacuum));
rvNakajimaPair(rvNakajimaMonomial[B]{(2,0,3)},rvNakajimaMonomial[B]{(2,23,3)});

// The first two middle colors have intersection -1 and squares -2.
let w = rvNakajimaMonomial[B]{(1,1,1),(1,2,1)};
rvNakajimaPair(w,w);
