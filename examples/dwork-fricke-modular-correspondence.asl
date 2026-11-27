rvDworkFrickeSeries[8];
rvDworkFrickeIdentity;

rvunram F[7,1]{0,1};
let t=rvUnramElem[F]{2};
let tt=rvDworkFrickeT(t);
assert tt==rvUnramElem[F]{4};
assert rvDworkFrickeT(tt)==t;
assert rvDworkFrickeZ(t)==rvDworkFrickeZ(tt);
assert rvDworkFrickeJ[1](t)==rvDworkFrickeJ[2](tt);
assert rvDworkFrickeJ[2](t)==rvDworkFrickeJ[1](tt);
rvDworkFrickeBase(rvDworkFrickeZ(t));
rvDworkFrickeFiber(t);

rvdworkfricke E,Q(t);
let I=rvDworkFrickeMap[E,Q](t);
let J=rvDworkFrickeDual[E,Q](t);
rvECMapDegree(I);
rvECMapCertify(I);
assert rvECMapCompose(J,I)==rvECMultMap[E,2];
assert rvECMapCompose(I,J)==rvECMultMap[Q,2];
let point=rvECPoint[E](rvUnramElem[F]{2},rvUnramElem[F]{3});
let image=rvECPoint[Q](rvUnramElem[F]{3},rvUnramElem[F]{1});
assert rvECMapApply(I,point)==image;
assert rvECMapApply(I,rvECPoint[E](rvUnramElem[F]{4},rvUnramElem[F]{}))==rvECInfinity[Q];
assert rvECMapApply(I,rvECInfinity[E])==rvECInfinity[Q];
rvECPointCount[E];
rvECPointCount[Q];

rvDworkFrickeBase(rvUnramElem[F]{3});
rvunram K[7,1]{1,4,1};
let r=rvUnramElem[K]{0,1};
assert rvDworkFrickeT(r)==rvUnramFrob[1](r);
assert rvDworkFrickeZ(r)==rvUnramElem[K]{3};
assert rvDworkFrickeJ[1](r)==rvUnramFrob[1](rvDworkFrickeJ[2](r));

rvDworkFrickeUnit[7,2,2];
rvDworkFrickeUnit[11,2,2];
rvDworkFrickeUnit[17,2,2];
rvDworkFrickeUnit[29,3,2];
rvDworkFrickeUnit[7,3,2];
