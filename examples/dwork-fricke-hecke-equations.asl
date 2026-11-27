rvDworkHeckeEquation[3];
rvDworkHeckeEquation[5];
rvDworkHeckeCoeff[3,3,1];
rvDworkHeckeCoeff[3,2,2];
rvDworkHeckeCoeff[3,1,1];
rvDworkHeckeCoeff[5,0,0];

rvunram F[41,1]{0,1};
let h=rvUnramElem[F]{};
let seven=rvUnramElem[F]{7};
let other=rvUnramElem[F]{26};
assert rvDworkHeckeEval[3](h,seven)==zero[RVUnram(F)];
assert rvDworkHeckeEval[3](h,other)==zero[RVUnram(F)];
assert rvDworkHeckeSpecCoeff[3,0](h)==rvUnramElem[F]{17};
rvDworkHeckeNeighbors[3](h);
rvDworkHeckeNeighbors[5](rvUnramElem[F]{1});

rvunram C[3,1]{-1,-1,0,1};
let w=rvUnramElem[C]{0,1};
assert rvDworkHeckeEval[3](w,rvUnramFrob[1](w))==zero[RVUnram(C)];
assert rvDworkHeckeEval[3](w,rvUnramFrob[2](w))==zero[RVUnram(C)];
rvDworkHeckeNeighbors[3](w);

// An actual level-three elliptic edge preserved by the degree-two quotient.
rvunram K[23,1]{0,1};
let t=rvUnramElem[K]{1};
let s=rvUnramElem[K]{12};
let ht=rvDworkHauptmodulT(t);
let hs=rvDworkHauptmodulT(s);
assert ht==rvUnramElem[K]{4};
assert hs==rvUnramElem[K]{17};
assert ht==rvDworkHauptmodulZ(rvDworkFrickeZ(t));
assert hs==rvDworkHauptmodulT(rvDworkFrickeT(s));
assert rvDworkHeckeEval[3](ht,hs)==zero[RVUnram(K)];
assert rvModularPolyEval[3](rvDworkFrickeJ[1](t),rvDworkFrickeJ[1](s))==zero[RVUnram(K)];
assert rvModularPolyEval[3](rvDworkFrickeJ[2](t),rvDworkFrickeJ[2](s))==zero[RVUnram(K)];
rvDworkHeckeNeighbors[3](ht);

rvdworkfricke E,V(t);
rvdworkfricke P,W(s);
let upper=rvDworkFrickeMap[E,V](t);
let lower=rvDworkFrickeMap[P,W](s);
let k=rvECPoint[E](rvUnramElem[K]{2},rvUnramElem[K]{6});
let kp=rvECPoint[V](rvUnramElem[K]{12},rvUnramElem[K]{14});
assert rvECMapApply(upper,k)==kp;
let kernel=rveckernel[E,3]{[21],[1]};
let kernelp=rveckernel[V,3]{[11],[1]};
rvecisogeny L[E,Q](kernel);
rvecisogeny R[V,U](kernelp);
let scale=rvUnramElem[K]{6};
let left=rvECMapCompose(rvECChart[Q,P](scale),rvECMap(rvECIso[L]));
let right=rvECMapCompose(rvECChart[U,W](scale),rvECMap(rvECIso[R]));
let first=rvECMapCompose(lower,left);
let second=rvECMapCompose(right,upper);
assert first==second;
rvECMapDegree(first);
rvECMapCertify(first);
