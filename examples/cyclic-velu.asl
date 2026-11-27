// The two curves have the same point count. Their rational kernel counts differ.
rvunram V[13,1]{0,1};
rvec C[V](rvUnramElem[V]{4},rvUnramElem[V]{1});
rvec D[V](rvUnramElem[V]{0},rvUnramElem[V]{2});
rvECPointCount[C];
rvECPointCount[D];
rvECIsogenies[C,3];
rvECIsogenies[D,3];

let kc = rveckernel[C,3]{[2],[1]};
rvecisogeny I[C,Q](kc);
let i = rvECIso[I];
assert rvECGetKernel(i) == kc;
assert rvECA[Q] == zero[RVUnram(V)];
assert rvECB[Q] == rvUnramElem[V]{11};
let p = rvECPoint[C](rvUnramElem[V]{0},rvUnramElem[V]{1});
let ip = rvECPoint[Q](rvUnramElem[V]{1},rvUnramElem[V]{8});
assert rvECApply(i,p) == ip;
assert rvECApply(i,rvECScale[2](p)) == add(ip,ip);
assert rvECApply(i,rvECInfinity[C]) == zero[RVECPoint(Q)];
show rvECCurve[Q];
rvECCertify[I];

// The degree-five isogeny uses an irreducible quadratic kernel over F7.
rvunram W[7,1]{0,1};
rvec E[W](rvUnramElem[W]{1},rvUnramElem[W]{1});
rvECPointCount[E];
rvECIsogenies[E,5];
let k = rveckernel[E,5]{[1],[3],[1]};
rvecisogeny J[E,T](k);
let j = rvECIso[J];
assert rvECA[T] == rvUnramElem[W]{2};
assert rvECB[T] == rvUnramElem[W]{1};
assert rvECXNum[4](j) == rvUnramElem[W]{6};
assert rvECXDen[2](j) == rvUnramElem[W]{4};
assert rvECDNum[1](j) == rvUnramElem[W]{6};
assert rvECDDen[3](j) == rvUnramElem[W]{3};
let r = rvECPoint[E](rvUnramElem[W]{2},rvUnramElem[W]{2});
let jr = rvECPoint[T](rvUnramElem[W]{1},rvUnramElem[W]{2});
assert rvECApply(j,r) == jr;
assert rvECApply(j,add(r,r)) == rvECScale[2](jr);
assert rvECApply(j,neg(r)) == neg(jr);
assert rvECScale[5](r) == rvECInfinity[E];
show j;
rvECCertify[J];

rvunram F[5,1]{2,0,1};
let z = rvUnramElem[F]{0,1};
rvec A[F](z,rvUnramElem[F]{1});
rvECPointCount[A];
let extension = rveckernel[A,3]{[1,1],[1]};
rvecisogeny L[A,B](extension);
assert rvECA[B] == z;
assert rvECB[B] == rvUnramElem[F]{4,2};
rvECCertify[L];
