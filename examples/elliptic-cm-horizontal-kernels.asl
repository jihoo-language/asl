rvcmclasspoly H3[-3];
rvcmclasspolylift H12[H3,2];
rvcmclasspolylift H48[H12,2];
rvunram F[13,1]{0,1};
rvec E[F](rvUnramElem[F]{1},rvUnramElem[F]{2});
rvecendorder M[E,H48];

let k0 = rvECCMKernel[M,7,0];
let k1 = rvECCMKernel[M,7,1];
assert k0 == rveckernel[E,7]{[8],[7],[8],[1]};
assert k1 == rveckernel[E,7]{[12],[6],[12],[1]};
rvECCMHorizontalData[M,7,0];
rvECCMHorizontalData[M,7,1];

rveccmhorizontal I[M,Q,7,0];
rveccmhorizontal J[M,R,7,1];
assert rvECGetKernel(rvECIso[I]) == k0;
assert rvECGetKernel(rvECIso[J]) == k1;
assert rvECEndJ[E] == rvUnramElem[F]{6};
assert rvECEndJ[Q] == rvUnramElem[F]{4};
assert rvECEndJ[R] == rvUnramElem[F]{4};
assert rvECA[Q] == rvUnramElem[F]{5};
assert rvECB[Q] == rvUnramElem[F]{8};
assert rvECA[R] == rvUnramElem[F]{2};
assert rvECB[R] == rvUnramElem[F]{5};

rvecendorder N[Q,H48];
rvECEndData[N];
let phi = rvECMap(rvECIso[I]);
let dual = rvECDual(phi);
assert rvECMapCompose(dual,phi) == rvECMultMap[E,7];
let back = rvECCMKernel[N,7,1];
assert rvECMapKCoeff[0](dual) == rvECKCoeff[0](back);
assert rvECMapKCoeff[1](dual) == rvECKCoeff[1](back);
assert rvECMapKCoeff[2](dual) == rvECKCoeff[2](back);
