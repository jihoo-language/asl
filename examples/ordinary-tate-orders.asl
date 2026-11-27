// Frobenius preserves both coordinates of multiplication by 5 on E[9].
rvunram F61[61,1]{0,1};
let z61 = rvUnramElem[F61]{0};
let b61 = rvUnramElem[F61]{2};
rvec E[F61](z61,b61);
assert rvEllFrobCoordinates[E,3,2,5] == (S.hi,S.hi);
assert rvEllFrobCoordinates[E,3,2,-5] == (S.hi,S.lo);
rvEllFrobCert[E,3,2,-5];
rvEllLocalOrder[E,3];
rvellipticcm K[E];
let omega = rvEllLocalGen[K,E,3];
assert add(rvQpQuad[K,5,0,1],rvQpQuadMul(rvQpQuad[K,9,0,1],omega)) == rvCMPi[K];
assert add(omega,rvQpQuadConj(omega)) == rvQpQuad[K,-1,0,1];
assert rvQpQuadMul(omega,rvQpQuadConj(omega)) == rvCMPi[K,0];

// These are actual degree-three Vélu isogenies over F13.
rvunram F13[13,1]{0,1};
let zero13 = rvUnramElem[F13]{0};
let one13 = rvUnramElem[F13]{1};
let two13 = rvUnramElem[F13]{2};
rvec Lower[F13](rvUnramElem[F13]{4},one13);
rvec Upper[F13](zero13,two13);
rvEllLocalOrder[Lower,3];
rvEllLocalOrder[Upper,3];
rvecisogeny Up[Lower,UpTarget](rveckernel[Lower,3]{[2],[1]});
rvecisogeny Flat[Upper,FlatTarget](rveckernel[Upper,3]{[0],[1]});
rvecisogeny DownA[Upper,DownATarget](rveckernel[Upper,3]{[6],[1]});
rvecisogeny DownB[Upper,DownBTarget](rveckernel[Upper,3]{[5],[1]});
rvecisogeny DownC[Upper,DownCTarget](rveckernel[Upper,3]{[2],[1]});
rvEllEdge[Up];
rvEllEdge[Flat];
rvEllEdge[DownA];
rvEllEdge[DownB];
rvEllEdge[DownC];
