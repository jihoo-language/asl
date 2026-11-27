// Two independent maps between the ordinary F41 curves of discriminant -20.
rvunram K[41,1]{0,1};
rvec E[K](rvUnramElem[K]{1},rvUnramElem[K]{17});
rvecisogeny I2[E,F](rveckernel[E,2]{[6],[1]});
rvecisogeny I3[E,F3](rveckernel[E,3]{[7],[1]});
let f=rvECMap(rvECIso[I2]);
let g3=rvECMap(rvECIso[I3]);
let g=rvECMapCompose(rvECChart[F3,F](rvUnramElem[K]{12}),g3);
rvcmhom H(f,g);

// Undo the target chart before applying the normalized dual of I3.
let fd=rvECDual(f);
let gd=rvECMapCompose(rvECDual(g3),rvECChart[F,F3](rvUnramElem[K]{24}));
rvcmhomdual HD[H](fd,gd);
let id=rvECMultMap[E,1];
let w=rvECSum(rvECMapCompose(fd,g),rvECNegate(id));
assert rvECMapCompose(w,w)==rvECMultMap[E,-5];
assert rvECMapCompose(gd,f)==rvECSum(id,rvECNegate(w));
rvcmhomframe R(id,w);

rvcmhomcompose P[HD,H,R];
rvCMHomProductData[P];
assert rvCMHomProductColumn[P,0]==rvCMHomCoord[R]{2,0};
assert rvCMHomProductColumn[P,1]==rvCMHomCoord[R]{1,1};
assert rvCMHomProductColumn[P,2]==rvCMHomCoord[R]{1,-1};
assert rvCMHomProductColumn[P,3]==rvCMHomCoord[R]{3,0};
assert rvCMHomCoordMap(rvCMHomProductBasis[P,0])==id;
assert add(rvCMHomProductColumn[P,3],neg(rvCMHomProductColumn[P,0]))==rvCMHomCoord[R]{1,0};
assert rvECSum(rvCMHomCoordMap(rvCMHomProductColumn[P,3]),rvECNegate(rvCMHomCoordMap(rvCMHomProductColumn[P,0])))==id;
assert rvCMHomCoordMap(rvCMHomProductBasis[P,1])==w;
rvCMHomCoordinates(rvCMHomCoords[R](w));

rvcmorder O[-20];
let ideal=rvCMHomClass[O,H];
assert ideal==rvCMHomClass[O,HD];
assert rvCMCompose(ideal,rvCMHomClass[O,HD])==rvCMPrincipal[O];
assert rvCMCompose(ideal,ideal)==rvCMPrincipal[O];
assert rvCMHomClass[O,R]==rvCMPrincipal[O];

let fc=rvCMHomCoords[H](f);
let gdc=rvCMHomCoords[HD](gd);
assert rvCMHomCoordCompose[R](gdc,fc)==rvCMHomCoord[R]{1,-1};
let negative=rvCMHomCoord[R]{-1,0};
assert rvCMHomCoordInverse[R](negative)==negative;
