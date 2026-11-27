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

let idF=rvECMultMap[F,1];
let wF=rvECSum(rvECMapCompose(g,fd),rvECNegate(idF));
rvcmhomframe RF(idF,wF);
rvcmabelian A[H,HD,R,RF];

// The first Hom generator gives the rational identification of the two factors.
let D=rvCMAbelianDiv[A,3,3](rvCMHomCoords[H](g));
let J=rvCMAbelianDiv[A,2,4](rvCMHomCoords[H](f));
let v=rvCMAbelianMukai[2,-1](D);
let u=rvCMAbelianMukai[3,5](J);
let M=rvCMAbelianMap[A](id,fd,g,idF);
rvCMAbelianMapData(M);
rvCMAbelianDet(D);
rvCMAbelianData(v);
rvCMAbelianData(rvCMAbelianPull(M,v));
rvCMAbelianData(rvCMAbelianFourier(rvCMAbelianLine(D)));
let unit=rvCMAbelianMukai[1,0](zero[RVCMAbelianDiv(A)]);
assert rvCMAbelianTensor(unit,D)==rvCMAbelianLine(D);
rvCMAbelianPair(v,u);
rvCMAbelianPair(rvCMAbelianFourier(v),rvCMAbelianFourier(u));
rvCMAbelianPair(rvCMAbelianTensor(v,J),rvCMAbelianTensor(u,J));
rvCMAbelianPair(rvCMAbelianPull(M,v),u);
rvCMAbelianPair(v,rvCMAbelianPush(M,u));
rvCMAbelianPair(rvCMAbelianPull(M,v),rvCMAbelianPull(M,u));
assert rvCMAbelianFourier(rvCMAbelianFourier(v))==v;
assert rvCMAbelianTensor(rvCMAbelianTensor(v,J),neg(J))==v;
assert rvCMAbelianPush(M,rvCMAbelianPull(M,v))==rvCMAbelianScale[5,1](v);
assert rvCMAbelianPull(M,rvCMAbelianPush(M,v))==rvCMAbelianScale[5,1](v);
let adjoint=rvCMAbelianDagger(M);
assert rvCMAbelianFourier(rvCMAbelianPull(M,v))==rvCMAbelianPush(adjoint,rvCMAbelianFourier(v));
assert rvCMAbelianBlock[1](adjoint)==gd;
assert rvCMAbelianBlock[2](adjoint)==f;
let determinantMap=rvCMAbelianCompose(M,rvCMAbelianAdjugate(M));
assert rvCMAbelianBlock[0](determinantMap)==rvECNegate(w);
assert rvCMAbelianBlock[3](determinantMap)==rvECNegate(wF);
let I=rvCMAbelianDiv[A,1,1](zero[RVCMHomCoordinate(H)]);
let su=rvCMAbelianFourier(rvCMAbelianTensor(u,I));
let ssu=rvCMAbelianFourier(rvCMAbelianTensor(su,I));
assert rvCMAbelianFourier(rvCMAbelianTensor(ssu,I))==u;
let graph=rvCMAbelianGraph(D);
assert graph==rvCMHomCycle[H]{4,6,0,-1};
assert rvCMAbelianFromGraph[A](graph)==D;
