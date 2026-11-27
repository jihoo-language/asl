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

let D=rvCMAbelianDiv[A,3,3](rvCMHomCoords[H](g));
rvcmthetakernel K6(D);
rvCMQuotientProfiles[K6];

// Inspect the actual polarized quotients of the nonprincipal F41 product.
rvCMQuotientTorelliProfiles[K6];
rvcmquotient Split(rvCMThetaLagrangian[K6,0]);
rvcmquotient Jacobian(rvCMThetaLagrangian[K6,4]);
rvCMQuotientNS[Split];
rvCMQuotientTrans[Split];
rvCMQuotientTorelli[Split];
rvCMQuotientFactorFrob[Split];
rvCMQuotientNS[Jacobian];
rvCMQuotientTrans[Jacobian];
rvCMQuotientTorelli[Jacobian];

let u=rvCMQuotientProjector[Split,0];
let v=rvCMQuotientProjector[Split,1];
assert rvCMQuotientEndCompose(u,u)==u;
assert rvCMQuotientEndCompose(v,v)==v;
assert rvCMQuotientRosati(u)==u;
assert rvCMQuotientRosati(v)==v;
let point=rvCMQuotientVec[Split]{1,2,3,4};
assert add(rvCMQuotientEndApply(u,point),rvCMQuotientEndApply(v,point))==point;
assert rvCMQuotientEndApply(rvCMQuotientEndCompose(u,v),point)==zero[RVCMQuotientLattice(Split)];
let ns=rvCMQuotientNSBasis[Split,0];
let tr=rvCMQuotientTransBasis[Split,0];
assert add(tr,neg(tr))==zero[RVCMQuotientWedge(Split)];
let symmetric=rvCMQuotientToRosati(ns);
assert rvCMQuotientRosati(symmetric)==symmetric;
rvCMQuotientWedgePair(ns,tr);
