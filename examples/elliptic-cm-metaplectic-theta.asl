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

// The non-diagonal divisor gives the six-dimensional theta space.
let D=rvCMAbelianDiv[A,3,3](rvCMHomCoords[H](g));
rvcmthetakernel K6(D);
rvCMMetaVerify[K6];
let S=rvCMMetaF[K6];
let Q6=rvCMMetaT[K6];
let Ga=rvCMMetaGauss[K6];
let one12=rvCycloRat[1](rvCyclo[12,0]);
let six12=rvCycloRat[1](rvCycloScale[6](rvCyclo[12,0]));
let gauss12=rvCycloRat[1](add(neg(rvCyclo[12,0]),add(rvCycloScale[2](rvCyclo[12,1]),add(rvCycloScale[2](rvCyclo[12,2]),neg(rvCyclo[12,3])))));
assert Ga==gauss12;
assert rvCycloQMul(Ga,rvCMMetaGaussBar[K6])==six12;
assert rvCMMetaCompose(S,rvCMMetaFInv[K6])==rvCMMetaIdentity[K6];
assert rvCMMetaCompose(Q6,rvCMMetaTInv[K6])==rvCMMetaIdentity[K6];
let F2=rvCMMetaCompose(S,S);
assert F2==rvCMMetaScale(rvCMMetaParity[K6],six12);
let FT=rvCMMetaCompose(S,Q6);
assert rvCMMetaCompose(rvCMMetaCompose(FT,FT),FT)==rvCMMetaScale(F2,Ga);
let sa=rvCMMetaFAuto[K6];
let ta=rvCMMetaTAuto[K6];
assert rvCMMetaIntertwine(sa)==S;
assert rvCMMetaIntertwine(ta)==Q6;
assert rvCMMetaCocycle(sa,sa)==six12;
assert rvCMMetaCoeff[5](rvCMMetaApply(S,rvCMMetaBasis[K6,0]))==one12;

// The non-product shear comes from the quadratic character x1*x2/2.
let off=zero[RVCMHomCoordinate(H)];
rvcmthetakernel K8(rvCMAbelianDiv[A,2,4](off));
rvCMMetaVerify[K8];
let X1=rvCMMetaLift[0](rvCMThetaElem[K8]{1,0,0,0});
let Z1=rvCMMetaLift[0](rvCMThetaElem[K8]{0,1,0,0});
let X2=rvCMMetaLift[0](rvCMThetaElem[K8]{0,0,1,0});
let Z2=rvCMMetaLift[0](rvCMThetaElem[K8]{0,0,0,1});
let Y1=rvCMMetaLift[0](rvCMThetaElem[K8]{1,0,0,2});
let Y2=rvCMMetaLift[0](rvCMThetaElem[K8]{0,1,1,0});
let M=rvCMMetaAuto[K8](Y1,Z1,Y2,Z2);
let U=rvCMMetaIntertwine(M);
let W=rvCMMetaF[K8];
let wa=rvCMMetaFAuto[K8];
let z8=rvCycloRat[1](rvCyclo[8,0]);
let n8=neg(z8);
assert rvCMMetaEntry[5,5](U)==n8;
assert rvCMMetaEntry[6,6](U)==z8;
assert rvCMMetaCompose(U,U)==rvCMMetaIdentity[K8];
assert rvCMMetaCompose(U,rvCMMetaHeis(X1))==rvCMMetaCompose(rvCMMetaHeis(Y1),U);
let wm=rvCMMetaAutoCompose(wa,M);
let scalar=rvCMMetaCocycle(wa,M);
assert rvCMMetaCompose(W,U)==rvCMMetaScale(rvCMMetaIntertwine(wm),scalar);
let ww=rvCMMetaAutoCompose(wa,wa);
let cww=rvCMMetaCocycle(wa,wa);
assert rvCycloQMul(cww,rvCMMetaCocycle(ww,M))==rvCycloQMul(scalar,rvCMMetaCocycle(wa,wm));

// A maximal isotropic subgroup of exponent two in type (2,4).
let Nonfree=rvCMMetaProjector[K8]{X1,rvCMMetaLiftCompose(X2,X2),rvCMMetaLiftCompose(Z2,Z2)};
assert rvCMMetaCompose(Nonfree,Nonfree)==Nonfree;
assert rvCMMetaTrace(Nonfree)==z8;
assert rvCMMetaEntry[0,0](Nonfree)==rvCycloRat[4](rvCyclo[8,0]);
show Ga;
show rvCMMetaCocycle(wa,wa);
