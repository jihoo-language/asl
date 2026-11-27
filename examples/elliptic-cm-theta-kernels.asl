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

// Use the positive Hermitian divisor with diagonal (3,3) and lower-left entry g.
let D=rvCMAbelianDiv[A,3,3](rvCMHomCoords[H](g));
rvcmthetakernel Theta(D);
rvCMThetaData[Theta];
rvCMThetaLagrangians[Theta];
let e=rvCMThetaElem[Theta]{0,0,1,0};
let ftheta=rvCMThetaElem[Theta]{0,0,0,1};
rvCMThetaPair(e,ftheta);
rvCMThetaPair(ftheta,e);
rvCMThetaOrder(e);
rvCMThetaLift(e);
rvCMThetaLift(ftheta);
assert rvCMThetaFromLattice[Theta]{0,0,1,0}==e;
assert rvCMThetaFromLattice[Theta]{-1,0,0,0}==ftheta;
let Sep=rvCMThetaSubgroup[Theta]{e};
assert rvCMThetaPerp(Sep)==Sep;
rvCMThetaSubgroupData(Sep);
let op=rvCMThetaOp[0](e);
let clock=rvCMThetaOp[0](ftheta);
let center=rvCMThetaOp[1](zero[RVCMTheta(Theta)]);
assert rvCMThetaMul(op,clock)==rvCMThetaMul(center,rvCMThetaMul(clock,op));
rvCMThetaAction[0,2](rvCMThetaMul(op,clock));

// Nonfree Lagrangians contribute to polarization type (2,4).
let off=zero[RVCMHomCoordinate(H)];
rvcmthetakernel Double(rvCMAbelianDiv[A,2,2](off));
rvcmthetakernel Mixed(rvCMAbelianDiv[A,2,4](off));
rvCMThetaLagrangians[Double];
rvCMThetaLagrangians[Mixed];

// The same four elliptic morphism blocks act on actual full E[3] and F[3].
rvunram L[41,1]{38,0,1};
rvhyperelliptic EE[K]{rvUnramElem[K]{17},rvUnramElem[K]{1},rvUnramElem[K]{0},rvUnramElem[K]{1}};
rvhyperelliptic FF[K]{rvUnramElem[K]{3},rvUnramElem[K]{30},rvUnramElem[K]{0},rvUnramElem[K]{1}};
rvjacbase EL[EE,L];
rvjacbase FL[FF,L];
let P=rvJacPoint[EL](rvUnramElem[L]{34},rvUnramElem[L]{6});
let Q=rvJacPoint[EL](rvUnramElem[L]{23},rvUnramElem[L]{0,16});
let R3=rvJacPoint[FL](rvUnramElem[L]{7},rvUnramElem[L]{8});
let T=rvJacPoint[FL](rvUnramElem[L]{31},rvUnramElem[L]{0,13});
rvjactorsionbasis BE[EL,3](rvJacWeil[3](P,Q)){P,Q};
rvjactorsionbasis BF[FL,3](rvJacWeil[3](R3,neg(T))){R3,neg(T)};
rvCMThetaTorsion[Theta,BE,BF];
