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

// Quotient by subgroups of orders two and six on a common isotropic line.
let e=rvCMThetaElem[K6]{0,0,1,0};
let e3=add(e,add(e,e));
rvcmquotient Small(rvCMThetaSubgroup[K6]{e3});
rvcmquotient Full(rvCMThetaSubgroup[K6]{e});
rvCMQuotientData[Small];
rvCMQuotientData[Full];
rvCMQuotientFrob[Small];
rvCMQuotientFrob[Full];
rvCMQuotientTransitionData[Small,Full];

let x=rvCMPolarVec[K6]{1,0,0,0};
let y=rvCMPolarVec[K6]{0,1,0,0};
let ys=rvCMQuotientInclude[Small](y);
let yf=rvCMQuotientInclude[Full](y);
assert rvCMQuotientTransition[Full](ys)==yf;
let dy=rvCMQuotientPolarize(y);
let dys=rvCMQuotientPolarize(ys);
let dyf=rvCMQuotientPolarize(yf);
assert rvCMQuotientPullback[Full](dyf)==dy;
assert rvCMQuotientTransitionPullback[Small](dyf)==dys;
rvCMQuotientPair(x,y);
rvCMQuotientPair(rvCMQuotientInclude[Full](x),yf);
show rvCMQuotientTorsionLift[Full](e);

// The full integral centralizer is larger than its diagonal quadratic center.
rvCMQuotientEndoData[Full];
let a=rvCMQuotientEndBasis[Full,0];
let b=rvCMQuotientEndBasis[Full,1];
assert rvCMQuotientRosati(rvCMQuotientRosati(a))==a;
assert rvCMQuotientRosati(rvCMQuotientEndCompose(a,b))==rvCMQuotientEndCompose(rvCMQuotientRosati(b),rvCMQuotientRosati(a));
let frob=rvCMQuotientFrobenius[Full];
let scalar41=rvCMQuotientEnd[Full]{41,0,0,0, 0,41,0,0, 0,0,41,0, 0,0,0,41};
assert rvCMQuotientEndCompose(rvCMQuotientRosati(frob),frob)==scalar41;
assert rvCMQuotientEndApply(scalar41,yf)==rvCMQuotientVec[Full]{0,41,0,0};
