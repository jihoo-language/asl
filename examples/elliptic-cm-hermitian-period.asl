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

// The central imaginary unit in the period is independent of CM conjugation.
let P=rvCMPeriod(J,D);
let SP=rvCMPeriodFourier(P);
rvCMPeriodData(P);
rvCMPeriodData(SP);
assert rvCMPeriodDet(P)==rvCMCharge[A,0,16,1];
assert rvCMPeriodFactor(SP)==rvCMCharge[A,0,16,1];
assert rvCMPeriodFourier(SP)==P;
let c=rvCMPeriodCharge(P,v);
let sc=rvCMPeriodCharge(SP,rvCMAbelianFourier(v));
assert c==rvCMCharge[A,17,-20,1];
assert sc==rvCMCharge[A,-20,-17,16];
assert rvCMChargeMul(rvCMPeriodFactor(SP),sc)==c;

let T=rvCMPeriodTensor(P,J);
assert rvCMPeriodCharge(T,rvCMAbelianTensor(v,J))==c;
let U=rvCMPeriodPull(P,M);
assert rvCMPeriodCharge(U,v)==rvCMPeriodCharge(P,rvCMAbelianPush(M,v));
assert rvCMPeriodDet(U)==rvCMCharge[A,0,80,1];
let I=rvCMAbelianDiv[A,1,1](zero[RVCMHomCoordinate(H)]);
let W1=rvCMPeriodFourier(rvCMPeriodTensor(P,I));
let W2=rvCMPeriodFourier(rvCMPeriodTensor(W1,I));
assert rvCMPeriodFactor(W2)==rvCMChargeMul(rvCMPeriodFactor(W1),rvCMPeriodDet(rvCMPeriodTensor(W1,I)));
assert rvCMPeriodFourier(rvCMPeriodTensor(W2,I))==P;

// Numerical alignment on the scalar slice Z=(b+i t) I.
let z0=zero[RVCMAbelianDiv(A)];
let a=rvCMAbelianMukai[1,-4](z0);
let b=rvCMAbelianLine(neg(I));
let wall=rvCMPeriodWall(I,z0,a,b);
rvCMPeriodWallData(wall);
rvCMPeriodWallAt[-5,2,3,2](wall);
let alignment=rvCMPeriod(rvCMAbelianScale[-5,2](I),rvCMAbelianScale[3,2](I));
assert rvCMPeriodCharge(alignment,a)==rvCMCharge[A,0,15,2];
assert rvCMPeriodCharge(alignment,b)==rvCMCharge[A,0,9,2];
rvCMChargeRay(rvCMPeriodCharge(alignment,a),rvCMPeriodCharge(alignment,b));
