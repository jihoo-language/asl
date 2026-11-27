// Alternating point and full boundaries give the order-three Fourier matrix.
rvgroup G[3];
rvdwcoeff P[G,3];
rvdouble D[P](rvPostnikovK[P]);
rvdoubleboundary O[D]{0};
rvdoubleboundary F[D]{0,1,2};
rvdoubledefects A[O,F];
rvdoubledefects B[F,O];
rvdoubledefects C[O,O];
rvdoubledefects E[F,F];
let x=rvDoubleDefectBasis[A,0];
let y=rvDoubleDefectBasis[B,0];
let z=rvDoubleDefectBasis[C,0];
let f=rvDoubleDefectSixJ[C,E](x,y,x,x);
rvDoubleDefectSixJData(f);
show rvDoubleDefectSixJNorm[0,0](f);
show rvDoubleDefectSixJNorm[1,0](f);
show rvDoubleDefectSixJEntry[0,0](f);
show rvDoubleDefectSixJEntry[1,1](f);
show rvDoubleDefectSixJInvEntry[1,1](f);
assert rvCycloQMul(rvDoubleDefectSixJEntry[1,1](f),rvDoubleDefectSixJInvEntry[1,1](f))==rvCycloRat[3](rvCyclo[9,0]);
rvDoubleDefectPentagon[C,E,C,A,B](x,y,x,y,z);

// The type-III cocycle on (C2)^3 is k8(a,b,c)=a1*b2*c3.
rvgroup G8{
  {0,1,2,3,4,5,6,7},
  {1,0,3,2,5,4,7,6},
  {2,3,0,1,6,7,4,5},
  {3,2,1,0,7,6,5,4},
  {4,5,6,7,0,1,2,3},
  {5,4,7,6,1,0,3,2},
  {6,7,4,5,2,3,0,1},
  {7,6,5,4,3,2,1,0},
};
rvdwcoeff P8[G8,2];
let k8=rvPostnikovCochain[P8,3]{
  0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,
};
rvdouble D8[P8](k8);

// Compose a plane-to-point sector with a point-to-plane sector.
rvdoubleboundary Plane[D8]{0,1,2,3};
rvdoubleboundary Point[D8]{0};
rvdoubledefects Y[Plane,Point];
rvdoubledefects X[Point,Plane];
rvdoubledefects PP[Plane,Plane];
rvdoubledefects QQ[Point,Point];
let y0=rvDoubleDefectBasis[Y,0];
let x4=rvDoubleDefectBasis[X,1];

let y4=rvDoubleDefectBasis[Y,1];
let x0=rvDoubleDefectBasis[X,0];
let mix=rvDoubleDefectSixJ[PP,QQ](y0,x4,y4,y0);
rvDoubleDefectSixJData(mix);
show rvDoubleDefectSixJEntry[0,0](mix);
show rvDoubleDefectSixJEntry[0,2](mix);
rvDoubleDefectPentagon[PP,QQ,PP,Y,X](y0,x4,y4,x0,rvDoubleDefectBasis[PP,0]);
