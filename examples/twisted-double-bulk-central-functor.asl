// The type-III cocycle on (C2)^3 is k(a,b,c)=a1*b2*c3.
rvgroup G{
  {0,1,2,3,4,5,6,7},
  {1,0,3,2,5,4,7,6},
  {2,3,0,1,6,7,4,5},
  {3,2,1,0,7,6,5,4},
  {4,5,6,7,0,1,2,3},
  {5,4,7,6,1,0,3,2},
  {6,7,4,5,2,3,0,1},
  {7,6,5,4,3,2,1,0},
};
rvdwcoeff P[G,2];
let k=rvPostnikovCochain[P,3]{
  0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,
};
rvdouble D[P](k);

// Compose a plane-to-point sector with a point-to-plane sector.
rvdoubleboundary Plane[D]{0,1,2,3};
rvdoubleboundary Point[D]{0};
rvdoubledefects Y[Plane,Point];
rvdoubledefects X[Point,Plane];
rvdoubledefects PP[Plane,Plane];
rvdoubledefects QQ[Point,Point];
let y0=rvDoubleDefectBasis[Y,0];
let x4=rvDoubleDefectBasis[X,1];

// Induce the bulk simple of flux 4 with its two-dimensional projective action.
let central=rvDoubleBulkBoundary[PP,14];
rvDoubleRelativeData(central);
rvDoubleRelativeAction[0,2,4,0](central);
rvDoubleBulkBoundaryValidate[PP];
rvDoubleBulkBoundaryValidate[QQ];
rvDoubleBulkBoundaryFusion[PP,14,14];

// Cross the distinct plane and point boundaries.
let f=rvDoubleCentralMove[PP,QQ,14](y0);
rvDoubleCentralData(f);
rvDoubleCentralImage[6,0](f);
rvDoubleCentralBasis[0,6,0](f);
show rvDoubleCentralEntry[6,1,0](f);
show rvDoubleCentralInvEntry[6,0,1](f);
rvDoubleCentralUnit[PP,14];
rvDoubleCentralWhisker[PP,QQ,PP,PP,14](y0,x4);

// On the point boundary of the twisted cyclic-three double, crossing a
// degree-one simple uses the inverse complex half-braiding.
rvgroup G3[3];
rvdwcoeff P3[G3,3];
rvdouble D3[P3](rvPostnikovCochain[P3,3]{0,1,1,1,0,2,2,2});
rvdoubleboundary Pnt[D3]{0};
rvdoubledefects E[Pnt,Pnt];
let a=rvDoubleDefectBasis[E,1];
let b=rvDoubleDefectBasis[E,2];
let f3=rvDoubleCentralMove[E,E,3](a);
rvDoubleCentralImage[2,0](f3);
show rvDoubleCentralEntry[2,0,0](f3);
show rvDoubleCentralInvEntry[2,0,0](f3);
rvDoubleCentralWhisker[E,E,E,E,3](a,b);
