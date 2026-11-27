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




// The plane has four invertibles and one dimension-two self defect.
rvdoubleboundary Plane[D]{0,1,2,3};
rvdoubledefects C[Plane,Plane];
rvDoubleDefectData[C];
rvDoubleDefectCoset[C,0];
rvDoubleDefectCoset[C,1];
rvDoubleDefectSimple[C,4];
let m=rvDoubleDefectBasis[C,4];
show rvDoubleDefectTheta[C,1,1,2];
rvDoubleDefectRep[1,1](m);
rvDoubleDefectRep[2,0](m);
assert rvDoubleDefectOpposite[C](m)==m;
rvDoubleDefectFuse(m,m);
rvDoubleDefectFusionProjector[0,0,3](m,m);

// The torsion difference between two Klein-four boundaries gives one 2D defect.
rvgroup K{{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}};
rvdwcoeff O[K,2];
rvdouble U[O](rvPostnikovK[O]);
rvdoubleboundary A[U]{0,1,2,3};
rvdoubleboundaryshift B[A]{0,(1,2),(1,2),0,0,0,0,(1,2),(1,2)};
rvdoubledefects AB[A,B];
rvdoubledefects BA[B,A];
rvDoubleDefectData[AB];
let v=rvDoubleDefectBasis[AB,0];
assert rvDoubleDefectOpposite[BA](v)==rvDoubleDefectBasis[BA,0];
