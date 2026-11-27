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

// The plane boundary has a nonpointed simple and a 40-dimensional tube.
rvdoubletube T[PP];
rvDoubleTubeData[T];
rvDoubleTubeAnnular[T,0,0];
rvDoubleTubeAnnular[T,0,1];
rvDoubleTubeAnnular[T,0,4];
rvDoubleTubeAnnular[T,4,4];
rvDoubleTubeBlockData[T,8];
rvDoubleTubeBasisData[T,0];
let x=rvDoubleTubeMatrixUnit[T,8,0,1];
let y=rvDoubleTubeMatrixUnit[T,8,1,0];
let xy=rvDoubleTubeProduct(x,y);
let yx=rvDoubleTubeProduct(y,x);
assert xy==rvDoubleTubeMatrixUnit[T,8,0,0];
assert yx==rvDoubleTubeMatrixUnit[T,8,1,1];
rvDoubleTubeDecompose(xy);
rvDoubleTubeDecompose(yx);
rvDoubleTubeStack[T,36,37];
rvDoubleTubeStack[T,5,8];

rvdoubletube U[QQ];
rvDoubleTubeData[U];
rvDoubleTubeAnnular[U,4,4];
rvDoubleTubeStack[U,33,34];
rvDoubleTubeStack[U,34,33];
let a=rvDoubleTubeBasis[U,33];
let b=rvDoubleTubeBasis[U,34];
assert rvDoubleTubeProduct(a,b)==rvDoubleTubeScale(neg(rvCycloRat[1](rvCyclo[4,0])),rvDoubleTubeProduct(b,a));

// The pointed cyclic-three tube retains the carry cocycle phase.
rvgroup G3[3];
rvdwcoeff P3[G3,3];
rvdouble D3[P3](rvPostnikovCochain[P3,3]{0,1,1,1,0,2,2,2});
rvdoubleboundary Triv[D3]{0};
rvdoubledefects C3[Triv,Triv];
rvdoubletube W[C3];
rvDoubleTubeData[W];
rvDoubleTubeStack[W,5,4];
let q=rvDoubleTubeBasis[W,4];
show rvDoubleTubeCoefficient[3](rvDoubleTubeProduct(rvDoubleTubeProduct(q,q),q));
