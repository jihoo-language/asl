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

// The same bulk modular action is recovered from the plane and point boundaries.
rvdoubletube T[PP];rvdoubletube U[QQ];
rvDoubleTubeModular[T];rvDoubleTubeModular[U];
assert rvDoubleTubeSEntry[T,0,14]==rvDoubleTubeSEntry[U,0,14];
assert rvDoubleTubeSEntry[T,0,14]==rvCycloRat[4](rvCyclo[4,0]);
let id=rvDoubleTubeCentral[T,14];
assert rvDoubleTubeWord{S,S,S,S}(id)==id;
let r=rvDoubleTubeRotate[T,4,2];
rvDoubleTubeRotationData(r);
rvDoubleTubeRotationData(rvDoubleTubeRotate[T,4,4]);

// Add the nontrivial cyclic term on the third coordinate.
let exterior=rvPostnikovCochain[P,3]{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1,
};
rvdouble J[P](add(k,exterior));
rvdoubleboundary Shifted[J]{0,1,2,3};
rvdoubledefects SS[Shifted,Shifted];
rvdoubletube V[SS];
rvDoubleTubeRotationData(rvDoubleTubeRotate[V,4,2]);
rvDoubleTubeRotationData(rvDoubleTubeRotate[V,4,4]);
rvDoubleTubeIndicatorPeriod[V,4];

// In the pointed twisted cyclic-three boundary the third indicator is complex.
rvgroup G3[3];rvdwcoeff P3[G3,3];
rvdouble D3[P3](rvPostnikovCochain[P3,3]{0,1,1,1,0,2,2,2});
rvdoubleboundary Triv[D3]{0};rvdoubledefects C3[Triv,Triv];
rvdoubletube W[C3];
rvDoubleTubeModular[W];
let rot=rvDoubleTubeRotate[W,1,3];
show rvDoubleTubeRotationEntry[0,0](rot);
show rvDoubleTubeLens[W,3];
