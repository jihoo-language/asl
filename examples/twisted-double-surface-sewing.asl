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

// Genus two includes the charged separating channels.
rvDoubleSurfaceDimensions[D,0];
rvDoubleSurfaceDimensions[D,1];
rvDoubleSurfaceDimensions[D,2];
rvDoubleSurfaceDimensions[D,3];
rvDoubleSurfaceSewing[D,1,1];
rvDoubleSurfaceCharge[D]{5,2};
rvDoubleSurfaceChain[D,Slide,0,1]{0,2,5,5};
let neutral=rvDoubleSurfaceBasis[D,0]{0,2,5,5};
let leg1=rvDoubleSurfaceBasis[D,5]{5,2};
let leg2=rvDoubleSurfaceBasis[D,5]{5,7};
let chargedcut=rvDoubleSurfaceSew(leg1,leg2);
let slid=rvDoubleSurfaceMove[Slide,0,1](neutral);
assert slid==rvDoubleSurfaceScale(neg(rvCycloRat[1](rvCyclo[4,0])),chargedcut);
assert rvDoubleSurfaceCut[1,5](slid)==slid;
assert rvDoubleSurfaceMove[SlideInv,0,1](slid)==neutral;
show slid;

// The torus Fourier transform intertwines the actual marking moves.
rvDoubleSurfaceModular[D];
rvdoubletube T[PP];
let v=rvDoubleSurfaceBasis[D,0]{4,4};
let x=rvDoubleSurfaceToTube[T](v);
assert rvDoubleSurfaceFromTube(x)==v;
assert rvDoubleSurfaceToTube[T](rvDoubleSurfaceMove[S,0](v))==rvDoubleTubeWord{S}(x);
assert rvDoubleSurfaceToTube[T](rvDoubleSurfaceMove[T,0](v))==rvDoubleTubeWord{T[1]}(x);

// A twisted cyclic-three slide carries a genuinely complex phase.
rvgroup G3[3];rvdwcoeff P3[G3,3];
rvdouble D3[P3](rvPostnikovCochain[P3,3]{0,1,1,1,0,2,2,2});
rvDoubleSurfaceDimensions[D3,2];
rvDoubleSurfaceModular[D3];
rvDoubleSurfaceChain[D3,Slide,0,1]{0,1,2,0};
let z=rvDoubleSurfaceMove[Slide,0,1](rvDoubleSurfaceBasis[D3,0]{0,1,2,0});
show rvDoubleSurfaceCoefficient{2,1,2,2}(z);
