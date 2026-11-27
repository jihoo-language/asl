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
let product=rvDoubleRelative[PP](y0,x4);
rvDoubleRelativeData(product);
rvDoubleRelativeHom[4](product);
rvDoubleRelativeInclusion[4,0,0,0](product);
rvDoubleRelativeInclusion[4,1,0,1](product);
rvDoubleRelativeBasis[4,0](product);
rvDoubleRelativeAction[0,1,4,0](product);
rvDoubleRelativeInducedAction[0,1,0,0](y0);
let y4=rvDoubleDefectBasis[Y,1];
rvDoubleRelativeAssociator[Y](y0,x4,y4);
rvDoubleRelativeAssocImage[Y,0,3](y0,x4,y4);
let left=rvDoubleRelative[Y](product,y4);
let right=rvDoubleRelative[Y](y0,rvDoubleRelative[QQ](x4,y4));
rvDoubleRelativeData(left);
rvDoubleRelativeData(right);

// The original cocycle is visible on the point boundary at grades 1,2,4.
let a=rvDoubleDefectBasis[QQ,1];
let b=rvDoubleDefectBasis[QQ,2];
let c=rvDoubleDefectBasis[QQ,4];
rvDoubleRelativeAssociator[QQ](a,b,c);
rvDoubleRelativeAssocImage[QQ,7,0](a,b,c);
