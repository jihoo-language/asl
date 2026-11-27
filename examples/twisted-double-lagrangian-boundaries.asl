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



// The span of the first two coordinates supports a subgroup boundary.
rvdoubleboundary Plane[D]{0,1,2,3};
rvDoubleBoundaryData[Plane];
show rvDoubleCondensate[Plane];

// The full type-III subgroup has an integral obstruction witness.
rvdoubleboundary Whole[D]{0,1,2,3,4,5,6,7};
rvDoubleBoundaryData[Whole];
rvDoubleBoundaryWitness[Whole];
rvDoubleBoundaryCycle[D]{0,1,2,3,4,5,6,7}{
  (1,2,4,1),(1,4,2,-1),(2,1,4,-1),
  (2,4,1,1),(4,1,2,1),(4,2,1,-1)
};

// On the untwisted Klein four group, the torsion shift is x1*y2/2.
rvgroup K{{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}};
rvdwcoeff C[K,2];
rvdouble U[C](rvPostnikovK[C]);
rvdoubleboundary Base[U]{0,1,2,3};
rvdoubleboundaryshift Torsion[Base]{0,(1,2),(1,2),0,0,0,0,(1,2),(1,2)};
show rvDoubleCondensate[Base];
show rvDoubleCondensate[Torsion];
show rvDoubleBoundaryMu[Torsion,1,2];
rvDoubleAnnulus[Base,Base];
rvDoubleAnnulus[Torsion,Torsion];
rvDoubleAnnulus[Base,Torsion];
