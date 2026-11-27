// Full cuspidal J+ Hecke data gives the zeta of the curve X0(2p).
rvmanin M23[23];
rvManinZeta[M23,5];
rvManinZeta[M23,5,2];
let P=rvManinFrobPoly[M23,5];
rvOreCoeff[10,0](P);

rvmanin M41[41];
rvManinZeta[M41,5];
rvManinZeta[M41,5,2];
rvManinZeta[M41,5,20];

rvunram F5[5,1]{0,1};
rvunram F25[5,1]{-2,0,1};
rvunram F49[7,1]{-3,0,1};
rvunram F121[11,1]{-2,0,1};

rvManinPointCount[M23,F5];
rvManinPointCount[M23,F25];
rvManinPointCount[M23,F49];
rvManinPointCount[M23,F121];
rvManinPointCount[M41,F5];
rvManinPointCount[M41,F25];
rvManinPointCount[M41,F49];
rvManinPointCount[M41,F121];

// A scalar Frobenius and a non-scalar repeated eigenvalue on E[23].
rvManinCurveLines[M23](rvUnramElem[F25]{0},rvUnramElem[F25]{1});
rvManinCurveLines[M23](rvUnramElem[F49]{1,1},rvUnramElem[F49]{4,1});

rvmanin M5[5];
rvunram F7[7,1]{0,1};
assert rvManinFrobPoly[M5,7,2]==rvOre[0,0];
rvManinZeta[M5,7,2];
rvManinPointCount[M5,F7];

// The same trace and discriminant give different residual Frobenius actions.
rvunram F31[31,1]{0,1};
rvManinCurveLines[M5](rvUnramElem[F31]{0},rvUnramElem[F31]{11});
rvManinCurveLines[M5](rvUnramElem[F31]{1},rvUnramElem[F31]{17});
rvManinPointCount[M5,F31];
