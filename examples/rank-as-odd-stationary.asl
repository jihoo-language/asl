rvasodd First[3] {(13,1),(11,2),(10,1),(1,1)};
rvasodd Second[3] {(13,1),(5,1),(2,1)};

// The full primitive character factors agree, coefficient by coefficient.
assert rvASOddCoeff[First,1,0] == rvASOddCoeff[Second,1,0];
assert rvASOddCoeff[First,1,1] == rvASOddCoeff[Second,1,1];
assert rvASOddCoeff[First,1,2] == rvASOddCoeff[Second,1,2];
assert rvASOddCoeff[First,1,3] == rvASOddCoeff[Second,1,3];
assert rvASOddCoeff[First,1,4] == rvASOddCoeff[Second,1,4];
assert rvASOddCoeff[First,1,5] == rvASOddCoeff[Second,1,5];
assert rvASOddCoeff[First,1,6] == rvASOddCoeff[Second,1,6];
assert rvASOddCoeff[First,1,7] == rvASOddCoeff[Second,1,7];
assert rvASOddCoeff[First,1,8] == rvASOddCoeff[Second,1,8];
assert rvASOddCoeff[First,1,9] == rvASOddCoeff[Second,1,9];
assert rvASOddCoeff[First,1,10] == rvASOddCoeff[Second,1,10];
assert rvASOddCoeff[First,1,11] == rvASOddCoeff[Second,1,11];
assert rvASOddCoeff[First,1,12] == rvASOddCoeff[Second,1,12];

rvASOddCharacterL[First,1];
rvASOddLPolynomial[First];
rvASOddLPolynomial[Second];
rvASOddPoints[First,12];
rvASOddPoints[Second,12];
rvASOddaNumber[First];
rvASOddaNumber[Second];
rvASOddEO[First];
rvASOddEO[Second];
