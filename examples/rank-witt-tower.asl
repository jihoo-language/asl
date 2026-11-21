let x = rvmon[[0]]{1};
let zero_coordinate = zero[RVPoly([0],0)];
rvwitttower CubicLength[[0]] {x,zero_coordinate,zero_coordinate};
rvWittTowerBreaks[CubicLength];
rvWittTowerGenus[CubicLength];
rvWittTowerCharacterL[CubicLength,3];
rvWittTowerLPolynomial[CubicLength];
rvWittTowerNewton[CubicLength];
rvWittTowerPoints[CubicLength,8];

// The genus-35 numerator needs extension degrees no larger than seven.
rvwitttower QuarticLength[[0]] {x,zero_coordinate,zero_coordinate,zero_coordinate};
rvWittTowerBreaks[QuarticLength];
rvWittTowerGenus[QuarticLength];
rvWittTowerLPolynomial[QuarticLength];
rvWittTowerNewton[QuarticLength];
rvWittTowerPoints[QuarticLength,1];
rvWittTowerPoints[QuarticLength,2];
rvWittTowerPoints[QuarticLength,3];
rvWittTowerPoints[QuarticLength,4];
rvWittTowerPoints[QuarticLength,5];
rvWittTowerPoints[QuarticLength,8];
