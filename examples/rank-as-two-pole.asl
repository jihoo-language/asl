rvastwo Kloosterman[5] {(1,1),(-1,1)};
rvASTwoGenus[Kloosterman];
rvASTwoPrank[Kloosterman];
rvASTwoCharacterL[Kloosterman,1];
rvASTwoLPolynomial[Kloosterman];
rvASTwoNewton[Kloosterman];
rvASTwoPoints[Kloosterman,1];
rvASTwoPoints[Kloosterman,3];

rvastwo Wild[5] {(4,1),(-2,1)};
rvastwo Swapped[5] {(-4,1),(2,1)};
rvastwo Rescaled[5] {(4,1),(-2,-1)};
// x -> 1/x exchanges the poles, and x -> 2x gives Rescaled.
assert rvASTwoCoeff[Wild,1,3] == rvASTwoCoeff[Swapped,1,3];
assert rvASTwoCoeff[Wild,1,3] == rvASTwoCoeff[Rescaled,1,3];
rvASTwoGenus[Wild];
rvASTwoPrank[Wild];
rvASTwoCharacterL[Wild,1];
rvASTwoCharacterNewton[Wild,1];
rvASTwoNewton[Wild];
rvASTwoPoints[Wild,6];

// All three genus-fourteen covers have pole orders 1 and 6.
rvastwo Pure[5] {(6,1),(-1,1)};
rvastwo Cubic[5] {(6,1),(3,1),(-1,1)};
rvastwo Quartic[5] {(6,1),(4,1),(-1,1)};
rvASTwoGenus[Pure];
rvASTwoGenus[Cubic];
rvASTwoGenus[Quartic];
rvASTwoPrank[Pure];
rvASTwoPrank[Cubic];
rvASTwoPrank[Quartic];
rvASTwoCharacterNewton[Pure,1];
rvASTwoCharacterNewton[Cubic,1];
rvASTwoCharacterNewton[Quartic,1];
rvASTwoPoints[Pure,7];
rvASTwoPoints[Cubic,7];
rvASTwoPoints[Quartic,7];
