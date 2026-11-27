// The square of the eventual standard symmetric-group representation.
let standard = rvStableSchurCharacter[[1]];
let square = rvStableSchurKronecker(standard,standard);
show square;
rvStableSchurCharacterCoeff[[2]](square);

// Straightening recovers the smaller symmetric groups from the same square.
show rvStableSchurSpecialize[2](square);
let at3 = rvStableSchurSpecialize[3](square);
show at3;
show rvStableSchurSpecialize[4](square);
let standard3 = rvSchurCharacter[[2,1]];
assert at3 == rvSchurCharacterKronecker(standard3,standard3);
let higher = rvStableSchurCharacter[[2,1]];
rvStableSchurCharacterCoeff[[2,1]](rvStableSchurKronecker(higher,higher));
