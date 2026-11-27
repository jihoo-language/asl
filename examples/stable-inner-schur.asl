// The permutation character is the sum of the trivial and standard tails.
let unit = rvStableSchurCharacter[[]];
let standard = rvStableSchurCharacter[[1]];
let v = add(unit,standard);
let symmetric = rvStableSchurSymmetric[3](v);
let mixed = rvStableInnerSchur[[2,1]](v);
let exterior = rvStableSchurExterior[3](v);
show symmetric;
show mixed;
show exterior;

// On a permutation with cycle lengths 3,2,1, its square fixes three letters.
rvSchurCharacterValue[[3,2,1]](rvStableSchurSpecialize[6](rvStableSchurAdams[2](v)));
rvSchurCharacterValue[[3,2,1]](rvStableSchurSpecialize[6](mixed));
