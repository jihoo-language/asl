// In each basis triple (a,b,c), the differential is x^a y^b w^c dx.
let x = rvmon[[0]]{1};
let x3 = rvmon[[0]]{3};
let empty = zero[RVPoly([0], 0)];
rvwittcover Elliptic[[0]] { x, empty };
rvwittcover GenusTwo[[0]] { x, x3 };
rvwittcover GenusSix[[0]] { x3, empty };

rvWittGaps[Elliptic];
assert rvWittCartier(rvWittOmega[Elliptic,3,0]) == zero[WittDiff(Elliptic,3)];
rvWittBasis[GenusTwo];
rvWittGaps[GenusSix];
rvWittBasis[GenusSix];
rvWittCartierRank[GenusSix, 0];
rvWittCartierRank[GenusSix, 1];
rvWittCartierRank[GenusSix, 2];
rvWittCartierRank[GenusSix, 3];
rvWittCartierRank[GenusSix, 4];
rvWittCartierRank[GenusSix, 5];
rvWittaNumber[GenusSix];

// For GenusTwo the pole from x^3 must be removed by w = z + x*y.
assert rvWittCartier(rvWittOmega[GenusTwo,3,1]) == rvWittOmega[GenusTwo,3,0];

// For GenusSix the chain is w -> x*y -> y -> x -> 1 -> 0.
let apply = \v : WittDiff(GenusSix, 3) => rvWittCartier(v);
assert apply(rvWittOmega[GenusSix,3,4]) == rvWittOmega[GenusSix,3,5];
assert apply(rvWittOmega[GenusSix,3,5]) == rvWittOmega[GenusSix,3,2];
assert apply(rvWittOmega[GenusSix,3,2]) == rvWittOmega[GenusSix,3,1];
assert apply(rvWittOmega[GenusSix,3,1]) == rvWittOmega[GenusSix,3,0];
assert apply(rvWittOmega[GenusSix,3,0]) == zero[WittDiff(GenusSix,3)];
assert apply(rvWittOmega[GenusSix,3,3]) == zero[WittDiff(GenusSix,3)];

// In F8 with modulus x^3+x+1, the square root of the coefficient 2 is 6.
assert rvWittCartier(rvWittDiffScale[2](rvWittOmega[GenusSix,11,4])) ==
          rvWittDiffScale[6](rvWittOmega[GenusSix,11,5]);
show rvWittCartier(rvWittDiffScale[2](rvWittOmega[GenusSix,11,4]));
