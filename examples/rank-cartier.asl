rvascover GenusTwo[[0,1]] {rvmon[[0,1]]{0,5},rvmon[[0,1]]{5,0}};
let a = rvASOmega[GenusTwo,11,0];
let b = rvASOmega[GenusTwo,11,1];
assert rvASCartier(b) == a;
assert rvASCartier(a) == zero[ASDiff(GenusTwo,11)];
// In F8 with X^3+X+1, the square root of the coefficient 2 is 6.
let scaled = rvASCartier(rvASDiffScale[2](b));
assert scaled == rvASDiffScale[6](a);
show scaled;
rvASCartierRank[GenusTwo,0];
rvASCartierRank[GenusTwo,1];
rvASCartierRank[GenusTwo,2];
rvASaNumber[GenusTwo];

rvascover GenusFour[[0,1]] {rvmon[[0,1]]{0,9},rvmon[[0,1]]{9,0}};
rvASCartierRank[GenusFour,0];
rvASCartierRank[GenusFour,1];
rvASCartierRank[GenusFour,2];
rvASCartierRank[GenusFour,3];
rvASaNumber[GenusFour];

rvascover Ordinary[[0,1]] {
    add(rvmon[[0,1]]{2,0},rvmon[[0,1]]{0,2}),rvmon[[0,1]]{1,1}
};
let ordinary = rvASOmega[Ordinary,3,0];
assert rvASCartier(ordinary) == ordinary;
rvASCartierRank[Ordinary,1000];
rvASaNumber[Ordinary];
