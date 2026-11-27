// Reduce the regularized quantum equation of Gr(2,5).
let transformed = rvOreFactorial[5](rvQGrAnnihilator[2,5]);
let cofactor = rvOreEulerReduce(transformed);
assert rvOreMul(rvOreEulerFactor(transformed),cofactor) == transformed;

// B = theta^2 - q(11*theta^2 + 11*theta + 3) - q^2(theta + 1)^2.
let theta = rvOre[0,1];
let one = rvOre[0,0];
let middle = add(rvOre[0,2,11],add(rvOre[0,1,11],rvOre[0,0,3]));
let shifted = add(theta,one);
let apery = add(rvOre[0,2],neg(add(rvOreMul(rvOre[1,0],middle),
    rvOreMul(rvOre[2,0],rvOreMul(shifted,shifted)))));
assert cofactor == apery;
show cofactor;
rvQGrRegularizedRhs[2,5];
rvQGrRegularized[2,5,5];

// The Laurent polynomial is (1+x)(1+y)(1+x+y)/(xy).
// In its n-th power the constant term is the x^n y^n coefficient
// of the numerator. Expand (1+x+y)^n in (1+x) and y to get the
// same binomial sum. The recurrence has leading coefficient n^2
// in degree n>0, so A(0)=1 determines the whole period.
rvperiod Mirror[2] {
    ([-1,-1],1), ([-1,0],2), ([-1,1],1), ([0,-1],2),
    ([0,0],3), ([0,1],1), ([1,-1],1), ([1,0],1)
};
// Verify both the Laurent declaration and the WZ certificate against the
// cofactor that was computed from the quantum connection.
rvCTWZCheck[Mirror,2,1,1,2](cofactor) {
    {(0,5,1),(1,4,6),(0,4,-5),(2,3,-11),(1,3,7)},
    {(1,0,0,1),(1,1,0,1),(1,1,-1,1)}
};
rvCTExact[Mirror,0];
rvCTExact[Mirror,1];
rvCTExact[Mirror,2];
rvCTExact[Mirror,3];

// The base-three digits of 7 are 1,2 from least to most significant.
let state7 = rvCTSection[2](rvCTSection[1](rvCTInitial[Mirror,3,2]));
let mod9 = rvCTRead(state7);
assert mod9 == rvCTMod[Mirror,3,2,7];
assert mod9 == rvCTResidue[3,2,6];
show mod9;

// The 38 base-three digits of 3^38-1 are all 2. Verify that the
// entire mod-27 state is fixed after the first such digit.
let state27 = rvCTSection[2](rvCTInitial[Mirror,3,3]);
assert rvCTSection[2](state27) == state27;
let mod27 = rvCTRead(state27);
assert mod27 == rvCTMod[Mirror,3,3,1350851717672992088];
assert mod27 == rvCTResidue[3,3,19];
show mod27;
