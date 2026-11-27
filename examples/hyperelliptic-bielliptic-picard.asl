rvunram F[7,1]{0,1};
let z = rvUnramElem[F]{0};
let o = rvUnramElem[F]{1};
let two = rvUnramElem[F]{2};
let three = rvUnramElem[F]{3};
let four = rvUnramElem[F]{4};
let six = rvUnramElem[F]{6};

// Use the branch point x=2 as the base point of the Jacobian.
rvhyperelliptic C[F]{six,z,four,z,six,z,o};
rvbielliptic C[Plus,Minus](o,z,z,o);
rvjacchart S[C](two);
rvJacCount[S];
rvJacCount[Plus];
rvJacCount[Minus];

let pushPlus = \d : RVJac(S) => rvBiellipticPush[Plus](d);
let pushMinus = \d : RVJac(S) => rvBiellipticPush[Minus](d);
let pullPlus = \p : RVJac(Plus) => rvBiellipticPull[S](p);
let pullMinus = \q : RVJac(Minus) => rvBiellipticPull[S](q);

assert (\p : RVJac(Plus) => pushPlus(pullPlus(p)))
    == (\p : RVJac(Plus) => rvJacScale[2](p));
assert (\q : RVJac(Minus) => pushMinus(pullMinus(q)))
    == (\q : RVJac(Minus) => rvJacScale[2](q));
assert (\p : RVJac(Plus) => pushMinus(pullPlus(p)))
    == (\p : RVJac(Plus) => zero[RVJac(Minus)]);
assert (\q : RVJac(Minus) => pushPlus(pullMinus(q)))
    == (\q : RVJac(Minus) => zero[RVJac(Plus)]);
assert (\d : RVJac(S) => add(pullPlus(pushPlus(d)),pullMinus(pushMinus(d))))
    == (\d : RVJac(S) => rvJacScale[2](d));

// One point above the old infinity is a pole for the plus quotient.
let pole = rvJacPoint[S](z,o);
assert pushPlus(pole) == rvJacPoint[Plus](four,z);
assert rvJacU[2](rvJacScale[2](pole)) == o;
assert pushPlus(rvJacScale[2](pole)) == zero[RVJac(Plus)];
assert pushPlus(add(pole,rvJacPoint[S](o,two))) == rvBiellipticPoint[Plus](three,two);
