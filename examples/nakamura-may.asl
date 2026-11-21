let u = cword{[0, 1]};
let b = cNakamura[1](u);
let c = cNakamura[2](b);
assert cobarD(b) == cNakamura[1](cobarD(u));
assert cobarD(c) == cNakamura[2](cobarD(b));

let b3 = mayClass[3, 0](b);
assert b3 == maygen[3, 2, 6, 0, 0];
show mayD(b3);

// The raw second square meets the E_6 condition. Correct it before reading d_7.
let c6 = mayClass[6, 0](c);
assert mayLiftTo[6, 0](c) == c;
let corrected = mayLiftTo[7, 0](c);
let correction = add(c, corrected);
assert mayLiftTo[0, 1](correction) == correction;
let c7 = mayClass[7, 0](corrected);
assert c7 == mayNext(c6);
assert mayD(c7) == maygen[7, 5, 12, 7, 0];
show mayD(c7);

// The same leading term over A(1) has a closed lift.
let finite = mayLiftTo[13, 0](resA[1](c));
assert cobarD(finite) == zero[CobarAn(1, 5, 12)];
assert extClass(finite) == extangen[1, 4, 12, 0];
show extClass(finite);
