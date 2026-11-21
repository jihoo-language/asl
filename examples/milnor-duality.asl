let xi2 = cword{[0, 1]};
let left = mul(sq[2], sq[1]);
let right = mul(sq[1], sq[2]);
let q1 = milnor{[0, 1]};
assert q1 == add(left, right);
show q1;
show milnorPair(left, xi2);
show milnorPair(right, xi2);

let q2 = milnor{[0, 0, 1]};
assert q2 == add(mul(sq[4], q1), mul(q1, sq[4]));
assert milnorPair(q2, cword{[0, 0, 1]}) == S.hi;
assert milnorPair(q2, cword{[7]}) == S.lo;
show q2;

let z = coact[0, 0](gen[One, 0, 0]);
assert z == cword{};
show milnorPair(milnor{[0]}, z);
