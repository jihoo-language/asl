space Two { a, b };
let first = zero[MapExt(Two, One, 1, 0)];
let second = zero[MapExt(One, One, 0, 0)];
let third = zero[MapExt(One, One, 1, 0)];
let fourth = zero[MapExt(One, Two, 0, 0)];

// The only contribution factors through the middle copy of One.
// On H_0 this gives exactly the two-by-two matrices of rank at most one.
let bracket = eBracket(fourth, third, second, first);
let collapse = emap(\x : Two => Two.a);
assert bracket(collapse) == S.hi;
assert bracket(eid[Two]) == S.lo;
assert eInBracket(fourth, third, second, first, collapse) == S.hi;
card bracket;
