let h0 = mapextgen[One, One, 1, 1, 0];
let h1 = mapextgen[One, One, 1, 2, 0];
let square = ecomp(h1, h1);
let triple = eBracket(h0, h1, h0);

// This algebraic triple Massey product on Adams E2 is the singleton h1 squared.
assert square == mapextgen[One, One, 2, 4, 0];
assert triple(square) == S.hi;
assert eInBracket(h0, h1, h0, zero[MapExt(One, One, 2, 4)]) == S.lo;
card triple;
