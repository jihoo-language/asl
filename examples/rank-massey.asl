// The five nonedges specify a flag Stanley-Reisner ring.
rvideal Flag[[1, 2, 3, 4, 5, 6]] {
    rvmon[[1, 2, 3, 4, 5, 6]]{1, 1, 0, 0, 0, 0},
    rvmon[[1, 2, 3, 4, 5, 6]]{0, 0, 1, 1, 0, 0},
    rvmon[[1, 2, 3, 4, 5, 6]]{0, 0, 0, 0, 1, 1},
    rvmon[[1, 2, 3, 4, 5, 6]]{1, 0, 1, 0, 0, 0},
    rvmon[[1, 2, 3, 4, 5, 6]]{0, 0, 0, 1, 1, 0}
};
let x1 = rvmon[[1, 2, 3, 4, 5, 6]]{1, 0, 0, 0, 0, 0};
let x3 = rvmon[[1, 2, 3, 4, 5, 6]]{0, 0, 1, 0, 0, 0};
let x5 = rvmon[[1, 2, 3, 4, 5, 6]]{0, 0, 0, 0, 1, 0};
let A = rvkterm[Flag, [2]](x1);
let B = rvkterm[Flag, [4]](x3);
let C = rvkterm[Flag, [6]](x5);
let a = rvTorClass(A);
let b = rvTorClass(B);
let c = rvTorClass(C);
assert rvkProd(A, B) == zero[RVKoszul(Flag, 2, 4)];
let bc = rvkProd(B, C);
let V = rvkterm[Flag, [3, 4, 6]](x5);
assert rvkD(V) == bc;
assert rvkD(rvkPrimitive(bc)) == bc;
let answer = rvkterm[Flag, [2, 3, 4, 6]](rvMul(x1, x5));
assert rvkProd(A, V) == answer;
let q = rvTorClass(answer);
assert q == rvtorgen[Flag, 4, 6, 0];
rvbetti[Flag, 3, 4];
rvbetti[Flag, 4, 6];
let bracket = rvTBracket(a, b, c);
card bracket;
assert bracket(q) == S.hi;
assert rvTInBracket(a, b, c, q) == S.hi;
assert rvTInBracket(a, b, c, zero[RVTor(Flag, 4, 6)]) == S.lo;
show q;
