// The five nonedges specify a flag Stanley-Reisner ring.
rvideal Flag[[1, 2, 3, 4, 5, 6]] {
    rvmon[[1, 2, 3, 4, 5, 6]]{1, 1, 0, 0, 0, 0},
    rvmon[[1, 2, 3, 4, 5, 6]]{0, 0, 1, 1, 0, 0},
    rvmon[[1, 2, 3, 4, 5, 6]]{0, 0, 0, 0, 1, 1},
    rvmon[[1, 2, 3, 4, 5, 6]]{1, 0, 1, 0, 0, 0},
    rvmon[[1, 2, 3, 4, 5, 6]]{0, 0, 0, 1, 1, 0}
};

rvhochster[Flag, [], -1];
rvhochster[Flag, [1, 2], 0];
rvhochster[Flag, [1, 2, 3, 4, 5, 6], 1];
rvhochster[Flag, [1, 2, 3, 4, 5, 6], 2];

// The empty face represents the augmented degree-minus-one cochain.
let empty = rvhface[Flag, [1, 2], []];
let v1 = rvhface[Flag, [1, 2], [1]];
let v2 = rvhface[Flag, [1, 2], [2]];
assert rvkD(empty) == add(v1, v2);
let a = rvTorClass(v1);
let b = rvTorClass(rvhface[Flag, [3, 4], [3]]);
let c = rvTorClass(rvhface[Flag, [5, 6], [5]]);

// The edge cochain pairs to one with the cycle 1–4–2–5–1.
let edge = rvhface[Flag, [1, 2, 3, 4, 5, 6], [1, 5]];
let observed = rvTorClass(edge);
assert observed == rvmtorgen[Flag, 4, [1, 1, 1, 1, 1, 1], 0];
assert rvkD(edge) == zero[RVMKoszul(Flag, 3, [1, 1, 1, 1, 1, 1])];
assert rvTHigher(a, b, c) == observed;
assert rvTInBracket(a, b, c, observed) == S.hi;
assert rvMForget(observed) == rvTHigher(rvMForget(a), rvMForget(b), rvMForget(c));
assert rvMForget(rvTLift(a,c,b)) == rvTLift(rvMForget(a),rvMForget(c),rvMForget(b));
show edge;
show observed;
