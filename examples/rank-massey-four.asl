// The relation x4*x5 is omitted from the triangular list.
rvideal Stair[[1, 2, 3, 4, 5, 6, 7, 8]] {
    rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{1, 0, 0, 0, 1, 0, 0, 0},
    rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 1, 0, 0, 1, 0, 0, 0},
    rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 1, 0, 0, 0, 1, 0, 0},
    rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 0, 1, 0, 1, 0, 0, 0},
    rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 0, 1, 0, 0, 1, 0, 0},
    rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 0, 1, 0, 0, 0, 1, 0},
    rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 0, 0, 1, 0, 1, 0, 0},
    rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 0, 0, 1, 0, 0, 1, 0},
    rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 0, 0, 1, 0, 0, 0, 1}
};

let A = rvkterm[Stair, [5]](rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{1, 0, 0, 0, 0, 0, 0, 0});
let a = rvTorClass(A);
let B = rvkterm[Stair, [6]](rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 1, 0, 0, 0, 0, 0, 0});
let b = rvTorClass(B);
let C = rvkterm[Stair, [7]](rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 0, 1, 0, 0, 0, 0, 0});
let c = rvTorClass(C);
let D = rvkterm[Stair, [8]](rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 0, 0, 1, 0, 0, 0, 0});
let d = rvTorClass(D);
// Every proper interval is supplied by the same formula.
let U = rvkterm[Stair, [1, 5, 6]](rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 1, 0, 0, 0, 0, 0, 0});
let V = rvkterm[Stair, [2, 6, 7]](rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 0, 1, 0, 0, 0, 0, 0});
let W = rvkterm[Stair, [3, 7, 8]](rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 0, 0, 1, 0, 0, 0, 0});
let X = rvkterm[Stair, [1, 2, 5, 6, 7]](rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 0, 1, 0, 0, 0, 0, 0});
let Y = rvkterm[Stair, [2, 3, 6, 7, 8]](rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 0, 0, 1, 0, 0, 0, 0});

let curvature = rvTCurvature(A, B, C, D){
    0 < 2 => U, 1 < 3 => V, 2 < 4 => W, 0 < 3 => X, 1 < 4 => Y
};
let value = rvTMassey(A, B, C, D){
    0 < 2 => U, 1 < 3 => V, 2 < 4 => W, 0 < 3 => X, 1 < 4 => Y
};
let omitted = rvkterm[Stair, [1, 2, 3, 5, 6, 7, 8]](rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 0, 0, 1, 0, 0, 0, 0});
let survivor = rvkterm[Stair, [1, 2, 3, 6, 7, 8]](rvmon[[1, 2, 3, 4, 5, 6, 7, 8]]{0, 0, 0, 1, 1, 0, 0, 0});
assert add(curvature, survivor) == rvkD(omitted);
assert value == rvTorClass(survivor);
assert value == rvTorClass(curvature);
assert a == rvtorgen[Stair, 1, 2, 0];
assert b == rvtorgen[Stair, 1, 2, 1];
assert c == rvtorgen[Stair, 1, 2, 3];
assert d == rvtorgen[Stair, 1, 2, 6];
assert value == rvtorgen[Stair, 6, 8, 0];
show a;
show b;
show c;
show d;
rvbetti[Stair, 6, 8];
show value;
