rvideal Squares[[0, 1, 2]] {
    rvmon[[0, 1, 2]]{2, 0, 0},
    rvmon[[0, 1, 2]]{0, 2, 0},
    rvmon[[0, 1, 2]]{0, 0, 2}
};
rvbarpath P[Squares] { [2, 0, 0], [0, 2, 0], [0, 0, 2] };
let a = rvmtorgen[Squares, 1, [2, 0, 0], 0];
let b = rvmtorgen[Squares, 1, [0, 2, 0], 0];
let c = rvmtorgen[Squares, 1, [0, 0, 2], 0];
let A = rvBCDesuspend(rvBarWord[P, [0, 1]](a));
let B = rvBCDesuspend(rvBarWord[P, [1, 2]](b));
let C = rvBCDesuspend(rvBarWord[P, [2, 3]](c));
let AB = rvBCDesuspend(rvBarWord[P, [0, 1, 2]](a, b));
let BC = rvBCDesuspend(rvBarWord[P, [1, 2, 3]](b, c));
let long = rvBCDesuspend(rvBarWord[P, [0, 1, 2, 3]](a, b, c));
let S1 = rvBCMul(A, BC);
let S2 = rvBCMul(AB, C);
let J1 = rvBCDesuspend(rvBarWord[P, [0, 2, 3]](rvTorProd(a,b), c));
let J2 = rvBCDesuspend(rvBarWord[P, [0, 1, 3]](a, rvTorProd(b,c)));

assert rvBCD(long) == add(add(S1,S2),add(J1,J2));
assert rvBCD(rvBCD(long)) == zero[RVBarCobar(P,0,3,3)];
assert rvBCD(S1) == rvBCMul(A,rvBCD(BC));
assert rvBCMul(rvBCMul(A,B),C) == rvBCMul(A,rvBCMul(B,C));
assert rvBCProject(rvBCMul(A,B)) == rvTorProd(a,b);
assert rvBCProject(A) == a;

assert rvBCHomotopy(rvBCD(long)) == long;
assert rvBCHomotopy(long) == zero[RVBarCobar(P,0,3,6)];
let mixed = add(S1,J2);
let lifted = rvBCHomotopy(mixed);
assert add(rvBCD(lifted),rvBCHomotopy(rvBCD(mixed))) == mixed;
assert rvBCProject(rvBCD(long)) == zero[RVMTor(Squares,4,[2,2,2])];
assert rvBCProject(lifted) == zero[RVMTor(Squares,5,[2,2,2])];
assert rvBCHomotopy(lifted) == zero[RVBarCobar(P,0,3,6)];
show rvBCD(long);
