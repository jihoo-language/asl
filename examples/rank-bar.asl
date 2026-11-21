rvideal Squares[[0, 1, 2]] {
    rvmon[[0, 1, 2]]{2, 0, 0},
    rvmon[[0, 1, 2]]{0, 2, 0},
    rvmon[[0, 1, 2]]{0, 0, 2}
};
rvbarpath P[Squares] { [2, 0, 0], [0, 2, 0], [0, 0, 2] };
let a = rvmtorgen[Squares, 1, [2, 0, 0], 0];
let b = rvmtorgen[Squares, 1, [0, 2, 0], 0];
let c = rvmtorgen[Squares, 1, [0, 0, 2], 0];
let x = rvBarWord[P, [0, 1, 2, 3]](a, b, c);
let left = rvBarWord[P, [0, 2, 3]](rvTorProd(a,b), c);
let right = rvBarWord[P, [0, 1, 3]](a, rvTorProd(b,c));
let combined = rvBarWord[P, [0, 3]](rvTorProd(rvTorProd(a,b),c));
assert rvBarD(x) == add(left,right);
assert rvBarD(left) == combined;
assert rvBarD(right) == combined;
assert rvBarD(rvBarD(x)) == zero[RVBar(P, 0, 3, 4)];
assert rvBarLength[2](rvBarD(x)) == rvBarD(x);
assert rvBarLength[3](rvBarD(x)) == zero[RVBar(P, 0, 3, 5)];
rvbardim[P, 0, 3, 6];
rvbardim[P, 0, 3, 5, 2];
rvbardim[P, 0, 3, 4];
show rvBarD(x);

// The binary part induces the first page differential.
let e0 = rvBarClass[0, 3](x);
assert rvBarD(e0) == zero[RVBarE(P, 0, 3, 0, 5, 3)];
let e1 = rvBarNext(e0);
assert rvBarPrev(e1) == e0;
assert rvBarRep(e1) == x;
let output = rvBarD(e1);
assert output == rvBarClass[1, 2](add(left,right));
assert rvBarNext(output) == zero[RVBarE(P, 0, 3, 2, 5, 2)];
