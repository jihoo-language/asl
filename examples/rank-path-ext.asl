rvideal Squares[[0, 1, 2]] {
    rvmon[[0, 1, 2]]{2, 0, 0},
    rvmon[[0, 1, 2]]{0, 2, 0},
    rvmon[[0, 1, 2]]{0, 0, 2}
};
rvbarpath P[Squares] { [2, 0, 0], [0, 2, 0], [0, 0, 2] };
let a = rvmtorgen[Squares, 1, [2, 0, 0], 0];
let b = rvmtorgen[Squares, 1, [0, 2, 0], 0];
let c = rvmtorgen[Squares, 1, [0, 0, 2], 0];
let word = rvBarWord[P, [0, 1, 2, 3]](a, b, c);
let left = rvBarWord[P, [0, 2, 3]](rvTorProd(a, b), c);
let right = rvBarWord[P, [0, 1, 3]](a, rvTorProd(b, c));
let top = rvBarWord[P, [0, 3]](rvTorProd(rvTorProd(a, b), c));
let f = rvBarDualize(top);
let coboundary = rvBarDualD(f);
assert coboundary == add(rvBarDualize(left), rvBarDualize(right));
assert rvBarPair(coboundary,left) == rvBarPair(f,rvBarD(left));
assert rvBarPair(coboundary,right) == rvBarPair(f,rvBarD(right));
assert rvBarPair(coboundary,left) == S.hi;
assert rvBarPair(coboundary,right) == S.hi;
assert rvBarDualD(coboundary) == zero[RVBarDual(P,0,3,6)];

// The contraction crossing vertex 1 vanishes in the pairing at that cut.
let first = rvBarDualize(rvBarWord[P,[0,1]](a));
let second = rvBarDualize(rvBarWord[P,[1,3]](rvTorProd(b,c)));
let cup = rvPathCup(first,second);
assert cup == rvBarDualize(right);
assert rvBarPair(cup,left) == S.lo;
assert rvBarPair(cup,rvBarD(word)) == S.hi;
let leibniz = add(rvPathCup(rvBarDualD(first),second),rvPathCup(first,rvBarDualD(second)));
assert rvBarDualD(cup) == leibniz;
assert leibniz == rvBarDualize(word);
assert rvPathClass(leibniz) == zero[RVPathExt(P,0,3,6)];

rvpathextdim[P,0,1,2];
rvpathextdim[P,0,3,4];
rvpathextdim[P,0,3,5];
rvpathextdim[P,0,3,6];
show coboundary;
