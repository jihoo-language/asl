rvunram K[5,1]{0,1};
let z = rvUnramElem[K]{};
let o = rvUnramElem[K]{1};
let two = rvUnramElem[K]{2};
let four = rvUnramElem[K]{4};
rvhyperelliptic E[K]{o,o,z,o};
rvhyperelliptic F[K]{o,z,o,o};
// The branch match x ↦ 1/x gives C: Y² = t⁶ + t² + 1.
rvellglue C[E,F](z,o,o,z);

let left = rvEllGlueLeft[C](two,two);
let right = rvEllGlueRight[C](two,two);
assert left == rvJacPoint[E](four,two);
assert right == rvJacPoint[F](four,four);
assert rvEllGlueLeft[C](neg(two),two) == left;
assert rvEllGlueRight[C](neg(two),two) == neg(right);
assert rvEllGlueRight[C](z,o) == zero[RVJac(F)];
show left;
show right;
rvHWWWeil[C];
rvHWWWeil[E];
rvHWWWeil[F];
