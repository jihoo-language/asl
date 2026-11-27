rvunram F[11,1]{0,1};
let zeroF = rvUnramElem[F]{};
let oneF = rvUnramElem[F]{1};
rvhyperelliptic C[F]{oneF,zeroF,zeroF,zeroF,zeroF,oneF};

let point = rvJacPoint[C](zeroF,oneF);
assert rvJacDivisible2(point) == S.hi;
rvJacHalfFiber[C](zeroF,oneF);
let first = rvJacHalfAt[C,0](zeroF,oneF);
let last = rvJacHalfAt[C,15](zeroF,oneF);
assert add(first,first) == point;
assert add(last,last) == point;
show first;
show last;
