// Genus two with y² = x(x−1)(x−2)(x−3)(x−4) over F7.
rvunram F[7,1]{0,1};
let z = rvUnramElem[F]{};
let one = rvUnramElem[F]{1};
let minusOne = rvUnramElem[F]{-1};
rvhyperelliptic C[F]{z,rvUnramElem[F]{3},minusOne,z,rvUnramElem[F]{4},one};
let t0 = rvJacPoint[C](z,z);
let t1 = rvJacPoint[C](one,z);
let t2 = rvJacPoint[C](rvUnramElem[F]{2},z);
let t01 = add(t0,t1);
let t12 = add(t1,t2);
let p = rvJacPoint[C](rvUnramElem[F]{5},one);

// The first pair shares infinity and the second shares the point at x=1.
assert rvJacWeil2(t0,t1) == minusOne;
assert rvJacWeil2(t01,t12) == minusOne;
assert rvJacWeil2(t01,t01) == one;
assert rvJacWeil2(t01,t2) == rvUnramMul(rvJacWeil2(t0,t2),rvJacWeil2(t1,t2));

// The point at x=0 uses the leading unit at a shared branch point.
assert rvJacTate2(t0,t0) == minusOne;
assert rvJacTate2(t0,p) == minusOne;
assert rvJacTate2(t0,add(t0,p)) == rvUnramMul(rvJacTate2(t0,t0),rvJacTate2(t0,p));
assert rvJacTate2(t0,add(p,p)) == one;
assert rvJacTate2(t01,p) == rvUnramMul(rvJacTate2(t0,p),rvJacTate2(t1,p));

show rvJacWeil2(t01,t12);
show rvJacTate2(t0,t0);
show rvJacTate2(t0,add(t0,p));
