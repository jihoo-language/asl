// The function y−1 has divisor 3[(0,1)−∞] on this elliptic curve.
rvunram F[7,1]{0,1};
let oneF = rvUnramElem[F]{1};
let zeroF = rvUnramElem[F]{};
rvhyperelliptic E[F]{oneF,zeroF,zeroF,oneF};
let t = rvJacPoint[E](zeroF,oneF);
let p = rvJacPoint[E](oneF,rvUnramElem[F]{3});
let third = rvJacTate[3](t,p);
assert third == rvUnramElem[F]{4};
assert rvJacTate[3](t,t) == rvUnramElem[F]{2};
assert rvJacTate[6](t,p) == third;
show third;

// On the genus-two curve y² = x⁵+(x+1)², the function y−x−1 has
// divisor 5[(0,1)−∞]. Its leading unit at (0,1) is 1/2 in parameter x.
rvunram W[11,1]{0,1};
let oneW = rvUnramElem[W]{1};
let zeroW = rvUnramElem[W]{};
rvhyperelliptic C[W]{oneW,rvUnramElem[W]{2},oneW,zeroW,zeroW,oneW};
let r = rvJacPoint[C](zeroW,oneW);
let q = rvJacPoint[C](oneW,rvUnramElem[W]{4});
let fifth = rvJacTate[5](r,q);
let shared = rvJacTate[5](r,r);
assert fifth == rvUnramElem[W]{4};
assert shared == rvUnramElem[W]{3};
assert rvJacTate[5](r,add(r,r)) == rvUnramElem[W]{9};
assert rvJacTate[5](r,add(r,q)) == rvUnramMul(shared,fifth);
assert rvJacTate[5](add(r,r),q) == rvUnramMul(fifth,fifth);
show fifth;
show shared;
