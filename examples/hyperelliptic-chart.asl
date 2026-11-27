// Send a rational branch point to infinity and compare divisor classes.
rvunram F[5,1]{0,1};
let z=rvUnramElem[F]{};
let o=rvUnramElem[F]{1};
let t=rvUnramElem[F]{2};
let three=rvUnramElem[F]{3};
let four=rvUnramElem[F]{4};
rvhyperelliptic C[F]{o,o,z,three,o,o};
rvjacchart N[C](t);
let branch=rvJacPoint[C](t,z);
let p=rvJacPoint[C](z,o);
let infinity=rvJacPoint[N](z,z);
let mapped=rvJacPoint[N](three,three);
assert rvJacRechart[N](branch)==infinity;
assert rvJacRechart[N](p)==add(mapped,infinity);
assert rvJacRechart[N](add(p,p))==add(mapped,mapped);
assert rvJacUnchart[N](rvJacRechart[N](p))==p;

// The even-degree model has nonsquare leading coefficient and a branch at x=2.
rvhyperelliptic E[F]{four,o,z,three,z,o,t};
rvjacchart M[E](t);
rvHWWWeil[E];
rvHWWWeil[M];
rvJacCount[M];
