// On y²=x³+2 over F7, y−3 and y−3x+1 have divisors 3T and 3U.
rvunram F[7,1]{0,1};
let z = rvUnramElem[F]{};
let one = rvUnramElem[F]{1};
let three = rvUnramElem[F]{3};
rvhyperelliptic E[F]{rvUnramElem[F]{2},z,z,one};
let t = rvJacPoint[E](z,three);
let u = rvJacPoint[E](three,one);

// The affine values are 5 and 4. The shared infinity contributes a minus sign.
let pairing = rvJacWeil[3](t,u);
assert pairing == rvUnramElem[F]{4};
assert rvJacWeil[3](u,t) == rvUnramElem[F]{2};
assert rvJacWeil[3](t,t) == one;
assert rvJacWeil[3](t,neg(t)) == one;
assert rvJacWeil[3](add(t,u),u) == pairing;
assert rvJacWeil[3](t,add(t,u)) == pairing;
assert rvJacWeil[6](t,u) == rvUnramMul(pairing,pairing);
show pairing;

// On this genus-two curve, f=x⁵+A²=(x−1)⁵+B² with A=2+9x² and B=4+9x+8x².
// The functions y−A and y−B give the Weil value −(B(1)−A(1))/(A(0)−B(0))=5.
rvunram W[11,1]{0,1};
let zeroW = rvUnramElem[W]{};
let oneW = rvUnramElem[W]{1};
let fourW = rvUnramElem[W]{4};
rvhyperelliptic C[W]{fourW,zeroW,rvUnramElem[W]{3},zeroW,fourW,oneW};
let p = rvJacPoint[C](zeroW,rvUnramElem[W]{2});
let q = rvJacPoint[C](oneW,rvUnramElem[W]{10});
let fifth = rvJacWeil[5](p,q);
assert fifth == rvUnramElem[W]{5};
let twice = add(p,p);
let both = add(p,q);
assert rvJacU[2](twice) == oneW;
assert rvJacWeil[5](twice,both) == rvUnramMul(fifth,fifth);
assert rvJacWeil[5](both,p) == rvUnramElem[W]{9};
assert rvJacWeil[5](twice,both) == rvJacWeil[5](both,neg(twice));
show fifth;
