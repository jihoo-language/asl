// In F3, C has the rational branch factor x and two quadratic branch factors.
rvunram F[3,1]{0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
let t = rvUnramElem[F]{2};
rvhyperelliptic C[F]{z,o,z,z,z,o};
let linear = rvJacPoint[C](z,z);
let quadratic = rvjac[C]{[2],[1],[1]}{};
let sum = add(linear,quadratic);

// The two support polynomials share x^2-x-1. That factor cancels in the sum.
rvJacKummer(linear);
rvJacKummer(quadratic);
rvJacKummer(sum);
assert rvJacKummerSupport[1](sum) == t;
assert rvJacDivisible2(sum) == S.lo;
assert rvJacDivisible2(add(sum,sum)) == S.hi;
rvJacKummer(add(sum,sum));
assert rvJacKummerSupport[0](add(sum,sum)) == o;
assert rvJacKummerSupport[1](add(sum,sum)) == z;

// Both branch orbits of M give nonsquare coordinates for (2,1).
rvhyperelliptic M[F]{o,z,z,o,t,o};
let point = rvJacPoint[M](t,o);
rvJacKummer(point);
assert rvJacDivisible2(point) == S.lo;
