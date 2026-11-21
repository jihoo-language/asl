let x = rvmon[[0, 1, 2]]{1, 0, 0};
let y = rvmon[[0, 1, 2]]{0, 1, 0};
let z = rvmon[[0, 1, 2]]{0, 0, 1};
let x2 = rvMul(x, x);
let y2 = rvMul(y, y);

// Remove the part of this homogeneous ideal supported only at the origin.
rvideal Thick[[0, 1, 2]] { x2, rvMul(x, y), rvMul(x, z), y2 };
rvprojective Projective[Thick];
rvcolon ColonX[Thick](x);
rvsaturate ChartZ[Thick](z);
rvgroebner[Projective];
rvgroebner[ColonX];
rvgroebner[ChartZ];
assert rvInIdeal[Thick](x) == S.lo;
assert rvInIdeal[Projective](x) == S.hi;
assert rvInIdeal[Projective](y) == S.lo;
assert rvInIdeal[Projective](y2) == S.hi;
rvhilbert[Thick, 1];
rvhilbert[Projective, 1];
rvhilbert[Projective, 1000];
let certificate = add(rvMul(x, rvMul(z, z)), rvMul(y2, z));
assert rvWitness[Projective, 0](certificate) == rvMul(z, z);
assert rvWitness[Projective, 1](certificate) == z;

rvideal Z[[0, 1, 2]] { z };
rvideal XY[[0, 1, 2]] { x, y };
rvintersect Corner[Z, XY];
rvgroebner[Corner];
assert rvInIdeal[Corner](rvMul(x, z)) == S.hi;
assert rvInIdeal[Corner](z) == S.lo;

// The Jacobian ideal of the cuspidal cubic cuts out a nonreduced
// singular scheme of length four at the single projective point [0:0:1].
let cubic = add(rvMul(y, y2), rvMul(x2, z));
rvideal Jacobian[[0, 1, 2]] {
    cubic, rvHasse[1, 0, 0](cubic), rvHasse[0, 1, 0](cubic), rvHasse[0, 0, 1](cubic)
};
rvprojective Singular[Jacobian];
rvgroebner[Singular];
rvhilbert[Singular, 1000];
assert rvInIdeal[Singular](x) == S.lo;
assert rvInIdeal[Singular](y) == S.lo;
assert rvInIdeal[Singular](x2) == S.hi;
assert rvInIdeal[Singular](y2) == S.hi;
