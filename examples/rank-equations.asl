// The active block of x0*Q0 + x1*Q1 + x2*Q2 on Cusp is
// [x1  0 x2; x0 x1  0;  0 x0 x1].
qmodule Cusp[0, 1, 2] {
    b1:0, b2:2, b3:4, a1:3, a2:5, a3:7
} {
    1(a1)=>[b1], 1(a2)=>[b2], 1(a3)=>[b3],
    0(a1)=>[b2], 0(a2)=>[b3], 2(a3)=>[b1]
};

let x0 = rvmon[[0, 1, 2]]{1, 0, 0};
let x1 = rvmon[[0, 1, 2]]{0, 1, 0};
let x2 = rvmon[[0, 1, 2]]{0, 0, 1};
let determinant = rvMinor[Cusp, [0, 1, 2], [(0,0), (2,0), (4,0)], [(3,0), (5,0), (7,0)]];
assert determinant == add(rvMul(rvMul(x1, x1), x1), rvMul(rvMul(x0, x0), x2));
show determinant;
rvminors[Cusp, [0, 1, 2], 3];
rvminors[Cusp, [0, 1, 2], 2];
assert rvInMinors[Cusp, [0, 1, 2], 3](determinant) == S.hi;
assert rvInMinors[Cusp, [0, 1, 2], 3](rvMul(rvMul(x1, x1), x1)) == S.lo;

// The first derivatives give x1^2 and x0^2, so the unique
// projective singular point of the cubic is [0:0:1].
assert rvHasse[1, 0, 0](determinant) == zero[RVPoly([0, 1, 2], 2)];
assert rvHasse[0, 1, 0](determinant) == rvMul(x1, x1);
assert rvHasse[0, 0, 1](determinant) == rvMul(x0, x0);
assert rvHasse[2, 0, 0](determinant) == x2;
show rvHasse[2, 0, 0](determinant);

// The F4 point [1:2:1] and F8 point [1:2:3] lie on the cubic.
assert rvVanish[7, [1, 2, 1]](determinant) == S.hi;
assert rvVanish[7, [1, 2, 2]](determinant) == S.lo;
assert rvVanish[11, [1, 2, 3]](determinant) == S.hi;
assert rvVanish[11, [1, 2, 1]](determinant) == S.lo;
rvdim[One, Cusp, [0, 1, 2], 7, [1, 2, 1]];
rvdim[One, Cusp, [0, 1, 2], 11, [1, 2, 3]];
