// Two quadrics whose first S-pair adds a new cubic to the basis.
let x = rvmon[[0, 1, 2]]{1, 0, 0};
let y = rvmon[[0, 1, 2]]{0, 1, 0};
let z = rvmon[[0, 1, 2]]{0, 0, 1};
let first = add(rvMul(x, y), rvMul(z, z));
let second = add(rvMul(x, x), rvMul(y, z));
rvideal Intersection[[0, 1, 2]] { first, second };
rvgroebner[Intersection];

let cubic = add(rvmon[[0, 1, 2]]{1, 0, 2}, rvmon[[0, 1, 2]]{0, 2, 1});
assert rvInIdeal[Intersection](cubic) == S.hi;
let u = rvWitness[Intersection, 0](cubic);
let v = rvWitness[Intersection, 1](cubic);
assert u == x;
assert v == y;
assert add(rvMul(u, first), rvMul(v, second)) == cubic;
assert rvNormal[Intersection](rvmon[[0, 1, 2]]{2, 1, 0}) == rvmon[[0, 1, 2]]{0, 2, 1};
rvhilbert[Intersection, 0];
rvhilbert[Intersection, 1];
rvhilbert[Intersection, 2];
rvhilbert[Intersection, 1000];

qmodule Cusp[0, 1, 2] {
    b1:0, b2:2, b3:4, a1:3, a2:5, a3:7
} {
    1(a1)=>[b1], 1(a2)=>[b2], 1(a3)=>[b3],
    0(a1)=>[b2], 0(a2)=>[b3], 2(a3)=>[b1]
};
rvfitting Cuspidal[Cusp, [0, 1, 2], 3];
rvgroebner[Cuspidal];
rvhilbert[Cuspidal, 1000];
let determinant = rvMinor[Cusp, [0, 1, 2], [(0,0), (2,0), (4,0)], [(3,0), (5,0), (7,0)]];
assert rvWitness[Cuspidal, 0](rvMul(z, determinant)) == z;
