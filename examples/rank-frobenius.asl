// The Ballico–Hefez parametrization with h=4.
rvcurve BH[[0,1]] {
    rvmon[[0,1]]{5,0},
    add(rvmon[[0,1]]{4,1},rvmon[[0,1]]{1,4}),
    rvmon[[0,1]]{0,5}
};
rvorders[BH];
rvFrobOrders[BH,1];
rvFrobOrders[BH,2];
rvFrobOrders[BH,4];

// Over F4 the generic Frobenius tangent skips order one.
// Its ordinary tangent contains the fourth-power point identically.
let x4 = rvmon[[0,1]]{20,0};
let y4 = add(rvmon[[0,1]]{16,4},rvmon[[0,1]]{4,16});
let z4 = rvmon[[0,1]]{0,20};
let tangentX = rvOscMinor[BH,1,[1,2]];
let tangentY = rvOscMinor[BH,1,[0,2]];
let tangentZ = rvOscMinor[BH,1,[0,1]];
assert add(add(rvMul(tangentX,x4),rvMul(tangentY,y4)),rvMul(tangentZ,z4)) == zero[RVPoly([0,1],28)];

let w = rvFrobWronskian[BH,2];
let f4 = add(rvmon[[0,1]]{1,4},rvmon[[0,1]]{4,1});
let f16 = add(rvmon[[0,1]]{1,16},rvmon[[0,1]]{16,1});
assert w == rvMul(f4,f16);
show w;

// Every F4 point, including infinity, is a double zero of w.
rvOrdAt[7,[1,0]](w);
rvOrdAt[7,[1,1]](w);
rvOrdAt[7,[1,2]](w);
rvOrdAt[7,[1,3]](w);
rvOrdAt[7,[0,1]](w);
// In F16 with modulus X^4+X+1, the point [1:2] is outside F4.
rvOrdAt[19,[1,2]](w);
