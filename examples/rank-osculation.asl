let a = rvmon[[0,1]]{5,0};
let b = rvmon[[0,1]]{4,1};
let c = rvmon[[0,1]]{1,4};
let d = rvmon[[0,1]]{0,5};
rvcurve Twist[[0,1]] { a,b,c,d };
rvorders[Twist];

// Tangent coordinates are in lexicographic order (01,02,03,12,13,23).
let middle = rvmon[[0,1]]{4,4};
let zeroTangent = zero[RVPoly([0,1],8)];
assert rvOscMinor[Twist,1,[0,1]] == rvmon[[0,1]]{8,0};
assert rvOscMinor[Twist,1,[0,2]] == zeroTangent;
assert rvOscMinor[Twist,1,[0,3]] == middle;
assert rvOscMinor[Twist,1,[1,2]] == middle;
assert rvOscMinor[Twist,1,[1,3]] == zeroTangent;
assert rvOscMinor[Twist,1,[2,3]] == rvmon[[0,1]]{0,8};
rvosculating Tangents[Twist,1,[2,3,4,5,6,7]];
rvgroebner[Tangents];

// The osculating planes have the original four forms in lexicographic order.
assert rvOscMinor[Twist,2,[0,1,2]] == a;
assert rvOscMinor[Twist,2,[0,1,3]] == b;
assert rvOscMinor[Twist,2,[0,2,3]] == c;
assert rvOscMinor[Twist,2,[1,2,3]] == d;
rvosculating Planes[Twist,2,[8,9,10,11]];
rvgroebner[Planes];
assert rvOscMinor[Twist,3,[0,1,2,3]] == rvmon[[0,1]]{0,0};

// The gap in the orders can be arbitrarily large. The tangent image is
// the same conic after dropping the common Frobenius exponent.
rvcurve Large[[0,1]] {
    rvmon[[0,1]]{1025,0}, rvmon[[0,1]]{1024,1}, rvmon[[0,1]]{1,1024}, rvmon[[0,1]]{0,1025}
};
rvorders[Large];
assert rvOscMinor[Large,1,[0,3]] == rvmon[[0,1]]{1024,1024};
assert rvOscMinor[Large,2,[0,2,3]] == rvmon[[0,1]]{1,1024};
rvosculating LargeTangents[Large,1,[2,3,4,5,6,7]];
rvgroebner[LargeTangents];
