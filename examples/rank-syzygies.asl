let one = rvmon[[0, 1, 2]]{0, 0, 0};
let x = rvmon[[0, 1, 2]]{1, 0, 0};
let y = rvmon[[0, 1, 2]]{0, 1, 0};
let z = rvmon[[0, 1, 2]]{0, 0, 1};
rvideal Quadrics[[0, 1, 2]] {
    add(rvMul(x, y), rvMul(z, z)),
    add(rvMul(x, x), rvMul(y, z))
};

rvbetti[Quadrics, 0, 0];
rvbetti[Quadrics, 1, 2];
rvbetti[Quadrics, 2, 4];
let u = add(rvkterm[Quadrics, [0]](y), rvkterm[Quadrics, [2]](z));
let v = add(rvkterm[Quadrics, [0]](x), rvkterm[Quadrics, [1]](z));
assert rvkterm[Quadrics, []](rvMul(x, x)) == rvkterm[Quadrics, []](rvMul(y, z));
assert rvkD(u) == zero[RVKoszul(Quadrics, 0, 2)];
assert rvkD(v) == zero[RVKoszul(Quadrics, 0, 2)];
let uv = rvTorProd(rvTorClass(u), rvTorClass(v));
assert uv == rvtorgen[Quadrics, 2, 4, 0];
assert rvTorClass(rvTorRep(uv)) == uv;
show uv;

let edge = rvkterm[Quadrics, [0, 1]](one);
let boundary = rvkD(edge);
assert rvkD(boundary) == zero[RVKoszul(Quadrics, 0, 2)];
assert rvTorClass(add(u, boundary)) == rvTorClass(u);
assert rvkProd(u, boundary) == rvkD(rvkProd(u, edge));
assert rvTorClass(rvkProd(u, boundary)) == zero[RVTor(Quadrics, 2, 4)];
show rvkProd(u, boundary);

rvideal DoublePoint[[0, 1]] {
    rvmon[[0, 1]]{2, 0}, rvmon[[0, 1]]{1, 1}, rvmon[[0, 1]]{0, 2}
};
rvbetti[DoublePoint, 0, 0];
rvbetti[DoublePoint, 1, 2];
rvbetti[DoublePoint, 2, 3];
let a = rvtorgen[DoublePoint, 1, 2, 0];
let b = rvtorgen[DoublePoint, 1, 2, 1];
show rvTorProd(a, b);
