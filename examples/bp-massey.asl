let a = bpt[3, 2, 1, 1, 1];
let v = bpv[3, 2, 1, 1];
let w = add(bpMul(v, a), bpPow[2](a));
assert bpD(w) == bpCup(a, a);

let m = bpMassey(a, a, a, w, w);
let beta = bpDelta(bpDelta(bpReduce[2](bpv[3, 2, 0, 2])));
show m;
assert bpD(m) == zero[BP(3, 2, 3)];
assert add(m, beta) == bpD(bpMul(bpPow[2](v), a));
