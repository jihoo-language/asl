// Brown–Peterson coefficients and the first two Hazewinkel coordinates.
let v1 = bpv[2, 2, 0, 1];
let v2 = bpv[2, 2, 0, 2];
let t1 = bpt[2, 2, 1, 1, 1];
let t2 = bpt[2, 2, 1, 1, 2];

show bpCoface[0](v1);
show bpCoface[0](v2);
show bpCoface[1](t2);
show bpD(t2);
assert bpD(t1) == zero[BP(2, 2, 2)];
assert bpD(bpD(v2)) == zero[BP(2, 2, 2)];
assert bpD(bpD(t2)) == zero[BP(2, 2, 3)];
