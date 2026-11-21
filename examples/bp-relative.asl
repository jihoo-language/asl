// Raise the height from one to two by retaining the second-level coordinates.
let t2 = bpt[2,2,1,1,2];
let v2 = bpv[2,2,0,2];
let q = bpRel[1](t2);
assert bpRelLift(q) == t2;
assert bpD(q) == zero[BPRel(2,1,2,2)];

let v = bpRel[1](v2);
assert bpD(v) == bpScale[2](q);
assert bpNull[6](q) == v;
assert bpRel[1](bpD(v2)) == bpD(v);

bprelext[2,1,2,1,6];
bporder[6](q);
show bpD(v);
