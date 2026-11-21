// Degree six: the relative class detects a boundary only visible at height two.
let t2 = bpt[2,2,1,1,2];
let q = bpRel[1](t2);
let x = bpRelClass[6](q);
assert x == bprelgen[2,1,2,1,6,0];
assert bpRelClass[6](bpRelRep(x)) == x;
let delta = bpHeightDelta(x);
assert delta == bpextgen[2,1,2,6,0];
assert bpRaise[2](delta) == zero[BPExt(2,2,2,6)];
let a = bpt[2,1,2,1,1];
let b = bpt[2,1,2,2,1];
let v = bpv[2,1,2,1];
let cocycle = add(bpMul(v,bpMul(a,b)),bpNeg(bpMul(a,bpPow[2](b))));
assert bpHeightDelta(q) == cocycle;
assert bpRaise[2](cocycle) == bpD(t2);
assert bpExtClass[6](cocycle) == delta;
bporder[6](bpExtRep(delta));
show delta;

// Degree eight: the lower class maps to twice the upper generator.
let alpha = bpExtClass[8](bpBockP[3](bpPow[4](bpv[2,1,0,1])));
let v1 = bpv[2,2,0,1];
let v2 = bpv[2,2,0,2];
let g = add(bpPow[4](v1),bpMul(bpint[2,2,0,8],bpMul(v1,v2)));
let sigma = bpExtClass[8](bpBockP[4](g));
let image = bpRel[1](sigma);
assert bpRaise[2](alpha) == bpScale[2](sigma);
assert image == bprelgen[2,1,2,1,8,0];
assert bpRel[1](bpRaise[2](alpha)) == zero[BPRelExt(2,1,2,1,8)];
assert bpHeightDelta(image) == zero[BPExt(2,1,2,8)];
assert bpRelClass[8](bpRel[1](bpExtRep(sigma))) == image;
bporder[8](bpRelRep(image));
show image;
