// The first-page differential agrees with the height connecting map.
let t2 = bpt[2,2,1,1,2];
let x = bpHeightClass[1,2,6](t2);
assert x == bpHeightgen[2,2,1,2,1,6,0];
let delta = bpHeightDelta(bpRel[1](t2));
assert bpHeightD(x) == bpHeightClass[1,1,6](bpRaise[2](delta));
show bpHeightD(x);

// A chosen page-one representative needs a lower-height correction.
let start = bpHeightgen[2,2,1,2,1,8,0];
let lift = bpHeightRep(start);
show bpD(lift);
let stable = bpHeightNext(start);
let cycle = bpHeightRep(stable);
assert bpD(cycle) == zero[BP(2,2,2)];
assert bpHeightNext(bpHeightPrev(stable)) == stable;
let correction = add(lift,bpNeg(cycle));
assert bpHeightClass[0,2,8](correction) == zero[BPHeight(2,2,0,2,1,8)];
assert bpD(correction) == bpD(lift);

// The two stable factors are the height grades of one cyclic group.
let v1 = bpv[2,2,0,1];
let v2 = bpv[2,2,0,2];
let g = add(bpPow[4](v1),bpMul(bpint[2,2,0,8],bpMul(v1,v2)));
let sigma = bpBockP[4](g);
let lower = add(bpMul(bpint[2,2,1,2],sigma),bpNeg(bpD(bpMul(v1,v2))));
let low = bpHeightClass[3,1,8](lower);
let high = bpHeightClass[3,2,8](sigma);
assert low == bpHeightgen[2,2,3,1,1,8,0];
assert high == bpHeightgen[2,2,3,2,1,8,0];
assert bpHeightClass[3,2,8](bpMul(bpint[2,2,1,2],sigma)) == zero[BPHeight(2,2,3,2,1,8)];
show bpScale[4](low);
assert bpScale[8](low) == zero[BPHeight(2,2,3,1,1,8)];
assert bpScale[2](high) == zero[BPHeight(2,2,3,2,1,8)];
assert bpHeightClass[3,1,8](bpHeightRep(low)) == low;
assert bpHeightClass[3,2,8](bpHeightRep(high)) == high;
show low;
show high;
bporder[8](bpHeightRep(low));
bporder[8](bpHeightRep(high));
