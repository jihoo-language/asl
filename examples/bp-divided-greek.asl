let v2 = bpReduce[1](bpv[2,2,0,2]);
let input = bpThick[2](bpPow[2](v2));
assert bpD(input) == zero[BPThick(2,2,1,2,2)];
let first = bpThickDelta(input);
let beta_2_2 = bpDelta(first);
assert bpD(beta_2_2) == zero[BP(2,2,3)];
show beta_2_2;
bporder[8](beta_2_2);

let odd = bpReduce[1](bpv[3,2,0,2]);
let beta_3_3 = bpDelta(bpThickDelta(bpThick[3](bpPow[3](odd))));
assert bpD(beta_3_3) == zero[BP(3,2,3)];
bporder[36](beta_3_3);
