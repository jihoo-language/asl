let v2 = bpReduce[1](bpv[2,2,0,2]);
let power = bpPow[2](v2);
let x = bpThickExtClass[12](bpThick[2](power));
let y = bpThickDelta(x);
let beta = bpDelta(y);
show x;
show y;
show beta;
assert y == add(bpIExtgen[2,2,1,1,8,0],bpIExtgen[2,2,1,1,8,2]);
assert beta == bpExtClass[8](bpDelta(bpThickDelta(bpThick[2](power))));
assert bpScale[2](beta) == zero[BPExt(2,2,2,8)];

let v2_odd = bpReduce[1](bpv[3,2,0,2]);
let odd = bpThickExtClass[48](bpThick[3](bpPow[3](v2_odd)));
show odd;
show bpThickDelta(odd);
show bpDelta(bpThickDelta(odd));
