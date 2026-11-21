let v = bpv[2, 1, 0, 1];
let fourth = bpPow[4](v);
let height_one = bpBockClass[3, 8, 0](fourth);
let d3 = bpBockD(height_one);
assert d3 == bpBockgen[2, 1, 3, 1, 8, 3, 0];
assert bpRaise[2](d3) == zero[BPBock(2, 2, 3, 1, 8, 3)];
let raised = bpRaise[2](height_one);
assert bpBockD(raised) == bpRaise[2](d3);

let v1 = bpRaise[2](v);
let v2 = bpv[2, 2, 0, 2];
let v1v2 = bpMul(v1, v2);
let g = add(bpRaise[2](fourth), bpMul(bpint[2, 2, 0, 8], v1v2));
let next = bpBockNext(raised);
assert next == bpBockClass[4, 8, 0](g);
show bpBockD(next);
assert bpBockD(next) == bpBockgen[2, 2, 4, 1, 8, 4, 0];

let alpha = bpBockP[3](fourth);
let sigma = bpBockP[4](g);
let twice_sigma = bpMul(bpint[2, 2, 1, 2], sigma);
assert add(bpRaise[2](alpha), bpD(v1v2)) == twice_sigma;
let higher_sigma = bpExtClass[8](sigma);
let map_alpha = bpRaise[2](bpExtClass[8](alpha));
assert map_alpha == bpScale[2](higher_sigma);
show map_alpha;
