let t = bpt[3, 1, 1, 1, 1];
let x = bpBockClass[0, 4, 0](t);
let y = bpBockClass[0, 8, 0](bpPow[2](t));
let product = bpCup(x, y);
let image = bpBockD(product);
assert image == pscale[2](bpCup(x, bpBockD(y)));
show image;

let v1 = bpv[2, 2, 0, 1];
let v2 = bpv[2, 2, 0, 2];
let g = add(bpPow[4](v1), bpMul(bpint[2, 2, 0, 8], bpMul(v1, v2)));
let cycle = bpBockClass[4, 8, 0](g);
let p = bpBockClass[4, 0, 1](bpint[2, 2, 0, 2]);
let shifted = bpBockD(bpCup(p, cycle));
assert shifted == bpCup(p, bpBockD(cycle));
assert shifted == bpBockgen[2, 2, 4, 1, 8, 5, 0];
show shifted;
