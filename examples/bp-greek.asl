// The two-step p = 3 connecting cocycle in weight 12.
let beta_mod_3 = bpDelta(bpReduce[2](bpv[3, 2, 0, 2]));
let beta = bpDelta(beta_mod_3);
show beta_mod_3;
show beta;
assert bpD(beta) == zero[BP(3, 2, 3)];
assert bpMul(bpint[3, 2, 2, 3], beta) == bpD(bpLift(beta_mod_3));

// At height one and p = 2, this representative detects order 8.
let fourth = bpPow[4](bpv[2, 1, 0, 1]);
let alpha_4_3 = bpBockP[3](fourth);
show alpha_4_3;
assert bpD(alpha_4_3) == zero[BP(2, 1, 2)];
assert bpMul(bpint[2, 1, 1, 8], alpha_4_3) == bpD(fourth);
