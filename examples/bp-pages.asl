// The height-two p-adic Bockstein survives one page beyond height one.
bpBockdim[2, 1, 3, 0, 8, 0];
bpBockdim[2, 1, 4, 0, 8, 0];
bpBockdim[2, 2, 4, 0, 8, 0];
bpBockdim[2, 2, 5, 0, 8, 0];

// Source and target of the invariant-ideal d_2 in weight 16.
bpNovdim[2, 2, 2, 1, 16, 1];
bpNovdim[2, 2, 2, 2, 16, 3];
bpNovdim[2, 2, 3, 1, 16, 1];
bpNovdim[2, 2, 3, 2, 16, 3];

let v1 = bpv[2, 2, 0, 1];
let v2 = bpv[2, 2, 0, 2];
let fourth = bpPow[4](v1);
let corrected = add(fourth, bpMul(bpint[2, 2, 0, 8], bpMul(v1, v2)));
let height_two_d4 = bpBockClass[4, 8, 0](corrected);
assert bpBockNext(bpBockClass[3, 8, 0](fourth)) == height_two_d4;
assert bpBockD(height_two_d4) == bpBockgen[2, 2, 4, 1, 8, 4, 0];
assert bpBockD(height_two_d4) == bpBockClass[4, 8, 4](bpD(corrected));
