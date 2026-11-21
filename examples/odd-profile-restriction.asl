// A corrected lift supports d_5 over the full algebra at p=3.
let source = pmaygen[3, 5, 3, 15, 6, 0];
let lift = mayRep(source);
assert mayD(source) == pmaygen[3, 5, 4, 15, 11, 0];

let restricted = resA[1](lift);
assert cobarD(restricted) == resA[1](cobarD(lift));
assert cobarD(restricted) == zero[PCobarAn(3, 1, 4, 15)];
let period = extClass(restricted);
assert period == pscale[2](pextangen[3, 1, 3, 15, 0]);
assert resA[0](period) == zero[PExtAn(3, 0, 3, 15)];

pextdim[3, 3, 15];
pextandim[3, 1, 3, 15];
pextandim[3, 0, 3, 15];
show period;
