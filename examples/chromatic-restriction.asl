// The full-algebra May d_7 target disappears upon restriction to A(1).
let lift = mayRep(maygen[7, 4, 12, 0, 0]);
let restricted = resA[1](lift);
assert cobarD(restricted) == zero[CobarAn(1, 5, 12)];
assert resA[1](cobarD(lift)) == cobarD(restricted);

let v1fourth = extClass(restricted);
assert v1fourth == extangen[1, 4, 12, 0];
assert resA[0](v1fourth) == zero[ExtAn(0, 4, 12)];

extdim[4, 12];
extandim[1, 4, 12];
extandim[0, 4, 12];
show v1fourth;
