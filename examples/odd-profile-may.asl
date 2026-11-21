// Restriction kills the full d_5 target and leaves the corrected lift closed.
let full = pmaygen[3, 5, 3, 15, 6, 0];
let finite = resA[1](full);
assert finite == pmayangen[3, 1, 5, 3, 15, 6, 0];
assert mayD(finite) == resA[1](mayD(full));
assert mayD(finite) == zero[PMayAn(3, 1, 5, 4, 15, 11)];
let survives = mayNext(finite);
assert survives == pmayangen[3, 1, 6, 3, 15, 6, 0];
assert extClass(mayRep(survives)) == pscale[2](pextangen[3, 1, 3, 15, 0]);
assert resA[0](survives) == zero[PMayAn(3, 0, 6, 3, 15, 6)];
show survives;
pmayandim[3, 1, 6, 3, 15, 6];
