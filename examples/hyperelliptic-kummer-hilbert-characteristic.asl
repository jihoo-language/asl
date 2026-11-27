// The branch-node action has both fixed points and cycles of length three.
rvunram F[3,1]{0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
let t = rvUnramElem[F]{2};
rvhyperelliptic C[F]{z,t,z,t,o,o};
rvKummerHilbChar[C,1,1,2,2,12];
rvKummerHilbChar[C,2,1,4,2,8];
rvKummerHilbChar[C,2,1,4,5,8];
assert rvKummerHilbCharCoeff[C,2,1,4,2,8] == rvPrimeElt[2,1];

// The presentation also permits reduction modulo the field characteristic.
rvKummerHilbChar[C,1,1,2,3,6];
assert rvKummerHilbCharCoeff[C,1,1,2,3,1] == rvPrimeElt[3,2];
