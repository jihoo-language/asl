// A chosen integral lift over F3. Coefficients go from constant to leading.
rvhyperelliptic C[3]{0,1,1,1,2,1};
rvHWStatus[C];
rvHWPoints[C,1];
rvHWPoints[C,2];
rvHWWeil[C];

let a = rvHW[C,3,3];
let b = rvHW[C,2,3];
let m = rvHWUnit[C,3];
assert a == rvHWMat[3,2,3]{11,2,11,22};
assert b == rvHWMat[3,2,3]{20,6,24,17};
assert m == rvHWMat[3,2,3]{16,4,13,11};
assert rvHWCompose(m,b) == a;
assert rvHWCompose(rvHWInv(b),a) == rvHWMat[3,2,3]{25,13,13,2};
assert rvHWCompose(b,rvHWInv(b)) == rvHW[C,0,3];
assert rvHWChar[0](m) == rvHWScalar[3,3,1];
assert rvHWChar[1](m) == zero[RVHWScalar(3,3)];
assert rvHWChar[2](m) == rvHWScalar[3,3,16];
rvHWVal(rvHWChar[2](m));
show a;
show b;
show m;
rvHWConverges[C,3];
rvHWCertify[C,3];
assert rvHWUnit[C,7] == rvHWMat[3,2,7]{1420,706,931,767};
rvHWConverges[C,7];
rvHWCertify[C,7];

// The three-by-three matrices also give distinct left and right ratios.
rvhyperelliptic G[3]{1,1,1,2,2,1,0,1};
let u = rvHWUnit[G,3];
assert u == rvHWMat[3,3,3]{7,22,12,25,8,1,10,23,0};
assert rvHWCompose(rvHWInv(rvHW[G,2,3]),rvHW[G,3,3])
    == rvHWMat[3,3,3]{13,4,3,10,26,19,19,14,3};
rvHWPoints[G,3];
rvHWWeil[G];
rvHWCertify[G,3];

// Cross-check against the signed Legendre Dwork unit over the prime field.
rvhyperelliptic E[7]{0,3,-4,1};
assert rvHWEntry[1,1](rvHWUnit[E,2]) == rvHWScalar[7,2,39];
rvHWWeil[E];
legendreDwork[7,3,2];
