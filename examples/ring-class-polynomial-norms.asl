rvcmclasspoly H3[-3];
rvcmclasspoly H4[-4];
rvcmclasspolylift H12[H3,2];
rvcmclasspolylift H16[H4,2];
rvcmclasspolylift H27[H3,3];
rvcmclasspolylift H48[H12,2];

// Reach conductor six in two orders and verify the agreement in the registry.
rvcmclasspolylift H108a[H12,3];
rvcmclasspolylift H108b[H27,2];
rvcmclasspolylift H432[H108a,2];

rvCMClassPolyData[H48];
rvCMClassPolyCoeff[H48,0];
rvCMClassPolyCoeff[H48,1];
rvCMClassPolyCoeff[H16,0];
rvCMClassPolyCoeff[H27,0];
rvCMClassPolyData[H108b];
rvCMClassPolyCoeff[H108a,0];
rvCMClassPolyCoeff[H108b,0];
rvCMClassPolyData[H432];

// H_-48 reduces to Y^2 in characteristic five.
rvunram F[5,1]{2,0,1};
let t = rvUnramElem[F]{0,1};
assert rvCMClassPolyModCoeff[H48,F,0] == zero[RVUnram(F)];
assert rvCMClassPolyModCoeff[H48,F,1] == zero[RVUnram(F)];
assert rvCMClassPolyModCoeff[H48,F,2] == rvUnramElem[F]{1};
assert rvCMClassPolyEval[H48](t) == rvUnramPow[2](t);
