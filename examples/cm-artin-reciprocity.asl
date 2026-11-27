rvcmorder O[-23];
let unit = rvCMPrincipal[O];
let alpha = rvCMClass[O,2,1,3];
let beta = rvCMInverse(alpha);
assert rvCMArtinClass[O,2,0] == alpha;
assert rvCMArtinClass[O,2,1] == beta;
assert rvCMArtinTranslate[O,2,0,2](unit) == beta;
assert rvCMArtinTranslate[O,2,1,-2](alpha) == unit;
rvCMArtinPrime[O,2,0];
rvCMArtinOrbit[O,2,0](unit);
rvCMArtinOrbit[O,2,1](unit);

rvcmclasspoly H3[-3];
rvcmclasspolylift H12[H3,2];
rvcmclasspolylift H48[H12,2];
rvunram F[7,1]{0,1};
rvunram L[7,1]{5,5,1};
rvunram W[7,3]{5,5,1};
rvCMArtinFactors[H48,F];
rvCMArtinFactors[H48,L];

rvec E[L](rvUnramElem[L]{6,4},rvUnramElem[L]{4,5});
rvecendorder M[E,H48];
rvCMArtinLift[M,W];
