// Frobenius over F25 with lambda = X and X^2 + 2 = 0.
rvunram W[5,5]{2,0,1};
let lambda = rvUnramElem[W]{0,1};
rvlegendrecm K(lambda);
let pi = rvCMPi[K];
let bar = rvQpQuadConj(pi);
assert rvQpQuadMul(pi,bar) == rvQpQuad[K,25,0,1];
let unit = rvQpQuadSpecialize[W](pi);
assert unit == rvLegendreUnit(lambda);
assert unit == rvUnramElem[W]{2823};
show unit;

rvCMOrder[K,3];
rvCMOrder[K,4];
rvCMOrder[K,6];
rvCMNorm[K,6];
rvCMIndex[K,3,6];
let g = rvCMTerm[K,6,0,1];
let h = rvCMInclude[3](g);
assert h == rvCMTerm[K,3,-15625,142];
assert rvCMInclude[1](h) == rvCMInclude[1](g);
assert rvCMInclude[1](g) == rvCMTerm[K,1,-8525,-2982];
assert rvCMLift(g) == rvCMPi[K,6];
assert rvCMLift(h) == rvQpQuadPow[6](pi);
let x = add(g,rvCMTerm[K,6,3,0]);
assert rvCMInclude[3](rvCMMul(x,g)) == rvCMMul(rvCMInclude[3](x),h);
assert rvCMLift(rvCMMul(x,g)) == rvQpQuadMul(rvCMLift(x),rvCMLift(g));

// A different precision and degree use the same split quadratic embedding.
rvunram P[5,4]{0,1};
assert rvQpQuadSpecialize[P](pi) == rvUnramElem[P]{323};
let y = rvUnramElem[P]{2};
rvlegendrecm Q(y);
let alpha = rvCMPi[Q];
assert rvQpQuadMul(alpha,rvQpQuadConj(alpha)) == rvQpQuad[Q,5,0,1];
assert rvQpQuadSpecialize[P](alpha) == rvLegendreUnit(y);
show rvQpQuadSpecialize[P](alpha);
