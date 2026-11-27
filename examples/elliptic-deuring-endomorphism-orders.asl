rvcmclasspoly H4[-4];
rvcmclasspolylift H16[H4,2];
rvunram F[13,1]{0,1};
rvec G[F](rvUnramElem[F]{4},rvUnramElem[F]{});
rvec N[F](rvUnramElem[F]{2},rvUnramElem[F]{1});

assert rvECEndJ[G] == rvUnramElem[F]{12};
assert rvECEndJ[N] == rvUnramElem[F]{1};
rvECEndCheck[G,H16];
rvECEndCheck[N,H4];
rvecendorder Maximal[G,H4];
rvecendorder Nonmaximal[N,H16];
rvECEndData[Maximal];
rvECEndData[Nonmaximal];

let pg = rvECEndFrobenius[Maximal];
let pn = rvECEndFrobenius[Nonmaximal];
assert pg == rvECEndElem[Maximal,3,2];
assert pn == rvECEndElem[Nonmaximal,3,1];
assert rvECEndMul(pg,rvECEndConj(pg)) == rvECEndElem[Maximal,13,0];
assert rvECEndMul(pn,rvECEndConj(pn)) == rvECEndElem[Nonmaximal,13,0];
assert add(add(rvECEndMul(pn,pn),rvECEndMul(rvECEndElem[Nonmaximal,-6,0],pn)),rvECEndElem[Nonmaximal,13,0]) == zero[RVECEnd(Nonmaximal)];
rvECEndTrace(pg);
rvECEndNorm(pn);

// The exceptional j=0 is a root of H_-3 in characteristic five,
// but this curve is supersingular.
rvcmclasspoly H3[-3];
rvunram S[5,1]{0,1};
rvec Sup[S](rvUnramElem[S]{},rvUnramElem[S]{1});
rvECEndCheck[Sup,H3];
