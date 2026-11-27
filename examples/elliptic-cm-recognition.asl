// Geometric Frobenius and the Lang map of an ordinary curve over F11.
rvunram F11[11,1]{0,1};
rvec E[F11](rvUnramElem[F11]{1},rvUnramElem[F11]{7});
rvellipticcm K[E];
let pi = rvECFrobenius[E];
let lang = rvECLang[E];
let kernel3 = rveckernel[E,3]{[7],[1]};
let kernel5 = rveckernel[E,5]{[7],[3],[1]};
rvecisogeny I[E,E3](kernel3);
let pushed5 = rvECPush(rvECIso[I],kernel5);
rvecisogeny J[E3,G](pushed5);
let composite = rvECMapCompose(rvECMap(rvECIso[J]),rvECMap(rvECIso[I]));
let loop = rvECMapCompose(rvECChart[G,E](rvUnramElem[F11]{1}),composite);
assert loop == lang;
assert (\point : RVECPoint(E) => rvECMapApply(loop,point)) == (\point : RVECPoint(E) => rvECInfinity[E]);
assert rvECCMRecognize[K](pi) == rvCMPi[K];
assert rvECCMRecognize[K](lang) == rvQpQuad[K,5,-1,2];
rvECCMCertify[K](loop);

// The same map over F121 uses a different Frobenius generator.
rvunram F121[11,1]{1,0,1};
rvec E2[F121](rvUnramElem[F121]{1,0},rvUnramElem[F121]{7,0});
rvellipticcm K2[E2];
let lifted = rvECExtend[E2,E2](loop);
assert rvECCMRecognize[K2](lifted) == rvQpQuad[K2,15,1,6];
rvECCMCertify[K2](lifted);
let pi_over_f11 = rvECExtend[E2,E2](pi);
let pi_over_f121 = rvECFrobenius[E2];
assert rvECMapCompose(pi_over_f11,pi_over_f11) == pi_over_f121;
assert rvECMapCompose(rvECMultMap[E2,3],lifted) == rvECSum(rvECMultMap[E2,14],pi_over_f121);

// A global automorphism of the j=0 curve lies outside Z[pi].
rvunram F13[13,1]{0,1};
rvec C[F13](rvUnramElem[F13]{0},rvUnramElem[F13]{2});
rvellipticcm L[C];
let alpha = rvECChart[C,C](rvUnramElem[F13]{3});
let conjugate = rvECChart[C,C](rvUnramElem[F13]{9});
let identity = rvECChart[C,C](rvUnramElem[F13]{1});
assert rvECCMRecognize[L](alpha) == rvQpQuad[L,-3,1,6];
assert rvECCMRecognize[L](conjugate) == rvQpQuad[L,-3,-1,6];
assert rvECSum(alpha,conjugate) == rvECNegate(identity);
assert rvECSum(rvECSum(alpha,alpha),alpha) == rvECSum(identity,rvECFrobenius[C]);
rvECCMCertify[L](alpha);
assert rvECMapCompose(alpha,conjugate) == identity;
rvEllLocalOrder[C,3];
let local_generator = rvEllLocalGen[L,C,3];
assert rvECCMRecognize[L](alpha) == add(rvCMPi[L,0],local_generator);
