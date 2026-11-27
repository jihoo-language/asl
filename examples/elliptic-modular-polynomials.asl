rvModularPolyDegree[2];
rvModularPolyCoeff[2,2,2];
rvModularPolyCoeff[2,0,0];
rvModularPolyCoeff[5,5,4];

// The degree-three Vélu quotient over F13 has j-invariants 3 and 0.
rvunram F[13,1]{0,1};
rvec E[F](rvUnramElem[F]{4},rvUnramElem[F]{1});
let kernel = rveckernel[E,3]{[2],[1]};
rvecisogeny I[E,Q](kernel);
assert rvECA[Q] == rvUnramElem[F]{0};
assert rvECB[Q] == rvUnramElem[F]{11};
assert rvModularPolyEval[3](rvUnramElem[F]{3},rvUnramElem[F]{0}) == zero[RVUnram(F)];

// In characteristic eleven the full integer polynomial is unnecessary.
rvunram W[11,1]{1,0,1};
let i = rvUnramElem[W]{0,1};
let one = rvUnramElem[W]{1};
assert rvModularPolyEval[11](i,one) == rvUnramElem[W]{2};
assert rvModularPolySpecCoeff[11,11](i) == i;
assert rvModularPolySpecCoeff[11,1](i) == neg(i);
show rvModularPolySpecCoeff[11,0](i);
