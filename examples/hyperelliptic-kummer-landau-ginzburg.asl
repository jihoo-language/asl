rvverbitskytate T{5};
rvmf P[T]{2};
rvMFObject[P];
rvMFHom[P,P,0];
rvMFHom[P,P,1];

let x=rvMFBasis[P,P,0,0,0,1];
let odd=rvMFBasis[P,P,1,0,0,0];
assert rvMFCompose(odd,odd)==neg(x);
let xodd=rvMFStableAction(rvVerbitskyTateTerm[T,stable]{1}{0}{0},odd);
assert xodd==rvMFCompose(x,odd);
rvMFTrace(xodd);
rvMFPair(x,odd);
assert rvMFStableAction(rvVerbitskyTateTerm[T,stable]{2}{0}{0},odd)==zero[RVMF(P,P,1)];
assert rvMFStableAction(rvVerbitskyTateTerm[T,stable]{0}{0}{-1000000},odd)==odd;

rvMFCone(x);
rvmfcone C(x);
rvMFObject[C];
rvmfshift Q[P];
rvMFObject[Q];
rvMFTrace(rvMFShift[Q,Q](xodd));

let h=rvMFTerm[P,P,1,0,0,0,1];
let boundary=rvMFDiff(h);
assert rvMFClass(boundary)==zero[RVMF(P,P,0)];
assert rvMFDiff(rvMFPrimitive(boundary))==boundary;
