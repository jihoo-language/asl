rvverbitskytate T{5};
rvmf P[T]{2};
rvknorrer K[P];
rvKnorrerObject[K];

let x=rvMFBasis[P,P,0,0,0,1];
let xi=rvKnorrerInclude[K,K](x);
let odd=rvMFBasis[P,P,1,0,0,0];
let oi=rvKnorrerInclude[K,K](odd);
assert rvKnorrerProject(oi)==odd;
assert rvKnorrerCompose(oi,oi)==neg(xi);
rvMFTrace(rvMFCompose(x,odd));
rvKnorrerTrace(rvKnorrerCompose(xi,oi));

let a=rvKnorrerTerm[K,K,1,0,1]{0,0,0};
let b=rvKnorrerTerm[K,K,1,1,0]{0,0,0};
assert rvKnorrerProject(a)==zero[RVMFCochain(P,P,1)];
assert rvKnorrerProject(b)==zero[RVMFCochain(P,P,1)];
assert rvKnorrerProject(rvKnorrerCompose(a,b))==rvMFTerm[P,P,0,0,0,0,0];

let f=rvKnorrerTerm[K,K,0,0,0,(3,2)]{1,1,2};
let h=rvKnorrerHomotopy(f);
assert h==rvKnorrerTerm[K,K,1,0,1,(3,2)]{1,0,2};
assert add(rvKnorrerDiff(h),rvKnorrerHomotopy(rvKnorrerDiff(f)))==f;
assert rvKnorrerHomotopy(h)==zero[RVKnorrerCochain(K,K,0)];
rvKnorrerWeights(f);
rvKnorrerWeights(h);

let g=rvKnorrerTerm[K,K,1,3,2]{1,2,1};
let boundary=rvKnorrerDiff(g);
let perturbation=add(rvKnorrerRepresentative(xi),boundary);
assert rvKnorrerClass(perturbation)==xi;
assert rvKnorrerDiff(rvKnorrerPrimitive(boundary))==boundary;
rvKnorrerCone(perturbation);
rvmfcone C(rvKnorrerProject(perturbation));
rvknorrer KC[C];
rvKnorrerObject[KC];

rvmfshift Pshift[P];
rvknorrer Kshift[Pshift];
assert rvKnorrerShift[Kshift,Kshift](oi)==rvKnorrerInclude[Kshift,Kshift](rvMFShift[Pshift,Pshift](odd));
assert rvKnorrerShift[Kshift,Kshift](boundary)==rvKnorrerDiff(rvKnorrerShift[Kshift,Kshift](g));
