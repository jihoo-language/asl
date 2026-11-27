rvDworkHeckeEquation[9];
rvDworkHeckeCoeff[9,0,0];

rvunram F[41,1]{0,1};
let h=rvUnramElem[F]{};
let b=rvUnramElem[F]{2};
rvDworkHeckeNeighbors[9](h);
rvDworkHeckeTower[3](h);
rvDworkHeckeTower[3,5](h);
assert rvDworkHeckeCyclicCoeff[5,5,30](h)==rvUnramElem[F]{1};
assert rvDworkHeckeCyclicCoeff[3,3,8](h)==rvDworkHeckeSpecCoeff[9,8](h);
let rhs=rvUnramMul(rvUnramPow[4](add(h,neg(b))),rvDworkHeckeEval[9](h,b));
assert rvDworkHeckeCompose[3](h,b)==rhs;
assert rvDworkHeckeCompose[3,5](h,b)==rvDworkHeckeCompose[5,3](h,b);

rvunram C[3,1]{-1,-1,0,1};
let c=rvUnramElem[C]{0,1};
rvDworkHeckeNeighbors[9](c);
rvDworkHeckeTower[3](c);

// These are three different Hauptmodul values on a cyclic 3-by-3 path.
let t=rvUnramElem[F]{7};
let mid=rvUnramElem[F]{24};
let end=rvUnramElem[F]{29};
let ha=rvDworkHauptmodulT(t);
let hb=rvDworkHauptmodulT(mid);
let hc=rvDworkHauptmodulT(end);
assert ha==rvUnramElem[F]{7};
assert hb==h;
assert hc==rvUnramElem[F]{26};
assert rvDworkHeckeEval[3](ha,hb)==zero[RVUnram(F)];
assert rvDworkHeckeEval[3](hb,hc)==zero[RVUnram(F)];
assert rvDworkHeckeEval[9](ha,hc)==zero[RVUnram(F)];

rvdworkfricke E,E2(t);
rvdworkfricke M,M2(mid);
rvdworkfricke N,N2(end);
let point=rvECPoint[E](rvUnramElem[F]{6},rvUnramElem[F]{9});
let three=rvECPoint[E](rvUnramElem[F]{7},rvUnramElem[F]{29});
assert rvECScale[3](point)==three;
assert rvECScale[9](point)==rvECInfinity[E];

let firstkernel=rveckernel[E,3]{[34],[1]};
rvecisogeny I[E,Q](firstkernel);
let first=rvECMap(rvECIso[I]);
let image=rvECPoint[Q](rvUnramElem[F]{6},rvUnramElem[F]{27});
assert rvECMapApply(first,point)==image;
let secondkernel=rveckernel[Q,3]{[35],[1]};
rvecisogeny J[Q,R](secondkernel);
let second=rvECMap(rvECIso[J]);
let scale=rvUnramElem[F]{14};
let normalizedmid=rvECChart[Q,M](scale);
let normalizedend=rvECChart[R,N](scale);
let cyclic=rvECMapCompose(normalizedend,rvECMapCompose(second,first));
assert rvECMapApply(cyclic,point)==rvECInfinity[N];
let normalizedimage=rvECPoint[M](rvUnramElem[F]{28},rvUnramElem[F]{1});
assert rvECMapApply(normalizedmid,image)==normalizedimage;
let normalizedkernel=rveckernel[M,3]{[13],[1]};
rvecisogeny Jmid[M,S](normalizedkernel);
let directend=rvECChart[S,N](rvUnramElem[F]{1});
let alternative=rvECMapCompose(directend,
    rvECMapCompose(rvECMap(rvECIso[Jmid]),rvECMapCompose(normalizedmid,first)));
assert alternative==cyclic;
rvECMapDegree(cyclic);
rvECMapCertify(cyclic);
