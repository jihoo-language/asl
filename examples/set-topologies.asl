space Diamond {r<a<t,r<b<t};
setomega O[Diamond];
settopology J[O] {a,b,t};
settopologyMap Recovered(setTopMap[J]);
assert setTopMap[Recovered] == setTopMap[J];
settopologysupport[Recovered];

setfamily F[Diamond] {r:1,a:2,b:2,t:2}
  {r<a:[0],r<b:[1],a<t:[0,1],b<t:[1,0]};
setfamily G[Diamond] {r:2,a:2,b:2,t:2}
  {r<a:[0,1],r<b:[0,1],a<t:[0,1],b<t:[0,1]};
setreflect R[J,F];
setreflectcheck[J,F];
setreflectcheck[J,R];
setreflectcheck[J,G];
setsections[F];
setsections[R];
elements RT[R];
card (\p : RT => S.hi);
setnat f[F,G] {r:[1],a:[1,0],b:[0,1],t:[1,0]};
let extension = setReflectExtend[R](f);
assert setComp(extension,setReflectUnit[R]) == f;
assert setReflectExtend[R](setComp(extension,setReflectUnit[R])) == extension;
assert setReflectExtend[R](setReflectUnit[R]) == setId[R];

let extend = \u : SetNat(F,G) => setReflectExtend[R](u);
let restrict = \u : SetNat(R,G) => setComp(u,setReflectUnit[R]);
assert (\u : SetNat(F,G) => restrict(extend(u))) == (\u : SetNat(F,G) => u);
assert (\u : SetNat(R,G) => extend(restrict(u))) == (\u : SetNat(R,G) => u);
card (\u : SetNat(F,G) => S.hi);
card (\u : SetNat(R,G) => S.hi);
setreflect RR[J,R];
assert setComp(setReflectInverse[RR],setReflectUnit[RR]) == setId[R];
assert setComp(setReflectUnit[RR],setReflectInverse[RR]) == setId[RR];

space Chain {x0<x1<x2<x3<x4<x5};
setomega P[Chain];
settopology E[P] {x0,x2,x4};
settopology U[P] {x1,x3,x5};
settopologymeet M[E,U];
settopologyjoin L[E,U];
assert setTopMap[M] == setId[P];
assert setTopMap[L] == setTop[P,P];
settopologyjoinsteps[E,U];
