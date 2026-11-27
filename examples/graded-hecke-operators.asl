let t = rvHeckeT[24,2];
let one = rvHeckeUnit[24];
rvHeckeActionTrace(t);
rvHeckeChar(t);

let relation = add(add(rvHeckeMul(t,t),rvHeckeScale[-1080](t)),rvHeckeScale[-20468736](one));
assert rvHeckeSameAction(relation,zero[RVHecke(24)]) == S.hi;
show relation;

let mixed = add(rvHeckeScale[3](one),add(rvHeckeScale[-2](t),rvHeckeT[24,6]));
rvHeckeChar(mixed);
