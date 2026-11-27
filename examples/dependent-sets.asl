// The current fiber over lo is empty, but a future choice at hi remains.
setfamily F[S] {lo:0, hi:1} {lo<hi:[]};
elements E[F];
setfamily G[E] {hi_0:2} {};
setsum Sum[E,G];
setproduct Pi[E,G];
elements SumTotal[Sum];
elements PiTotal[Pi];
show setProjection[SumTotal];
show setProjection[PiTotal];

setfamily H[S] {lo:2, hi:2} {lo<hi:[1,0]};
let p = setProjection[E];
setpullback PH[H](p);
setpullback PP[Pi](p);
setnat choose[PH,G] {hi_0:[0,1]};
let curry = setTranspose[Pi](choose);
let evaluation = setEval[Pi,PP];
assert setComp(evaluation, setPullNat[PH,PP](curry)) == choose;
assert setUntranspose[PH](curry) == choose;
assert setTranspose[Pi](setUntranspose[PH](curry)) == curry;
show curry;

// The dependent sum has its own adjunction with pullback.
setnat tag[G,PH] {hi_0:[1,0]};
let folded = setSumLift[Sum](tag);
assert setSumRestrict[PH](folded) == tag;
assert setSumLift[Sum](setSumRestrict[PH](folded)) == folded;
