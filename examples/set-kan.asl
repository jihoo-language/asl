// The two elements at a identify choices from both branches in the colimit.
space Fork {a,b,c,a<b,a<c};
setfamily F[Fork] {a:2,b:2,c:3} {a<b:[0,0],a<c:[0,1]};
let f = tab Fork -> S {a=>S.lo,b=>S.lo,c=>S.hi};
setlan L[F](f);
setran R[F](f);
setfamily H[S] {lo:2,hi:2} {lo<hi:[1,0]};
setpullback PH[H](f);
setpullback PL[L](f);
setpullback PR[R](f);

setnat u[F,PH] {a:[0,0],b:[0,1],c:[1,1,0]};
setnat v[PH,F] {a:[0,1],b:[0,0],c:[1,0]};
let l = setLanLift[L](u);
let r = setRanTranspose[R](v);
assert setLanRestrict[PH](l) == u;
assert setRanUntranspose[PH](r) == v;
assert setLanLift[L](setLanRestrict[PH](l)) == l;
assert setRanTranspose[R](setRanUntranspose[PH](r)) == r;
assert setComp(setPullNat[PL,PH](l),setLanUnit[L,PL]) == u;
assert setComp(setRanCounit[R,PR],setPullNat[PH,PR](r)) == v;
show l;
show r;

setlan LH[PH](f);
setran RH[PH](f);
assert setLanRestrict[PH](setLanCounit[LH]) == setId[PH];
assert setRanUntranspose[PH](setRanUnit[RH]) == setId[PH];
assert setLanLift[L](setLanUnit[L,PL]) == setId[L];
assert setRanTranspose[R](setRanCounit[R,PR]) == setId[R];

// Two strict pullback squares with empty upper-left corner.
space Empty {};
setfamily Two[One] {pt:2} {};
let first = tab One -> S {pt=>S.lo};
let second = tab One -> S {pt=>S.hi};
let e = tab Empty -> One {};
setpullback E[Two](e);
setlan A[Two](first);
setran B[Two](second);
setlan AE[E](e);
setran BE[E](e);
setpullback PA[A](second);
setpullback PB[B](first);
show setLanBC[AE,PA];
show setId[PA];
show setRanBC[BE,PB];
show setId[BE];
