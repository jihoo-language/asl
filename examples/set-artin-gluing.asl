// Relabeled strata. The open arrow forces compatibility of the section coordinates.
space X {v,c,d,u,c<d,c<u,d<v,u<v};
space C {dc,cc,cc<dc};
space U {vu,uu,uu<vu};
let closed = tab C -> X {dc=>X.d,cc=>X.c};
let open = tab U -> X {vu=>X.v,uu=>X.u};
setrecollement R(closed,open);

setfamily F[X] {c:2,d:3,u:2,v:3} {
    c<d:[1,2],c<u:[1,0],d<v:[0,2,2],u<v:[2,2]
};
setrestrict A[R,F,closed];
setrestrict B[R,F,open];
setgluelimit T[R,B];
let alpha = setAttach[A,T];
show alpha;
setglue G[R](alpha);
let recover = setGlueRecover[G,F];
let forget = setGlueForget[G,F];
assert setComp(forget,recover) == setId[F];
assert setComp(recover,forget) == setId[G];
setrecollement Alias(closed,open);
setglue Identical[Alias](alpha);
assert setComp(setGlueRecover[Identical,F],forget) == setGlueMap[G,Identical](setId[A],setId[B]);

setfamily A2[C] {dc:2,cc:1} {cc<dc:[1]};
setfamily B2[U] {vu:2,uu:1} {uu<vu:[1]};
setgluelimit T2[R,B2];
setnat alpha2[A2,T2] {dc:[0,1],cc:[0]};
setglue G2[R](alpha2);
setnat a[A,A2] {dc:[0,1,1],cc:[0,0]};
setnat b[B,B2] {vu:[0,0,1],uu:[0,0]};
assert setComp(setGlueLimitMap[T,T2](b),alpha) == setComp(alpha2,a);
let m = setGlueMap[G,G2](a,b);
show m;
assert setGlueClosed[G,G2](m) == a;
assert setGlueOpen[G,G2](m) == b;
assert setGlueMap[G,G2](setGlueClosed[G,G2](m),setGlueOpen[G,G2](m)) == m;
assert setGlueMap[G,G](setId[A],setId[B]) == setId[G];
let throughStrata = \m:SetNat(G,G2) => setGlueMap[G,G2](setGlueClosed[G,G2](m),setGlueOpen[G,G2](m));
let identity = \m:SetNat(G,G2) => m;
assert throughStrata == identity;

setran RB[B](open);
setpullback CRB[RB](closed);
setnat bridge[T,CRB] {dc:[0,1,2],cc:[0,1]};
setnat back[CRB,T] {dc:[0,1,2],cc:[0,1]};
assert setComp(back,bridge) == setId[T];
assert setComp(bridge,back) == setId[CRB];

// Regroup three strata by gluing the two open layers first.
let openLower = tab One -> U {pt=>U.uu};
let openUpper = tab One -> U {pt=>U.vu};
setrecollement Inside(openLower,openUpper);
setrestrict L[Inside,B,closed];
setrestrict H[Inside,B,open];
setgluelimit ULimit[Inside,H];
setglue BAgain[Inside](setAttach[L,ULimit]);
let rebuildB = setGlueRecover[BAgain,B];
let forgetB = setGlueForget[BAgain,B];
setgluelimit TAgain[R,BAgain];
let attachedAgain = setComp(setGlueLimitMap[T,TAgain](rebuildB),alpha);
setglue GAgain[R](attachedAgain);
let regroup = setGlueMap[G,GAgain](setId[A],rebuildB);
let ungroup = setGlueMap[GAgain,G](setId[A],forgetB);
assert setComp(ungroup,regroup) == setId[G];
assert setComp(regroup,ungroup) == setId[GAgain];
