// Full geometric level 9 on an ordinary genus-two Jacobian.
rvunram K[37,1]{0,1};
rvhyperelliptic E[K]{rvUnramElem[K]{0},rvUnramElem[K]{1},rvUnramElem[K]{0},rvUnramElem[K]{1}};
rvhyperelliptic F[K]{rvUnramElem[K]{0},rvUnramElem[K]{1},rvUnramElem[K]{0},rvUnramElem[K]{1}};
rvellglue C[E,F](rvUnramElem[K]{1},rvUnramElem[K]{6},rvUnramElem[K]{18},rvUnramElem[K]{1});
rvjacchart S[C](rvUnramElem[K]{9});
rvunram L[37,1]{35,0,0,1};
rvjacbase EL[E,L];
rvjacbase FL[F,L];
rvellglue CL[EL,FL](rvUnramElem[L]{1},rvUnramElem[L]{6},rvUnramElem[L]{18},rvUnramElem[L]{1});
rvjacchart SL[CL](rvUnramElem[L]{9});
let pE=rvJacPoint[EL](rvUnramElem[L]{9,16,7},rvUnramElem[L]{7,3,34});
let qE=rvJacPoint[EL](rvUnramElem[L]{7,10,26},rvUnramElem[L]{23,33,31});
let pF=rvJacPoint[FL](rvUnramElem[L]{9,16,7},rvUnramElem[L]{7,3,34});
let qF=rvJacPoint[FL](rvUnramElem[L]{7,10,26},rvUnramElem[L]{23,33,31});
let e1=rvEllGlueMap[SL](pE,zero[RVJac(FL)]);
let e2=rvEllGlueMap[SL](zero[RVJac(EL)],pF);
let f1=rvEllGlueMap[SL](qE,zero[RVJac(FL)]);
let f2=rvEllGlueMap[SL](zero[RVJac(EL)],qF);
assert rvJacWeil[9](pE,qE)==rvUnramElem[L]{12};
assert rvJacWeil[9](e1,f1)==rvUnramElem[L]{33};
rvjacpowerbasis B[SL,3,2](rvUnramElem[L]{33}){e1,e2,f1,f2};


rvkummerbrauerpairing P[B,S,1];
rvkummerbrauerpairing P2[B,S,2];
rvkummerbrauerpairing P3[B,S,3];
rvkummerbrauerpairing P4[B,S,4];
rvkummerbrauerpairing P6[B,S,6];
rvkummerbrauerquotient Q[P,1];
rvkummerbrauerquotient Q2[P2,1];
rvkummerbrauerquotient Q3[P3,1,P];
rvkummerbrauerquotient Q4[P4,1];
rvKummerBrauerQuotientData[Q];
rvKummerBrauerQuotientData[Q3];

let v=rvKummerBrauerFixed[P,1]{1,0,0,0,0,0};
let w=rvKummerBrauerFixed[P,1]{0,0,1,0,0,0};
let divisor=rvKummerBrauerFixed[P,1]{0,1,0,0,0,0};
assert rvKummerBrauerAlgebraic[Q](divisor)==S.hi;
assert rvKummerBrauerAlgebraic[Q](v)==S.lo;
let x=rvKummerBrauerClass[Q](v);
assert rvKummerBrauerClass[Q](add(v,divisor))==x;
assert x==rvKummerBrauerClass[Q](rvKummerBrauerLift(x));
assert add(x,neg(x))==zero[RVKummerBrauerClass(Q)];
assert rvKummerBrauerScale[3](x)==zero[RVKummerBrauerClass(Q)];

let x2=rvKummerBrauerRestrict[Q2](x);
assert x2==rvKummerBrauerClass[Q2](rvKummerBrauerRestrict[P2](v));
assert rvKummerBrauerCorestrict[Q](x2)==rvKummerBrauerScale[2](x);
assert rvKummerBrauerRestrict[Q2](rvKummerBrauerCorestrict[Q](x2))==rvKummerBrauerScale[2](x2);
let x4=rvKummerBrauerRestrict[Q4](x);
assert rvKummerBrauerRestrict[Q4](x2)==x4;
assert rvKummerBrauerCorestrict[Q](rvKummerBrauerCorestrict[Q2](x4))==rvKummerBrauerCorestrict[Q](x4);
assert rvKummerBrauerCorestrict[Q](x4)==x;

let w2=rvKummerBrauerRestrict[P2](w);
rvKummerBrauerLink(v,w);
rvKummerBrauerLink(rvKummerBrauerRestrict[P2](v),w2);
rvKummerBrauerLink(v,rvKummerBrauerCorestrict[P](w2));
let w2class=rvKummerBrauerClass[Q2](w2);
rvKummerBrauerLink(rvKummerBrauerLift(x2),rvKummerBrauerLift(w2class));
rvKummerBrauerLink(rvKummerBrauerLift(x),rvKummerBrauerLift(rvKummerBrauerCorestrict[Q](w2class)));

let x3=rvKummerBrauerRestrict[Q3](x);
assert rvKummerBrauerAlgebraic[Q3](rvKummerBrauerLift(x3))==S.lo;
assert rvKummerBrauerCorestrict[Q](x3)==zero[RVKummerBrauerClass(Q)];
let w3=rvKummerBrauerRestrict[P3](w);
assert rvKummerBrauerCorestrict[P](w3)==zero[RVKummerBrauerFixed(P,1)];
rvKummerBrauerLink(rvKummerBrauerLift(x3),w3);

// The level-9 fixed lifts can be transferred without a level-9 quotient.
let z3=rvKummerBrauerFixed[P3,2]{1,0,0,0,0,0};
let z6=rvKummerBrauerRestrict[P6](z3);
assert rvKummerBrauerCorestrict[P3](z6)==rvKummerBrauerScale[2](z3);
