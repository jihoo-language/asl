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

// The same elliptic basis is used in both factors, then transported to B.
rvjacpowerbasis EB[EL,3,2](rvUnramElem[L]{12}){pE,qE};
rvkummerbrauerpairing P[B,S,1];
rvkummerbrauerquotient Q[P,1];
rvkummercmcycles A[P,EB](rvUnramElem[K]{6});

// The normalized Mukai Frobenius sees the chosen B=v/3 modulo 3.
let v=rvKummerBrauerFixed[P,1]{1,0,0,0,0,0};
rvbrauermukai M[Q,1](v);
let u=rvBrauerMukai[M]{1,0,0,0,0,0,0,0};
let t=rvBrauerMukai[M]{0,0,0,0,0,0,0,1};
let fu=rvBrauerMukaiFrob[1](u);
assert fu==rvBrauerMukai[M]{1,1,0,1,2,0,0,2};
assert rvBrauerMukaiFrob[2](fu)==rvBrauerMukaiFrob[3](u);
assert rvBrauerMukaiFrob[3](u)==u;
assert rvBrauerMukaiFrob[0](u)==u;
assert add(fu,neg(fu))==zero[RVBrauerMukai(M)];
rvBrauerMukaiPair(u,t);
rvBrauerMukaiPair(fu,rvBrauerMukaiFrob[1](t));
rvBrauerMukaiData[M];
let witness=rvBrauerMukaiWitness[M];
assert witness==rvBrauerMukai[M]{0,2,0,0,0,0,0,0};
assert rvBrauerMukaiFrob[1](witness)==witness;

// The cubic restriction still has exact order three.
rvkummerbrauerpairing P3[B,S,3];
rvkummerbrauerquotient Q3[P3,1,P];
rvbrauermukaiextend M3[M,Q3];
let u3=rvBrauerMukaiTransport[M3](u);
assert rvBrauerMukaiFrob[1](u3)==rvBrauerMukaiTransport[M3](rvBrauerMukaiFrob[3](u));
rvBrauerMukaiData[M3];

let graph=rvKummerCMCycle[A]{0,0,0,1};
rvbrauermukaialgebraic V[M](graph);
assert rvBrauerMukaiFrob[1](rvBrauerMukaiTransport[V](u))==rvBrauerMukaiTransport[V](fu);
assert rvBrauerMukaiTransport[M](rvBrauerMukaiTransport[V](u))==u;
rvBrauerMukaiData[V];
rvbrauermukai G[Q,1](rvKummerCMCycleLift[1](graph));
rvBrauerMukaiData[G];
let g=rvBrauerMukaiWitness[G];
assert rvBrauerMukaiFrob[1](g)==g;
