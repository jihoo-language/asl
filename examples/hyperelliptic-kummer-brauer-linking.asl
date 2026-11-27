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
let v=rvKummerBrauerFixed[P,1]{1,0,0,0,0,0};
let w=rvKummerBrauerFixed[P,1]{0,0,1,0,0,0};
rvKummerBrauerLink(v,w);
rvKummerBrauerLink(w,v);
rvKummerBrauerLink(v,v);
assert rvKummerBrauerLinkNonzero(v,w)==S.hi;
assert rvKummerBrauerLinkNonzero(v,v)==S.lo;
assert add(v,neg(v))==zero[RVKummerBrauerFixed(P,1)];
assert rvKummerBrauerScale[-2](v)==v;
rvKummerBrauerLink(add(v,v),w);
rvKummerBrauerLinkForm[P];

// Change the primitive root from 33 to 33^2, and change the f-basis with it.
rvjacpowerbasis Btw[SL,3,2](rvUnramElem[L]{16}){
    e1,e2,rvJacScale[2](f1),rvJacScale[2](f2)
};
rvkummerbrauerpairing T[Btw,S,1];
let vt=rvKummerBrauerChange[T](v);
let wt=rvKummerBrauerChange[T](w);
assert vt==rvKummerBrauerFixed[T,1]{2,0,0,0,0,0};
assert wt==rvKummerBrauerFixed[T,1]{0,0,1,0,0,0};
assert rvKummerBrauerChange[P](vt)==v;
rvKummerBrauerLink(vt,wt);

// Over the cubic base extension the order-three subgroup is isotropic.
rvkummerbrauerpairing P3[B,S,3];
let v3=rvKummerBrauerFixed[P3,1]{1,0,0,0,0,0};
let w3=rvKummerBrauerFixed[P3,1]{0,0,1,0,0,0};
rvKummerBrauerLink(v3,w3);
rvKummerBrauerLinkForm[P3];
rvKummerK3BrauerArtinTate[B,S,3];
