// An ordinary Kummer surface with varying 3-primary Brauer group.
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
rvJacPowerFrob[B,S,1];
rvjactorsionbasis B3[SL,3](rvUnramElem[L]{10}){
    rvJacScale[3](e1),rvJacScale[3](e2),rvJacScale[3](f1),rvJacScale[3](f2)
};
rvKummerK3Brauer[B3,S,1];
rvKummerK3Brauer[B,S,1];
rvKummerK3Brauer[B,S,3];
rvKummerK3BrauerArtinTate[B,S,3];
rvKummerK3BrauerArtinTate[B,S,9];
rvKummerK3ArtinTate[S,1];

// A superspecial Kummer surface with trivial 3-primary Brauer group over F25.
rvunram K5[5,1]{0,1};
rvhyperelliptic X[K5]{rvUnramElem[K5]{0},rvUnramElem[K5]{1},rvUnramElem[K5]{0},
                    rvUnramElem[K5]{0},rvUnramElem[K5]{0},rvUnramElem[K5]{1}};
rvunram W[5,1]{3,0,1};
rvjacbase XW[X,W];
let a=rvUnramElem[W]{0,1};
let u=rvJacScale[2](rvJacPoint[XW](rvUnramElem[W]{1},a));
let v=add(u,rvJacScale[2](rvJacPoint[XW](rvUnramElem[W]{1,2},a)));
let w=rvJacScale[4](rvJacPoint[XW](rvUnramElem[W]{1,1},a));
let t=rvJacScale[4](rvJacPoint[XW](rvUnramElem[W]{1,3},a));
rvjacpowerbasis Z[XW,3,1](rvUnramElem[W]{2,2}){u,v,w,t};
rvKummerK3Brauer[Z,X,2];
