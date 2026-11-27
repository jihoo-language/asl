// These quivers have identical nearby monodromy and local cohomology ranks.
shdisk P[5,2,2] [1,0,0,0] [0,1,0,0];
shdisk Q[5,2,2] [1,0,0,0] [0,0,0,1];
shdiskmat[P,T];
shdiskmat[Q,T];
shdiskmat[P,S];
shdiskmat[Q,S];
shdiskranks[P];
shdiskranks[Q];
shdiskhomdim[P,Q];
shdiskmap f[P,Q] [0,0,1,2] [0,3,0,1];
shdiskmap g[Q,Q] [2,0,0,2] [2,0,0,2];
assert shDiskComp(g,f) == pscale[2](f);
assert add(f,pscale[4](f)) == zero[ShDiskHom(P,Q)];

shdiskdual DP[P];
shdiskdual DQ[Q];
shdiskdual DDP[DP];
shdiskdual DDQ[DQ];
let df = shDiskDualMap[DQ,DP](f);
let ddf = shDiskDualMap[DDP,DDQ](df);
let bp = shDiskBidual[P,DDP];
let bq = shDiskBidual[Q,DDQ];
assert shDiskComp(ddf,bp) == shDiskComp(bq,f);
assert shDiskComp(shDiskBidualInv[P,DDP],bp) == shDiskId[P];
assert shDiskComp(bp,shDiskBidualInv[P,DDP]) == shDiskId[DDP];
assert shDiskDualMap[DQ,DP](shDiskComp(g,f))
    == shDiskComp(df,shDiskDualMap[DQ,DQ](g));
shdiskmapmat(df);

// Rank-one Picard–Lefschetz quivers over the same nearby vector space.
shdisk Left[5,2,1] [0,1] [-1,0];
shdisk Right[5,2,1] [1,0] [0,1];
shdiskmat[Left,T];
shdiskmat[Right,T];
sheaf Near[One,5] {pt:2} {};
shdiskoperator TLeft[Left,T,Near,Near];
shdiskoperator TRight[Right,T,Near,Near];
let l = shMapC[TLeft];
let r = shMapC[TRight];
assert shExtCompose(l,shExtCompose(r,l)) == shExtCompose(r,shExtCompose(l,r));

shdiskshriek Lower[Left];
shdiskstar Upper[Left];
shdiskmiddle Mid[Left];
let epi = shDiskCanonical[Lower,Mid];
let mono = shDiskCanonical[Mid,Upper];
assert shDiskComp(mono,epi) == shDiskCanonical[Lower,Upper];
assert shDiskComp(shDiskCanonical[Left,Upper],shDiskCanonical[Lower,Left])
    == shDiskCanonical[Lower,Upper];
let twice = pscale[2](shDiskId[Left]);
assert shDiskComp(shDiskExtMap[Mid,Mid](twice),epi)
    == shDiskComp(epi,shDiskExtMap[Lower,Lower](twice));
shdiskranks[Mid];

// The labels on the open circle and the cone are independent.
space Circle {d,a,c,b,a<c,a<d,b<c,b<d};
space Cone {cd,cb,tip,cc,ca,tip<ca,tip<cb,ca<cc,ca<cd,cb<cc,cb<cd};
let j = tab Circle -> Cone {a=>Cone.ca,b=>Cone.cb,c=>Cone.cc,d=>Cone.cd};
let i = tab One -> Cone {pt=>Cone.tip};
shdiskrealize Real[P,j,tip,a,b,c,d];
shdiskrealize RealLeft[Left,j,tip,a,b,c,d];
shextracomplex AtTip[i,Real];
shextracomplex AtLeft[i,RealLeft];
shlocaldim[AtTip,pt,relative,1];
shlocaldim[AtTip,pt,relative,2];
shlocaldim[AtLeft,pt,relative,1];
shlocaldim[AtLeft,pt,relative,2];

let w = shTotEmbed[Real__costalk,1](shTerm[Real__costalk__1,[pt],1]);
assert shD(shDiskCostalkLift[Real](w)) == shDiskCostalkLift[Real](shD(w));
let wleft = shTotEmbed[RealLeft__costalk,1](shTerm[RealLeft__costalk__1,[pt],0]);
assert shD(shDiskCostalkLift[RealLeft](wleft)) == shDiskCostalkLift[RealLeft](shD(wleft));
let closedw = shTotEmbed[Real__costalk,1](shTerm[Real__costalk__1,[pt],0]);
assert shLocalClass[AtTip,pt,relative](shDiskCostalkLift[Real](closedw))
    == shlocalgen[AtTip,pt,relative,1,0];
let tail = shTotEmbed[Real__costalk,2](shTerm[Real__costalk__2,[pt],1]);
assert shLocalClass[AtTip,pt,relative](shDiskCostalkLift[Real](tail))
    == shScale[-1](shlocalgen[AtTip,pt,relative,2,0]);

shdiskstalk Stalk[P];
shdiskcostalk Costalk[P];
shhyperdim[Stalk,0];
shhyperdim[Stalk,1];
shhyperdim[Costalk,1];
shhyperdim[Costalk,2];
