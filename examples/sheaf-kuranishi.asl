space X {a,l,p,q,b,c,a<l<b<c,a<p<q<b};
sheaf F[X,5] {a:1,l:1,p:1,q:1,b:1,c:1}
  {a<l:[0],l<b:[0],a<p:[0],p<q:[0],q<b:[0],b<c:[0]};
let x0 = shClass(shHomTerm[F,F,[a,l],0,0]);
let x1 = shClass(shHomTerm[F,F,[l,b],0,0]);
let y0 = shClass(shHomTerm[F,F,[a,p],0,0]);
let y1 = shClass(shHomTerm[F,F,[p,q],0,0]);
let y2 = shClass(shHomTerm[F,F,[q,b],0,0]);
let tail = shClass(shHomTerm[F,F,[b,c],0,0]);
shkuranishi K[F,4] {x0,x1,y0,y1,y2,tail};
shkuranishiDims[K];
shkuranishiEq[K,0];

let r2 = shKurR[K,[1,1,0,0,0,0]];
let kappa2 = shKurKappa[K,[1,1,0,0,0,0]];
assert shClass(r2) == kappa2;
let r3 = shKurR[K,[1,1,0,0,0,1]];
let a3 = shKurA[K,[1,1,0,0,0,1]];
let q3 = add(r3,shScale[-1](shD(a3)));
let e = shKurA[K,[0,0,0,0,0,1]];
assert shD(q3) == add(shExtCompose(e,r2),shScale[-1](shExtCompose(r2,e)));
show shD(q3);

let realized = shKurArc[K,[1,1,1,1,1,1],[1,2,1,1,1,1],4];
assert shKurArcKappa[K,[1,1,1,1,1,1],[1,2,1,1,1,1],3] == zero[ShExt(F,F,2)];
assert shDefTerm[4](realized) == shHomTerm[F,F,[a,c],0,0];
assert shDefObstruction(realized) == zero[ShExt(F,F,2)];
show shDefTerm[4](realized);

shkuranishi Diagonal[F,2] {add(x0,x1)};
shkuranishiEq[Diagonal,0];
