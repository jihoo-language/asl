// Four punctures with scalar monodromy −1 over F5.
shdisk P[5,1,1] [1] [2];
shdiskseq Sphere[P,P,P,P];
shdiskprod[Sphere];
shparrank[Sphere];

// Gauge a4=0. These cocycles have coordinates (a1,a2)=(1,0),(0,1).
let xr = shParRaw[Sphere]{[1],[0],[-1],[0]};
let yr = shParRaw[Sphere]{[0],[1],[1],[0]};
let x = shParClass(xr);
let y = shParClass(yr);
assert shParSum(xr) == zero[ShParVec(Sphere)];
assert shParSum(yr) == zero[ShParVec(Sphere)];
let nearby = shParVec[Sphere][-1];
assert shParPrimitive(shParBoundary(nearby)) == nearby;
assert shParClass(shParBoundary(nearby)) == zero[ShPar(Sphere)];
assert shParClass(shParRep(x)) == x;
assert shParClass(shParRep(y)) == y;

shhurwitz H1[Sphere] [1];
shhurwitz H2[Sphere] [2];
shdiskmonosame[Sphere,H1];
shdiskseqsame[Sphere,H1];
let s1 = \z:ShPar(Sphere) => shParIdentify[Sphere](shParBraid[H1](z));
let s2 = \z:ShPar(Sphere) => shParIdentify[Sphere](shParBraid[H2](z));
assert s1(x) == pscale[4](y);
assert s1(y) == add(x,pscale[2](y));
assert s2(x) == add(x,pscale[4](y));
assert s2(y) == y;
shparrows(shParRep(s1(s2(x))));
shparrows(shParRep(s2(s1(x))));

// The nearby basis is unchanged. The moved vanishing basis changes sign.
shdiskpoint Moved[H1,1];
shdiskmap neg[Moved,P] [1] [-1];
shdiskseqmap fix[H1,Sphere]{shDiskId[P],neg,shDiskId[P],shDiskId[P]};
shdiskseqmapat[fix,1];
assert shParMap[fix](shParBraid[H1](x)) == s1(x);
assert shParMap[fix](shParBraid[H1](y)) == s1(y);

shhurwitz A[Sphere] [1,2,1];
shhurwitz B[Sphere] [2,1,2];
shdiskseqsame[A,B];
assert shParIdentify[Sphere](shParBraid[A](xr)) == shParIdentify[Sphere](shParBraid[B](xr));
assert shParIdentify[Sphere](shParBraid[A](yr)) == shParIdentify[Sphere](shParBraid[B](yr));
shhurwitz Full[Sphere] [1,2,1,2,1,2];
shhurwitz Order1[Sphere] [1,1,1,1,1];
shhurwitz Order2[Sphere] [2,2,2,2,2];
let full = \z:ShPar(Sphere) => shParIdentify[Sphere](shParBraid[Full](z));
let o1 = \z:ShPar(Sphere) => shParIdentify[Sphere](shParBraid[Order1](z));
let o2 = \z:ShPar(Sphere) => shParIdentify[Sphere](shParBraid[Order2](z));
let id = \z:ShPar(Sphere) => z;
let minus = \z:ShPar(Sphere) => pscale[4](z);
assert full == minus;
assert o1 == id;
assert o2 == id;
