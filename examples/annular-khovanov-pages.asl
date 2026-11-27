// The E1 differential is the annular-degree -2 action on an integral cycle.
khbraid Cross[2,[1]];
let high = khABasis[Cross,0,3,2,0];
let pHigh = khPageClass[1,2](khMod[3,2](high));
let odd = khPageD(pHigh);
assert odd == khPageClass[1,0](khMod[3,0](khVMinus(high)));
assert odd == khPageGen[Cross,3,1,1,3,0,0];
assert khPageD(odd) == zero[KhPage(Cross,3,1,2,3,-2)];
assert khD(khPageRep(khPageBoundary(odd))) == khMod[3,2](khPageRep(odd));
khPageDim[Cross,3,1,0,3,2];
khPageDim[Cross,3,2,0,3,2];

// The exact Smith null detects two adjacent UCT contributions at p=2.
khbraid Trefoil[2,[1,1,1]];
khPageDim[Trefoil,2,1,2,7,0];
khPageDim[Trefoil,2,1,3,7,0];
khPageDim[Trefoil,3,1,2,7,0];
khPageDim[Trefoil,3,1,3,7,0];

// A class of the transverse cycle is hit on page two for A.
khbraid A[4,[3,-2,-2,3,3,2,-3,-1,2,1,1]];
let a = khMod[2,-4](khPsi[A]);
let a1 = khPageClass[1,-4](a);
let a2 = khPageNext(a1);
assert a2 == khPageGen[A,2,2,0,-1,-4,0];
assert khPagePrev(a2) == a1;
let hitA = khPageBoundary(a2);
assert khPageD(hitA) == a2;
assert khD(khPageRep(hitA)) == khMod[2,0](a);
assert khPageNext(a2) == zero[KhPage(A,2,3,0,-1,-4)];
let below = khKappaDual[A];
assert khD(below) == zero[KhFD(A,-1,-1,-2)];
assert khPair(below,khToF2[-2](khPsi[A])) == S.hi;
show hitA;

// For B the hit is already on page one. Positive stabilization changes it.
khbraid B[4,[3,-2,-2,3,3,2,-3,1,1,2,-1]];
let b = khMod[2,-4](khPsi[B]);
let b1 = khPageClass[1,-4](b);
let hitB = khPageBoundary(b1);
assert khPageD(hitB) == b1;
assert khD(khPageRep(hitB)) == khMod[2,-2](b);
assert khPageNext(b1) == zero[KhPage(B,2,2,0,-1,-4)];
show hitB;
khbraid BPlus[5,[3,-2,-2,3,3,2,-3,1,1,2,-1,4]];
let bplus = khMod[2,-5](khPsi[BPlus]);
let plus2 = khPageNext(khPageClass[1,-5](bplus));
let hitPlus = khPageBoundary(plus2);
assert khPageD(hitPlus) == plus2;
assert khD(khPageRep(hitPlus)) == khMod[2,-1](bplus);
assert khPageNext(plus2) == zero[KhPage(BPlus,2,3,0,-1,-5)];
show hitPlus;
