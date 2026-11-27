// A homogeneous cell with both annular degrees of the integral differential.
khbraid Mixed[3,[2,1,-2]];
let x = khAnnular[1](khTerm[Mixed,[],[1]]);
let d0 = khD0(x);
let dm = khDMinus2(x);
assert d0 == khAnnular[1](khTerm[Mixed,[0],[1,2]]);
assert dm == khAnnular[-1](add(khTerm[Mixed,[1],[0]],khTerm[Mixed,[2],[1,2]]));
assert khD(x) == add(khForgetK(d0),khForgetK(dm));
assert khD0(d0) == zero[KhAC(Mixed,1,0,1)];
assert khDMinus2(dm) == zero[KhAC(Mixed,1,0,-3)];
let mixed = khDMinus2(d0);
assert mixed == khAnnular[-1](khTerm[Mixed,[0,1],[0,1],-1]);
assert add(mixed,khD0(dm)) == zero[KhAC(Mixed,1,0,-1)];
khAnnularCheck[Mixed];

// The same unknot in S^3 can have different annular embeddings.
khbraid U[1,[]];
khbraid P[3,[1,2]];
khbraid N[2,[-1]];
khEuler[U];
khEuler[P];
khAnnularEuler[U];
khAnnularEuler[P];
khAnnularGroup[P,0,5,3];
khAnnularGroup[U,0,5,3];
khKappa[U];
khKappa[P];
khKappa[N];
let trivial = khKappaNull[N];
assert trivial == khToF2[0](khTerm[N,[],[0]]);
assert khD(trivial) == khToF2[0](khPsi[N]);

// The four-braids A and B have different filtered boundary thresholds.
khbraid A[4,[3,-2,-2,3,3,2,-3,-1,2,1,1]];
khbraid B[4,[3,-2,-2,3,3,2,-3,1,1,2,-1]];
khKappa[A];
khKappa[B];
let a = khPsi[A];
let nullA = khKappaNull[A];
let belowA = khKappaDual[A];
assert khD(nullA) == khToF2[0](a);
assert khD(belowA) == zero[KhFD(A,-1,-1,-2)];
assert khPair(belowA,khToF2[-2](a)) == S.hi;
let b = khPsi[B];
let nullB = khKappaNull[B];
let belowB = khKappaDual[B];
assert khD(nullB) == khToF2[-2](b);
assert khD(belowB) == zero[KhFD(B,-1,-1,-4)];
assert khPair(belowB,b) == S.hi;

// Positive stabilization can raise kappa, as it does for B.
khbraid BPlus[5,[3,-2,-2,3,3,2,-3,1,1,2,-1,4]];
khKappa[BPlus];
let stabilized = khPsi[BPlus];
let nullPlus = khKappaNull[BPlus];
let belowPlus = khKappaDual[BPlus];
assert khD(nullPlus) == khToF2[-1](stabilized);
assert khD(belowPlus) == zero[KhFD(BPlus,-1,-1,-3)];
assert khPair(belowPlus,khToF2[-3](stabilized)) == S.hi;
