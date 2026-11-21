qmodule X[0] { x0:0, x1:0 } {};

// Frobenius exchanges the two distinct eigenvalues of the norm.
let a = add(rvgen[X,X,[0],7,[1],1], rvScale[2](rvgen[X,X,[0],7,[1],2]));
let n = rvFrobNorm[3,1](a);
assert n == add(rvScale[2](rvgen[X,X,[0],7,[1],0]), rvScale[3](rvgen[X,X,[0],7,[1],3]));
let c = rvShintani[3,1](a);
let expected = add(rvgen[X,X,[0],3,[1],1], add(rvgen[X,X,[0],3,[1],2], rvgen[X,X,[0],3,[1],3]));
assert c == expected;
let frame = rvShintaniFrame[3,1](a);
assert ecomp(n, frame) == ecomp(frame, rvExtend[7,1](c));
assert ecomp(frame, rvInverse(frame)) == rvId[X,[0],7,[1]];
rvInvariantFactors[3,1](a);
show c;

// A sigma-conjugate has the same canonical descent.
let s = add(rvId[X,[0],7,[1]], rvScale[3](rvgen[X,X,[0],7,[1],2]));
let a2 = ecomp(rvInverse(s), ecomp(a, rvFrob[1](s)));
assert rvFrobNorm[3,1](a2) == ecomp(rvInverse(s), ecomp(n,s));
assert rvShintani[3,1](a2) == c;

// These two norms have characteristic polynomial T^2+1.
let u = rvId[X,[0],7,[1]];
let unipotent = add(u, rvScale[2](rvgen[X,X,[0],7,[1],2]));
rvInvariantFactors[3,1](u);
rvInvariantFactors[3,1](unipotent);
let unipotentNorm = rvFrobNorm[3,1](unipotent);
assert unipotentNorm == add(u, rvgen[X,X,[0],7,[1],2]);
let unipotentFrame = rvShintaniFrame[3,1](unipotent);
assert ecomp(unipotentNorm, unipotentFrame) ==
    ecomp(unipotentFrame, rvExtend[7,1](rvShintani[3,1](unipotent)));

// Over F8/F2 this ordered norm has noncommuting factors.
let a00 = rvgen[X,X,[0],11,[1],0];
let a10 = rvgen[X,X,[0],11,[1],1];
let a01 = rvgen[X,X,[0],11,[1],2];
let a11 = rvgen[X,X,[0],11,[1],3];
let b = add(add(rvScale[2](a00), rvScale[4](a10)), add(rvScale[3](a01), a11));
assert ecomp(b, rvFrob[1](b)) == add(add(rvScale[2](a00), rvScale[2](a01)), rvScale[3](a11));
assert ecomp(rvFrob[1](b), b) == add(add(a00, rvScale[3](a10)), rvScale[2](a01));
let cubicNorm = rvFrobNorm[3,1](b);
assert cubicNorm == add(add(rvScale[3](a00), rvScale[6](a10)), add(rvScale[7](a01), rvScale[3](a11)));
assert cubicNorm == ecomp(b, ecomp(rvFrob[1](b), rvFrob[2](b)));
let cubicFrame = rvShintaniFrame[3,1](b);
assert ecomp(cubicNorm, cubicFrame) == ecomp(cubicFrame, rvExtend[11,1](rvShintani[3,1](b)));
assert ecomp(rvInverse(cubicFrame), cubicFrame) == rvId[X,[0],11,[1]];
