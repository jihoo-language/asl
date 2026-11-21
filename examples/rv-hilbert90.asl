qmodule X[0] { x0:0, x1:0 } {};
qmodule Y[0] { y0:0 } {};

// Over F4/F2 this unipotent matrix has norm one.
let id = rvId[X,[0],7,[1]];
let upper = rvgen[X,X,[0],7,[1],2];
let a = add(id, upper);
let frame = rvHilbert90[3,1](a);
assert frame == add(id, rvScale[3](upper));
assert ecomp(a, rvFrob[1](frame)) == frame;
assert ecomp(rvInverse(frame), frame) == id;
show frame;

// Descent respects composition even for maps between different ranks.
let trivial = rvId[Y,[0],7,[1]];
let f = rvgen[X,Y,[0],3,[1],1];
let g = rvgen[Y,X,[0],3,[1],1];
let largeF = rvTwistedExtend[3,1](trivial, f, a);
let largeG = rvTwistedExtend[3,1](a, g, trivial);
assert ecomp(trivial, rvFrob[1](largeF)) == ecomp(largeF, a);
assert ecomp(a, rvFrob[1](largeG)) == ecomp(largeG, trivial);
assert rvTwistedDescend[3,1](trivial, largeF, a) == f;
assert rvTwistedDescend[3,1](a, largeG, trivial) == g;
assert ecomp(largeG, largeF) == rvTwistedExtend[3,1](a, ecomp(g, f), a);
assert rvTwistedDescend[3,1](a, ecomp(largeG, largeF), a) == ecomp(g, f);
assert rvTwistedDescend[3,1](a, rvId[X,[0],7,[1]], a) == rvId[X,[0],3,[1]];
show rvTwistedDescend[3,1](trivial, largeF, a);

// Over F8/F2 the first two cocycle factors do not commute.
let x00 = rvgen[X,X,[0],11,[1],0];
let x10 = rvgen[X,X,[0],11,[1],1];
let x01 = rvgen[X,X,[0],11,[1],2];
let x11 = rvgen[X,X,[0],11,[1],3];
let cubic = add(add(rvScale[7](x00), rvScale[3](x10)), add(rvScale[5](x01), rvScale[5](x11)));
let firstTwo = ecomp(cubic, rvFrob[1](cubic));
assert firstTwo == add(add(rvScale[5](x00), rvScale[2](x10)), add(rvScale[5](x01), rvScale[4](x11)));
assert ecomp(rvFrob[1](cubic), cubic) == add(add(rvScale[4](x10), rvScale[2](x01)), x11);
assert ecomp(firstTwo, rvFrob[2](cubic)) == rvId[X,[0],11,[1]];
let cubicFrame = rvHilbert90[3,1](cubic);
assert ecomp(cubic, rvFrob[1](cubicFrame)) == cubicFrame;
assert ecomp(cubicFrame, rvInverse(cubicFrame)) == rvId[X,[0],11,[1]];
show cubicFrame;
