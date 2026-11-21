qmodule X[0] { x0:0,x1:0 } {};
let one = rvU[0](rvId[X,[0],7,[1]]);

// The fixed F2 classes form 12 of the 80 F4 vertices at radius three.
rvTreeShell[3,1,3](one);
let x10 = rvgen[X,X,[0],7,[1],1];
let x01 = rvgen[X,X,[0],7,[1],2];
let a = add(rvU[0](x10),rvU[1](x01));

// The two-cycle exchanges the endpoints of the standard edge.
rvTreeShell[3,1,0](a);
rvTreeShell[3,1,3](a);
rvTreeBall[3,1,2](a);
rvTreeCount[3,1,2,[-1,2]](a);
rvTreeClosureCount[3,1,2,[-1,2]](a);

// A non-Frobenius-fixed basis of the standard lattice permutes the shell.
let basis = add(one,rvU[1](rvScale[2](x01)));
let changed = ecomp(basis,ecomp(a,rvFrob[1](basis)));
rvTreeShell[3,1,3](changed);
let twisted = ecomp(rvU[-2](rvId[X,[0],7,[1]]),a);
rvTreeShell[3,1,3](twisted);

// Relative Frobenius now has order three.
rvTreeShell[3,1,2](rvU[0](rvId[X,[0],11,[1]]));
