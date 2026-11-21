qmodule X[0] { x0:0, x1:0 } {};
qmodule Y[0] { y0:0, y1:0, y2:0 } {};

// The ordinary isocrystal has eigenvalues 1 and u.
let x00 = rvgen[X,X,[0],3,[1],0];
let x11 = rvgen[X,X,[0],3,[1],3];
let ordinary = add(rvU[0](x00), rvU[1](x11));
rvNewton[3,1](ordinary);
let extended = add(rvU[0](rvExtend[7,1](x00)), rvU[1](rvExtend[7,1](x11)));
rvNewton[3,1](extended);

// Relative Frobenius over F4/F2 exchanges the coefficients omega and omega^2.
let a = add(rvU[0](rvgen[X,X,[0],7,[1],1]), rvU[1](rvScale[2](rvgen[X,X,[0],7,[1],2])));
let norm = add(rvU[1](rvScale[2](rvgen[X,X,[0],7,[1],0])), rvU[1](rvScale[3](rvgen[X,X,[0],7,[1],3])));
assert rvUNorm[3,1](a) == norm;
assert ecomp(a,rvFrob[1](a)) == norm;
let unit = rvId[One,[0],3,[1]];
assert rvNewtonCoeff[3,1,0](a) == rvU[2](unit);
assert rvNewtonCoeff[3,1,1](a) == rvU[1](unit);
assert rvNewtonCoeff[3,1,2](a) == rvU[0](unit);
rvNewton[3,1](a);

// This Laurent shear is its own inverse. Sigma-conjugacy preserves the polynomial.
let s = add(rvU[0](rvId[X,[0],7,[1]]), rvU[-1](rvgen[X,X,[0],7,[1],2]));
let changed = ecomp(s,ecomp(a,rvFrob[1](s)));
assert rvUNorm[3,1](changed) == ecomp(s,ecomp(norm,s));
assert rvNewtonCoeff[3,1,0](changed) == rvU[2](unit);
assert rvNewtonCoeff[3,1,1](changed) == rvU[1](unit);

// The constant terms in det([[1,1],[1,1+u]]) cancel.
let allOnes = add(add(x00,rvgen[X,X,[0],3,[1],1]), add(rvgen[X,X,[0],3,[1],2],x11));
let cancel = add(rvU[0](allOnes), rvU[1](x11));
assert rvNewtonCoeff[3,1,0](cancel) == rvU[1](unit);
assert rvNewtonCoeff[3,1,1](cancel) == rvU[1](unit);
rvNewton[3,1](cancel);

// The three-cycle has characteristic polynomial T^3+u.
let cycle = add(rvU[0](add(rvgen[Y,Y,[0],3,[1],1],rvgen[Y,Y,[0],3,[1],5])), rvU[1](rvgen[Y,Y,[0],3,[1],6]));
assert rvNewtonCoeff[3,1,0](cycle) == rvU[1](unit);
assert rvNewtonCoeff[3,1,1](cycle) == zero[RVLaurent(One,One,[0],3,[1])];
rvNewton[3,1](cycle);
let twisted = ecomp(rvU[-2](rvId[Y,[0],3,[1]]),cycle);
assert rvNewtonCoeff[3,1,0](twisted) == rvU[-5](unit);
rvNewton[3,1](twisted);
