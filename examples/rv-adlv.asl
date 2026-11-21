qmodule X[0] { x0:0,x1:0 } {};
let x00 = rvgen[X,X,[0],3,[1],0];
let x10 = rvgen[X,X,[0],3,[1],1];
let x01 = rvgen[X,X,[0],3,[1],2];
let x11 = rvgen[X,X,[0],3,[1],3];
let one = rvU[0](rvId[X,[0],3,[1]]);
let a = add(rvU[0](x10),rvU[1](x01));
let thick = add(rvU[2](x00),rvU[0](x11));

// Newton slopes belong to the isocrystal. Hodge exponents depend on its lattice.
rvNewton[3,1](a);
rvHodge(a);
rvPosition(one,thick);
rvMuAt[3,1](a,one);
rvMuAt[3,1](a,thick);
assert rvInADLV[3,1,[0,1]](a,one) == S.hi;
assert rvInADLV[3,1,[0,1]](a,thick) == S.lo;
assert rvInADLV[3,1,[-1,2]](a,thick) == S.hi;
assert rvInADLVClosure[3,1,[-1,2]](a,one) == S.hi;
assert rvInADLVClosure[3,1,[-1,2]](a,thick) == S.hi;
assert rvInADLVClosure[3,1,[0,1]](a,thick) == S.lo;
assert rvInADLVClosure[3,1,[-2,2]](a,one) == S.lo;

// A basis change of determinant 1+u is a unit over the valuation ring.
let unitFrame = add(one,rvU[1](x00));
rvPosition(unitFrame,one);
let changedBasis = ecomp(thick,unitFrame);
rvMuAt[3,1](a,changedBasis);
assert rvInADLV[3,1,[-1,2]](a,changedBasis) == S.hi;

// The minimum valuation of an entry is zero. The determinant is u by cancellation.
let allOnes = add(add(x00,x10),add(x01,x11));
let cancel = add(rvU[0](allOnes),rvU[1](x11));
rvHodge(cancel);
