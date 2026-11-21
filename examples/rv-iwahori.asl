qmodule X[0] { x0:0,x1:0 } {};
let unit = rvU[0](rvId[X,[0],3,[1]]);
let e10 = rvgen[X,X,[0],3,[1],1];
let e01 = rvgen[X,X,[0],3,[1],2];
let finite = rvU[0](add(e10,e01));
let affine = add(rvU[-1](e10),rvU[1](e01));

// Both matrices preserve the standard lattice. Only the identity preserves its flag.
rvHodge(unit);
rvHodge(finite);
rvIw(unit);
rvIw(finite);
rvIwLength(finite);
rvIw(affine);
rvIwLength(affine);
assert rvInIwADLV[3,1,[1,0]](finite,unit) == S.hi;
assert rvInIwADLV[3,1,[0,1]](finite,unit) == S.lo;
assert rvInIwADLVClosure[3,1,[1,0]](unit,unit) == S.hi;
assert rvInIwADLVClosure[3,1,[-1,2]](unit,unit) == S.hi;
assert rvInIwADLVClosure[3,1,[-1,2]](finite,unit) == S.lo;
assert rvInIwADLVClosure[3,1,[1,0]](affine,unit) == S.lo;

// Lower triangular modulo u gives the stabilizer of this descending chain.
let lower = add(unit,rvU[0](e10));
let upper = add(unit,rvU[0](e01));
rvIw(lower);
rvIw(upper);
let changedBasis = add(unit,rvU[1](rvgen[X,X,[0],3,[1],0]));
rvIwPosition(changedBasis,finite);

// Over F4 the coefficient Frobenius changes an upper triangular unit frame.
let largeUnit = rvU[0](rvId[X,[0],7,[1]]);
let x = add(largeUnit,rvU[0](rvScale[2](rvgen[X,X,[0],7,[1],2])));
rvIwAt[3,1](largeUnit,x);
assert rvInIwADLV[3,1,[1,0]](largeUnit,x) == S.hi;
assert rvInIwADLVClosure[3,1,[0,1]](largeUnit,x) == S.lo;
