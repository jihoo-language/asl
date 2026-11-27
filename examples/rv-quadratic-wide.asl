// The individual squared terms exceed i128, while the norm is 7^31.
rvqpadic K[7,2,3];
let x = rvQpQuad[K,48912917532184290021,-34586655674624740457,1];
let expectedInverse = rvQpQuad[K,48912917532184290021,34586655674624740457,157775382034845806615042743];
let one = rvQpQuad[K,1,0,1];
assert x == rvQpQuadPow[31](rvQpQuad[K,3,-1,1]);
assert rvQpQuadInv(x) == expectedInverse;
assert rvQpQuadMul(x,expectedInverse) == one;
rvQpQuadVal(x);
rvQpQuadVal(rvQpQuadConj(x));
show rvResSpecializeQuad[K,2,[0,31]]{[one,x]};
rvQpQuadResidue[32](x);

// The numerator and denominator have equal valuation here.
// Only the requested precision is needed in the resulting unit.
let y = rvQpQuad[K,48912917532184290021,-34586655674624740457,157775382034845806615042743];
rvQpQuadResidue[31](y);
rvQpQuadResidue[32](y);
let scaled = rvQpQuad[K,48912917532184290021,-34586655674624740457,147];
assert rvQpQuadMul(scaled,rvQpQuadInv(scaled)) == one;
rvQpQuadVal(scaled);
show rvResSpecializeQuad[K,2,[0,29]]{[one,scaled]};
rvqpadic Other[7,2,4];
let other = rvQpQuad[Other,48912917532184290021,-34586655674624740457,147];
rvQpQuadVal(other);
show rvResSpecializeQuad[Other,2,[0,-2]]{[rvQpQuad[Other,1,0,1],other]};

// This noncoordinate flag requires a wide exact pivot over K.
let p = rvQpQuad[K,7,0,1];
let z = zero[RVQpQuad(K)];
rvphodgequad Flag[K] {
    {x,rvQpQuadMul(p,x)}, {(0,1,one)},
    {(33,[one,x]),(30,[one,z])}
};
rvPHodgeData[Flag,[0]];
rvPHodgeWeak[Flag];
