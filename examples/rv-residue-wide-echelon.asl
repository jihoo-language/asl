// The exact ratio 1/(2x) does not fit as a stored field element.
rvqpadic K[7,2,3];
let z = zero[RVQpQuad(K)];
let one = rvQpQuad[K,1,0,1];
let x = rvQpQuadMul(rvQpQuad[K,2,0,1],rvQpQuadPow[46](rvQpQuad[K,3,-1,1]));
let original = rvResSpecializeQuad[K,3,[46,0,0]]{[x,one,z],[z,z,one]};
let dependent = rvResSpecializeQuad[K,3,[46,0,0]]{[z,z,z],[z,z,one],[x,one,z],[x,one,z]};
assert original == dependent;
assert original == rvResSpace[7,3]{[1,4,0],[0,0,1]};
show original;
rvResSupport(original);

// In the other embedding x is a unit. The chart has the same angular ratio.
rvqpadic C[7,2,4];
let c0 = zero[RVQpQuad(C)];
let c1 = rvQpQuad[C,1,0,1];
let y = rvQpQuadMul(rvQpQuad[C,2,0,1],rvQpQuadPow[46](rvQpQuad[C,3,-1,1]));
let other = rvResSpecializeQuad[C,3,[0,0,0]]{[y,c1,c0],[c0,c0,c1]};
assert other == original;
show other;
