// A large determinant-one change of projective coordinates.
let standard = rvStableLine[7]{[1,0],[0,1],[-1,1],[1,1]};
let moved = rvStableLine[7]{
    [10000000000000000000000000,9999999999999999999999999],
    [10000000000000000000000001,10000000000000000000000000],
    [1,1],
    [20000000000000000000000001,19999999999999999999999999]};
assert standard == moved;
rvStableComponents(moved);

rvqpadic K[7,2,3];
let one = rvQpQuad[K,1,0,1];
let z = zero[RVQpQuad(K)];
let two = rvQpQuad[K,2,0,1];
let h = rvQpQuad[K,10000000000000000000000000,0,1];
let up = add(h,one);
let down = add(h,neg(one));
let movedQuad = rvStableLineQuad[K]{
    [h,down],[up,h],[one,one],[add(h,up),add(h,down)]};
assert movedQuad == standard;
rvStableComponents(movedQuad);

// The inverse of x has an oversized canonical denominator.
// The weighted line only needs its angular component.
let x = rvQpQuadPow[46](rvQpQuad[K,3,-1,1]);
show rvResSpecializeQuad[K,2,[0,0]]{[x,x]};
show rvResSpecializeQuad[K,2,[46,0]]{[x,two]};
show rvResSpecializeQuad[K,3,[0,0,0]]{[one,z,z],[z,x,rvQpQuadMul(two,x)]};

// A 2-by-2 chart has determinant one after the signed products cancel.
assert rvResSpecializeQuad[K,2,[0,0]]{[h,up],[down,h]}
    == rvResSpace[7,2]{[1,0],[0,1]};
