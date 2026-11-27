rvqpadic K[5,-1,2];
let zeroK = zero[RVQpQuad(K)];
let one = rvQpQuad[K,1,0,1];
let beta = rvQpQuad[K,-1,-2,1];
let delta = rvQpQuadPow[19](beta);
let tilted = add(one,delta);
rvQpQuadResidue[19](delta);
rvQpQuadResidue[20](delta);

let M = rvValMatQuad[K,4]{[one,zeroK,one,one],[zeroK,one,one,tilted]};
rvValPlucker[[2,3]](M);
rvValInitial[[0,0,0,0]](M);
rvValInitial[[0,0,19,19]](M);
rvValExchange[[0,1],[2,3],0](M);
rvValTie[[0],[1,2,3]](M);
rvValInitial[[0,0,-19,-19]](rvValDual(M));
assert rvValDual(rvValDual(M)) == M;

// Replacing a row by a nonzero multiple changes every minor by the same factor.
let scaled = rvValMatQuad[K,4]{[beta,zeroK,beta,beta],[zeroK,one,one,tilted]};
assert scaled == M;
rvphodgequad D[K] {
    {one,rvQpQuad[K,2,0,1],rvQpQuad[K,5,0,1],rvQpQuad[K,10,0,1]},
    {},
    {(1,[one,zeroK,one,one]),(1,[zeroK,one,one,tilted]),
     (0,[zeroK,zeroK,one,zeroK]),(0,[zeroK,zeroK,zeroK,one])}
};
assert rvValFil[D,1] == M;
rvValSupport(rvValFil[D,2]);
rvValSupport(rvValFil[D,0]);
