// The characteristic polynomial of Legendre lambda=2 over F5 is T^2+2T+5.
// Choose i=2 mod 5 and lift along that branch of Q(i).
rvqpadic K[5,-1,2];
let zeroK = zero[RVQpQuad(K)];
let one = rvQpQuad[K,1,0,1];
let alpha = rvQpQuad[K,-1,2,1];
let beta = rvQpQuadConj(alpha);
assert rvQpQuadMul(alpha,beta) == rvQpQuad[K,5,0,1];
assert rvQpQuadInv(beta) == rvQpQuad[K,-1,2,5];
rvQpQuadVal(alpha);
rvQpQuadVal(beta);
rvQpQuadResidue[8](alpha);
legendreDwork[5,2,8];
let lift = rvQpQuadPow[8](beta);
rvQpQuadVal(lift);
rvQpQuadResidue[8](lift);
rvQpQuadResidue[9](lift);

// These are chosen filtrations on the same Frobenius, with N=0.
rvphodgequad Axis[K] {
    {alpha,beta}, {}, {(1,[one,zeroK]),(0,[zeroK,one])}
};
rvphodgequad Tilted[K] {
    {alpha,beta}, {}, {(1,[one,lift]),(0,[zeroK,one])}
};
rvPHodgeWeak[Axis];
rvPHodgeHN[Axis];
rvPHodgeWeak[Tilted];
rvPHodgeHN[Tilted];

// The same rational Frobenius roots at the conjugate embedding.
rvqpadic Conjugate[5,-1,3];
rvQpQuadVal(rvQpQuad[Conjugate,-1,2,1]);
rvQpQuadVal(rvQpQuad[Conjugate,-1,-2,1]);
