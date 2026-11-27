// The input factors are ascending and are certified against the actual operator.
rvunram F23[23,1]{-5,0,1};
rvfrickess B23[F23];
rvfrickeheckeorder A[B23,0]{[-4,1],[0,1],[-5,0,1]};
rvFrickeOrderData[A];
rvFrickeOrderCongruence[A,0,2,11];
rvFrickeOrderCongruence[A,0,2,5];

let one=rvFrickeOrderBasis[A,0,0];
assert rvFrickeOrderIntegral(rvFrickeOrderProject[A,0](one))==one;
assert rvFrickeOrderIntegral(rvFrickeOrderProject[A,1](one))==zero[RVFrickeSSDiv(B23)];
show rvFrickeOrderProject[A,0](rvFrickeSSBasis[B23,0]);

rvunram F41[41,1]{-3,0,1};
rvfrickess B41[F41];
rvfrickeheckeorder B[B41,0]{[-4,1],[2,1],[2,-4,0,1]};
rvFrickeOrderData[B];
rvFrickeOrderCongruence[B,0,2,5];
rvFrickeOrderCongruence[B,1,2,2];

let anti=rvFrickeOrderBasis[B,1,0];
assert rvFrickeSSFrob(anti)==neg(anti);
assert rvFrickeOrderIntegral(rvFrickeOrderProject[B,1](anti))==anti;
show rvFrickeOrderProject[B,1](rvFrickeSSBasis[B41,1]);

// B3+2B5 separates these blocks, while B3+B5 has a repeated eigenvalue.
rvfrickeheckeorder D[B23,2]{[-16,1],[-8,1],[-1,4,1]};
rvFrickeOrderData[D];
rvfrickeheckeorder C[B23,1]{[-10,1],[-4,1],[1,1]};
rvFrickeOrderData[C];
