rvunram L[5,1]{3,0,1};
let z = rvUnramElem[L]{};
let one = rvUnramElem[L]{1};
let alpha = rvUnramElem[L]{0,1};
rvhyperelliptic C[L]{z,one,z,z,z,one};
rvrichelot D[C]{[2],[],[1]}{[-2],[],[1]}{[],[1]};
rvjacchart N[D](z);

// The first two classes and the last two classes form the two halves of a level-3 basis.
let gen0 = rvJacScale[2](rvJacPoint[C](one,alpha));
let gen1 = add(gen0,rvJacScale[2](rvJacPoint[C](rvUnramElem[L]{1,2},alpha)));
let gen2 = rvJacScale[4](rvJacPoint[C](rvUnramElem[L]{1,1},alpha));
let gen3 = rvJacScale[4](rvJacPoint[C](rvUnramElem[L]{1,3},alpha));
let zeta = rvUnramElem[L]{2,2};
rvjactorsionbasis B3[C,3](zeta){gen0,gen1,gen2,gen3};
let im0 = rvRichelotPush[N](gen0);
let im1 = rvRichelotPush[N](gen1);
let im2 = rvRichelotPush[N](gen2);
let im3 = rvRichelotPush[N](gen3);

// In the image basis, the primitive root is zeta squared and the matrix is the identity.
rvjactorsionbasis H3[N,3](rvUnramPow[2](zeta)){im0,im1,im2,im3};
rvJacRichelotMatrix[B3,H3];
// With the same primitive root on both curves, the matrix has diagonal (1,1,2,2).
rvjactorsionbasis K3[N,3](zeta){im0,im1,rvJacScale[2](im2),rvJacScale[2](im3)};
rvJacRichelotMatrix[B3,K3];

// The branch points with x=0, alpha, -alpha, and 2 alpha give level-2 bases.
let twiceAlpha = rvUnramElem[L]{0,2};
let t0 = rvJacPoint[C](z,z);
let t1 = rvJacPoint[C](alpha,z);
let t2 = rvJacPoint[C](neg(alpha),z);
let t3 = rvJacPoint[C](twiceAlpha,z);
let q0 = rvJacPoint[N](z,z);
let q1 = rvJacPoint[N](alpha,z);
let q2 = rvJacPoint[N](neg(alpha),z);
let q3 = rvJacPoint[N](twiceAlpha,z);
let e2 = add(add(t0,t1),t2);
rvjactorsionbasis B2[C,2](neg(one)){t0,e2,t1,add(add(t0,t1),t3)};
rvjactorsionbasis H2[N,2](neg(one)){q0,add(add(q0,q1),q2),q1,add(add(q0,q1),q3)};
assert rvRichelotPush[N](t0) == zero[RVJac(N)];
assert rvRichelotPush[N](e2) == zero[RVJac(N)];
assert rvRichelotPush[N](add(t0,e2)) == zero[RVJac(N)];
rvJacRichelotMatrix[B2,H2];
