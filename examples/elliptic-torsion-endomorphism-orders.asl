rvunram F[13,1]{0,1};
rvec E[F](rvUnramElem[F]{4},rvUnramElem[F]{0});
rvec V[F](rvUnramElem[F]{2},rvUnramElem[F]{1});
rvunram L[13,1]{11,0,1};
rvec VL[L](rvUnramElem[L]{2},rvUnramElem[L]{1});

rvecplainlevel A[E,E,2](
    rvECPoint[E](rvUnramElem[F]{0},rvUnramElem[F]{0}),
    rvECPoint[E](rvUnramElem[F]{3},rvUnramElem[F]{0}));
rvecplainlevel B[V,VL,2](
    rvECPoint[VL](rvUnramElem[L]{2},rvUnramElem[L]{0}),
    rvECPoint[VL](rvUnramElem[L]{12,2},rvUnramElem[L]{0}));
assert rvECLevelFrob[A]==rvECLevelMat[A,1,0,0,1];
assert rvECLevelFrob[B]==rvECLevelMat[B,1,0,1,1];
rvECEndTorsionCheck[E,A];
rvECEndTorsionCheck[V,B];
rvecendtorsion M[E,A];
rvecendtorsion N[V,B];
rvECEndData[M];
rvECEndData[N];

rvunram G[37,1]{0,1};
rvec U[G](rvUnramElem[G]{1},rvUnramElem[G]{0});
let p=rvECPoint[U](rvUnramElem[G]{2},rvUnramElem[G]{11});
let q=rvECPoint[U](rvUnramElem[G]{11},rvUnramElem[G]{11});
rvecplainlevel U2[U,U,2](rvECScale[3](p),rvECScale[3](q));
rvecplainlevel U3[U,U,3](rvECScale[2](p),rvECScale[2](q));
rvecplainlevel U6[U,U,6](p,q);
rvECEndTorsionCheck[U,U2];
rvECEndTorsionCheck[U,U3];
rvECEndTorsionCheck[U,U2,U3];
rvECPlainLevelScalar[U6,2];
rvecendtorsion R[U,U6];
rvECEndData[R];
