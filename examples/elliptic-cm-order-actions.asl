rvcmclasspoly H8[-8];
rvunram F[17,1]{0,1};
rvec E[F](rvUnramElem[F]{1},rvUnramElem[F]{7});
rvecendorder M[E,H8];
rveccmaction FrobeniusOnly[M];
rvECCMActionData[FrobeniusOnly];
rveccmcycle K[M,3,1];
rveccmaction Full[M,K,K];
rvECCMActionData[Full];

let w=rvECCMActionGenerator[Full];
assert w==rvECEndElem[M,0,1];
rvECMapDegree(rvECCMActionMap[Full](w));
assert rvECCMActionGenerator[FrobeniusOnly]==rvECEndElem[M,0,2];
assert (\p:RVECPoint(E)=>rvECCMAct[Full](w,rvECCMAct[Full](w,p)))
    ==(\p:RVECPoint(E)=>rvECScale[-2](p));

rvunram L[17,1]{14,0,1};
rvec EL[L](rvUnramElem[L]{1},rvUnramElem[L]{7});
rveccmactionbase OL[FrobeniusOnly,EL];
rveccmactionbase AL[Full,EL];
rvECCMActionData[AL];
let p=rvECPoint[EL](rvUnramElem[L]{0,1},rvUnramElem[L]{2,1});
let wp=rvECPoint[EL](rvUnramElem[L]{13,6},rvUnramElem[L]{7,9});
let fp=rvECPoint[EL](rvUnramElem[L]{0,16},rvUnramElem[L]{2,16});
assert rvECCMAct[AL](w,p)==wp;
assert rvECCMAct[AL](w,wp)==rvECScale[-2](p);
assert rvECCMAct[OL](rvECEndFrobenius[M],p)==fp;
assert rvECCMAct[AL](rvECEndFrobenius[M],p)==fp;
