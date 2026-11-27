// A genus-two ordinary curve over F27 in its named Witt coefficient ring.
rvunram W[3,3]{1,-1,0,1};
let f0 = rvUnramElem[W]{0,1,2};
let f1 = rvUnramElem[W]{1,2,0};
let f2 = rvUnramElem[W]{1,1,0};
let f3 = rvUnramElem[W]{2,2,1};
let f4 = rvUnramElem[W]{2};
let f5 = rvUnramElem[W]{1};
rvhyperelliptic C[W]{f0,f1,f2,f3,f4,f5};
let A = rvHWW[C,3];
let B = rvHWW[C,2];
let M = rvHWWRatio[C];
assert M == rvHWWCompose(A,rvHWWInv(rvHWWFrob[1](B)));
assert rvHWWEntry[1,1](M) == rvUnramElem[W]{16,4,3};
assert rvHWWEntry[1,1](rvHWWCompose(A,rvHWWInv(B))) == rvUnramElem[W]{26,13,6};

let M1 = rvHWWFrob[1](M);
let M2 = rvHWWFrob[2](M);
let P = rvHWWUnit[C];
assert P == rvHWWCompose(rvHWWCompose(M,M1),M2);
assert rvHWWFrob[1](P) == rvHWWCompose(rvHWWCompose(rvHWWInv(M),P),M);
assert rvHWWEntry[1,1](P) == rvUnramElem[W]{3,3,7};
assert rvHWWChar[1](P) == rvUnramElem[W]{25};
assert rvHWWChar[2](P) == rvUnramElem[W]{7};
let reverse = rvHWWCompose(rvHWWCompose(M2,M1),M);
assert rvHWWChar[1](reverse) == rvUnramElem[W]{10};
assert rvHWWChar[2](reverse) == rvUnramElem[W]{7};
show P;
rvHWWPoints[C,1];
rvHWWPoints[C,2];
rvHWWWeil[C];
rvHWWVerify[C](P);

// Change an integral lift while keeping the same residue curve.
rvhyperelliptic Lift[W]{add(f0,rvUnramElem[W]{3}),f1,f2,f3,f4,f5};
let changed = rvHWWUnit[Lift];
assert rvHWWEntry[1,1](rvHWWRatio[Lift]) == rvUnramElem[W]{25,4,3};
assert rvHWWEntry[1,1](changed) == rvUnramElem[W]{21,3,16};
assert rvHWWChar[1](changed) == rvHWWChar[1](P);
assert rvHWWChar[2](changed) == rvHWWChar[2](P);
rvHWWVerify[Lift];
