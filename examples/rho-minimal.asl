let x6 = rhobockgen[6,1,2,-3,0,0];
let z6 = rhoRep(x6);
let p6 = rhoMinP(z6);
let d6 = rhoMinD(p6);
assert d6 == rhoMinP(cobarD(z6));
assert rhoMinP(rhoMinI(p6)) == p6;
assert rhoClass[6,0](rhoMinI(p6)) == x6;
assert rhoClass[6,6](rhoMinI(d6)) == rhoD(x6);
show d6;
show rhoD(x6);

let x7 = rhobockgen[7,2,4,-2,0,0];
let z7 = rhoRep(x7);
let p7 = rhoMinP(z7);
let d7 = rhoMinD(p7);
assert d7 == rhoMinP(cobarD(z7));
assert rhoMinD(d7) == zero[RhoMin(4,4,-2)];
assert rhoClass[7,7](rhoMinI(d7)) == rhoD(x7);
show rhoD(x7);

// A correction at filtration seven cancels the component at eight.
let a = rhoMin[3](mextgen[3,6,-1,0]);
let b = rhoMin[7](mextgen[3,10,3,0]);
let lo = rhoMin[8](mextgen[4,11,4,0]);
let hi = rhoMin[10](mextgen[4,13,6,0]);
assert rhoMinD(a) == add(lo,hi);
assert rhoMinD(b) == lo;
let corrected = add(a,b);
assert rhoMinD(corrected) == hi;
let source = rhoClass[7,3](rhoMinI(corrected));
assert rhoD(source) == rhoClass[7,10](rhoMinI(hi));
assert rhoD(source) == rhobockgen[7,4,3,-4,10,0];
show rhoD(source);

let tau32 = rhoMin[0](mextgen[0,0,-32,0]);
let arrow32 = rhoMinD(tau32);
let raw32 = rmcword[0,0]{([],[32])};
assert rhoMinD(arrow32) == zero[RhoMin(2,0,-32)];
assert rhoMinI(arrow32) == cobarD(rmcword[0,32]{});
assert rhoMinI(arrow32) == rhoShift[32](raw32);
assert rhoMinP(rhoMinI(arrow32)) == arrow32;
rhoorder(raw32);
show arrow32;
show rhoMinCoef[32](arrow32);
