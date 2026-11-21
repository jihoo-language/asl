// In each tridegree the tau power supports its first nonzero rho differential.
let x0 = rhoClass[1,0](rmcword[0,1]{});
let y0 = rhoD(x0);
assert y0 == rhobockgen[1,1,0,-1,1,0];
assert rhoClass[0,1](cobarD(rhoRep(x0))) == rhoClass[0,1](rmcword[1,0]{([],[1])});
assert rhoNext(y0) == zero[RhoBock(2,1,0,-1,1)];
assert rhoD(y0) == zero[RhoBock(1,2,0,-1,2)];
show y0;

let x1 = rhoClass[2,0](rmcword[0,2]{});
let y1 = rhoD(x1);
assert y1 == rhobockgen[2,1,0,-2,2,0];
assert rhoClass[0,2](cobarD(rhoRep(x1))) == rhoClass[0,2](rmcword[2,1]{([1],[])});
assert rhoNext(y1) == zero[RhoBock(3,1,0,-2,2)];
assert rhoD(y1) == zero[RhoBock(2,2,0,-2,4)];
show y1;

let before1 = rhoClass[1,0](rmcword[0,2]{});
assert rhoD(before1) == zero[RhoBock(1,1,0,-2,1)];
assert rhoNext(before1) == x1;
assert rhoNext(rhoPrev(x1)) == x1;

let x2 = rhoClass[4,0](rmcword[0,4]{});
let y2 = rhoD(x2);
assert y2 == rhobockgen[4,1,0,-4,4,0];
assert rhoClass[0,4](cobarD(rhoRep(x2))) == rhoClass[0,4](rmcword[4,2]{([2],[])});
assert rhoNext(y2) == zero[RhoBock(5,1,0,-4,4)];
assert rhoD(y2) == zero[RhoBock(4,2,0,-4,8)];
show y2;

let before2 = rhoClass[3,0](rmcword[0,4]{});
assert rhoD(before2) == zero[RhoBock(3,1,0,-4,3)];
assert rhoNext(before2) == x2;
assert rhoNext(rhoPrev(x2)) == x2;

let x3 = rhoClass[8,0](rmcword[0,8]{});
let y3 = rhoD(x3);
assert y3 == rhobockgen[8,1,0,-8,8,0];
assert rhoClass[0,8](cobarD(rhoRep(x3))) == rhoClass[0,8](rmcword[8,4]{([4],[])});
assert rhoNext(y3) == zero[RhoBock(9,1,0,-8,8)];
assert rhoD(y3) == zero[RhoBock(8,2,0,-8,16)];
show y3;

let before3 = rhoClass[7,0](rmcword[0,8]{});
assert rhoD(before3) == zero[RhoBock(7,1,0,-8,7)];
assert rhoNext(before3) == x3;
assert rhoNext(rhoPrev(x3)) == x3;

let x4 = rhoClass[16,0](rmcword[0,16]{});
let y4 = rhoD(x4);
assert y4 == rhobockgen[16,1,0,-16,16,0];
assert rhoClass[0,16](cobarD(rhoRep(x4))) == rhoClass[0,16](rmcword[16,8]{([8],[])});
assert rhoNext(y4) == zero[RhoBock(17,1,0,-16,16)];
assert rhoD(y4) == zero[RhoBock(16,2,0,-16,32)];
show y4;

let before4 = rhoClass[15,0](rmcword[0,16]{});
assert rhoD(before4) == zero[RhoBock(15,1,0,-16,15)];
assert rhoNext(before4) == x4;
assert rhoNext(rhoPrev(x4)) == x4;

// A surviving E1 class can need a corrected cochain representative.
let lead = rmcword[2,3]{([1],[])};
let x = rhoClass[1,2](lead);
let z = add(lead,add(rmcword[3,2]{([1],[1])},rmcword[3,2]{([],[0,1])}));
assert rhoD(x) == zero[RhoBock(1,2,0,-4,3)];
assert rhoNext(x) == rhoClass[2,2](z);
assert rhoRep(rhoNext(x)) == z;
assert rhoNext(rhoPrev(rhoNext(x))) == rhoNext(x);
