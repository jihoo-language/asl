let h0 = extClass(rmcword[0,0]{([],[1])});
let t2 = extClass(rmcword[0,0]{([],[2])});
let t4 = extClass(rmcword[0,0]{([],[4])});
let h1 = extClass(rmcword[0,0]{([1],[])});

let product = yoneda(h0,t2);
let coordinates = rmbarCoords(product);
assert rmbarClass(coordinates) == product;
assert extClass(rmbarLift(coordinates)) == product;
show coordinates;
rhodepth(h0);
rhodepth(t2);
rhodepth(coordinates);
rhoorder(product);
rhohidden(h0,t2);

let quotient = rhoDivide[1](product);
assert rhoShift[1](quotient) == product;
let raw = concat(extRep(h0),extRep(t2));
let reconstructed = rhoShift[1](rmbarLift(rmbarCoords(quotient)));
let witness = rhoProductNull[1](h0,t2);
assert witness == add(rmcword[0,1]{([1],[1])},rmcword[0,1]{([],[0,1])});
assert cobarD(witness) == add(raw,reconstructed);

let a = rmcword[0,0]{([],[1]),([],[0,1])};
let b = rmcword[0,0]{([1],[1]),([],[1])};
let c = rmcword[0,0]{([],[0,1]),([],[1])};
assert rhoShift[1](extClass(add(a,add(b,c)))) == product;

// Use a common stable page to compare the graded and total products.
let leading0 = rhoEinf[4,0](h0);
let leading2 = rhoEinf[4,0](t2);
assert rhoMul(leading0,leading2) == rhoEinf[4,0](product);
assert rhoEinf[4,0](product) == zero[RhoBock(4,2,3,0,0)];
let leading_product = rhoEinf[1](product);
assert leading_product == rhobockgen[4,2,3,0,1,0];
assert extClass(rhoRep(leading_product)) == product;
show leading_product;

// Two further depth-two products and their checked divisors.
rhohidden(h1,t4);
rhohidden(t2,t4);
let p1 = yoneda(h1,t4);
assert rhoShift[2](rhoDivide[2](p1)) == p1;
let p = yoneda(t2,t4);
let q = rhoDivide[2](p);
assert rhoShift[2](q) == p;
assert cobarD(rhoProductNull[2](t2,t4)) == add(concat(extRep(t2),extRep(t4)),rhoShift[2](rmbarLift(rmbarCoords(q))));
