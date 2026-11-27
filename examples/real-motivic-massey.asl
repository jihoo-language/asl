let rho = rmcword[1,0]{};
let h0 = rmcword[0,0]{([],[1])};
let h1 = rmcword[0,0]{([1],[])};
let u = rmMasseyNull(rho,h0);
let v = rmMasseyNull(h0,h1);
assert cobarD(u) == concat(rho,h0);
assert cobarD(v) == concat(h0,h1);
let z = rmMassey(rho,h0,h1);
assert z == add(concat(u,h1),concat(rho,v));
assert cobarD(z) == zero[RMCobar(2,2,0)];

// The full coset contains two nonzero classes.
let bracket = rmTriple(extClass(rho),h0,extClass(h1));
let base = rmcword[0,0]{([],[2])};
let direction = rhoShift[2](rmcword[0,0]{([2],[])});
assert rmTripleDir[0](bracket) == extClass(direction);
rmTripleDim(bracket);
rmTripleContains(bracket,z);
rmTripleContains(bracket,base);
rmTripleContains(bracket,add(base,direction));
rmTripleContains(bracket,zero[RMExt(1,2,0)]);

// The adjacent nulls are both tau^4. Their scalar cohomology vanishes.
let x = rmcword[0,0]{([],[4])};
let r = rmcword[4,0]{};
let bounded = rhoShift[4](rmcword[0,0]{([],[8])});
assert rmMasseyNull(x,r) == rmcword[0,4]{};
assert rmMasseyNull(r,x) == rmcword[0,4]{};
assert rmMassey(x,r,x) == bounded;
let singleton = rmTriple(x,r,x);
rmTripleDim(singleton);
assert rmTripleBase(singleton) == extClass(bounded);
rhoorder(bounded);

// The sparse defining cycle still computes at exponent 32.
let x32 = rmcword[0,0]{([],[32])};
let r32 = rmcword[32,0]{};
let z32 = rmMassey(x32,r32,x32);
assert z32 == rhoShift[32](rmcword[0,0]{([],[64])});
assert rmMasseyNull(x32,r32) == rmcword[0,32]{};
assert cobarD(rhoNull[32](z32)) == rhoShift[32](z32);
rhoorder(z32);
