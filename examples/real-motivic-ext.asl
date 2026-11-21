// h0 and rho h1 lie in the same real-motivic grade.
let h0 = extClass(rmcword[0,0]{([],[1])});
let h1 = extClass(rmcword[0,0]{([1],[])});
let shifted = rhoShift[1](h1);
assert add(h0,shifted) == add(rmextgen[1,1,0,0],rmextgen[1,1,0,1]);
rmextdim[1,1,0];
rhoorder(h0);
rhoorder(add(h0,shifted));
assert rhoShift[1](h0) == zero[RMExt(1,0,-1)];
assert rhoProject(h0) == zero[Ext(1,1)];
assert rhoProject(shifted) == extgen[1,1,0];
assert realComplex(h0) == mextgen[1,1,0,0];
assert realComplex(shifted) == zero[MExt(1,1,0)];

// An incoming-only solve certifies the non-scalar torsion class h0^2.
let square = yoneda(h0,h0);
rhoorder(square);
let null = rhoNull[1](square);
assert null == rmcword[0,1]{([],[1])};
assert cobarD(null) == rhoShift[1](extRep(square));

// A direct scalar primitive and its absence at positive topological
// degree certify the exact order of the large raw class.
let large = rmcword[0,0]{([],[32])};
rhoorder(large);
let witness = rhoNull[32](large);
assert witness == rmcword[0,32]{};
assert cobarD(witness) == rhoShift[32](large);
show witness;
