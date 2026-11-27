let rho = extClass(rmcword[1,0]{});
let h0 = extClass(rmcword[0,0]{([],[1])});
let h1 = extClass(rmcword[0,0]{([1],[])});
let square = yoneda(h0,h0);
let unit = extClass(rmcword[0,0]{});
assert ainfty(h0) == zero[RMExt(2,1,0)];
assert alift(h0) == extRep(h0);
assert ainfty(rho,h1) == yoneda(rho,h1);
assert ainfty(h0,unit,h1) == zero[RMExt(1,3,1)];
assert alift(h0,unit,h1) == zero[RMCobar(0,3,1)];

// With both adjacent products exact, the chosen value lies in the full coset.
let triple = ainfty(rho,h0,h1);
rmTripleContains(rmTriple(rho,h0,h1),triple);

// Both binary lifts are tau^4, and the triple coset has no indeterminacy.
let x = extClass(rmcword[0,0]{([],[4])});
let r = extClass(rmcword[4,0]{});
assert alift(x,r) == rmcword[0,4]{};
assert alift(r,x) == rmcword[0,4]{};
let tower = ainfty(x,r,x);
assert tower == extClass(rhoShift[4](rmcword[0,0]{([],[8])}));
assert tower == rmTripleBase(rmTriple(x,r,x));
rhoorder(tower);

// These nonzero higher operations use the chosen real-motivic contraction.
show ainfty(rho,h0,rho,square);
show ainfty(rho,h0,h0,h1,square);
