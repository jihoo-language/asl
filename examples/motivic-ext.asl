// The classical generators in topological degrees one, two, and four have these motivic lifts.
let h0 = extClass(mcword[0]{([],[1])});
let h1 = extClass(mcword[0]{([1],[])});
let h2 = extClass(mcword[0]{([2],[])});
assert yoneda(h0,h1) == zero[MExt(2,3,1)];
assert yoneda(h1,h0) == zero[MExt(2,3,1)];

let square = yoneda(h1,h1);
let cube = yoneda(h1,square);
let cubic = yoneda(yoneda(h0,h0),h2);
assert tauShift[1](cube) == cubic;
assert cubic == mextgen[3,6,2,0];
show cubic;

let fourth = yoneda(square,square);
assert fourth == mextgen[4,8,4,0];
assert extClass(extRep(fourth)) == fourth;
assert tauShift[1](fourth) == zero[MExt(4,8,3)];
assert realize(fourth) == zero[Ext(4,8)];
show fourth;

// Both a computed primitive and this five-term primitive certify tau h1^4=0.
let killed = tauShift[1](extRep(fourth));
assert cobarD(extPrimitive(killed)) == killed;
let a = mcword[0]{([],[0,1]),([],[0,1]),([1],[])};
let b = mcword[0]{([1],[]),([],[1,1]),([1],[])};
let c = mcword[0]{([1],[0,1]),([],[1]),([1],[])};
let d = mcword[0]{([2],[]),([],[1]),([],[0,1])};
let e = mcword[0]{([2],[]),([],[1]),([1],[1])};
let witness = add(add(a,b),add(c,add(d,e)));
assert cobarD(witness) == killed;
mextdim[4,8,4];
mextdim[4,8,3];
