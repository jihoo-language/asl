// Each factor gives the xi powers and then the tau_i powers.
let tau0 = mcword[0]{([],[1])};
let tau1 = mcword[0]{([],[0,1])};
let xi1 = mcword[0]{([1],[])};
let xi2 = mcword[0]{([0,1],[])};

// Normalizing a square contributes a global coefficient tau.
assert mcword[0]{([],[2])} == tauShift[1](xi1);
assert mcword[0]{([],[0,2])} == tauShift[1](xi2);
assert mcword[0]{([],[2]),([],[0,2])} == tauShift[2](concat(xi1,xi2));
assert cobarD(tau1) == concat(xi1,tau0);

// A reduced coproduct can also introduce a global tau.
let both = mcword[0]{([],[1,1])};
let a = mcword[0]{([1],[1]),([],[1])};
let b = tauShift[1](concat(xi1,xi1));
let c = concat(tau0,tau1);
let d = concat(tau1,tau0);
assert cobarD(both) == add(add(a,b),add(c,d));
assert cobarD(cobarD(both)) == zero[MCobar(3,4,1)];
assert realize(cobarD(both)) == cobarD(realize(both));
show cobarD(both);
show realize(cobarD(both));

// This normalized degree-six chain bounds the motivic cubic relation.
let lift = add(add(mcword[0]{([],[0,1]),([],[0,1])},mcword[0]{([],[1]),([2],[1])}),
               add(add(mcword[0]{([1],[]),([],[1,1])},mcword[0]{([1],[0,1]),([],[1])}),
                   mcword[0]{([2],[1]),([],[1])}));
let h1cube = concat(xi1,concat(xi1,xi1));
let h0sqh2 = concat(tau0,concat(tau0,mcword[0]{([2],[])}));
assert cobarD(lift) == add(tauShift[1](h1cube),h0sqh2);
