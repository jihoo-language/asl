// The two units of the real-motivic Hopf algebroid differ on tau.
let tau = rmcword[0,1]{};
let rho = rmcword[1,0]{};
let t0 = rmcword[0,0]{([],[1])};
assert cobarD(tau) == rhoShift[1](t0);
assert cobarD(rho) == zero[RMCobar(1,-1,-1)];
show cobarD(tau);

// The coefficient on the right travels through the preceding arrow.
let square = rmcword[0,0]{([],[2])};
assert concat(t0,tau) == add(concat(tau,t0),rhoShift[1](square));
assert realComplex(square) == mcword[1]{([1],[])};

// Frobenius gives rho^2 [tau_0^2], expanded as three normal monomials.
let d = cobarD(rmcword[0,2]{});
assert d == rhoShift[2](square);
assert cobarD(d) == zero[RMCobar(2,0,-2)];
assert realComplex(d) == zero[MCobar(1,0,-2)];
assert realComplex(rho) == zero[MCobar(0,-1,-1)];
show d;
