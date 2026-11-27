// The first Bockstein of tau^32 modulo rho^17 is nonzero.
let a = extClass(rhoReduce[17](rmcword[0,32]{}));
let b = rhoBeta(a);
show b;

// The second Bockstein is zero in Ext. Its canonical lift gives a nonzero
// cochain differential, and rhoBetaNull solves for an explicit primitive.
assert rhoBeta(b) == zero[RMRhoExt(17,2,34,2)];
let raw = rhoBetaRep(b);
let witness = rhoBetaNull(b);
assert cobarD(witness) == raw;
assert extPrimitive(raw) == witness;
show raw;
show witness;
