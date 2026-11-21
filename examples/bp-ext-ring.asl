// The first cobar generator has a nonzero square of order two.
let eta = bpExtClass[2](bpt[2,2,1,1,1]);
let square = bpCup(eta,eta);
assert square == bpextgen[2,2,2,4,0];
assert add(square,square) == zero[BPExt(2,2,2,4)];
assert bpExtClass[4](bpExtRep(square)) == square;
let one = bpextgen[2,2,0,0,0];
assert bpCup(one,eta) == eta;
assert bpCup(eta,one) == eta;
show square;

// Two invariant factors give four distinct classes in this bidegree.
bpext[2,2,2,8];
let x = bpextgen[2,2,2,8,0];
let y = bpextgen[2,2,2,8,1];
let xy = add(x,y);
assert bpExtClass[8](bpExtRep(xy)) == xy;
show xy;
card(\z : BPExt(2,2,2,8) => S.hi);

// Passing to Ext compares the two signed representatives.
let t = bpt[3,2,1,1,1];
let v = bpv[3,2,1,1];
let w = add(bpMul(v,t),bpPow[2](t));
let m = bpExtClass[12](bpMassey(t,t,t,w,w));
let beta = bpExtClass[12](bpDelta(bpDelta(bpReduce[2](bpv[3,2,0,2]))));
assert m == bpNeg(beta);
assert bpScale[3](beta) == zero[BPExt(3,2,2,12)];
bporder[12](bpExtRep(beta));
show add(m,beta);
