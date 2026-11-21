let alpha = bpextgen[3, 2, 1, 4, 0];
let bracket = bpTriple(alpha, alpha, alpha);
let beta = bpExtClass[12](bpDelta(bpDelta(bpReduce[2](bpv[3,2,0,2]))));
bpTripleDefined(alpha, alpha, alpha);
bpTripleContains(bracket, bpNeg(beta));
bpTripleContains(bracket,beta);
assert bpTripleBase(bracket) == bpNeg(beta);
show bracket;
bpTripleSize(bracket);

let eta = bpextgen[2,2,1,2,0];
let empty_middle = zero[BPExt(2,2,1,4)];
let double = bpTriple(eta,empty_middle,eta);
let target0 = bpextgen[2,2,2,8,0];
let target1 = bpextgen[2,2,2,8,1];
bpTripleContains(double,target0);
bpTripleContains(double,target1);
bpTripleContains(double,add(target0,target1));
bpTripleSize(double);
assert double == bpTriple(zero[BPExt(2,2,1,2)],empty_middle,eta);
bpTripleDefined(eta, eta, eta);

// All inputs are nonzero. The coset is the other translate of the same subgroup.
let alpha_2 = bpextgen[2,2,1,4,0];
let affine = bpTriple(eta,alpha_2,eta);
assert bpTripleBase(affine) == target1;
show affine;
bpTripleContains(affine,target0);
bpTripleContains(affine,target1);
bpTripleContains(affine,add(target0,target1));
