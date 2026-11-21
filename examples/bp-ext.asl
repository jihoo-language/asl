// The first-height class has order 8.
let v = bpv[2,1,0,1];
let t = bpt[2,1,1,1,1];
let c = bpv[2,1,1,1];
let alpha = add(bpMul(bpPow[3](c),t),
    add(bpMul(bpint[2,1,1,3],bpMul(bpPow[2](c),bpPow[2](t))),
    add(bpMul(bpint[2,1,1,4],bpMul(c,bpPow[3](t))),
        bpMul(bpint[2,1,1,2],bpPow[4](t)))));
bpext[2,1,1,8];
bporder[8](alpha);
assert bpNull[8](alpha) == bpPow[4](v);

// The corrected class in height two has order 16.
let v1 = bpv[2,2,1,1];
let v2 = bpv[2,2,1,2];
let t1 = bpt[2,2,1,1,1];
let t2 = bpt[2,2,1,1,2];
let base_terms = add(bpMul(v1,t2),add(bpMul(v2,t1),bpNeg(bpMul(bpPow[3](v1),t1))));
let lower_terms = add(bpMul(bpint[2,2,1,-5],bpMul(v1,bpPow[3](t1))),
                     bpMul(bpint[2,2,1,-4],bpMul(bpPow[2](v1),bpPow[2](t1))));
let arrow_terms = add(bpMul(bpint[2,2,1,2],bpMul(t1,t2)),
                     bpMul(bpint[2,2,1,-3],bpPow[4](t1)));
let sigma = add(base_terms,add(lower_terms,arrow_terms));
let x = bpv[2,2,0,1];
let z = bpv[2,2,0,2];
let correction = add(bpPow[4](x),bpMul(bpint[2,2,0,8],bpMul(x,z)));
bpext[2,2,1,8];
bporder[8](sigma);
assert bpNull[8](sigma) == correction;
assert bpD(correction) == bpMul(bpint[2,2,1,16],sigma);

// A two-cocycle of order 3 at the odd prime.
let y = bpv[3,2,2,1];
let a = bpt[3,2,2,1,1];
let b = bpt[3,2,2,2,1];
let beta = add(bpMul(bpint[3,2,2,4],bpMul(y,bpMul(a,b))),
    add(bpMul(bpint[3,2,2,5],bpMul(bpPow[2](a),b)),bpNeg(bpMul(a,bpPow[2](b)))));
bpext[3,2,1,12];
bpext[3,2,2,12];
bporder[12](beta);
assert bpD(bpNull[12](beta)) == bpMul(bpint[3,2,2,3],beta);
