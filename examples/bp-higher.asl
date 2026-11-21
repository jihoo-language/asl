// An explicit integral fivefold defining system, constant along each interval length.
let t = bpt[5,2,1,1,1];
let v = bpv[5,2,1,1];
let x = bpExtClass[8](t);
let w2 = add(bpNeg(bpMul(v,t)),bpMul(bpint[5,2,1,-3],bpPow[2](t)));
let w3 = add(bpMul(bpint[5,2,1,-5],bpMul(bpPow[2](v),t)),
             add(bpMul(bpint[5,2,1,-24],bpMul(v,bpPow[2](t))),
                 bpMul(bpint[5,2,1,-39],bpPow[3](t))));
let w4 = add(add(bpMul(bpint[5,2,1,-77],bpMul(bpPow[3](v),t)),
                 bpMul(bpint[5,2,1,-573],bpMul(bpPow[2](v),bpPow[2](t)))),
             add(bpMul(bpint[5,2,1,-1903],bpMul(v,bpPow[3](t))),
                 bpMul(bpint[5,2,1,-2369],bpPow[4](t))));
let m5 = bpHigherSystem[5](x,x,x,x,x,w2,w2,w2,w2,w3,w3,w3,w4,w4);
let beta = bpExtClass[40](bpDelta(bpDelta(bpReduce[2](bpv[5,2,0,2]))));
assert m5 == bpNeg(beta);
assert bpChosenHigher(x,x,x,x,x) == m5;
assert bpScale[5](m5) == zero[BPExt(5,2,2,40)];
show m5;
