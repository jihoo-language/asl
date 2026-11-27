rvhodgewick C[20]{0,1}{[0,1],[1,0]};

let chain=rvVerbitskyCyclicTerm[C,1,hh,3,4,0]{0,1,2,3};
assert rvVerbitskyAdams[3](rvVerbitskyAdams[2](chain))==rvVerbitskyAdams[6](chain);
let p1=rvVerbitskyHodgeProject[1](chain);
let p2=rvVerbitskyHodgeProject[2](chain);
let p3=rvVerbitskyHodgeProject[3](chain);
assert add(add(p1,p2),p3)==chain;

rvVerbitskyHodgeDim[C,1,hh,2,2,2];
rvVerbitskyHodgeDim[C,1,hh,2,3,1];
rvVerbitskyHodgeDim[C,1,hh,4,6,2];
rvMonomialCICotangentDim[4,2,6]{[2,1],[2,1]};
rvMonomialCICotangentDiff{[2,1],[2,1]}{0,0}{2,1}{};

let p=rvVerbitskyHHClass(rvVerbitskyHHTerm[C,1,2,0,3]{2,1});
show rvVerbitskyHodgeProject[2](p);

rvhodgewick B[21]{0,1}{[2]};
let tower=rvVerbitskyCyclicBasis[B,1,hc,2,0,0];
assert rvVerbitskyHodgeProject[1](tower)==tower;
let down=rvVerbitskyCyclicS(tower);
assert rvVerbitskyHodgeProject[0](down)==down;
let x=rvVerbitskyCyclicBasis[B,1,hc,0,1,0];
assert rvVerbitskyHodgeProject[0](x)==x;
let bx=rvVerbitskyCyclicB(x);
assert rvVerbitskyHodgeProject[1](bx)==bx;
show bx;
