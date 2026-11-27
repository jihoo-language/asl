rvhodgewick B[21]{0,1}{[2]};
rvVerbitskyCyclicDims[B,1,hc,2,0];
rvVerbitskyCyclicDims[B,1,hc,2,4];
rvVerbitskyCyclicDims[B,1,hc,2,5];
rvVerbitskySBI[B,1,2,0];
rvVerbitskySBI[B,1,2,1];

let x=rvVerbitskyCyclicTerm[B,1,hh,0,1,0]{1};
let bx=rvVerbitskyChainConnes(x);
assert bx==rvVerbitskyCyclicTerm[B,1,hh,1,1,0]{0,1};
let bxclass=rvVerbitskyCyclicB(rvVerbitskyCyclicI(rvVerbitskyCyclicClass(x)));
assert bxclass==rvVerbitskyCyclicClass(bx);
show bxclass;

let lift=rvVerbitskyNegativeLift(bxclass);
assert rvVerbitskyNegativeEval(lift)==bxclass;
assert rvVerbitskyNegativeU(lift)==zero[RVVerbitskyCyclic(B,1,hn,-1,1)];
assert rvVerbitskyNegativePeriodic(lift)==zero[RVVerbitskyCyclic(B,1,hp,1,1)];
show lift;
let raw=rvVerbitskyNegativePeriodic(rvVerbitskyCyclicRepresentative(lift));
let primitive=rvVerbitskyCyclicPrimitive(raw);
assert primitive==rvVerbitskyCyclicTerm[B,1,hp,2,1,-1]{1};
assert rvVerbitskyCyclicDiff(primitive)==raw;
show primitive;
