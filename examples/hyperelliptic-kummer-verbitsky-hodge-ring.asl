rvhodgewick E[11]{-1,-1,1}{[2]};
rvVerbitskyRank[E,2,0];
rvVerbitskyRank[E,2,2];
rvVerbitskyRank[E,2,4];
rvVerbitskyRank[E,2,6];
rvVerbitskyRank[E,2,8];
rvVerbitskyRelations[E,2,6];

let d=rvVerbitskyClass(rvHodgeD[E,2,0]);
let q=rvVerbitskyClass(rvHodgeQ[E,2]{1});
let theta=rvVerbitskyClass(rvHodgeQ[E,2]{0,1});
let center=add(rvVerbitskyScale[2,1](theta),neg(q));
let d3=rvVerbitskyCup(d,rvVerbitskyCup(d,d));
assert rvVerbitskyScale[3,1](rvVerbitskyCup(d,q))==rvVerbitskyScale[11,1](d3);
assert rvVerbitskyScale[6,1](rvVerbitskyCup(d,theta))==rvVerbitskyScale[11,1](d3);
assert rvVerbitskyCup(d,center)==zero[RVVerbitsky(E,2,6)];
rvVerbitskyPair(center,center);
rvVerbitskySignature[4](d);
