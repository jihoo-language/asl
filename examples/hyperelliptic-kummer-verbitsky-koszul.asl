rvhodgewick B[21]{0,1}{[2,0],[0,-2]};
rvVerbitskyTorTable[B,2];
rvVerbitskyTorBetti[B,2,1,4];
rvVerbitskyTorBetti[B,2,2,4];

let x=rvVerbitskyBasis[B,2,2,0];
let y=rvVerbitskyBasis[B,2,2,1];
let z=rvVerbitskyClass(rvHodgeQ[B,2]{1});
let x2=rvVerbitskyCup(x,x);
let y2=rvVerbitskyCup(y,y);
let r=rvVerbitskyTorClass(rvVerbitskyKoszulTerm{0}(add(x2,rvVerbitskyScale[3,1](y2))));
rvVerbitskyMinimalRelation(r);

let x3=rvVerbitskyCup(x,x2);
let a=rvVerbitskyTorClass(add(rvVerbitskyKoszulTerm{2}(rvVerbitskyScale[4,1](z)),rvVerbitskyKoszulTerm{0}(rvVerbitskyScale[-161,1](x3))));
let b=rvVerbitskyTorClass(rvVerbitskyKoszulTerm{0,1}(add(rvVerbitskyScale[8,1](z),add(rvVerbitskyScale[-21,1](x2),rvVerbitskyScale[21,1](y2)))));
assert rvVerbitskyTorCup(a,b)==rvVerbitskyTorScale[-2100,1](rvVerbitskyTorBasis[B,2,3,8,0]);
