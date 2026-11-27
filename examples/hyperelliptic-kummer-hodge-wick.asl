rvhodgewick E[11]{-1,-1,1}{[2]};
rvHodgeTrace[E]{0,1};
rvHodgeTrace[E]{(1,2),(1,3)};

let d=rvHodgeD[E,2,0];
let d2=rvHodgeCup(d,d);
let q=rvHodgeQ[E,2]{1};
let theta=rvHodgeQ[E,2]{0,1};
rvHodgePair(d2,d2);
rvHodgePair(d2,q);
rvHodgePair(d2,theta);
rvHodgePair(q,q);
rvHodgePair(q,theta);
rvHodgePair(theta,theta);

let central=rvHodgeQ[E,3]{-1,2};
rvHodgePair(rvHodgeQ[E,3]{1},rvHodgeCup(central,central));
let bb=rvHodgeInverseBB[E,3];
rvHodgeIntegral(rvHodgeCup(bb,rvHodgeCup(bb,bb)));
