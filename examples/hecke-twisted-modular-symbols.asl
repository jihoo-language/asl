rvheckepacket Delta(rvHeckeT[12,2]);
rvheckeeigenperiod P[Delta,plus,2];
rvheckeeigenperiod Q[Delta,minus,1];

let a = rvHeckeSymbol[P,-3,7,2,5];
let b = rvHeckeSymbol[P,2,5,1,0];
assert add(a,b) == rvHeckeSymbol[P,-3,7,1,0];
show rvHeckeSymbolCoeff[10](rvHeckeSymbol[P,0,1,1,0]);

show rvHeckeTwistRatio[P,5,3,3];
show rvHeckeTwistRatio[Q,-3,3,2];
show rvHeckeTwistRatio[P,-4,6,3];

rvheckepacket E(rvHeckeT[24,2]);
rvheckeeigenperiod Raised[E,plus,2];
show rvHeckeTwistRatio[Raised,5,3,3];
rvheckeeigenperiod RaisedOdd[E,minus,1];
show rvHeckeTwistRatio[RaisedOdd,-3,3,2];
