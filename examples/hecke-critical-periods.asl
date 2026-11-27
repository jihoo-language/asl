rvheckepacket Delta(rvHeckeT[12,2]);
rvheckeeigenperiod Plus[Delta,plus,2];
rvheckeeigenperiod Minus[Delta,minus,1];
show rvHeckePeriodCoeff[Plus,0];
show rvHeckePeriodCoeff[Plus,10];
show rvHeckeCriticalRatio[Plus,1,3];
show rvHeckeCriticalRatio[Minus,2,4];

rvheckepacket E(rvHeckeT[24,2]);
rvheckeeigenperiod P[E,plus,2];
rvheckeeigenperiod Q[E,plus,4];
assert rvHeckeCriticalRatio[P,1,3] == rvHeckeCriticalRatio[Q,1,3];
show rvHeckeCriticalRatio[P,1,3];
rvheckeeigenperiod Odd[E,minus,1];
show rvHeckeCriticalRatio[Odd,2,4];
