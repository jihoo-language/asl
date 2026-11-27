rvheckepacket Delta(rvHeckeT[12,2]);
rvheckeeigenperiod P[Delta,plus,2];

let theta0 = rvMazurTateTheta[P,3,0,2];
let theta1 = rvMazurTateTheta[P,3,1,2];
let theta2 = rvMazurTateTheta[P,3,2,2];
let rhs = add(rvMazurTateScale(rvHeckeEigenCoeff[Delta,3],theta1),
              neg(rvMazurTateScale(rvHeckeEigenScalar[Delta,59049,1],rvMazurTateTransfer(theta0))));
assert rvMazurTateTrace(theta2) == rhs;
show rvMazurTateCoeff[1](theta1);
show rvMazurTateCoeff[1](rvMazurTateTrace(theta2));

rvmazurordinary O[P,11,12];
rvMazurTateOrdinaryData[O];
show rvMazurTateMeasure[O,1,2,2];
