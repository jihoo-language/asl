let t = rvHeckeT[24,2];
rvheckepacket E(t);
rvHeckeEigenData[E];

let a2 = rvHeckeEigenCoeff[E,2];
let a3 = rvHeckeEigenCoeff[E,3];
let a4 = rvHeckeEigenCoeff[E,4];
let a6 = rvHeckeEigenCoeff[E,6];
assert a2 == rvHeckeEigenGenerator[E];
assert rvHeckeEigenMul(a2,a3) == a6;
assert a4 == add(rvHeckeEigenMul(a2,a2),rvHeckeEigenScalar[E,-8388608,1]);
show a3;
rvHeckeEigenTrace(a3);
rvHeckeEigenNorm(a3);

let shifted = add(rvHeckeScale[5](t),rvHeckeScale[7](rvHeckeUnit[24]));
rvheckepacket F(shifted);
let b3 = rvHeckeEigenCoeff[F,3];
show b3;
rvHeckeEigenCoord[1](b3);
assert rvHeckeEigenCoeff[F,1] == rvHeckeEigenScalar[F,1,1];
