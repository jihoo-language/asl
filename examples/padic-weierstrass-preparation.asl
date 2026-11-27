let unramified = rvTatePoly[3,4](1,3,1,12);
rvTatePrepare(unramified);
rvTateFactors(unramified);
rvTateUnramified[0](unramified);

let ramified = rvTatePoly[3,4](-3,-9,1,3);
rvTatePrepare(ramified);
rvTateNewton(ramified);
rvTateUnramified[0](ramified);

rvheckepacket Delta(rvHeckeT[12,2]);
rvheckeeigenperiod P[Delta,plus,2];
rvmazurordinary O[P,11,4];
let disc = rvMazurDisc[0,1](rvMazurAmice[O,3,0,0]);
rvTatePrepare(disc);
rvTateUnramified[0](disc);
rvTatePrepare(rvMazurDiscRestrict[5,1](disc));
