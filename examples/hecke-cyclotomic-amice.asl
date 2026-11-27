rvheckepacket Delta(rvHeckeT[12,2]);
rvheckeeigenperiod P[Delta,plus,2];
rvmazurordinary O[P,11,4];

let levelOne = rvMazurAmice[O,1,0,0];
let levelTwo = rvMazurAmice[O,2,0,0];
let levelThree = rvMazurAmice[O,3,0,0];
assert rvMazurAmiceProject(levelThree) == levelTwo;
assert rvMazurAmiceProject(levelTwo) == levelOne;

rvMazurAmiceCoeff[0](levelTwo);
rvMazurAmiceCoeff[1](levelTwo);
rvMazurAmiceCoeff[1](levelThree);
rvMazurAmiceIwasawa[4](levelTwo);
rvMazurAmiceIwasawa[4](rvMazurAmice[O,2,0,2]);
rvMazurAmiceCharacter[1](levelThree);
