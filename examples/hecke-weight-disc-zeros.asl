rvheckepacket Delta(rvHeckeT[12,2]);
rvheckeeigenperiod P[Delta,plus,2];
rvmazurordinary O[P,11,4];

let low = rvMazurDisc[0,1](rvMazurAmice[O,2,0,0]);
let amice = rvMazurAmice[O,3,0,0];
let high = rvMazurDisc[0,1](amice);
rvMazurDiscZeros(low);
show rvMazurDiscRoot[5](low);
rvMazurDiscZeros(high);
show rvMazurDiscRoot[5](high);
let near = rvMazurDiscRestrict[5,1](high);
assert near == rvMazurDisc[55,2](amice);
rvMazurDiscZeros(near);
show rvMazurDiscRoot[10](near);
rvMazurDiscZeros(rvMazurDiscRestrict[4,1](high));
rvMazurDiscZeros(rvMazurDisc[0,1](rvMazurAmiceProject(rvMazurAmiceProject(amice))));
