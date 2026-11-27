// The carry cocycle on C3 changes the order of a unit flux from three to nine.
rvgroup G[3];
rvdwcoeff P[G,3];
let k=rvPostnikovCochain[P,3]{0,1,1,1,0,2,2,2};
rvdouble D[P](k);
rvDoubleData[D];

let x=rvDoubleSimple[D,1]{0,1,2};
let y=rvDoubleFuse(x,x);
let cube=rvDoubleFuse(x,y);
let charge=rvDoubleSimple[D,0]{0,6,3};
assert x==rvDoubleBasis[D,3];
assert cube==charge;
assert rvDoubleFuse(y,cube)==rvDoubleFuse(x,rvDoubleFuse(y,y));
assert rvDoubleFuse(cube,rvDoubleFuse(cube,cube))==rvDoubleBasis[D,0];
assert rvDoubleFuse(x,rvDoubleDual(x))==rvDoubleBasis[D,0];
rvDoubleSector[D,3];
show cube;
show rvDoubleLens[D,3];
let expected=rvCycloRat[3](add(rvCyclo[9,0],rvCycloScale[2](rvCyclo[9,3])));
assert rvDoubleLens[D,3]==expected;
assert rvDoubleLens[D,-3]==rvCycloQGalois[-1](expected);
