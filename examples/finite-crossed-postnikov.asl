rvgroup H[4];
rvgroup E[4];
rvcrossed C[H,E]{0,2,0,2}{{0,1,2,3},{0,3,2,1},{0,1,2,3},{0,3,2,1}};
rvCrossedData[C];
assert rvCrossedBoundary[C](rvGroupElem[H,1])==rvGroupElem[E,2];
assert rvCrossedAssoc[C,1,1,1]==rvGroupElem[H,2];

rvpostnikov P[C,2];
let k=rvPostnikovK[P];
rvPostnikovData[P];
rvPostnikovClass(k);
rvPostnikovCoeff[1,1,1](k);
let id=rvPostnikovCell[P,1,0];
let assoc=rvPostnikovAssoc[1,1,1](k);
assert assoc==rvPostnikovCell[P,1,1];
// The two paths around the pentagon at (1,1,1,1).
let left=rvPostnikovVertical(rvPostnikovAssoc[0,1,1](k),rvPostnikovAssoc[1,1,0](k));
let right=rvPostnikovVertical(rvPostnikovHorizontal(id,assoc),rvPostnikovHorizontal(assoc,id));
assert left==right;

// A genuinely twisted kernel action with a gauge-dependent representative.
rvgroup H3[3];
rvgroup E2{{0,1},{1,0}};
rvcrossed Twist[H3,E2]{0,0,0}{{0,1,2},{0,2,1}};
rvpostnikov Q[Twist,1];
let canonical=rvPostnikovK[Q];
let change=rvPostnikovCochain[Q,2]{1};
let changed=rvPostnikovGauge(canonical,change);
rvPostnikovData[Q];
rvPostnikovCoeff[1,1,1](changed);
rvPostnikovClass(changed);
assert rvPostnikovGauge(canonical,rvPostnikovGaugeWitness(canonical,changed))==changed;
assert rvPostnikovHorizontal(rvPostnikovCell[Q,1,1],rvPostnikovCell[Q,0,1])==rvPostnikovCell[Q,1,0];
