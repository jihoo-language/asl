// SU(3) at level two. The common cyclotomic order is D=60.
let u = rvWzwBasis[3,2,[0,0,0]];
let f = rvWzwBasis[3,2,[0,0,1]];
let fbar = rvWzwBasis[3,2,[0,1,1]];
assert rvWzwCharge(f) == fbar;
assert rvWzwCharge(fbar) == f;
let v = add(u,rvWzwScale(rvCyclo[60,4],f));
assert rvWzwT[-3](rvWzwT[3](v)) == v;
assert rvWzwS(rvWzwS(v)) == rvWzwScale(rvCycloScale[-75](rvCyclo[60,0]),rvWzwCharge(v));
assert rvWzwS(rvWzwSAdjoint(v)) == rvWzwScale(rvCycloScale[75](rvCyclo[60,0]),v);
// Read the framing scalar from the vacuum and check its square.
let stup = rvWzwS(rvWzwT[1](u));
let sttwo = rvWzwS(rvWzwT[1](stup));
let stthree = rvWzwS(rvWzwT[1](sttwo));
let b = rvWzwCoord[[0,0,0]](stthree);
assert b == rvCycloScale[375](add(rvCyclo[60,4],rvCyclo[60,14]));
assert rvCycloMul(b,b) == rvCycloScale[-421875](rvCyclo[60,48]);
assert stthree == rvWzwScale(b,u);
show b;
// The raw S entry is the vacuum amplitude times the conjugate character.
let vacuum = rvWzwSEntry[3,2,[0,0,0],[0,0,1]];
let character = rvVerlinde[[0,0,1]](rvFusionBasis[2,[0,0,1]]);
let actual = rvWzwSEntry[3,2,[0,0,1],[0,0,1]];
assert actual == rvCycloMul(vacuum,rvCycloConj(rvCycloEmbed[60](character)));
assert actual == rvWzwCoord[[0,0,1]](rvWzwS(f));
show actual;
show rvWzwTEntry[3,2,[0,0,1]];
