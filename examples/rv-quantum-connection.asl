// The operators use left coefficients q^a before powers of theta.
let theta = rvOre[0,1];
let q = rvOre[1,0];
assert rvOreMul(theta,q) == add(rvOre[1,1],q);
let u = rvQGrOne[2,4];
let x = rvQGrShift[1](u);
assert rvOreApply(rvOreMul(theta,q),u) == rvOreApply(theta,x);
assert rvQGrNabla(u) == rvQGrBasis[2,4,[1]];
let l = rvQGrAnnihilator[2,4];
assert l == add(rvOre[0,5],add(rvOre[1,1,-4],rvOre[1,0,-2]));
assert rvOreApply(l,u) == zero[RVQGr(2,4)];
show l;
rvQGrPeriod[2,4,4];
// The vector recurrence fixes the resonant coefficients of Gr(2,5).
let resonance = rvQGrAnnihilator[2,5];
assert rvOreApply(resonance,rvQGrOne[2,5]) == zero[RVQGr(2,5)];
rvOreCoeff[0,7](resonance);
rvOreCoeff[0,10](resonance);
rvQGrPeriod[2,5,3];
// The unit generates a rank-14 cyclic module inside the rank-20 connection.
let cyclic = rvQGrAnnihilator[3,6];
assert rvOreApply(cyclic,rvQGrOne[3,6]) == zero[RVQGr(3,6)];
rvOreCoeff[0,14](cyclic);
rvQGrPeriod[1,5,4];
