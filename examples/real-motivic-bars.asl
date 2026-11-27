// The first module has two finite intervals and one free interval.
rmbars[2,4];
rmbarRank[2,4,4,1];
rmbarRank[2,4,4,2];
let free = rmbarGen[2,4,2];
assert rhoProject(free) == extgen[2,4,0];
rhoorder(free);

// This bar is born in topological degree eight and has length six.
rmbars[2,5];
let z6 = rmbarRep[2,5,1];
let g6 = rmbarGen[2,5,1];
let u6 = rmbarNull[2,5,1];
assert extClass(z6) == g6;
rhoorder(g6);
assert cobarD(u6) == rhoShift[6](z6);
let d6 = rhoD(rhoClass[6,0](u6));
assert d6 == rhoClass[6,6](rhoShift[6](z6));
assert d6 == rhoD(rhobockgen[6,1,2,-3,0,0]);
assert rhoNext(d6) == zero[RhoBock(7,2,2,-3,6)];
show d6;

// The length-seven bar gives a separate non-dyadic differential.
rmbars[3,6];
let z7 = rmbarRep[3,6,2];
let g7 = rmbarGen[3,6,2];
let u7 = rmbarNull[3,6,2];
assert extClass(z7) == g7;
rhoorder(g7);
assert cobarD(u7) == rhoShift[7](z7);
let d7 = rhoD(rhoClass[7,0](u7));
assert d7 == rhoClass[7,7](rhoShift[7](z7));
assert d7 == rhoD(rhobockgen[7,2,4,-2,0,0]);
assert rhoNext(d7) == zero[RhoBock(8,3,4,-2,7)];
show d7;
