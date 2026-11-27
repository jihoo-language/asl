// Both 14-dimensional lattices use the basis e_1,...,e_7,f_1,...,f_7.
rvcrystal Base[3,7] {
    (1,0,1), (2,1,1), (3,2,1), (4,3,1), (5,4,1), (6,5,1), (7,6,1),
    (8,7,3), (9,8,3), (10,9,3), (11,10,3), (12,11,3), (13,12,3), (0,13,-3)
};
rvcrystal Tilted[3,7] {
    (1,0,1), (2,1,1), (3,2,1), (4,3,1), (5,4,1), (6,5,1), (7,6,1),
    (8,7,3), (9,8,3), (10,9,3), (11,10,3), (12,11,3), (13,12,3), (0,13,-3),
    (0,6,27)
};

assert rvCrystalCongruent[Base,Tilted,3] == S.hi;
assert rvCrystalCongruent[Base,Tilted,4] == S.lo;
let first = rvCrystalEOType[Base];
assert first == eoType{0,1,2,3,4,5,6};
assert rvCrystalEOType[Tilted] == first;
show first;
rvCrystalChar[Base];
rvCrystalChar[Tilted];
rvCrystalNewton[Base];
rvCrystalNewton[Tilted];
