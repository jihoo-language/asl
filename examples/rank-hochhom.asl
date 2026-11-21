rvideal Exterior[[0]] { rvmon[[0]]{1} };
let one = rvmtorgen[Exterior, 0, [0], 0];
let e = rvmtorgen[Exterior, 1, [1], 0];
let x = rvHHClass(rvHHWord(e,e,e));
let y = rvHHClass(rvHHWord(one,e,e,e));
assert x == rvHHgen[Exterior,5,[3],0];
assert y == rvHHgen[Exterior,6,[3],0];
assert rvHHRep(x) == rvHHWord(e,e,e);
assert rvConnesB(x) == y;
assert rvConnesB(y) == zero[RVHH(Exterior,7,[3])];
assert rvConnesB(rvHHClass(rvHHWord(e,e))) == zero[RVHH(Exterior,4,[2])];

rvideal Flag[[1,2,3,4,5,6]] {
    rvmon[[1,2,3,4,5,6]]{1,1,0,0,0,0},
    rvmon[[1,2,3,4,5,6]]{0,0,1,1,0,0},
    rvmon[[1,2,3,4,5,6]]{0,0,0,0,1,1},
    rvmon[[1,2,3,4,5,6]]{1,0,1,0,0,0},
    rvmon[[1,2,3,4,5,6]]{0,0,0,1,1,0}
};
let unit = rvmtorgen[Flag,0,[0,0,0,0,0,0],0];
let a = rvmtorgen[Flag,1,[1,1,0,0,0,0],0];
let b = rvmtorgen[Flag,1,[0,0,1,1,0,0],0];
let c = rvmtorgen[Flag,1,[0,0,0,0,1,1],0];
let top = rvmtorgen[Flag,4,[1,1,1,1,1,1],0];
assert rvTHigher(a,b,c) == top;
let left = rvHHWord(unit,top);
let ends = add(rvHHWord(a,b,c),rvHHWord(c,a,b));
let w = rvHHWord(unit,a,b,c);
assert rvHHd(w) == add(left,ends);
assert rvHHClass(rvHHd(w)) == zero[RVHH(Flag,5,[1,1,1,1,1,1])];
assert rvHHClass(left) == rvHHClass(ends);
let coefficient = rvHHClass(rvHHWord(top));
assert coefficient == rvHHgen[Flag,4,[1,1,1,1,1,1],0];
assert rvConnesB(coefficient) == rvHHClass(left);
assert rvHHClass(left) == rvHHgen[Flag,5,[1,1,1,1,1,1],0];
assert rvHHClass(rvHHRep(rvConnesB(coefficient))) == rvConnesB(coefficient);
rvHHdim[Flag,4,[1,1,1,1,1,1]];
rvHHdim[Flag,5,[1,1,1,1,1,1]];
rvHHdim[Flag,6,[1,1,1,1,1,1]];
show rvConnesB(coefficient);
