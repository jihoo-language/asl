rvideal Exterior[[0]] { rvmon[[0]]{1} };
let one = rvmtorgen[Exterior, 0, [0], 0];
let e = rvmtorgen[Exterior, 1, [1], 0];

let odd = rvHHWord(e,e,e);
let even = rvHHWord(e,e);
assert rvConnesB(odd) == rvHHWord(one,e,e,e);
assert rvConnesB(even) == zero[RVHHChain(Exterior,4,[2])];
assert rvHHd(odd) == zero[RVHHChain(Exterior,4,[3])];
assert rvHHd(rvHHWord(one)) == zero[RVHHChain(Exterior,-1,[0])];
assert rvConnesB(rvHHWord(one)) == zero[RVHHChain(Exterior,1,[0])];
assert rvConnesB(rvConnesB(odd)) == zero[RVHHChain(Exterior,7,[3])];
assert rvHHLength[2](odd) == odd;
assert rvHHLength[3](odd) == zero[RVHHChain(Exterior,5,[3])];
card \x : RVHHChain(Exterior, 6, [3]) => S.hi;
show rvConnesB(odd);
