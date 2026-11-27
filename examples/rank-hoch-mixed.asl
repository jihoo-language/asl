rvideal Exterior[[0]] { rvmon[[0]]{1} };
rvHHMixed Odd[Exterior,[3]];
rvHHMixed Even[Exterior,[2]];
let unit = rvmtorgen[Exterior,0,[0],0];
let e = rvmtorgen[Exterior,1,[1],0];
let oddLow = rvHHWord(e,e,e);
let oddHigh = rvHHWord(unit,e,e,e);
let raw = rvHHToMixed[Odd](oddLow);
assert rvMixedToHH(raw) == oddLow;
assert mixB(raw) == rvHHToMixed[Odd](oddHigh);
assert mixb(raw) == zero[MixedChain(Odd,4)];

let x = rvHHToMixed[Odd,4](rvHHClass(oddLow));
assert x == mixgen[Odd,1,-3,4,0];
assert rvMixedToHH(x) == rvHHClass(oddLow);
assert mixD(x) == rvHHToMixed[Odd,5](rvConnesB(rvHHClass(oddLow)));
assert rvMixedToHH(mixD(x)) == rvHHClass(oddHigh);
assert mixNext(mixD(x)) == zero[MixedE(Odd,2,-4,5)];

let evenLow = rvHHClass(rvHHWord(e,e));
let y = rvHHToMixed[Even](evenLow);
assert rvMixedToHH(y) == evenLow;
assert mixD(y) == zero[MixedE(Even,1,2,1)];
assert mixNext(y) == mixgen[Even,2,3,0,0];
show x;
show mixD(x);
show mixNext(y);
