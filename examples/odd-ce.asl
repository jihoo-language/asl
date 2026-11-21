let a = resA[1](pword[3]{([], [0])});
let x = resA[1](pword[3]{([1], [])});
let v = resA[1](pword[3]{([], [1])});
assert ceD(ceClass[2, 0](v)) == ceMul(ceClass[2, 1](x), ceClass[2, 1](a));

let source = ceClass[2, 1](concat(v, concat(v, x)));
let third = ceNext(source);
let beta = add(pword[3]{([1], []), ([2], [])}, pword[3]{([2], []), ([1], [])});
let a3 = ceClass[3, 1](a);
let target = ceMul(ceMul(a3, a3), ceClass[3, 2](resA[1](beta)));
assert ceD(third) == pscale[2](target);
assert ceNext(ceD(third)) == zero[PCEAn(3, 1, 4, 4, 14, 4)];
show ceD(third);
pceandim[3, 1, 4, 3, 14, 1];
pceandim[3, 1, 4, 4, 14, 4];

let vp = ceClass[2, 0](concat(v, concat(v, v)));
let permanent = ceNext(ceNext(ceNext(vp)));
assert extClass(ceRep(permanent)) == pscale[2](pextangen[3, 1, 3, 15, 0]);
