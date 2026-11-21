mixedcomplex Stair[3] {
    x0:0, x1:2, x2:4, y1:1, y2:3, y3:5
} {
    x1 => [(y1,-1)], x2 => [(y2,-1)]
} {
    x0 => [(y1,1)], x1 => [(y2,1)], x2 => [(y3,1)]
};
let a = mixU[0](mixchain[Stair,0][(x0,1)]);
let b = mixU[1](mixchain[Stair,2][(x1,1)]);
let c = mixU[2](mixchain[Stair,4][(x2,1)]);
let z = add(add(a,b),c);
let y = mixU[3](mixchain[Stair,5][(y3,1)]);

let e0 = mixgen[Stair,0,0,0,0];
let e1 = mixNext(e0);
let e2 = mixNext(e1);
let e3 = mixNext(e2);
assert e3 == mixgen[Stair,3,0,0,0];
assert mixPrev(e3) == e2;
assert mixPageRep(e3) == z;
assert mixD(z) == y;
assert mixD(e0) == zero[MixedE(Stair,0,-1,0)];
assert mixD(e1) == zero[MixedE(Stair,1,-1,1)];
assert mixD(e2) == zero[MixedE(Stair,2,-1,2)];
assert mixD(e3) == mixPageClass[3,3](y);
assert mixD(pscale[2](e3)) == pscale[2](mixgen[Stair,3,-1,3,0]);
assert mixNext(mixD(e3)) == zero[MixedE(Stair,4,-1,3)];
mixdim[Stair,4,0,0];
show mixPageRep(e3);
show mixD(e3);
