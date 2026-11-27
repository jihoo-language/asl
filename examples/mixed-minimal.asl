mixedcomplex Stair[3] {
    x0:0, x1:2, x2:4, y1:1, y2:3, y3:5
} {
    x1 => [(y1,-1)], x2 => [(y2,-1)]
} {
    x0 => [(y1,1)], x1 => [(y2,1)], x2 => [(y3,1)]
};
mixminimal Min[Stair];
let x = mixMinGen[Min,0,0];
let y = mixMinGen[Min,5,0];
let X = mixU[0](x);
let Y = mixU[3](y);
let lift = add(add(mixU[0](mixchain[Stair,0][(x0,1)]),mixU[1](mixchain[Stair,2][(x1,1)])),
               mixU[2](mixchain[Stair,4][(x2,1)]));
assert mixDelta[0](x) == zero[MixedMin(Min,-1)];
assert mixDelta[1](x) == zero[MixedMin(Min,1)];
assert mixDelta[2](x) == zero[MixedMin(Min,3)];
assert mixDelta[3](x) == y;
assert mixDelta[4](x) == zero[MixedMin(Min,7)];
assert mixD(X) == Y;
assert mixD(Y) == zero[MixedMinU(Min,-2)];
assert mixMinI(X) == lift;
assert mixMinI(X) == mixPageRep(mixgen[Stair,3,0,0,0]);
assert mixMinP[Min](lift) == X;
assert mixMinH[Min](lift) == zero[MixedU(Stair,1)];
assert mixD(lift) == mixMinI(Y);
assert mixMinI(mixU[5](X)) == mixU[5](lift);

let z = mixU[0](mixchain[Stair,1][(y1,1)]);
assert mixMinP[Min](z) == mixU[2](y);
let h = mixMinH[Min](z);
assert h == add(mixU[0](mixchain[Stair,2][(x1,-1)]),mixU[1](mixchain[Stair,4][(x2,-1)]));
assert mixD(h) == add(z,pscale[2](mixMinI(mixMinP[Min](z))));
assert mixMinP[Min](h) == zero[MixedMinU(Min,2)];
assert mixMinH[Min](h) == zero[MixedU(Stair,3)];
mixmindim[Min,0];
mixmindim[Min,2];
mixmindim[Min,5];
show mixD(X);
show mixMinI(X);
