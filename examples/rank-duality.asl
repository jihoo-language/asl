let x2 = rvmon[[0,1,2]]{2,0,0};
let y2 = rvmon[[0,1,2]]{0,2,0};
let z2 = rvmon[[0,1,2]]{0,0,2};
let cusp = add(rvmon[[0,1,2]]{0,3,0},rvmon[[0,1,2]]{2,0,1});
rvgauss CuspGraph[[3,4,5]](cusp);
rvdual CuspDual[[3,4,5]](cusp);
let cuspLine = rvGenerator[CuspDual,0];
assert cuspLine == rvmon[[3,4,5]]{1,0,0};
rvdual CuspBidual[[6,7,8]](cuspLine);
assert rvInIdeal[CuspBidual](rvmon[[6,7,8]]{1,0,0}) == S.lo;

// On [x:y:z]=[s^3:s^2*t:t^3], the Gauss map extends as [0:t^2:s^2].
// In particular z*W=y*V in the saturated graph.
let zW = rvmon[[0,1,2,3,4,5]]{0,0,1,0,0,1};
let yV = rvmon[[0,1,2,3,4,5]]{0,1,0,0,1,0};
assert rvInIdeal[CuspGraph](add(zW,yV)) == S.hi;
rvgroebner[CuspGraph];
rvgroebner[CuspDual];
rvgroebner[CuspBidual];

let conic = add(rvmon[[0,1,2]]{1,0,1},y2);
rvdual ConicDual[[3,4,5]](conic);
let conicLine = rvGenerator[ConicDual,0];
assert conicLine == rvmon[[3,4,5]]{0,1,0};
rvdual ConicBidual[[6,7,8]](conicLine);
assert rvInIdeal[ConicBidual](rvmon[[6,7,8]]{0,1,0}) == S.lo;
rvgroebner[ConicDual];
rvgroebner[ConicBidual];

let fermat = add(add(rvmon[[0,1,2]]{3,0,0},rvmon[[0,1,2]]{0,3,0}),rvmon[[0,1,2]]{0,0,3});
assert rvHasse[1,0,0](fermat) == x2;
assert rvHasse[0,1,0](fermat) == y2;
assert rvHasse[0,0,1](fermat) == z2;
rvgauss FermatGraph[[3,4,5]](fermat);
rvdual FermatDual[[3,4,5]](fermat);
// The three minors identify the dual coordinates with the Frobenius image.
assert rvInIdeal[FermatGraph](add(rvmon[[0,1,2,3,4,5]]{2,0,0,0,1,0},rvmon[[0,1,2,3,4,5]]{0,2,0,1,0,0})) == S.hi;
assert rvInIdeal[FermatGraph](add(rvmon[[0,1,2,3,4,5]]{2,0,0,0,0,1},rvmon[[0,1,2,3,4,5]]{0,0,2,1,0,0})) == S.hi;
assert rvInIdeal[FermatGraph](add(rvmon[[0,1,2,3,4,5]]{0,2,0,0,0,1},rvmon[[0,1,2,3,4,5]]{0,0,2,0,1,0})) == S.hi;
rvgroebner[FermatDual];
