let x = rvmon[[0,1,2]]{1,0,0};
let y = rvmon[[0,1,2]]{0,1,0};
let z = rvmon[[0,1,2]]{0,0,1};
let cubic = add(rvmon[[0,1,2]]{0,3,0}, rvmon[[0,1,2]]{2,0,1});
rvideal Cusp[[0,1,2]] { cubic };
rvrees Blow[Cusp,[3,4]] { x,y };
rvnormalcone Cone[Cusp,[3,4]] { x,y };

// U and V represent x*T and y*T, respectively.
let U = rvmon[[0,1,2,3,4]]{0,0,0,1,0};
let V = rvmon[[0,1,2,3,4]]{0,0,0,0,1};
let X = rvmon[[0,1,2,3,4]]{1,0,0,0,0};
let Y = rvmon[[0,1,2,3,4]]{0,1,0,0,0};
let Z = rvmon[[0,1,2,3,4]]{0,0,1,0,0};
let u2 = rvMul(U,U);
let v2 = rvMul(V,V);
let y2 = rvMul(Y,Y);
assert rvInIdeal[Blow](add(rvMul(X,V),rvMul(Y,U))) == S.hi;
assert rvInIdeal[Blow](add(rvMul(y2,V),rvMul(rvMul(X,Z),U))) == S.hi;
assert rvInIdeal[Blow](add(rvMul(Y,v2),rvMul(Z,u2))) == S.hi;
assert rvInIdeal[Blow](add(rvMul(Y,y2),rvMul(rvMul(X,X),Z))) == S.hi;
assert rvInIdeal[Blow](U) == S.lo;
rvgroebner[Blow];

// In the chart z != 0, the exceptional normal direction U is doubled.
rvsaturate Chart[Cone](Z);
assert rvInIdeal[Cone](u2) == S.lo;
assert rvInIdeal[Cone](rvMul(Z,u2)) == S.hi;
assert rvInIdeal[Chart](u2) == S.hi;
assert rvInIdeal[Chart](U) == S.lo;
rvgroebner[Chart];
rvhilbert[Chart,1000];

rvideal Plane[[0,1]] {};
rvrees Quadratic[Plane,[2,3,4]] {
    rvmon[[0,1]]{2,0}, rvmon[[0,1]]{1,1}, rvmon[[0,1]]{0,2}
};
// A, B, C represent x^2*T, x*y*T, y^2*T. At x=y=0 the
// two linear syzygies vanish and the exceptional fiber is AC+B^2=0.
let acb2 = add(rvmon[[0,1,2,3,4]]{0,0,1,0,1},rvmon[[0,1,2,3,4]]{0,0,0,2,0});
assert rvInIdeal[Quadratic](acb2) == S.hi;
rvgroebner[Quadratic];
