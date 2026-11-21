let u = rvSphBasis[[0,1]];
let q2 = rvSphQ[2];
let i2 = rvSphId[2];
let square = rvSphMul(u,u);
let expected = add(rvSphBasis[[0,2]],rvSphMul(add(q2,i2),rvSphBasis[[1,1]]));
assert square == expected;
rvSphCoeff[[0,2]](square);
rvSphCoeff[[1,1]](square);
rvSphAt[4,[1,1]](square);

let a = rvSphBasis[[0,1,1]];
let b = rvSphBasis[[0,0,1]];
let together = rvSphMul(a,b);
assert together == rvSphMul(b,a);
rvSphCoeff[[1,1,1]](together);
rvSphAt[2,[1,1,1]](together);

// Iwahori volume one gives an extra factor [GL3(O):I].
let iq = rvHeckeQ[3];
let ii = rvHeckeId[3];
let poincare = rvHeckeMul(add(ii,iq),add(add(ii,iq),rvHeckeMul(iq,iq)));
assert rvSphIwahori(rvSphId[3]) == rvHeckeKL[[2,1,0]];
assert rvHeckeMul(rvSphIwahori(a),rvSphIwahori(b)) ==
       rvHeckeMul(poincare,rvSphIwahori(together));

let c = rvSphBasis[[0,0,2]];
let repeat = rvSphMul(c,c);
rvSphCoeff[[0,2,2]](repeat);
rvSphAt[8,[0,2,2]](repeat);

let central = rvSphBasis[[-4,-4]];
assert rvSphMul(central,u) == rvSphBasis[[-4,-3]];
assert rvSphMul(u,central) == rvSphBasis[[-4,-3]];
assert add(neg(together),together) == zero[RVSphHecke(3)];
