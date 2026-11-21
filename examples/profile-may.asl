// A(0) has a single cobar word in each diagonal degree.
let a0h0 = mayangen[0, 3, 1, 1, 0, 0];
assert mayMul(a0h0, a0h0) == mayangen[0, 3, 2, 2, 0, 0];
assert mayD(a0h0) == zero[MayAn(0, 3, 2, 1, 3)];
mayandim[0, 8, 4, 4, 0];

// A(1) retains the first differential of xi_2.
let xi2 = mayClass[1, 0](cword{[0, 1]});
let atOneXi2 = mayClass[1, 0](resA[1](cword{[0, 1]}));
assert atOneXi2 == resA[1](xi2);
assert mayD(atOneXi2) == mayangen[1, 1, 2, 3, 1, 0];
assert mayD(atOneXi2) == resA[1](mayD(xi2));
show mayD(atOneXi2);

// The full d_7 on b_{2,0}^2 is nonzero. Its target vanishes in A(1).
let full = maygen[7, 4, 12, 0, 0];
let fullD = mayD(full);
assert fullD == maygen[7, 5, 12, 7, 0];
let restricted = resA[1](full);
assert restricted == mayangen[1, 7, 4, 12, 0, 0];
assert mayD(restricted) == resA[1](fullD);
assert mayD(restricted) == zero[MayAn(1, 7, 5, 12, 7)];

// The A(1) lift survives to v_1^4 and its restriction to A(0) is zero.
let atEight = mayClass[8, 0](mayRep(restricted));
assert atEight == mayangen[1, 8, 4, 12, 0, 0];
assert mayClass[20, 0](mayRep(atEight)) == mayangen[1, 20, 4, 12, 0, 0];
let v1four = extClass(mayRep(atEight));
assert v1four == extangen[1, 4, 12, 0];
assert resA[0](restricted) == resA[0](full);
assert resA[0](full) == zero[MayAn(0, 7, 4, 12, 0)];
assert resA[0](v1four) == zero[ExtAn(0, 4, 12)];
show fullD;
show mayD(restricted);
show atEight;
show v1four;
mayandim[1, 20, 4, 12, 0];
