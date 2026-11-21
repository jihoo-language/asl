// The transition chooses the higher-filtration corrections automatically.
let b1 = mayClass[1, 0](cword{[0, 1], [0, 1]});
let b2 = mayNext(b1);
let b3 = mayNext(b2);
assert b3 == maygen[3, 2, 6, 0, 0];
assert mayNext(mayPrev(b3)) == b3;
show b3;

// b_{2,0} supports a d_3, and its image vanishes on E_4.
let relation = mayD(b3);
show relation;
assert mayNext(relation) == zero[May(4, 3, 6, 3)];
show mayNext(relation);

// Squaring h_1 commutes with the transition.
let h = mayClass[3, 1](cword{[2]});
let h4 = mayNext(h);
assert mayNext(mayMul(h, h)) == mayMul(h4, h4);
assert mayNext(mayPrev(h4)) == h4;
