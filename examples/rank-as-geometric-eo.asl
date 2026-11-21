let constant = rvmon[[0,1]]{10,0};
let top = rvmon[[0,1]]{0,10};
let middle = rvmon[[0,1]]{4,6};
rvascover A[[0,1]] {
    add(constant,add(top,middle)),
    rvmon[[0,1]]{9,1}
};
rvascover B[[0,1]] {
    add(add(constant,top),add(middle,add(rvmon[[0,1]]{6,4},rvmon[[0,1]]{8,2}))),
    rvmon[[0,1]]{5,5}
};
rvASLPolynomial[A];
rvASLPolynomial[B];
rvASPrank[A];
rvASPrank[B];
rvASaNumber[A];
rvASaNumber[B];
rvASPoleOrders[A];
rvASPoleOrders[B];
rvASEO[A];
rvASEO[B];
