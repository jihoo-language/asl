rvascover C[[0,1]] {
    add(rvmon[[0,1]]{0,7},rvmon[[0,1]]{4,3}),
    rvmon[[0,1]]{7,0}
};
let eta1 = rvASDeRhamGen[C,11,3];
let eta2 = rvASDeRhamGen[C,11,4];
let x = rvASDeRhamScale[2](eta1);
show rvASF(x);
show rvASV(x);
show rvASDeRhamPair(eta1,eta2);
show rvASDeRhamPair(rvASF(x),eta1);
show rvASDeRhamPair(x,rvASV(eta1));
assert rvASDeRhamPair(rvASF(x),eta1) == rvFieldElt[11,4];
assert rvASDeRhamPair(x,rvASV(eta1)) == rvFieldElt[11,2];
assert rvASDeRhamPair(x,x) == zero[BinaryScalar(11)];
assert rvASF(add(x,eta2)) == add(rvASF(x),rvASF(eta2));
assert rvASV(add(x,eta2)) == add(rvASV(x),rvASV(eta2));
assert add(rvFieldElt[11,2],rvFieldElt[11,4]) == rvFieldElt[11,6];
let omega = rvASDiffScale[2](rvASOmega[C,11,1]);
assert rvASV(rvASDeRhamIncl(omega)) == rvASDeRhamIncl(rvASCartier(omega));
rvASEO[C];
