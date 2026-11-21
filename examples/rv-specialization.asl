qmodule Cusp[0,1,2] {
    b1:0, b2:2, b3:4, a1:3, a2:5, a3:7
} {
    1(a1)=>[b1], 1(a2)=>[b2], 1(a3)=>[b3],
    0(a1)=>[b2], 0(a2)=>[b3], 2(a3)=>[b1]
};

// Over F4 the point (1,3,1) lies on the cusp.
let rawB1 = qterm[One,0,0,Cusp,0,0,[0,1,2]]{0,0,0};
let rawB2 = qterm[One,0,0,Cusp,2,0,[0,1,2]]{0,0,0};
let rawB3 = qterm[One,0,0,Cusp,4,0,[0,1,2]]{0,0,0};
let b1 = rvSpec[7,[1,3,1]](qClass(rawB1));
let b2 = rvSpec[7,[1,3,1]](qClass(rawB2));
let b3 = rvSpec[7,[1,3,1]](rawB3);
assert b1 == rvgen[One,Cusp,[0,1,2],7,[1,3,1],0];
assert b2 == rvScale[3](b1);
assert b3 == rvScale[2](b1);
assert rvSpec[7,[1,3,1]](qid[Cusp,[0,1,2]]) == rvId[Cusp,[0,1,2],7,[1,3,1]];
assert ecomp(rvId[Cusp,[0,1,2],7,[1,3,1]], b3) == b3;
show b1;
show b2;
show b3;

// The second basis class is represented by 2*a1 + 3*a2 + a3.
let rawFromA1 = qterm[Cusp,3,0,One,0,0,[0,1,2]]{0,0,0};
let fromA1 = rvSpec[7,[1,3,1]](rawFromA1);
assert fromA1 == rvScale[2](rvgen[Cusp,One,[0,1,2],7,[1,3,1],1]);
assert rvSpec[7,[1,3,1]](ecomp(rawB1, rawFromA1)) == ecomp(b1, fromA1);
assert ecomp(fromA1, b1) == zero[RVHom(One,One,[0,1,2],7,[1,3,1])];

// A Koszul boundary becomes zero. At the off-curve point all homology vanishes.
let a1 = qterm[One,0,0,Cusp,3,0,[0,1,2]]{0,0,0};
assert rvSpec[7,[1,3,1]](qD(a1)) == zero[RVHom(One,Cusp,[0,1,2],7,[1,3,1])];
let off = zero[RVHom(One,Cusp,[0,1,2],7,[2,3,1])];
assert rvSpec[7,[2,3,1]](rawB1) == off;
assert rvSpec[7,[2,3,1]](rawB2) == off;
assert rvSpec[7,[2,3,1]](rawB3) == off;
assert rvSpec[7,[2,3,1]](qid[Cusp,[0,1,2]]) == rvId[Cusp,[0,1,2],7,[2,3,1]];
show rvId[Cusp,[0,1,2],7,[2,3,1]];

// Laurent powers evaluate at the same affine point as polynomial powers.
let unit = vLocal(qid[One,[1]]);
let atThree = rvId[One,[1],7,[3]];
assert rvSpec[7,[3]](unit) == atThree;
assert rvSpec[7,[3]](vShift[1](unit)) == rvScale[3](atThree);
assert rvSpec[7,[3]](vShift[-1](unit)) == rvScale[2](atThree);
assert rvSpec[7,[3]](qRep(vShift[-2](unit))) == rvScale[3](atThree);
