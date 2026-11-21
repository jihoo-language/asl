qmodule Cusp[0,1,2] {
    b1:0, b2:2, b3:4, a1:3, a2:5, a3:7
} {
    1(a1)=>[b1], 1(a2)=>[b2], 1(a3)=>[b3],
    0(a1)=>[b2], 0(a2)=>[b3], 2(a3)=>[b1]
};

let rawB1 = qterm[One,0,0,Cusp,0,0,[0,1,2]]{0,0,0};
let rawB2 = qterm[One,0,0,Cusp,2,0,[0,1,2]]{0,0,0};
let b1 = rvSpec[7,[1,3,1]](qClass(rawB1));
let b2 = rvSpec[7,[1,3,1]](qClass(rawB2));
assert b2 == rvScale[3](b1);

// Squaring sends the affine point (1,3,1) to (1,2,1).
let conjugate = rvFrob[1](b2);
assert conjugate == rvScale[2](rvFrob[1](b1));
assert conjugate == rvSpec[7,[1,2,1]](qClass(rawB2));
assert rvFrob[2](b2) == b2;
show conjugate;

// The root 6 in F16 sends the encoded F4 scalar 3 to 7.
let extended = rvExtend[19,6](b2);
assert extended == rvSpec[19,[1,7,1]](rawB2);
assert rvDescend[7,6](extended) == b2;
assert rvExtend[19,7](b2) == rvFrob[1](extended);
assert rvExtend[19,6](conjugate) == rvFrob[1](extended);
show extended;

let unit = rvId[One,[0],19,[1]];
let relative = rvTrace[7,6](rvScale[2](unit));
assert relative == rvId[One,[0],7,[1]];
assert rvTrace[7,6](unit) == zero[RVHom(One,One,[0],7,[1])];
assert rvTrace[3,1](rvScale[2](relative)) == rvId[One,[0],3,[1]];
show relative;

// Both ways down the degree-four tower give the same nonzero trace.
let cubic = rvScale[8](unit);
assert rvTrace[3,1](cubic) == rvTrace[3,1](rvTrace[7,6](cubic));
assert rvTrace[3,1](cubic) == rvId[One,[0],3,[1]];

// Frobenius and field extension also commute with Laurent specialization.
let inverse = vShift[-1](vLocal(qid[One,[1]]));
let atThree = rvSpec[7,[3]](inverse);
assert rvFrob[1](atThree) == rvSpec[7,[2]](inverse);
assert rvExtend[19,6](atThree) == rvSpec[19,[7]](inverse);
