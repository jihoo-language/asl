// A finite graded module over E(Q_0,Q_1).
qmodule Flash[0, 1] {
    a1:2, b1:1, a2:4, b2:3, a3:6, b3:5, a4:8, b4:7
} {
    0(a1)=>[b1], 0(a2)=>[b2], 0(a3)=>[b3], 0(a4)=>[b4],
    1(a2)=>[b1], 1(a3)=>[b2], 1(a4)=>[b3]
};

// The intermediate terms cancel, leaving that fourth power of v_0.
let witness = add(add(qterm[One,0,0,Flash,2,0,[0,1]]{0,3}, qterm[One,0,0,Flash,4,0,[0,1]]{1,2}),
                  add(qterm[One,0,0,Flash,6,0,[0,1]]{2,1}, qterm[One,0,0,Flash,8,0,[0,1]]{3,0}));
let fourth = qterm[One,0,0,Flash,7,0,[0,1]]{4,0};
assert qD(witness) == fourth;
assert qClass(fourth) == zero[QExt(One, Flash, [0,1], 4, 11)];

// The third power is nonzero, so the annihilator has exact order four.
let third = qClass(qterm[One,0,0,Flash,7,0,[0,1]]{3,0});
assert third == qextgen[One, Flash, [0,1], 3, 10, 0];
assert ecomp(qid[Flash,[0,1]], third) == third;
assert ecomp(third, qid[One,[0,1]]) == third;
assert ecomp(third, qClass(qRes[0,1](cword{[1]}))) == qClass(fourth);
show third;
qextdim[One, Flash, [0,1], 4, 11];
