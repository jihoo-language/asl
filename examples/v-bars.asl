qmodule Flash[0, 1] {
    a1:2, b1:1, a2:4, b2:3, a3:6, b3:5, a4:8, b4:7
} {
    0(a1)=>[b1], 0(a2)=>[b2], 0(a3)=>[b3], 0(a4)=>[b4],
    1(a2)=>[b1], 1(a3)=>[b2], 1(a4)=>[b3]
};

let cycle = vbarRep[One, Flash, [0, 1], 7, 4, 0];
let null = vbarNull[One, Flash, [0, 1], 7, 4, 0];
let generator = qClass(cycle);
assert qD(cycle) == zero[VKoszul(One, Flash, [0, 1], 1, 7)];
assert qD(null) == vRatio[4](cycle);
assert vRatio[3](generator) == vextgen[One, Flash, [0, 1], 0, 1, 0];
assert vRatio[4](generator) == zero[VExt(One, Flash, [0, 1], 0, -1)];
vbars[One, Flash, [0, 1]];
vbardim[One, Flash, [0, 1], 0, 1];
vbardim[One, Flash, [0, 1], 0, -1];
