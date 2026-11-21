qmodule Cusp[0, 1, 2] {
    b1:0, b2:2, b3:4, a1:3, a2:5, a3:7
} {
    1(a1)=>[b1], 1(a2)=>[b2], 1(a3)=>[b3],
    0(a1)=>[b2], 0(a2)=>[b3], 2(a3)=>[b1]
};

rvpath Axis[[0,1,2], 7, 7, [[],[0,1],[1]]];
rvpath Curve[[0,1,2], 7, 7, [[0,0,0,1],[0,0,1],[1]]];

let corrected = add(jetterm[One,0,0,Cusp,3,0,Axis,0],
                    jetterm[One,0,0,Cusp,7,0,Axis,1]);
let d2 = jetterm[One,0,0,Cusp,4,0,Axis,2];
assert jetD(corrected) == d2;
assert jetD(d2) == zero[JetHom(One,Cusp,Axis,1)];

let cycle = add(jetterm[One,0,0,Cusp,3,0,Curve,0],
                add(jetterm[One,0,0,Cusp,5,0,Curve,1],
                    jetterm[One,0,0,Cusp,7,0,Curve,2]));
assert jetD(cycle) == zero[JetHom(One,Cusp,Curve,0)];

show jetD(corrected);
rvjet[One,Cusp,Axis];
rvjet[One,Cusp,Curve];
