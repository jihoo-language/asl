qmodule Cusp[0, 1, 2] {
    b1:0, b2:2, b3:4, a1:3, a2:5, a3:7
} {
    1(a1)=>[b1], 1(a2)=>[b2], 1(a3)=>[b3],
    0(a1)=>[b2], 0(a2)=>[b3], 2(a3)=>[b1]
};

rvpath Axis[[0,1,2], 7, 7, [[],[0,1],[1]]];
rvpath Transverse[[0,1,2], 7, 7, [[0,1],[],[1]]];
rvpath Curve[[0,1,2], 7, 7, [[0,0,0,1],[0,0,1],[1]]];

let a1 = jetterm[One,0,0,Cusp,3,0,Axis,0];
let a2 = jetterm[One,0,0,Cusp,5,0,Axis,0];
let ub2 = jetterm[One,0,0,Cusp,2,0,Axis,1];
let u2b3 = jetterm[One,0,0,Cusp,4,0,Axis,2];
let a1e1 = jetClass[1,0](a1);
let a2e1 = jetClass[1,0](a2);
assert jetD(a1e1) == zero[JetPage(One,Cusp,Axis,1,0,1)];
assert jetD(a2e1) == jetClass[1,1](ub2);
assert jetD(a2e1) == jetgen[One,Cusp,Axis,1,0,1,0];

let a1e2 = jetNext(a1e1);
assert jetPrev(a1e2) == a1e1;
assert jetRep(a1e2) == add(a1,jetterm[One,0,0,Cusp,7,0,Axis,1]);
assert jetClass[2,0](jetRep(a1e2)) == a1e2;
let d2 = jetD(a1e2);
assert d2 == jetClass[2,2](u2b3);
assert d2 == jetgen[One,Cusp,Axis,2,0,2,0];
assert jetNext(d2) == zero[JetPage(One,Cusp,Axis,3,0,2)];
show jetRep(a1e2);
show d2;

let ta = jetClass[1,0](jetterm[One,0,0,Cusp,3,0,Transverse,0]);
let tb = jetClass[1,0](jetterm[One,0,0,Cusp,5,0,Transverse,0]);
assert jetD(ta) == jetClass[1,1](jetterm[One,0,0,Cusp,2,0,Transverse,1]);
assert jetD(tb) == jetClass[1,1](jetterm[One,0,0,Cusp,4,0,Transverse,1]);
jetdim[One,Cusp,Transverse,2,1,0];
jetdim[One,Cusp,Transverse,3,1,0];

let c = add(jetterm[One,0,0,Cusp,3,0,Curve,0],
            add(jetterm[One,0,0,Cusp,5,0,Curve,1],
                jetterm[One,0,0,Cusp,7,0,Curve,2]));
let permanent = jetClass[7,0](c);
assert permanent == jetgen[One,Cusp,Curve,7,1,0,0];
assert jetD(jetRep(permanent)) == zero[JetHom(One,Cusp,Curve,0)];
show permanent;
