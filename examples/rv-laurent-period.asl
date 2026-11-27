// (1+x)(1+y)(1+x+y)/(xy), the Apéry period for zeta(2).
rvperiod Apery[2] {
    ([-1,-1],1), ([-1,0],2), ([-1,1],1), ([0,-1],2),
    ([0,0],3), ([0,1],1), ([1,-1],1), ([1,0],1)
};
rvCTExact[Apery,0];
rvCTExact[Apery,1];
rvCTExact[Apery,2];
rvCTExact[Apery,3];
let start = rvCTInitial[Apery,3,2];
let slide = rvCTSection[1](start);
let seven = rvCTRead(rvCTSection[2](slide));
assert seven == rvCTMod[Apery,3,2,7];
assert seven == rvCTResidue[3,2,6];
show seven;
show rvCTMod[Apery,3,3,1350851717672992088];
show rvCTMod[Apery,2,4,38685626227668133590597631];

// Here even the scalar Lucas digit rule over F2 loses the state.
rvperiod Carry[1] {([0],1), ([2],1), ([-1],1)};
let digit = \x : RVCTState(Carry,2,1) => rvCTSection[1](x);
let mid = digit(rvCTInitial[Carry,2,1]);
show mid;
assert rvCTRead(mid) == rvCTResidue[2,1,1];
assert rvCTRead(digit(mid)) == rvCTResidue[2,1,0];
rvCTExact[Carry,3];
