rvverbitskytate T{5};
rvmf P[T]{2};
rvlgdefect D[T]{0,1};
rvLGDefect[D];
rvlgdefectpush A[D,P];
rvLGDefectPush[A];
rvlgdefectresult Q[A];
rvMFObject[Q];

let s=add(rvCyclo[5,0],rvCyclo[5,1]);
let si=neg(add(rvCyclo[5,1],rvCyclo[5,3]));
assert rvLGDefectCoeff[D,0,1]==neg(s);
assert rvLGDefectMatrix[A,0,0,3,2]==rvCycloRat[1](s);
assert rvLGDefectMatrix[A,0,2,3,3]==rvCycloRat[1](si);
let y=rvMFBasis[P,P,0,0,0,1];
assert rvLGDefectMap[A,A,1,0,2](y)==rvCycloRat[1](si);
assert rvLGDefectMap[A,A,0,1,0](y)==rvCycloRat[1](neg(si));
assert rvLGDefectMap[A,A,1,1,1](y)==rvCycloRat[1](si);
assert rvLGDefectMap[A,A,1,1,2](rvMFIdentity[P])==zero[RVCycloQ(5)];

assert rvLGDefectBulk[D,0]==add(rvCyclo[5,0],rvCyclo[5,-1]);
assert rvLGDefectBulk[D,2]==add(rvCyclo[5,0],rvCyclo[5,-3]);
assert rvLGDefectReverse[D,0]==s;
show rvLGDefectBulk[D,0];
rvLGDefectFusion[D,D];

rvverbitskytate S{4};
rvmf Ptwo[S]{2};
rvlgdefect Even[S]{0,2};
rvlgdefectpush E[Even,Ptwo];
rvLGDefectPush[E];
assert rvLGDefectBulk[Even,0]==zero[RVCyclo(4)];
assert rvLGDefectBulk[Even,1]==rvCycloScale[2](rvCyclo[4,0]);
assert rvLGDefectMap[E,E,1,0,0](rvMFBasis[Ptwo,Ptwo,0,0,0,1])==rvCycloRat[1](rvCyclo[4,0]);
