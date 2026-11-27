// The figure-eight improves on the raw oriented Lee cycle.
khbraid Figure[3,[1,-2,1,-2]];
leeS[Figure];
leeDtProfile[Figure];
leeDt[Figure,0,1];
leeDt[Figure,1,2];
leeDt[Figure,1,1];
leeDt[Figure,3,2];
leeDt[Figure,2,1];

let orientation = leeOrientation[Figure];
let optimal = leeDtRep[Figure,1,2];
let primitive = leeDtNull[Figure,1,2];
let obstruction = leeDtDual[Figure,1,2];
let one = rvCycloRat[1](rvCyclo[1,0]);
assert leeD(orientation) == zero[LeeC(Figure,1,-13)];
assert leeD(optimal) == zero[LeeC(Figure,1,-13)];
assert add(orientation,leeD(primitive)) == optimal;
assert leeDtDualD(obstruction) == zero[LeeDual(Figure,-1)];
assert leeDtPair(obstruction,orientation) == one;
assert leeDtPair(obstruction,optimal) == one;
leeDtWeight[1,2](orientation);
leeDtWeight[1,2](optimal);
leeDtWeight[1,2](obstruction);

// The positive trefoil has d_t = 3 - 2 |1-t|.
khbraid Trefoil[2,[1,1,1]];
khbraid Mirror[2,[-1,-1,-1]];
leeS[Trefoil];
leeS[Mirror];
leeDtProfile[Trefoil];
leeDtProfile[Mirror];

// This four-braid has genuine rational breakpoints on both sides of t=1.
khbraid A[4,[3,-2,-2,3,3,2,-3,-1,2,1,1]];
leeS[A];
leeDtProfile[A];
leeDt[A,3,5];
leeDt[A,2,3];
let oA = leeOrientation[A];
let repA = leeDtRep[A,2,3];
let primitiveA = leeDtNull[A,2,3];
let obstructionA = leeDtDual[A,2,3];
assert leeD(repA) == zero[LeeC(A,1,-27)];
assert add(oA,leeD(primitiveA)) == repA;
assert leeDtDualD(obstructionA) == zero[LeeDual(A,-1)];
assert leeDtPair(obstructionA,oA) == rvCycloRat[1](rvCyclo[1,0]);
leeDtWeight[2,3](repA);
leeDtWeight[2,3](obstructionA);
