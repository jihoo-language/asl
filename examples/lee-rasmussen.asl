khbraid Trefoil[2,[1,1,1]];
khbraid Mirror[2,[-1,-1,-1]];
khbraid FigureEight[3,[1,-2,1,-2]];
khbraid Long[3,[1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2]];
leeS[Trefoil];
leeS[Mirror];
leeS[FigureEight];
leeLevels[Long,0];
leeRank[Long,0,16];
leeRank[Long,0,14];
leeRank[Long,0,13];
leeS[Long];

let hi = leeRep[Long,0,15,0];
let lo = leeRep[Long,0,13,0];
assert leeD(hi) == zero[LeeC(Long,1,15)];
assert leeD(lo) == zero[LeeC(Long,1,13)];
leeFiltration(hi);
leeFiltration(lo);
leeIndependent(hi,lo);

// The low-quantum cochain X tensor X can be moved to a higher level.
// A specified rational multiple of a boundary has an exact preimage.
khbraid Neg[2,[-1]];
let x = leeLift(khTerm[Neg,[],[0]]);
let raw = leeLift(khTerm[Neg,[0],[0,1]]);
let raised = leeRep[Neg,0,1,0];
assert add(raw,leeFilter[-3](raised)) == leeD(x);
leeFiltration(raw);
let half = leeScale[1,2](x);
let boundary = leeD(half);
let primitive = leeNull(boundary);
assert primitive == leeFilter[-6](half);
assert leeD(primitive) == leeFilter[-6](boundary);
leeFiltration(boundary);
