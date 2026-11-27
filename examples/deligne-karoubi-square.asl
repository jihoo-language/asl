// Seven primitive projectors split the twofold permutation tensor power.
rvkaroubi A(rvDeligneSquareProjector[0]);
rvkaroubi B(rvDeligneSquareProjector[1]);
rvkaroubi C(rvDeligneSquareProjector[2]);
rvkaroubi G(rvDeligneSquareProjector[3]);
rvkaroubi H(rvDeligneSquareProjector[4]);
rvkaroubi K(rvDeligneSquareProjector[5]);
rvkaroubi L(rvDeligneSquareProjector[6]);
rvKaroubiDimension[A];
rvKaroubiDimension[H];
rvKaroubiDimension[K];
rvKaroubiDimension[L];

// The off-diagonal units of the trivial and standard matrix blocks.
let ag = rvKaroubiMap[A,G](rvDeligneSquareUnit[3,0]);
let ga = rvKaroubiMap[G,A](rvDeligneSquareUnit[0,3]);
assert rvKaroubiCompose(ga,ag) == rvKaroubiIdentity[A];
assert rvKaroubiCompose(ag,ga) == rvKaroubiIdentity[G];
let bc = rvKaroubiMap[B,C](rvDeligneSquareUnit[2,1]);
let cb = rvKaroubiMap[C,B](rvDeligneSquareUnit[1,2]);
let bh = rvKaroubiMap[B,H](rvDeligneSquareUnit[4,1]);
let hb = rvKaroubiMap[H,B](rvDeligneSquareUnit[1,4]);
assert rvKaroubiCompose(cb,bc) == rvKaroubiIdentity[B];
assert rvKaroubiCompose(bh,hb) == rvKaroubiIdentity[H];
assert rvKaroubiCompose(bh,cb) == rvKaroubiMap[C,H](rvDeligneSquareUnit[4,2]);
rvKaroubiTrace(rvKaroubiCompose(hb,bh));

// At N=3 the whole symmetric residual projector specializes to zero.
rvPartitionZeroAt[3](rvKaroubiUnderlying(rvKaroubiIdentity[K]));
rvPartitionTraceAt[4](rvDeligneSquareProjector[5]);
rvPartitionTraceAt[5](rvDeligneSquareProjector[5]);
