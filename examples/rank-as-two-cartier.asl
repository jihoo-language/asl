rvastwo K[5] {(1,1),(-1,1)};
rvastwo Twist[5] {(1,1),(-1,1),(0,1)};
rvASTwoBasis[K];
show rvASTwoCartier(rvASTwoOmega[K,-1,2]);
show rvASTwoCartier(rvASTwoOmega[K,-1,3]);
rvASTwoNewton[K];
rvASTwoNewton[Twist];
rvASTwoOrdinaryRanks[K];
rvASTwoOrdinaryRanks[Twist];
rvASTwoHasseWitt[K];
rvASTwoHasseWitt[Twist];

rvastwo Wild[5] {(4,1),(-2,1)};
rvASTwoCartierRank[Wild,0];
rvASTwoCartierRank[Wild,1];
rvASTwoCartierRank[Wild,2];
rvASTwoStableRank[Wild];
rvASTwoaNumber[Wild];
rvASTwoOrdinaryRanks[Wild];
rvASTwoHasseWitt[Wild];
show rvASTwoOrdinaryGen[Wild,2];

// This measured cyclic vector for the ordinary deck action gives dx/x.
let v = rvASTwoCartierIter[12](rvASTwoOmega[Wild,-1,3]);
let d1 = add(rvASTwoTranslate[1](v),pscale[4](v));
let d2 = add(rvASTwoTranslate[1](d1),pscale[4](d1));
let d3 = add(rvASTwoTranslate[1](d2),pscale[4](d2));
assert d3 == rvASTwoOmega[Wild,-1,0];
