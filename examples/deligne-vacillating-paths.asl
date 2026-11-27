// The same final shape (2,1) is reached through two different paths.
let throughTwo = rvDelignePathProjector[3]{[],[],[1],[1],[2],[2],[2,1]};
let throughColumn = rvDelignePathProjector[3]{[],[],[1],[1],[1,1],[1,1],[2,1]};
rvDelignePathEigenvalues[3,6]{[],[],[1],[1],[2],[2],[2,1]};
rvDelignePathEigenvalues[3,6]{[],[],[1],[1],[1,1],[1,1],[2,1]};
rvPartitionTraceAt[6](throughTwo);
rvPartitionTraceAt[6](throughColumn);
rvkaroubi Mixed(throughTwo);
rvKaroubiDimension[Mixed];
assert rvPartitionCompose(throughTwo,throughTwo) == throughTwo;
assert rvPartitionCompose(throughTwo,throughColumn) == zero[RVPartitionLocal(3)];

// A half-step sum of transpositions leaves the distinguished color fixed.
let h = rvDeligneTowerOperator[3,5];
rvPartitionEntry[6,[1,0,2],[0,1,2]](h);
rvPartitionEntry[6,[2,1,2],[0,1,2]](h);
rvDelignePathCounts[3];
rvDelignePathCounts[4];
