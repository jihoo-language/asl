// The diagram with isolated top and bottom vertices is the all-ones matrix.
let j = rvPartitionDiagram[1]{0,1};
let diagonal = rvPartitionIdentity[1];
let jj = rvPartitionCompose(j,j);
assert jj == rvPartitionScale[1,1](j);
rvPartitionTrace(jj);
rvPartitionEntry[3,[0],[1]](jj);

// Exact equality patterns are obtained by Möbius inversion.
let unequal = rvPartitionOrbit[1]{7,9};
assert unequal == add(j,neg(diagonal));
rvPartitionEntry[3,[0],[1]](unequal);
rvPartitionEntry[3,[1],[1]](unequal);

// The swap moves the all-ones factor from one tensor coordinate to the other.
let swap = rvPartitionDiagram[2]{0,1,1,0};
let a = rvPartitionTensor(diagonal,j);
show rvPartitionCompose(swap,a);
show rvPartitionCompose(a,swap);
rvPartitionDimensions[2,2];
rvPartitionCentralizerCheck[2,2];
