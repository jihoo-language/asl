// Both hives in the adjoint tensor-square slice
rvHiveCount[[0,1,2],[0,1,2],[1,2,3]];
let low = rvHiveWitness[[0,1,2],[0,1,2],[1,2,3]];
let high = rvHive[[0,1,2],[0,1,2],[1,2,3]]{5};
assert low == rvHive[[0,1,2],[0,1,2],[1,2,3]]{4};
show low;
rvHiveAt[1,1](high);
card \h : RVHive([0,1,2],[0,1,2],[1,2,3]) => S.hi;
// Scaling all three highest weights by seven gives eight lattice points.
rvHiveCount[[0,7,14],[0,7,14],[7,14,21]];
// Independent determinant shifts by -3 and 2 preserve multiplicity.
let shifted = rvHiveWitness[[-3,-2,-1],[2,3,4],[0,1,2]];
rvHiveAt[1,1](shifted);
rvHiveCount[[-3,-2,-1],[2,3,4],[0,1,2]];
rvHorn[[-3,-2,-1],[2,3,4],[0,1,2]];
// Matching totals alone do not ensure that the slice is inhabited.
rvHorn[[0,0],[0,2],[1,1]];
rvHorn[[0,0],[0,2],[1,2]];
card \h : RVHive([0,0],[0,2],[1,1]) => S.hi;
