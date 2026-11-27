// Three persistent edges and a moving cross-ratio on one residue component.
let A = rvStableLine[7]{[0,1],[11398895185373143,1],[22797790370746286,1],[9387491736542939678792,1],
                          [1,1],[1977326744,1],[2,1],[3,1],[1,0]};
let B = rvStableLine[7]{[0,1],[11398895185373143,1],[34196685556119429,1],[9387491736542939678792,1],
                          [1,1],[1977326744,1],[2,1],[3,1],[1,0]};
let T = rvStableTree(A);
assert rvStableTree(B) == T;
rvStableSplits(T);
rvStableDirections[1](T);
rvStableResidue[1,8](A);
rvStableResidue[1,8](B);
rvStableComponents(A);

let P = rvValMat[7,9]{[0,11398895185373143,22797790370746286,9387491736542939678792,1,1977326744,2,3,1],
                      [1,1,1,1,1,1,1,1,0]};
assert rvStableTree(P) == T;
