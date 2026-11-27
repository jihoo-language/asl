rvunram B[3,7]{0,1};
rvcyclotomiclocal K[B,1];
rvcubicdiv D[K,2,5];

rvCubicTreeSphere[D,1];
rvCubicTreeSphere[D,2];
rvCubicTreeSphere[D,3];
let o = rvCubicTreeRoot[D];
let x = rvCubicTreeAt[D,1,3];
let y = rvCubicTreeChild[11](x);
let z = rvCubicTreeChild[9](x);
assert rvCubicTreeAdjacent(x,y) == S.hi;
assert rvCubicTreeLca(y,z) == x;
rvCubicTreeLevel(y);
rvCubicTreeDistance(y,z);
rvCubicTreePath(y);

let theta = rvCubicDivTheta[D];
let pi = rvCubicDivPi[D];
let one = rvCubicDivElem[D]([1]);
let g = rvCubicTreeMatrix(one,theta,pi,one);
let sigma = rvCubicTreePi[D];
let w = rvCubicTreeApply(rvCubicTreeCompose(g,sigma),y);
assert w == rvCubicTreeApply(g,rvCubicTreeApply(sigma,y));
assert rvCubicTreeApply(rvCubicTreeInverse(rvCubicTreeCompose(g,sigma)),w) == y;
assert rvCubicTreeFixes(sigma,o) == S.hi;
rvCubicTreeDistance(w,rvCubicTreeApply(rvCubicTreeCompose(g,sigma),z));

let target = rvCubicTreeTranslate(o);
assert rvCubicTreeAdjacent(o,target) == S.hi;
assert rvCubicTreeTranslateBack(rvCubicTreeTranslate(y)) == y;
rvCubicTreeHecke[1](x,x);
rvCubicTreeHecke[1](y,z);
rvCubicTreeHecke[2](x,y);
