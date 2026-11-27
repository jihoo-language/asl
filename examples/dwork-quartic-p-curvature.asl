rvDworkPCurvature[3,2];
rvDworkPCurvature[11,3];

let regular=rvDworkPCurvatureAt[11,3,1];
assert regular == rvHWMat[11,3,1]{7,4,0,2,3,4,1,3,1};
show regular;
rvDworkPCurvatureFiber[11,3,1];

let double=rvDworkPCurvatureAt[11,3,-1];
assert double == rvHWMat[11,3,1]{0,0,0,9,0,0,9,9,0};
let second=rvDworkPCurvatureAt[11,3,-1,2];
assert second == rvHWCompose(double,double);
assert second == rvHWMat[11,3,1]{0,0,0,0,0,0,4,0,0};
assert rvDworkPCurvatureAt[11,3,-1,1000000] == rvHWMat[11,3,1]{0,0,0,0,0,0,0,0,0};
rvDworkPCurvatureFiber[11,3,-1];
rvDworkPCurvatureFiber[13,3,9];

rvunram K[3,1]{1,0,1};
let x=rvUnramElem[K]{0,1};
let E=rvDworkPCurvatureClosed[3](x);
assert rvHWWEntry[1,2](E) == rvUnramElem[K]{2,2};
assert rvDworkPCurvatureClosed[3,2](x) == rvHWWCompose(E,E);
show E;
rvDworkPCurvatureClosedFiber[3](x);
