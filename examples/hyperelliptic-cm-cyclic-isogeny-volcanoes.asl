rvcmorder Top[-23];
rvcmorder Mid[-92];
rvcmorder Bottom[-368];

let source = rvCMPrincipal[Top];
let primeClass = rvCMClass[Top,2,1,3];
let horizontal = rvCMNeighbor[Top,2,0](source);
assert horizontal == primeClass;
assert rvCMNeighbor[Top,2,1](source) == rvCMInverse(primeClass);

let first = rvCMNeighbor[Mid,2,2](source);
let second = rvCMNeighbor[Bottom,2,2](first);
assert first == rvCMPrincipal[Mid];
assert second == rvCMPrincipal[Bottom];
assert rvCMNeighbor[Top,2,1](first) == source;
assert rvCMNeighbor[Mid,2,0](second) == first;

rvCMNeighbors[2](source);
rvCMNeighborData[2,1](first);
rvCMNeighbors[2](first);
rvCMNeighborCounts[2](second);

rvcmorder Gaussian[-4];
rvcmorder Eisenstein[-3];
rvCMNeighbors[2](rvCMPrincipal[Gaussian]);
rvCMNeighbors[2](rvCMPrincipal[Eisenstein]);
