rvhodgewick B[21]{0,1}{[2,0],[0,-2]};
rvVerbitskyDeformationTable[B,2];
rvVerbitskyDeformationDims[B,2,-2];

let g=add(rvVerbitskyTangentBasis[B,2,-1,0],rvVerbitskyTangentBasis[B,2,-1,3]);
let t=rvVerbitskyDeformationClass(g);
assert rvVerbitskyDeformationClass(rvVerbitskyDeformationRepresentative(t))==t;
show g;
rvVerbitskyObstruction(g);

rvhodgewick C[22]{0,1}{[2]};
rvVerbitskyDeformationTable[C,2];
rvVerbitskyObstruction(rvVerbitskyTangentBasis[C,2,-3,0]);
