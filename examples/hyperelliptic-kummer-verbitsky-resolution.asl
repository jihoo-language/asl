rvhodgewick B[21]{0,1}{[2,0],[0,-2]};
rvVerbitskyResolution[B,2];
rvVerbitskyResolutionEntry[B,2,1,0,0];
rvVerbitskyResolutionEntry[B,2,2,0,0];
rvVerbitskyResolutionEntry[B,2,2,4,0];
rvVerbitskyPfaffian[B,2];

let g=rvVerbitskyResolutionBasis[B,2,2,0];
let y0=rvVerbitskyResolutionBasis[B,2,1,0]{0,1};
let y1=rvVerbitskyResolutionBasis[B,2,1,1]{0,1};
let x2=rvVerbitskyResolutionBasis[B,2,1,2]{1};
let x3=rvVerbitskyResolutionBasis[B,2,1,3]{1};
assert rvVerbitskyResolutionDiff(g)==add(add(y0,neg(y1)),add(neg(x2),x3));
assert rvVerbitskyResolutionDiff(rvVerbitskyResolutionDiff(g))==zero[RVVerbitskyFree(B,2,0,4)];
let top=rvVerbitskyResolutionBasis[B,2,3,0];
assert rvVerbitskyResolutionDiff(rvVerbitskyResolutionDiff(top))==zero[RVVerbitskyFree(B,2,1,8)];
