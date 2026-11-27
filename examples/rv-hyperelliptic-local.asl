// Two elliptic normalizations and a one-dimensional torus.
let split = rvHyperArithmetic[7,1]{0,49,98,1,50,2,3,4};
let twist = rvHyperArithmetic[7,3]{0,49,98,1,50,2,3,4};
rvHyperLocal(split);
rvHyperLocalFactors(split);
rvHyperLocal(twist);
rvHyperLocalFactors(twist);
rvHyperNodalPoints[1](split);
rvHyperNodalPoints[2](split);
rvHyperNodalPoints[3](split);
rvHyperNodalPoints[4](split);
rvHyperNodalPoints[5](split);
rvHyperNodalPoints[1](twist);
rvHyperNodalPoints[2](twist);
rvHyperNodalPoints[3](twist);
rvHyperNodalPoints[4](twist);
rvHyperNodalPoints[5](twist);

// Changing the sign on both loops exchanges their eigenvalues.
let mixed = rvHyperArithmetic[7,1]{0,1,2,345,3,16810};
let mixedTwist = rvHyperArithmetic[7,3]{0,1,2,345,3,16810};
rvHyperLocal(mixed);
rvHyperLocal(mixedTwist);
rvHyperFixed(mixed);
rvHyperFixed(mixedTwist);
