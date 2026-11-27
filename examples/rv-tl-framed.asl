// Colors and framings follow components ordered by their least strand.
// The level-five Hopf link colors its two components by 1 and 2.
let hopf = rvTLFramed[5,2,[1,1],[1,2],[0,0]];
let d1 = neg(add(rvCyclo[28,2],rvCyclo[28,-2]));
let d3 = add(d1,neg(add(rvCyclo[28,6],rvCyclo[28,-6])));
// The fusion channels are 1 and 3. Their full twist ratios are A^-8 and A^4.
let expected = add(rvCycloMul(d1,rvCyclo[28,-8]),rvCycloMul(d3,rvCyclo[28,4]));
assert hopf == rvCycloRat[1](expected);
show hopf;
assert rvTLFramed[5,2,[-1,-1],[1,2],[0,0]] == rvCycloQGalois[-1](hopf);

// The two last letters stabilize the second component separately.
assert rvTLFramed[5,3,[1,1,2],[1,2],[0,0]] == hopf;
assert rvTLFramed[5,3,[1,1,-2],[1,2],[0,0]] == hopf;
rvTLLinking[3,[1,-2,1],[8,9]];
