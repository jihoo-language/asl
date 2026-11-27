// Project the second and third elementary strands onto color two.
let a = rvTLPath[5,[0,1,0,1]];
let projected = rvTLProject[2,2](a);
assert rvTLProject[2,2](projected) == projected;
assert rvTLE[2](projected) == zero[RVTL(5,3,1)];
let delta = rvCycloRat[1](neg(add(rvCyclo[28,2],rvCyclo[28,-2])));
assert rvTLScale(delta,projected) == add(rvTLScale(delta,a),neg(rvTLE[2](a)));
show projected;

// Color two on the two-strand trefoil at level five has channels 0, 2, 4.
let t = rvTLColoredJones[5,2,2,[1,1,1]];
let d2 = add(rvCyclo[28,4],add(rvCyclo[28,0],rvCyclo[28,-4]));
let d4 = add(d2,add(rvCyclo[28,8],rvCyclo[28,-8]));
let channels = add(rvCyclo[28,-24],
    add(neg(rvCycloMul(d2,rvCyclo[28,-12])),rvCycloMul(d4,rvCyclo[28,12])));
let weighted = rvCycloMul(rvCyclo[28,-24],channels);
assert rvCycloQMul(rvCycloRat[1](d2),t) == rvCycloRat[1](weighted);
show t;

// Six elementary strands give the color-two figure-eight. Check both
// Markov stabilizations by cabling on eight elementary strands.
let f = rvTLColoredJones[5,3,2,[1,-2,1,-2]];
assert f == rvTLColoredJones[5,4,2,[1,-2,1,-2,3]];
assert f == rvTLColoredJones[5,4,2,[1,-2,1,-2,-3]];
show f;
