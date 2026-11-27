// Three close pairs give an even split of depth 2, 3, and 5.
let points = rvStableLine[7]{[0,1],[49,1],[1,1],[344,1],[2,1],[16809,1]};
let theta = rvHyperCover(points);
rvHyperGraph(theta);
rvHyperHarmonic(theta);
rvHyperGenus(theta);
rvHyperCycles(theta);
rvHyperMonodromy(theta);
rvHyperComponent(theta);
rvHyperComponent[2](theta);

// Move one mark to infinity by an integral projective change of coordinates.
let inverted = rvStableLine[7]{[1,0],[1,49],[1,1],[1,344],[1,2],[1,16809]};
assert rvHyperCover(inverted) == theta;

// A ramified bridge joins two genus-one components. The even split
// produces one cycle; its pairing doubles with the ramification multiplier.
let mixed = rvHyperCover(rvStableLine[7]{[0,1],[7,1],[14,1],[1,1],[50,1],[2,1],[3,1],[4,1]});
rvHyperGraph(mixed);
rvHyperHarmonic(mixed);
rvHyperGenus(mixed);
rvHyperCycles(mixed);
rvHyperMetric[2](mixed);
rvHyperMonodromy[2](mixed);
rvHyperComponent[2](mixed);
