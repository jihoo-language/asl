rvhodgewick B[20]{0,1}{[0,1],[1,0]};
let tu=add(
    add(rvVerbitskyHHTerm[B,1,2,-2,0]{1,1},rvVerbitskyHHTerm[B,1,2,-2,2]{1,3}),
    rvVerbitskyHHTerm[B,1,2,-2,2]{3,1}
);
let tv=add(
    add(rvVerbitskyHHTerm[B,1,2,-2,0]{2,2},rvVerbitskyHHTerm[B,1,2,-2,1]{2,3}),
    rvVerbitskyHHTerm[B,1,2,-2,1]{3,2}
);
let p=rvVerbitskyHHTerm[B,1,2,0,3]{2,1};

rvverbitskyassoc Mix[3]{tu,tv,p};
rvVerbitskyHHDims[B,1,3,-2];
rvVerbitskyAssocReport[Mix];
rvVerbitskyAssocDirection[Mix]{1,0,1};
rvVerbitskyAssocDirection[Mix]{0,1,1};
rvVerbitskyAssocDirection[Mix]{1,1,0};
assert rvVerbitskyAssocObstruction[Mix]{1,0,1}==rvVerbitskyHHClass(rvVerbitskyHHBracket(p,tu));
assert rvVerbitskyAssocObstruction[Mix]{0,1,1}==rvVerbitskyHHClass(rvVerbitskyHHBracket(p,tv));
rvVerbitskyHHCoeff[2]{2,1,1}(rvVerbitskyAssocResidual[Mix]{1,0,1});
rvVerbitskyHHCoeff[1]{2,2,1}(rvVerbitskyAssocResidual[Mix]{0,1,1});
rvVerbitskyAssocProduct[Mix,1,1];
rvVerbitskyAssocProduct[Mix,2,1];
rvVerbitskyAssocProduct[Mix,1,2];

rvverbitskyassoc Carry[4]{tu,tv};
rvVerbitskyAssocReport[Carry];
rvVerbitskyAssocProduct[Carry,3,3];
assert rvVerbitskyAssocCoefficient[Carry]{1,1}==rvVerbitskyHHTerm[B,1,2,-4,0]{3,3};
rvverbitskyassoc Quantum[4]{p};
rvVerbitskyAssocReport[Quantum];
rvVerbitskyAssocProduct[Quantum,2,1];
