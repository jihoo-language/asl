// The monodromy chain makes the high plane weakly admissible.
rvphodge Good[3] {
    {(1,1),(3,1),(9,1)},
    {(0,1,1,1),(1,2,1,1)},
    {(-1,[1,0,0]),(2,[0,1,0]),(2,[0,0,1])}
};
rvPHodgeWeak[Good];
rvPHodgeData[Good,[0,1]];
rvPHodgeHN[Good];

rvphodge NoN[3] {
    {(1,1),(3,1),(9,1)},
    {},
    {(-1,[1,0,0]),(2,[0,1,0]),(2,[0,0,1])}
};
rvPHodgeWeak[NoN];

rvphodge Moved[3] {
    {(1,1),(3,1),(9,1)},
    {(0,1,1,1),(1,2,1,1)},
    {(2,[1,0,0]),(2,[0,1,0]),(-1,[0,0,1])}
};
rvPHodgeWeak[Moved];
rvPHodgeHN[Moved];

// This noncoordinate flag has a fractional first slope.
rvphodge Fractional[3] {
    {(1,1),(3,1),(9,1)},
    {(0,1,1,1),(1,2,1,1)},
    {(0,[1,0,0]),(1,[0,1,1]),(2,[1,1,0])}
};
rvPHodgeData[Fractional,[0,1]];
rvPHodgeHN[Fractional];
