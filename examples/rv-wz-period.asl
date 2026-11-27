// The Apéry numbers for zeta(2) satisfy the indicated Ore equation.
let apery = add(rvOre[0,2],add(rvOre[1,2,-11],add(rvOre[1,1,-11],
    add(rvOre[1,0,-3],add(rvOre[2,2,-1],add(rvOre[2,1,-2],rvOre[2,0,-1]))))));
rvWZCheck[2,1,1,2](apery) {
    {(0,5,1),(1,4,6),(0,4,-5),(2,3,-11),(1,3,7)},
    {(1,0,0,1),(1,1,0,1),(1,1,-1,1)}
};

// Independent recurrence for the Franel cube sums.
let franel = add(rvOre[0,2],add(rvOre[1,2,-7],add(rvOre[1,1,-7],
    add(rvOre[1,0,-2],add(rvOre[2,2,-8],add(rvOre[2,1,-16],rvOre[2,0,-8]))))));
rvWZCheck[3,0,1,2](franel) {
    {(3,3,-14),(2,4,27),(2,3,10),(1,5,-18),(1,4,-15),(0,6,4),(0,5,6)},
    {(1,0,0,3),(1,0,-1,1)}
};

// The signed central binomial sum vanishes in odd degrees.
let alternating = add(rvOre[0,1],add(rvOre[2,1,4],rvOre[2,0,4]));
rvWZCheck[2,0,-1,2](alternating) {
    {(0,4,-2),(1,3,6),(0,3,-2),(2,2,-5),(1,2,3)},
    {(1,0,0,2),(1,0,-1,1)}
};
