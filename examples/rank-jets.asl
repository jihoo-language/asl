qmodule Cusp[0, 1, 2] {
    b1:0, b2:2, b3:4, a1:3, a2:5, a3:7
} {
    1(a1)=>[b1], 1(a2)=>[b2], 1(a3)=>[b3],
    0(a1)=>[b2], 0(a2)=>[b3], 2(a3)=>[b1]
};

// All three jets have residue (0,0,1) over F4.
rvdim[One, Cusp, [0,1,2], 7, [0,0,1]];
rvjet[One, Cusp, [0,1,2], 7, 7, [[0,1],[],[1]]];
rvjet[One, Cusp, [0,1,2], 7, 7, [[],[0,1],[1]]];
rvjet[One, Cusp, [0,1,2], 7, 7, [[0,0,0,1],[0,0,1],[1]]];

// Use the commuting primitive actions at both ends.
rvjet[Cusp, Cusp, [0,1,2], 7, 4, [[],[0,1],[1]]];
rvjet[Cusp, Cusp, [0,1,2], 7, 4, [[0,0,0,1],[0,0,1],[1]]];
