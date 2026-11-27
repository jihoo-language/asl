rvcmclasspoly H3[-3];
rvcmclasspoly H4[-4];
rvcmclasspoly H7[-7];
rvcmclasspoly H8[-8];
rvcmclasspolylift H12[H3,2];
rvcmclasspolylift H16[H4,2];
rvcmclasspolylift H28[H7,2];
rvcmclasspolylift H48[H12,2];
rvcmclasspolylift H63[H7,3];

rvfrickecmnorm N3[H3];
rvfrickecmnorm N4[H4];
rvfrickecmnorm N7[H7];
rvfrickecmnorm N12[H12];
rvFrickeCMPolyData[N3];
rvFrickeCMPolyData[N4];
rvFrickeCMPolyData[N7];
rvFrickeCMPolyData[N12];

rvfrickecmedge E312[H3,H12];
rvfrickecmedge E416[H4,H16];
rvfrickecmedge E728[H7,H28];
rvfrickecmedge E1248[H12,H48];
rvFrickeCMPolyData[E312];
rvFrickeCMPolyData[E416];
rvFrickeCMPolyData[E728];
rvFrickeCMPolyData[E1248];
rvFrickeCMPolyCoeff[E1248,0];

rvfrickecmedge L4[H4,H4];
rvfrickecmedge L7[H7,H7];
rvfrickecmedge L8[H8,H8];
rvfrickecmedge L63[H63,H63];
rvFrickeCMPolyData[L4];
rvFrickeCMPolyData[L7];
rvFrickeCMPolyData[L8];
rvFrickeCMPolyData[L63];

rvunram U[5,3]{0,1};
let u=rvUnramElem[U]{1};
assert rvFrickeCMPolyEval[N3](u)==rvUnramElem[U]{124};
assert rvFrickeCMPolyModCoeff[E312,U,0]==rvUnramElem[U]{123};
assert rvFrickeCMPairCoeff[0](u)==rvUnramElem[U]{124};

rvunram F[13,1]{0,1};
rvFrickeCMEdgeEndpoints[E1248](rvUnramElem[F]{2});
rvunram G[11,1]{0,1};
rvFrickeCMEdgeEndpoints[E1248](rvUnramElem[G]{6});
rvunram K[7,1]{0,1};
rvFrickeCMEndpoints(rvUnramElem[K]{6});
rvFrickeCMEdgeEndpoints[L7](rvUnramElem[K]{5});
rvunram V[17,1]{0,1};
rvFrickeCMPolyFactors[E1248,V];

rvfrickecmedge Empty[H3,H7];
rvFrickeCMPolyData[Empty];
rvunram T[5,1]{0,1};
let collision=rvUnramElem[T]{2};
assert rvFrickeCMPolyEval[N3](collision)==zero[RVUnram(T)];
assert rvFrickeCMPolyEval[N7](collision)==zero[RVUnram(T)];
assert rvFrickeCMPolyEval[Empty](collision)==rvUnramElem[T]{1};
rvFrickeCMPolyFactors[N7,T];
