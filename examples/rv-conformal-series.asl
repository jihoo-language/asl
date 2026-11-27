// The handle gives a rational vacuum series in the genus variable t.
rvConformalSeries[3,2,[]];
rvConformalSeries[3,3,[]];
rvConformalMod[3,2,1000000000000000000,1000000007,[]];
rvConformalMod[3,3,1000000000000000000,1000000007,[]];
// A fundamental paired with its dual cancels a pole at level three.
rvConformalNum[3,3,[[0,0,1],[0,1,1]]];
rvConformalDen[3,3,[[0,0,1],[0,1,1]]];
// A single fundamental has nonzero central charge and gives the zero series.
rvConformalSeries[3,3,[[0,0,1]]];
rvConformalMod[3,3,1000000000000000000,1000000000,[[0,0,1]]];
// The positive handle operator gives a dimension beyond signed 128 bits.
rvConformalExact[3,2,100,[]];
