// The determinant for the combined action is lambda_1^3 + lambda_0^2 lambda_2.
qmodule Cusp[0, 1, 2] {
    b1:0, b2:2, b3:4, a1:3, a2:5, a3:7
} {
    1(a1)=>[b1], 1(a2)=>[b2], 1(a3)=>[b3],
    0(a1)=>[b2], 0(a2)=>[b3], 2(a3)=>[b1]
};

// With lambda_2=1 the curve is (lambda_0,lambda_1)=(r^3,r^2).
// Each pair below has the same nonzero coordinate positions.
rvdim[One, Cusp, [0,1,2], 7, [1,3,1]];
rvdim[One, Cusp, [0,1,2], 7, [2,3,1]];
rvspectrum[One, Cusp, [0,1,2], 7];

rvdim[One, Cusp, [0,1,2], 11, [3,4,1]];
rvdim[One, Cusp, [0,1,2], 11, [2,4,1]];
rvspectrum[One, Cusp, [0,1,2], 11];

// The affine cusp and the point at infinity have different ranks.
rvdim[One, Cusp, [0,1,2], 7, [0,0,1]];
rvdim[One, Cusp, [0,1,2], 7, [1,0,0]];
