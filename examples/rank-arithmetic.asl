qmodule TwistedCusp[0, 1, 2] {
    b1a:0, b1b:0, b2a:2, b2b:2, b3a:4, b3b:4,
    a1a:3, a1b:3, a2a:5, a2b:5, a3a:7, a3b:7
} {
    0(a1a)=>[b2a], 0(a1b)=>[b2b], 0(a2a)=>[b3a], 0(a2b)=>[b3b],
    1(a1a)=>[b1a], 1(a1b)=>[b1b], 1(a2a)=>[b2a], 1(a2b)=>[b2b],
    1(a3a)=>[b3a], 1(a3b)=>[b3b],
    2(a3a)=>[b1b], 2(a3b)=>[b1a, b1b]
};

// The two conjugate components acquire nonzero-torus points over F4.
rvorbit[7, [1, 1, 2]];
rvdim[One, TwistedCusp, [0, 1, 2], 7, [1, 1, 2]];
rvclosed[One, TwistedCusp, [0, 1, 2], 1];
rvclosed[One, TwistedCusp, [0, 1, 2], 6];
rvzeta[One, TwistedCusp, [0, 1, 2], 6];
