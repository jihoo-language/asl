rvverbitskytate T{5};
rvlgdefect I[T]{0,1,2};
rvlgdefect J[T]{1,2,3,4};
rvlgdefectconvolution F[I,J];
rvLGConvolution[F];
rvLGConvolutionSplit[F];

let one=rvCycloRat[1](rvCyclo[5,0]);
assert rvLGConvolutionTerm[F,0,4,0,0,0]==one;
assert rvLGConvolutionTerm[F,1,4,0,0,0]==neg(one);
assert rvLGConvolutionTerm[F,2,4,0,3,0]==one;
assert rvLGConvolutionTerm[F,3,0,0,0,0]==one;
assert rvLGConvolutionTerm[F,4,0,0,0,0]==one;
assert rvLGConvolutionTerm[F,5,0,4,0,5]==one;
rvlgconvolutionpart K[F,2];
rvLGDefect[K];
assert rvLGDefectCoeff[K,0,1]==neg(add(rvCyclo[5,3],rvCyclo[5,4]));

rvverbitskytate S{4};
rvlgdefect E[S]{0,2};
rvlgdefect O[S]{1,3};
rvlgdefectconvolution Q[E,O];
rvLGConvolution[Q];

let unit=rvCycloRat[1](rvCyclo[4,0]);
assert rvLGConvolutionTerm[Q,0,0,2,0,2]==neg(unit);
assert rvLGConvolutionTerm[Q,1,0,2,0,2]==unit;
assert rvLGConvolutionTerm[Q,2,0,0,1,2]==neg(unit);
assert rvLGConvolutionTerm[Q,2,0,2,2,2]==neg(unit);
