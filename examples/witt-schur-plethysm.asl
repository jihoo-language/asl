// Stable Schur identities in degree four.
let s2 = rvSchurCharacter[[2]];
let s11 = rvSchurCharacter[[1,1]];
let square = rvSchurPlethysm(s2,s2);
assert square == add(rvSchurCharacter[[4]],rvSchurCharacter[[2,2]]);
show square;
rvSchurCharacterCoeff[[2,2]](square);

// The virtual inner alphabet participates in every Adams power.
let g = add(s2,neg(s11));
let f = rvSchurPlethysm(s2,g);
show f;
let x = rvWittSpectrum{[1,-3,2]};
let composed = rvSchurCharacterEvaluate(f,x);
assert composed == rvSchurCharacterEvaluate(s2,rvSchurCharacterEvaluate(g,x));
rvWittSpectrumGhost[1](composed);
