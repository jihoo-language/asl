// The major-index residues distinguish a primitive eigenvalue from its full
// rational cyclotomic block.
rvWittSchurMajorIndex[[3,2,1]];
rvWittCyclicSchurMultiplicity[[3,2,1],6];
rvWittLieSchurMultiplicity[[3,2,1]];

// The degree-six free Lie algebra on an even two-dimensional module.
let x = rvWittSpectrum{[1,-2,1]};
let a = rvWittSpectrumSchur[[5,1]](x);
let b = rvWittSpectrumSchur[[4,2]](x);
let c = rvWittSpectrumSchur[[3,3]](x);
let lie = add(a,add(b,c));
assert lie == rvWittSpectrumLie[6](x);
assert add(lie,lie) == rvWittSpectrumCyclic[6,6](x);
rvWittSpectrumSchurGhost[[4,2],1](x);
rvWittSpectrumGhost[1](lie);
