// The signed flip splits the tensor square into its two rational factors.
let x = add(rvWittSpectrumLine[2],neg(rvWittSpectrumLine[1]));
let invariant = rvWittSpectrumCyclic[2,1](x);
let alternating = rvWittSpectrumCyclic[2,2](x);
assert add(invariant,alternating) == rvWittSpectrumTensor(x,x);
assert alternating == rvWittSpectrumLie[2](x);
rvWittSpectrumGhost[1](invariant);
rvWittSpectrumGhost[1](alternating);

// The primitive cyclic part for n=6 contains two copies of Lie_6.
let binary = rvWittSpectrum{[1,-2,1]};
let lie = rvWittSpectrumLie[6](binary);
assert rvWittSpectrumCyclic[6,6](binary) == add(lie,lie);
rvWittSpectrumCyclicGhost[6,6,1](binary);
rvWittSpectrumLieGhost[6,1](binary);
rvWittSpectrumLieGhost[14,1](binary);
