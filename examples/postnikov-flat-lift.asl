// Use the face poset of the antipodal quotient of the boundary of the
// 4-cross-polytope. In a signed face, p and m denote plus and minus.
// Fix the first sign of every representative to plus.
space RP3 {
  ap, bp, cp, dp, apbp, apbm, apcp, apcm,
  apdp, apdm, bpcp, bpcm, bpdp, bpdm, cpdp, cpdm,
  apbpcp, apbpcm, apbmcp, apbmcm, apbpdp, apbpdm, apbmdp, apbmdm,
  apcpdp, apcpdm, apcmdp, apcmdm, bpcpdp, bpcpdm, bpcmdp, bpcmdm,
  apbpcpdp, apbpcpdm, apbpcmdp, apbpcmdm, apbmcpdp, apbmcpdm, apbmcmdp, apbmcmdm,
  ap<apbp, ap<apbm, ap<apcp, ap<apcm, ap<apdp, ap<apdm,
  bp<apbp, bp<apbm, bp<bpcp, bp<bpcm, bp<bpdp, bp<bpdm,
  cp<apcp, cp<apcm, cp<bpcp, cp<bpcm, cp<cpdp, cp<cpdm,
  dp<apdp, dp<apdm, dp<bpdp, dp<bpdm, dp<cpdp, dp<cpdm,
  apbp<apbpcp, apbp<apbpcm, apbp<apbpdp, apbp<apbpdm, apbm<apbmcp, apbm<apbmcm,
  apbm<apbmdp, apbm<apbmdm, apcp<apbpcp, apcp<apbmcp, apcp<apcpdp, apcp<apcpdm,
  apcm<apbpcm, apcm<apbmcm, apcm<apcmdp, apcm<apcmdm, apdp<apbpdp, apdp<apbmdp,
  apdp<apcpdp, apdp<apcmdp, apdm<apbpdm, apdm<apbmdm, apdm<apcpdm, apdm<apcmdm,
  bpcp<apbpcp, bpcp<apbmcm, bpcp<bpcpdp, bpcp<bpcpdm, bpcm<apbpcm, bpcm<apbmcp,
  bpcm<bpcmdp, bpcm<bpcmdm, bpdp<apbpdp, bpdp<apbmdm, bpdp<bpcpdp, bpdp<bpcmdp,
  bpdm<apbpdm, bpdm<apbmdp, bpdm<bpcpdm, bpdm<bpcmdm, cpdp<apcpdp, cpdp<apcmdm,
  cpdp<bpcpdp, cpdp<bpcmdm, cpdm<apcpdm, cpdm<apcmdp, cpdm<bpcpdm, cpdm<bpcmdp,
  apbpcp<apbpcpdp, apbpcp<apbpcpdm, apbpcm<apbpcmdp, apbpcm<apbpcmdm, apbmcp<apbmcpdp, apbmcp<apbmcpdm,
  apbmcm<apbmcmdp, apbmcm<apbmcmdm, apbpdp<apbpcpdp, apbpdp<apbpcmdp, apbpdm<apbpcpdm, apbpdm<apbpcmdm,
  apbmdp<apbmcpdp, apbmdp<apbmcmdp, apbmdm<apbmcpdm, apbmdm<apbmcmdm, apcpdp<apbpcpdp, apcpdp<apbmcpdp,
  apcpdm<apbpcpdm, apcpdm<apbmcpdm, apcmdp<apbpcmdp, apcmdp<apbmcmdp, apcmdm<apbpcmdm, apcmdm<apbmcmdm,
  bpcpdp<apbpcpdp, bpcpdp<apbmcmdm, bpcpdm<apbpcpdm, bpcpdm<apbmcmdp, bpcmdp<apbpcmdp, bpcmdp<apbmcpdm,
  bpcmdm<apbpcmdm, bpcmdm<apbmcpdp,
};

rvgroup H[4];
rvgroup E[4];
rvcrossed C[H,E]{0,2,0,2}{{0,1,2,3},{0,3,2,1},{0,1,2,3},{0,3,2,1}};
rvpostnikov P[C,2];
let k = rvPostnikovK[P];

// An edge label records whether taking the boundary face reverses all signs.
rvflat Cover[RP3,P]{
  ap<apbp:0, ap<apbm:0, ap<apcp:0, ap<apcm:0, ap<apdp:0, ap<apdm:0,
  bp<apbp:0, bp<apbm:1, bp<bpcp:0, bp<bpcm:0, bp<bpdp:0, bp<bpdm:0,
  cp<apcp:0, cp<apcm:1, cp<bpcp:0, cp<bpcm:1, cp<cpdp:0, cp<cpdm:0,
  dp<apdp:0, dp<apdm:1, dp<bpdp:0, dp<bpdm:1, dp<cpdp:0, dp<cpdm:1,
  apbp<apbpcp:0, apbp<apbpcm:0, apbp<apbpdp:0, apbp<apbpdm:0, apbm<apbmcp:0, apbm<apbmcm:0,
  apbm<apbmdp:0, apbm<apbmdm:0, apcp<apbpcp:0, apcp<apbmcp:0, apcp<apcpdp:0, apcp<apcpdm:0,
  apcm<apbpcm:0, apcm<apbmcm:0, apcm<apcmdp:0, apcm<apcmdm:0, apdp<apbpdp:0, apdp<apbmdp:0,
  apdp<apcpdp:0, apdp<apcmdp:0, apdm<apbpdm:0, apdm<apbmdm:0, apdm<apcpdm:0, apdm<apcmdm:0,
  bpcp<apbpcp:0, bpcp<apbmcm:1, bpcp<bpcpdp:0, bpcp<bpcpdm:0, bpcm<apbpcm:0, bpcm<apbmcp:1,
  bpcm<bpcmdp:0, bpcm<bpcmdm:0, bpdp<apbpdp:0, bpdp<apbmdm:1, bpdp<bpcpdp:0, bpdp<bpcmdp:0,
  bpdm<apbpdm:0, bpdm<apbmdp:1, bpdm<bpcpdm:0, bpdm<bpcmdm:0, cpdp<apcpdp:0, cpdp<apcmdm:1,
  cpdp<bpcpdp:0, cpdp<bpcmdm:1, cpdm<apcpdm:0, cpdm<apcmdp:1, cpdm<bpcpdm:0, cpdm<bpcmdp:1,
  apbpcp<apbpcpdp:0, apbpcp<apbpcpdm:0, apbpcm<apbpcmdp:0, apbpcm<apbpcmdm:0, apbmcp<apbmcpdp:0, apbmcp<apbmcpdm:0,
  apbmcm<apbmcmdp:0, apbmcm<apbmcmdm:0, apbpdp<apbpcpdp:0, apbpdp<apbpcmdp:0, apbpdm<apbpcpdm:0, apbpdm<apbpcmdm:0,
  apbmdp<apbmcpdp:0, apbmdp<apbmcmdp:0, apbmdm<apbmcpdm:0, apbmdm<apbmcmdm:0, apcpdp<apbpcpdp:0, apcpdp<apbmcpdp:0,
  apcpdm<apbpcpdm:0, apcpdm<apbmcpdm:0, apcmdp<apbpcmdp:0, apcmdp<apbmcmdp:0, apcmdm<apbpcmdm:0, apcmdm<apbmcmdm:0,
  bpcpdp<apbpcpdp:0, bpcpdp<apbmcmdm:1, bpcpdm<apbpcpdm:0, bpcpdm<apbmcmdp:1, bpcmdp<apbpcmdp:0, bpcmdp<apbmcpdm:1,
  bpcmdm<apbpcmdm:0, bpcmdm<apbmcpdp:1,
};
rvFlatData[Cover];
rvFlatObstruction[Cover](k);
let pulled = rvFlatPullback[Cover](k);
rvFlatCoeff[dp,cpdm,bpcmdp,apbmcpdm](pulled);
assert rvFlatCoboundary(pulled)==zero[RVFlatCochain(Cover,4)];

// Sum all 192 tetrahedra. Each triangle occurs twice in its boundary.
rvFlatCyclePair(pulled){
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
};

let split = zero[RVPostnikovCochain(P,3)];
rvFlatObstruction[Cover](split);
let b = rvFlatLift[Cover](split);
let b2 = add(b,rvFlatH2Basis[Cover,0]);
rvFlatLiftCheck[Cover](split,b2);
assert rvFlatLiftDifference[Cover](split,b,b2)==rvFlatH2Basis[Cover,0];

rvflat Trivial[RP3,P]{
  ap<apbp:0, ap<apbm:0, ap<apcp:0, ap<apcm:0, ap<apdp:0, ap<apdm:0,
  bp<apbp:0, bp<apbm:0, bp<bpcp:0, bp<bpcm:0, bp<bpdp:0, bp<bpdm:0,
  cp<apcp:0, cp<apcm:0, cp<bpcp:0, cp<bpcm:0, cp<cpdp:0, cp<cpdm:0,
  dp<apdp:0, dp<apdm:0, dp<bpdp:0, dp<bpdm:0, dp<cpdp:0, dp<cpdm:0,
  apbp<apbpcp:0, apbp<apbpcm:0, apbp<apbpdp:0, apbp<apbpdm:0, apbm<apbmcp:0, apbm<apbmcm:0,
  apbm<apbmdp:0, apbm<apbmdm:0, apcp<apbpcp:0, apcp<apbmcp:0, apcp<apcpdp:0, apcp<apcpdm:0,
  apcm<apbpcm:0, apcm<apbmcm:0, apcm<apcmdp:0, apcm<apcmdm:0, apdp<apbpdp:0, apdp<apbmdp:0,
  apdp<apcpdp:0, apdp<apcmdp:0, apdm<apbpdm:0, apdm<apbmdm:0, apdm<apcpdm:0, apdm<apcmdm:0,
  bpcp<apbpcp:0, bpcp<apbmcm:0, bpcp<bpcpdp:0, bpcp<bpcpdm:0, bpcm<apbpcm:0, bpcm<apbmcp:0,
  bpcm<bpcmdp:0, bpcm<bpcmdm:0, bpdp<apbpdp:0, bpdp<apbmdm:0, bpdp<bpcpdp:0, bpdp<bpcmdp:0,
  bpdm<apbpdm:0, bpdm<apbmdp:0, bpdm<bpcpdm:0, bpdm<bpcmdm:0, cpdp<apcpdp:0, cpdp<apcmdm:0,
  cpdp<bpcpdp:0, cpdp<bpcmdm:0, cpdm<apcpdm:0, cpdm<apcmdp:0, cpdm<bpcpdm:0, cpdm<bpcmdp:0,
  apbpcp<apbpcpdp:0, apbpcp<apbpcpdm:0, apbpcm<apbpcmdp:0, apbpcm<apbpcmdm:0, apbmcp<apbmcpdp:0, apbmcp<apbmcpdm:0,
  apbmcm<apbmcmdp:0, apbmcm<apbmcmdm:0, apbpdp<apbpcpdp:0, apbpdp<apbpcmdp:0, apbpdm<apbpcpdm:0, apbpdm<apbpcmdm:0,
  apbmdp<apbmcpdp:0, apbmdp<apbmcmdp:0, apbmdm<apbmcpdm:0, apbmdm<apbmcmdm:0, apcpdp<apbpcpdp:0, apcpdp<apbmcpdp:0,
  apcpdm<apbpcpdm:0, apcpdm<apbmcpdm:0, apcmdp<apbpcmdp:0, apcmdp<apbmcmdp:0, apcmdm<apbpcmdm:0, apcmdm<apbmcmdm:0,
  bpcpdp<apbpcpdp:0, bpcpdp<apbmcmdm:0, bpcpdm<apbpcpdm:0, bpcpdm<apbmcmdp:0, bpcmdp<apbpcmdp:0, bpcmdp<apbmcpdm:0,
  bpcmdm<apbpcmdm:0, bpcmdm<apbmcpdp:0,
};
rvFlatObstruction[Trivial](k);
rvFlatLiftCheck[Trivial](k,rvFlatLift[Trivial](k));
