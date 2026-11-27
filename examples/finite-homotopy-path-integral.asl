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

// No background is supplied. The normalized group presentation
// has one generator satisfying its square relation.
rvhomotopy Nonsplit[RP3,P](k);
rvHomotopyCardinality[Nonsplit];
rvHomotopyComponent[Nonsplit,0];
rvHomotopyClass[Nonsplit,0,0];
rvHomotopyClass[Nonsplit,0,1];

let split=zero[RVPostnikovCochain(P,3)];
rvhomotopy Split[RP3,P](split);
rvHomotopyCardinality[Split];

space Circle {a,b,c,d,a<c,a<d,b<c,b<d};
rvhomotopy Loops[Circle,P](k);
rvHomotopyCardinality[Loops];
rvHomotopyComponent[Loops,0];
