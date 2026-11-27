// SU(2) level two. WRT is normalized to give S^3 the value one.
let one = rvCycloRat[1](rvCyclo[16,0]);
assert rvTLSurgery[2,0,[],[]] == one;
assert rvTLSurgery[2,1,[],[1]] == one;
assert rvTLSurgery[2,1,[],[-1]] == one;
assert rvTLSurgery[2,2,[1,1],[0,0]] == one;

// Slide the first component of the (2,3)-framed unlink over the second.
// The new framings are (5,3) and their linking number is three.
let unlink = rvTLSurgery[2,2,[],[2,3]];
let slide = rvTLSurgery[2,2,[1,1,1,1,1,1],[5,3]];
assert slide == unlink;
let expected = add(rvCyclo[16,1],neg(rvCyclo[16,3]));
assert slide == rvCycloRat[1](expected);
show slide;
assert rvTLSurgery[2,2,[-1,-1,-1,-1,-1,-1],[-5,-3]] == rvCycloQGalois[-1](slide);

// +1 surgery on the trefoil differs from +1 surgery on the unknot.
let trefoil = rvTLSurgery[2,2,[1,1,1],[1]];
assert trefoil == neg(one);
show trefoil;
assert rvTLSurgery[2,3,[1,1,1],[1,-1]] == trefoil;

// The doubled invariant uses the same normalization as rvLensTV.
// The Galois exponent nine sends A to -A in the level-two root field.
let doubled = rvTLSurgeryTV[2,1,[],[2]];
assert doubled == rvCycloQGalois[9](rvLensTV[2,2,2,1]);
show doubled;
assert rvTLSurgeryTV[2,1,[],[0]] == one;
assert rvTLSurgeryTV[2,2,[1,1],[1,1]] == one;
