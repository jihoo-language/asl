let h0 = extgen[1, 1, 0];
let h1 = extgen[1, 2, 0];
let h2 = extgen[1, 4, 0];

// The transferred ternary product of h0, h1, h0 gives h1 squared.
let triple = ainfty(h0, h1, h0);
assert triple == yoneda(h1, h1);
show triple;

// The product h2 h0 obstructs a Massey defining system for these inputs.
// The transferred operation still uses the chosen cobar contraction.
assert yoneda(h2, h0) == extgen[2, 5, 0];
let fourth = ainfty(h0, h1, h2, h0);
assert fourth == yoneda(h2, h2);
show fourth;
show alift(h0, h1, h2, h0);
