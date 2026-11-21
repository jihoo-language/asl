// The classical mod-2 Adams E2 page for the sphere.
let one = extClass(cword{});
let h0 = extClass(cword{[1]});
let h1 = extClass(cword{[2]});
let h2 = extClass(cword{[4]});

assert yoneda(one, h1) == h1;
let h0Square = yoneda(h0, h0);
let h1Square = yoneda(h1, h1);
assert h0Square == extgen[2, 2, 0];
assert h1Square == extgen[2, 4, 0];

// Both concatenations of h0 and h1 are cobar boundaries.
assert yoneda(h0, h1) == zero[Ext(2, 3)];
assert yoneda(h1, h0) == zero[Ext(2, 3)];
assert cobarD(extPrimitive(cword{[1], [2]})) == cword{[1], [2]};

let h1Cube = yoneda(h1Square, h1);
assert h1Cube == yoneda(h0Square, h2);
assert h1Cube == extgen[3, 6, 0];

extdim[2, 3];
extdim[3, 6];
show h1Cube;
