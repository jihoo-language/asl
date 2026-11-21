let h0 = extgen[1, 1, 0];
let h1 = extgen[1, 2, 0];
let h2 = extgen[1, 4, 0];
let h3 = extgen[1, 8, 0];
let x = yoneda(h0, h2);

assert extSq[0](h0) == h1;
assert extSq[0](h2) == h3;
assert extSq[2](x) == yoneda(x, x);
let middle = extSq[1](x);
assert middle == yoneda(yoneda(h0, h0), h3);
assert middle == extgen[3, 10, 0];
show middle;
