let h0 = extgen[1, 1, 0];
let h1 = extgen[1, 2, 0];
let h2 = extgen[1, 4, 0];
let h0square = yoneda(h0, h0);
let h0fourth = yoneda(h0square, h0square);
let h1cube = yoneda(yoneda(h1, h1), h1);

assert resHigher[1](h0) == resA[1](h0);
assert resA[1](h2) == zero[ExtAn(1, 1, 4)];
assert resA[1](h1cube) == zero[ExtAn(1, 3, 6)];

let r2 = resHigher[1](h0, h1cube);
let v1fourth = resHigher[1](h2, h0square, h1cube);
assert r2 == extangen[1, 3, 7, 0];
assert v1fourth == extangen[1, 4, 12, 0];
assert resHigher[1](h0fourth, h1, h2, h1) == v1fourth;
assert resHigher[1](resA[2](h2), resA[2](h0square), resA[2](h1cube)) == v1fourth;
show r2;
show v1fourth;
