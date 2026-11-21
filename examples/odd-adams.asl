// Tau_0 has stem zero and supports a polynomial tower in odd-prime Ext.
let tau0 = pword[3]{([], [0])};
let xi1 = pword[3]{([1], [])};
let a0 = extClass(tau0);
let x = extClass(xi1);
let tower = yoneda(yoneda(a0, a0), a0);
assert tower == pextgen[3, 3, 3, 0];
show tower;

// The inverse of 2 in F_3 gives a primitive for [xi_1 | xi_1].
let square = concat(xi1, xi1);
let xi1square = pword[3]{([2], [])};
assert cobarD(xi1square) == pscale[2](square);
assert extPrimitive(square) == pscale[2](xi1square);
assert yoneda(x, x) == zero[PExt(3, 2, 8)];

// The coproducts of tau_1 and xi_1 tau_0 kill the two mixed products.
assert cobarD(pword[3]{([], [1])}) == concat(xi1, tau0);
assert cobarD(pword[3]{([1], [0])}) == add(concat(xi1, tau0), pscale[2](concat(tau0, xi1)));
assert yoneda(x, a0) == zero[PExt(3, 2, 5)];
assert yoneda(a0, x) == zero[PExt(3, 2, 5)];

// The divided binomial coproduct represents b_{1,0} at p=3 and p=5.
let b3 = extClass(add(pword[3]{([1], []), ([2], [])}, pword[3]{([2], []), ([1], [])}));
assert b3 == pextgen[3, 2, 12, 0];
assert extClass(extRep(b3)) == b3;
show b3;

let b5lift = add(add(pword[5]{([1], []), ([4], [])}, pword[5]{([4], []), ([1], [])}),
                 pscale[2](add(pword[5]{([2], []), ([3], [])}, pword[5]{([3], []), ([2], [])})));
let b5 = extClass(b5lift);
assert b5 == pextgen[5, 2, 40, 0];
assert extClass(extRep(b5)) == b5;
show b5;

// Both x and y have odd stem. Their products are nonzero and anticommute.
let xi1ninth = pword[3]{([9], [])};
let y = extClass(xi1ninth);
let xy = yoneda(x, y);
let yx = yoneda(y, x);
assert xy == pextgen[3, 2, 40, 0];
assert yx == pscale[2](xy);
assert cobarD(pword[3]{([10], [])}) == add(concat(xi1, xi1ninth), concat(xi1ninth, xi1));
show xy;
show yx;
pextdim[3, 2, 40];
