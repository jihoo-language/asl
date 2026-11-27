rvunram F[7,1]{0,1};
let z = rvUnramElem[F]{0};
let o = rvUnramElem[F]{1};
let two = rvUnramElem[F]{2};
let three = rvUnramElem[F]{3};
let four = rvUnramElem[F]{4};
let six = rvUnramElem[F]{6};

// (x²+1)(x²+2)(x²+3)
rvhyperelliptic C[F]{six,z,four,z,six,z,o};
rvbielliptic C[Plus,Minus](o,z,z,o);

let p = rvBiellipticPoint[Plus](three,two);
let q = rvBiellipticPoint[Minus](three,two);
assert p == rvJacPoint[Plus](two,two);
assert q == rvJacPoint[Minus](three,two);
assert rvBiellipticPoint[Plus](neg(three),two) == p;
assert rvBiellipticPoint[Minus](neg(three),neg(two)) == q;
show p;
show q;

rvHWWWeil[C];
rvHWWWeil[Plus];
rvHWWWeil[Minus];
