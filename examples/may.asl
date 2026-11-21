// E_0 uses the component of the cobar differential that preserves filtration.
let xi2 = cword{[0, 1]};
let e0 = mayClass[0, 0](xi2);
assert mayD(e0) == zero[May(0, 2, 3, 0)];

// The next component gives d_1(h_{2,0}) = h_1 h_0.
let e1 = mayClass[1, 0](xi2);
let h0e1 = mayClass[1, 0](cword{[1]});
let h1e1 = mayClass[1, 1](cword{[2]});
assert mayD(e1) == mayMul(h1e1, h0e1);
assert mayD(e1) == maygen[1, 2, 3, 1, 0];
show mayD(e1);

// Correct the square [xi_2 | xi_2] in higher filtration to lift b_{2,0} to E_3.
let b20lift = add(add(cword{[0, 1], [0, 1]}, cword{[2, 1], [1]}),
                  add(cword{[2], [1, 1]}, cword{[5], [1]}));
assert cobarD(b20lift) == add(cword{[1], [4], [1]}, cword{[2], [2], [2]});
let b20 = mayClass[3, 0](b20lift);
assert b20 == maygen[3, 2, 6, 0, 0];
assert mayClass[3, 0](mayRep(b20)) == b20;

// Its d_3 identifies h_0^2 h_2 with h_1^3 on E_4.
let h0 = mayClass[3, 0](cword{[1]});
let h1 = mayClass[3, 1](cword{[2]});
let h2 = mayClass[3, 3](cword{[4]});
let h0h2h0 = mayMul(mayMul(h0, h2), h0);
let h1cubed = mayMul(mayMul(h1, h1), h1);
assert mayD(b20) == add(h0h2h0, h1cubed);
show mayD(b20);
maydim[3, 2, 6, 0];
maydim[4, 2, 6, 0];
maydim[3, 3, 6, 3];
maydim[4, 3, 6, 3];
