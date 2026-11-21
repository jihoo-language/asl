let xi1 = cword{[1]};
let xi2 = cword{[0, 1]};
let xi3 = cword{[0, 0, 1]};
let one = cword{};

assert cobarD(xi1) == zero[Cobar(2, 1)];
assert cobarD(xi2) == cword{[2], [1]};
assert cobarD(xi3) == add(cword{[0, 2], [1]}, cword{[4], [0, 1]});
assert cobarD(cobarD(xi3)) == zero[Cobar(3, 7)];
assert concat(one, xi3) == xi3;
assert cobarD(concat(xi2, xi3)) ==
    add(concat(cobarD(xi2), xi3), concat(xi2, cobarD(xi3)));

let boundary = \x : Cobar(1, 3) => cobarD(x);
show boundary(xi2);
show cobarD(xi3);
betti Cobar(1, 3);
