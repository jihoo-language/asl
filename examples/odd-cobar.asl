let xi = pword[3]{([1], [])};
let tau0 = pword[3]{([], [0])};
let tau1 = pword[3]{([], [1])};
let square = pword[3]{([2], [])};

assert cobarD(square) == pscale[2](concat(xi, xi));
assert cobarD(pword[3]{([3], [])}) == zero[PCobar(3, 2, 12)];
assert cobarD(tau1) == concat(xi, tau0);

let signed = cobarD(pword[3]{([], [0, 1])});
assert signed == add(
    pscale[2](pword[3]{([1], [0]), ([], [0])}),
    add(pscale[2](concat(tau0, tau1)), concat(tau1, tau0))
);
assert cobarD(signed) == zero[PCobar(3, 3, 6)];
assert cobarD(concat(xi, tau1)) == pscale[2](concat(xi, concat(xi, tau0)));
show cobarD(square);
show signed;
betti PCobar(3, 1, 16);
