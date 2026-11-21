// Three Milnor primitives give three commuting Koszul variables.
let v0 = qClass(qRes[0, 1, 2](cword{[1]}));
let v1 = qClass(qRes[0, 1, 2](cword{[0, 1]}));
let v2 = qClass(qRes[0, 1, 2](cword{[0, 0, 1]}));
assert ecomp(v0, v1) == ecomp(v1, v0);
assert qRes[0, 1, 2](cword{[1], [0, 1]}) == qRes[0, 1, 2](cword{[0, 1], [1]});
let mixed = cword{[1, 1]};
assert qRes[0, 1, 2](cobarD(mixed)) == qD(qRes[0, 1, 2](mixed));
let a = ecomp(v2, ecomp(v0, v0));
let b = ecomp(v1, ecomp(v1, v1));
assert add(a, b) == add(qextgen[One, One, [0, 1, 2], 3, 9, 0], qextgen[One, One, [0, 1, 2], 3, 9, 1]);
assert qClass(qRep(a)) == a;
assert qRes[1](a) == zero[QExt(One, One, [1], 3, 9)];
assert qRes[1](b) == qextgen[One, One, [1], 3, 9, 0];
qextdim[One, One, [0, 1, 2], 3, 9];

// The Ext class v_1^4 over A(1) survives restriction to Q_1.
let periodic = extangen[1, 4, 12, 0];
let restricted = qRes[1](periodic);
assert restricted == qextgen[One, One, [1], 4, 12, 0];
assert restricted == qClass(qRes[1](extRep(periodic)));
assert restricted == qRes[1](qRes[0, 1](periodic));
assert qRes[0](periodic) == zero[QExt(One, One, [0], 4, 12)];
show restricted;

// Restriction commutes with composition of endpoint maps.
space Two { a, b };
let into = emap(\pt : One => Two.a);
let out = emap(\pt : Two => One.pt);
assert ecomp(qRes[0, 1](out), qRes[0, 1](into)) == qRes[0, 1](ecomp(out, into));
assert ecomp(qRes[](out), qRes[](into)) == qRes[](eid[One]);
qextdim[Two, Two, [], 0, 0];
