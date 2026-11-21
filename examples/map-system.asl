let h0 = mapextgen[One, One, 1, 1, 0];
let h1 = mapextgen[One, One, 1, 2, 0];
let u = mword[One, 0, 0, One, 0, 0](add(cword{[3]}, cword{[0, 1]}));
let v = mword[One, 0, 0, One, 0, 0](cword{[0, 1]});

// The vertices follow the path from source to target.
let omega = eCurvature(h0, eRep(h1), h0){0 < 2 => u, 1 < 3 => v};
assert eD(omega) == zero[MapCobar(One, One, 3, 4)];
let chosen = eMassey(h0, h1, h0){0 < 2 => u, 1 < 3 => v};
assert chosen == eClass(omega);
assert chosen == ecomp(h1, h1);
show chosen;
