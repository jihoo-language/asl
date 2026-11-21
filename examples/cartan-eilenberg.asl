// The central extension by the exterior Milnor primitive Q_1 in A(1).
let v = ceClass[2, 0](resA[1](cword{[0, 1]}));
let h0 = ceClass[2, 1](resA[1](cword{[1]}));
let h1 = ceClass[2, 1](resA[1](cword{[2]}));
assert ceD(v) == ceMul(h0, h1);
assert ceNext(ceD(v)) == zero[CEAn(1, 3, 2, 3, 2)];
show ceD(v);

let square = ceNext(ceMul(v, v));
let h = ceNext(h1);
assert ceD(square) == ceMul(h, ceMul(h, h));
assert ceNext(ceD(square)) == zero[CEAn(1, 4, 3, 6, 3)];
show ceD(square);

// In cobar length five, filtration six is empty. This lift is closed.
let fourth = ceMul(square, square);
let survivor = ceNext(ceNext(ceNext(fourth)));
let periodic = extClass(ceRep(survivor));
assert periodic == extangen[1, 4, 12, 0];
show periodic;
