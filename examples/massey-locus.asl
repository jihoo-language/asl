// H¹ of the graph with four parallel edges has dimension three.
space G { p < e0, q < e0, p < e1, q < e1, p < e2, q < e2, p < e3, q < e3 };
let l = \v : G * G => fst(v);
let r = \v : G * G => snd(v);
let u = gen[G, 1, 0];
let v = gen[G, 1, 1];
let w = gen[G, 1, 2];
let ax = cup(pull(l, u), pull(r, u));
let by = cup(pull(l, v), pull(r, v));
let cz = cup(pull(l, w), pull(r, w));

// Every adjacent input is zero. The fourfold product still has 344 values.
let z = zero[G * G, 1];
let locus = bracket(z, z, z, z);
let rankTwo = add(ax, by);
let rankThree = add(rankTwo, cz);
card locus;
show (locus(rankTwo), locus(cz));
show locus(rankThree);

// On a torus the empty bracket and the singleton containing zero differ.
space C { p < e, q < e, p < f, q < f };
let t = gen[C, 1, 0];
let a = pull((\v : C * C => fst(v)), t);
let b = pull((\v : C * C => snd(v)), t);
let o = zero[C * C, 1];
card bracket(a, b, a);
card bracket(o, a, o);
card bracket(a, o, b);
