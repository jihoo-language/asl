// In H*(BC8; F2), x is the nonzero degree-one character.
let x = bar[8]{1, 3, 5, 7};
let b2 = bar[8]{2, 3, 6, 7};
assert x == binom[8, 1];
assert b2 == binom[8, 2];
let b3 = binom[8, 3];
let b4 = binom[8, 4];
let b5 = binom[8, 5];
let b6 = binom[8, 6];
let b7 = binom[8, 7];

// The full four-input defining system spells out the five inner witnesses.
let four = massey(x, x, x, x) {
    0 < 2 => b2, 1 < 3 => b2, 2 < 4 => b2,
    0 < 3 => b3, 1 < 4 => b3
};
assert four == powerdef(x, b2, b3);

// One output per arity, from 2 to 8. The last output is the carry class.
show powerdef(x);
show powerdef(x, b2);
show four;
show powerdef(x, b2, b3, b4);
show powerdef(x, b2, b3, b4, b5);
show powerdef(x, b2, b3, b4, b5, b6);
show powerdef(x, b2, b3, b4, b5, b6, b7);

// The same triple on C2 gives two different classes as one witness varies.
let y = bar[2]{1};
let z = zero[Bar(2)];
let triple = \w : Bar(2) => massey(y, z, y) { 0 < 2 => z, 1 < 3 => w };
show triple(z);
show triple(y);
