// Multiplication in the mod-two Steenrod algebra is composition, right to left.
let afterTwo = \theta : A(3) => mul(sq[2], theta);
show afterTwo(sq[3]);
assert afterTwo(sq[3]) == add(sq[5], mul(sq[4], sq[1]));
assert mul(sq[1], sq[1]) == zero[A(2)];

// Squaring one operation in degree 128 gives seven admissible words.
show mul(sq[128], sq[128]);
