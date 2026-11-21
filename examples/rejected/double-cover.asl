// This connected double cover exchanges its two sheets on the edge ca.
space Circle { a < ab, b < ab, b < bc, c < bc, c < ca, a < ca };
space Cover {
  a0 < ab0, b0 < ab0, b0 < bc0, c0 < bc0, c0 < ca0, a1 < ca0,
  a1 < ab1, b1 < ab1, b1 < bc1, c1 < bc1, c1 < ca1, a0 < ca1
};
let project = tab Cover -> Circle {
  a0 => Circle.a, a1 => Circle.a, b0 => Circle.b, b1 => Circle.b,
  c0 => Circle.c, c1 => Circle.c, ab0 => Circle.ab, ab1 => Circle.ab,
  bc0 => Circle.bc, bc1 => Circle.bc, ca0 => Circle.ca, ca1 => Circle.ca
};

let atA = tab Star(Circle, a) -> Cover { a => Cover.a0, ab => Cover.ab0, ca => Cover.ca1 };
let atB = tab Star(Circle, b) -> Cover { b => Cover.b0, ab => Cover.ab0, bc => Cover.bc0 };
let atC = tab Star(Circle, c) -> Cover { c => Cover.c0, bc => Cover.bc0, ca => Cover.ca0 };

let id = \x : Circle => x;
assert (\x : Star(Circle, a) => project(atA(x))) == res[a](id);
assert (\x : Star(Circle, b) => project(atB(x))) == res[b](id);
assert (\x : Star(Circle, c) => project(atC(x))) == res[c](id);

// Every local section checked out. At ca the first and third lifts disagree.
let section = sew Circle -> Cover { a => atA, b => atB, c => atC };
