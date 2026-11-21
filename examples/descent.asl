// Rotate the boundary of a triangle by specifying the map on three open stars.
space Circle { a < ab, b < ab, b < bc, c < bc, c < ca, a < ca };

let atA = tab Star(Circle, a) -> Circle {
  a => Circle.b, ab => Circle.bc, ca => Circle.ab
};
let atB = tab Star(Circle, b) -> Circle {
  b => Circle.c, ab => Circle.bc, bc => Circle.ca
};
let atC = tab Star(Circle, c) -> Circle {
  c => Circle.a, bc => Circle.ca, ca => Circle.ab
};

let turn = sew Circle -> Circle { a => atA, b => atB, c => atC };
let id = \x : Circle => x;
assert res[a](turn) == atA;
assert (\x : Circle => turn(turn(turn(x)))) == id;

show turn(Circle.ab);
betti Circle;
betti Circle * Circle;
homotopy turn ~ id;
