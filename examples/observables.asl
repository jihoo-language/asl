// Each map Fork -> S observes an open subset of Fork.
space Fork { root < left, root < right };
let leftOpen = tab Fork -> S { root => S.lo, left => S.hi, right => S.lo };
let rightOpen = tab Fork -> S { root => S.lo, left => S.lo, right => S.hi };

// A point acts on its observables by evaluation.
let eta = \x : Fork => \u : Fork -> S => u(x);
show eta(Fork.left)(leftOpen);
show eta(Fork.root)(leftOpen);

// Intersection lifts the meet on S to the function space Fork -> S.
let meet = \u : S => tab S -> S { lo => S.lo, hi => u };
let cap = \u : Fork -> S => \v : Fork -> S => \x : Fork => meet(u(x))(v(x));
let empty = \x : Fork => S.lo;
assert cap(leftOpen)(rightOpen) == empty;
assert cap(leftOpen)(leftOpen) == leftOpen;

// One -> Fork is isomorphic to Fork. The shortest zigzag goes through root.
let leftPoint = \p : One => Fork.left;
let rightPoint = \p : One => Fork.right;
homotopy leftPoint ~ rightPoint;
