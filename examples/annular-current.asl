// The five essential circles carry alternating standard and dual actions.
khbraid Five[5,[]];
khIrreps[Five,0,0];
let top = khHighest[Five,0,0,5,0];
assert top == khAnnular[5](khTerm[Five,[],[]]);
let rung = khF(top);
assert khE(rung) == khScale[5](top);
assert khE(khF(rung)) == khScale[8](rung);
let bottom = khF(khF(khF(khF(rung))));
assert bottom == khAnnular[-5](khTerm[Five,[],[0,1,2,3,4],120]);
assert khF(bottom) == zero[KhAC(Five,0,-7,-7)];
let other = khHighest[Five,0,0,3,3];
assert khE(other) == zero[KhAC(Five,0,5,5)];
assert khD0(other) == zero[KhAC(Five,1,3,3)];
show khClass(other);
khActionCheck[Five];

// The lowering odd operator detects a nonzero class for the positive crossing.
khbraid Crossing[2,[1]];
let high = khAGen[Crossing,0,3,2,0];
let odd = khVMinus(high);
assert odd == khAGen[Crossing,1,3,0,0];
assert khVMinus(odd) == zero[KhAH(Crossing,2,3,-2)];
show odd;

// The same even action on actual integral groups retains torsion.
khbraid Twisted[3,[1,1,1]];
let torsion = khAGen[Twisted,3,8,1,0];
assert khScale[2](torsion) == zero[KhAH(Twisted,3,8,1)];
assert khE(khF(torsion)) == torsion;
khAnnularGroup[Twisted,3,8,1];

// Phi_0 is an explicit null for the nonzero odd defects on this cycle.
khbraid Figure[3,[1,-2,1,-2]];
let z = khAnnular[1](khTerm[Figure,[],[1,2]]);
assert khD0(z) == zero[KhAC(Figure,-1,-3,1)];
let defect = add(khVMinus(khWPlus(z)), khWPlus(khVMinus(z)));
let witness = khScale[2](khLee0(z));
assert khForgetK(defect) == add(
   add(khTerm[Figure,[0,1],[],2],khTerm[Figure,[0,3],[],2]),
   add(khTerm[Figure,[1,2],[],-2],khTerm[Figure,[2,3],[],2])
);
assert khD0(witness) == defect;
assert khClass(defect) == zero[KhAH(Figure,0,1,1)];
assert khVZero(khVZero(z)) == neg(khD0(khLee0(z)));
assert khWPlus(z) == khScale[-2](khLeePlus(z));
khActionCheck[Figure];
