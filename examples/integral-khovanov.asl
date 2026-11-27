khbraid Trefoil[2,[1,1,1]];
khbraid Mirror[2,[-1,-1,-1]];
khbraid FigureEight[3,[1,-2,1,-2]];
khbraid Plus[3,[1,1,1,2]];
khbraid Minus[3,[1,1,1,-2]];

khGroup[Trefoil,0,1];
khGroup[Trefoil,0,3];
khGroup[Trefoil,2,5];
khGroup[Trefoil,3,7];
khGroup[Trefoil,3,9];
khGroup[Mirror,-2,-7];
khGroup[Plus,3,7];
khGroup[Minus,3,7];

let t = khGen[Trefoil,3,7,0];
let z = khRep(t);
let w = khNull[Trefoil,3,7,0];
assert khD(z) == zero[KhC(Trefoil,4,7)];
assert khClass(z) == t;
assert khD(w) == khScale[2](z);
assert khClass(khD(w)) == zero[KhH(Trefoil,3,7)];
assert khScale[2](t) == zero[KhH(Trefoil,3,7)];
show t;
show z;
show w;
show khD(w);

let c = khTerm[Trefoil,[],[]];
let first = add(khTerm[Trefoil,[0],[]],add(khTerm[Trefoil,[1],[]],khTerm[Trefoil,[2],[]]));
assert khD(c) == first;
assert khD(first) == zero[KhC(Trefoil,2,5)];

khEuler[Trefoil];
khEuler[FigureEight];
let delta = neg(rvCycloRat[1](add(rvCyclo[20,2],rvCyclo[20,-2])));
assert khEulerAt[Trefoil,3] == rvCycloQMul(delta,rvTLJones[3,2,[1,1,1]]);
assert khEulerAt[FigureEight,3] == rvCycloQMul(delta,rvTLJones[3,3,[1,-2,1,-2]]);
assert khEulerAt[Plus,3] == khEulerAt[Trefoil,3];
assert khEulerAt[Minus,3] == khEulerAt[Trefoil,3];
