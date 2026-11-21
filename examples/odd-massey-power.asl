// At p=3, 1/2! = 2. The compact and explicit defining systems agree.
let x3 = pword[3]{([1], [])};
let b2at3 = pscale[2](pword[3]{([2], [])});
let value3 = pPowerDef(x3, b2at3);
let positive3 = add(pword[3]{([1], []), ([2], [])}, pword[3]{([2], []), ([1], [])});
assert cobarD(positive3) == zero[PCobar(3, 3, 12)];
assert value3 == pscale[2](extClass(positive3));
assert value3 == pMassey(extClass(x3), x3, x3){0<2 => b2at3, 1<3 => b2at3};
assert pCurvature(x3, x3, x3){0<2 => b2at3, 1<3 => b2at3} == pscale[2](positive3);
show value3;

// At p=5 the coefficients of 1/2!, 1/3!, and 1/4! are 3, 1, and 4.
let x5 = pword[5]{([1], [])};
let b2at5 = pscale[3](pword[5]{([2], [])});
let b3at5 = pword[5]{([3], [])};
let b4at5 = pscale[4](pword[5]{([4], [])});
let value5 = pPowerDef(x5, b2at5, b3at5, b4at5);
let positive5 = add(
  add(pword[5]{([1], []), ([4], [])}, pscale[2](pword[5]{([2], []), ([3], [])})),
  add(pscale[2](pword[5]{([3], []), ([2], [])}), pword[5]{([4], []), ([1], [])})
);
assert cobarD(positive5) == zero[PCobar(5, 3, 40)];
assert value5 == pscale[4](extClass(positive5));
show value5;
