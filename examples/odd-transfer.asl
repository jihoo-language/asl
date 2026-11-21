let x = extClass(pword[3]{([1], [])});
let a = extClass(pword[3]{([], [0])});
let triple = ainfty(x, x, x);
assert triple == pscale[2](pextgen[3, 2, 12, 0]);
assert alift(x, x) == pword[3]{([2], [])};
assert ainfty(x, a, x) == pscale[2](pextgen[3, 2, 9, 0]);
show triple;

let raw = pword[5]{([1], [])};
let y = extClass(raw);
let fifth = ainfty(y, y, y, y, y);
assert fifth == pextgen[5, 2, 40, 0];
assert ainfty(pscale[2](y), pscale[3](y), y, y, y) == fifth;
assert alift(pscale[2](y), pscale[4](y), y, y) == pscale[3](alift(y, y, y, y));
let b2 = pscale[3](pword[5]{([2], [])});
let b3 = pword[5]{([3], [])};
let b4 = pscale[4](pword[5]{([4], [])});
assert pPowerDef(raw, b2, b3, b4) == pscale[4](fifth);
show fifth;
