let e = rvHeckeId[3];
let q = rvHeckeQ[3];
let s = rvHeckeT[[1,0,2]];
let t = rvHeckeT[[0,2,1]];
let square = rvHeckeMul(s,s);

// For complete GL3 flags these are the counts with middle position s twice.
rvHeckeCoeff[[0,1,2]](square);
rvHeckeCoeff[[1,0,2]](square);
rvHeckeAt[2,[0,1,2]](square);
rvHeckeAt[2,[1,0,2]](square);
rvHeckeAt[4,[0,1,2]](square);
rvHeckeAt[4,[1,0,2]](square);
let fourth = rvHeckeMul(square,square);
rvHeckeCoeff[[0,1,2]](fourth);
rvHeckeAt[2,[0,1,2]](fourth);
rvHeckeAt[4,[0,1,2]](fourth);

let sts = rvHeckeMul(s,rvHeckeMul(t,s));
assert sts == rvHeckeMul(t,rvHeckeMul(s,t));
assert sts == rvHeckeT[[2,1,0]];
let bs = rvHeckeKL[[1,0,2]];
let bt = rvHeckeKL[[0,2,1]];
assert bs == add(s,e);
assert rvHeckeMul(bs,bs) == rvHeckeMul(add(q,e),bs);
assert rvHeckeMul(bs,rvHeckeMul(bt,bs)) ==
       add(rvHeckeKL[[2,1,0]],rvHeckeMul(q,bs));

// Right multiplication by a rotation shifts the arguments of the window.
let rotate = rvHeckeT[[1,2,3]];
assert rvHeckeMul(s,rotate) == rvHeckeT[[0,2,4]];
assert rvHeckeMul(rotate,s) == rvHeckeT[[2,1,3]];
assert rvHeckeMul(rotate,rvHeckeT[[-1,0,1]]) == e;
assert rvHeckeKL[[1,2,3]] == rotate;
let central = rvHeckeT[[3,4,5]];
assert rvHeckeMul(s,central) == rvHeckeMul(central,s);
let mixed = add(bs,neg(t));
assert rvHeckeMul(rvHeckeMul(mixed,rotate),bt) ==
       rvHeckeMul(mixed,rvHeckeMul(rotate,bt));
