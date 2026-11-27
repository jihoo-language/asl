space Circle { a < ab, b < ab, b < bc, c < bc, c < ca, a < ca };

setfamily Twisted[Circle] {a:3, b:3, c:3, ab:3, bc:3, ca:3}
  {a<ab:[1,2,0], b<ab:[0,1,2], b<bc:[0,1,2],
   c<bc:[0,1,2], c<ca:[0,1,2], a<ca:[0,1,2]};
setfamily Trivial[Circle] {a:3, b:3, c:3, ab:3, bc:3, ca:3}
  {a<ab:[0,1,2], b<ab:[0,1,2], b<bc:[0,1,2],
   c<bc:[0,1,2], c<ca:[0,1,2], a<ca:[0,1,2]};
elements Cover[Twisted];
elements Three[Trivial];

betti Cover;
betti Three;
setsections[Twisted];
setsections[Trivial];

setnat rotate[Twisted,Twisted]
  {a:[1,2,0], b:[1,2,0], c:[1,2,0], ab:[1,2,0], bc:[1,2,0], ca:[1,2,0]};
let deck = setTotal[Cover,Cover](rotate);
let p = setProjection[Cover];
assert deck(Cover.a_0) == Cover.a_1;
assert (\x : Cover => p(deck(x))) == p;
assert setComp(rotate,setComp(rotate,rotate)) == setId[Twisted];
let q = setProjection[Three];
let section = setsect[Three,2];
assert (\x : Circle => q(section(x))) == (\x : Circle => x);
