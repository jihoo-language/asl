// The target is a monic quintic over F9 and the image splits into two points.
rvunram U[3,1]{1,0,1};
let zeroU = rvUnramElem[U]{};
let oneU = rvUnramElem[U]{1};
let a = rvUnramElem[U]{0,1};
let b = rvUnramElem[U]{1,1};
let c = rvUnramElem[U]{2,1};
rvhyperelliptic C[U]{zeroU,b,zeroU,c,zeroU,oneU};
rvrichelot D[C]{[1],[],[1]}{[1,1],[],[1]}{[],[1]};
let image = rvRichelotPoint[D](b,a);
let x1 = c;
let y1 = rvUnramElem[U]{2,2};
let x2 = rvUnramElem[U]{1,2};
let y2 = x2;
assert image == add(rvJacPoint[D](x1,y1),rvJacPoint[D](x2,y2));
assert rvRichelotPoint[D](b,neg(a)) == neg(image);

rvrichelotdual E[D];
let dual = add(rvRichelotPoint[E](x1,y1),rvRichelotPoint[E](x2,y2));
let doubled = rvJacScale[2](rvJacPoint[E](b,add(a,a)));
assert dual == doubled;
show image;
rvHWWWeil[C];
rvHWWWeil[D];

// A sextic target needs a branch chart. The full image here lies at old infinity.
rvunram F[7,1]{0,1};
let zeroF = rvUnramElem[F]{};
let oneF = rvUnramElem[F]{1};
let twoF = rvUnramElem[F]{2};
let fourF = rvUnramElem[F]{4};
let fiveF = rvUnramElem[F]{5};
rvhyperelliptic H[F]{zeroF,fourF,twoF,fiveF,twoF,oneF};
rvrichelot R[H]{[4],[2],[1]}{[],[1]}{[1],[],[1]};
rvjacchart N[R](twoF);
let charted = rvRichelotPoint[N](twoF,oneF);
assert charted == rvJacScale[2](rvJacPoint[N](zeroF,oneF));
show charted;
