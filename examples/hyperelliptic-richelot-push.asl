// A nonsplit support and a tangent divisor over F7.
rvunram F[7,1]{0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
let t = rvUnramElem[F]{2};
let h = rvUnramElem[F]{4};
let e = rvUnramElem[F]{5};
rvhyperelliptic C[F]{z,h,t,e,t,o};
rvrichelot D[C]{[4],[2],[1]}{[],[1]}{[1],[],[1]};
rvjacchart N[D](t);
let closed = rvjac[C]{[5],[3],[1]}{[6],[1]};
let image = rvRichelotPush[N](closed);
assert image == rvjac[N]{[2],[3],[1]}{[],[3]};
show image;

let tangent = rvjac[C]{[4],[3],[1]}{[2],[3]};
let p = rvJacPoint[C](t,o);
assert tangent == add(p,p);
let tangentImage = rvRichelotPush[N](tangent);
assert tangentImage == rvJacScale[2](rvRichelotPoint[N](t,o));
show tangentImage;

// Check the first push after adjoining the two conjugate points.
rvunram W[7,1]{1,0,1};
rvjacbase L[C,W];
rvrichelot J[L]{[4],[2],[1]}{[],[1]}{[1],[],[1]};
rvjacchart M[J](rvUnramElem[W]{2});
let x = rvUnramElem[W]{2,1};
let y = rvUnramElem[W]{1,1};
let xbar = rvUnramFrob[1](x);
let ybar = rvUnramFrob[1](y);
assert rvJacExtend[L](closed) == add(rvJacPoint[L](x,y),rvJacPoint[L](xbar,ybar));
assert rvJacExtend[M](image) == add(rvRichelotPoint[M](x,y),rvRichelotPoint[M](xbar,ybar));
