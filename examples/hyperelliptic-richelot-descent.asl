// Frobenius swaps x²+x+1+α and x²+x+1−α and fixes the pair {0,∞}.
rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let zeroF = rvUnramElem[F]{};
let oneF = rvUnramElem[F]{1};
let twoF = rvUnramElem[F]{2};
rvhyperelliptic C[F]{zeroF,twoF,twoF,zeroF,twoF,oneF};
rvjacbase CL[C,W];
rvrichelot DL[CL]{[1,1],[1],[1]}{[1,2],[1],[1]}{[],[1]};
rvrichelotdescend D[C,DL];
rvRichelotDescent[D];
rvHWWWeil[C];
rvHWWWeil[D];
rvJacCount[C];
rvJacCount[D];

let kernel = rvJacPoint[C](zeroF,zeroF);
assert add(kernel,kernel) == zero[RVJac(C)];
assert rvRichelotDescPush[D](kernel) == zero[RVJac(D)];
let p = rvJacPoint[C](oneF,oneF);
let image = rvRichelotDescPush[D](p);
assert image == rvjac[D]{[2],[],[1]}{[2],[1]};
assert rvRichelotDescPush[D](add(p,kernel)) == image;
assert rvJacExtend[DL](image) == rvRichelotPush[DL](rvJacExtend[CL](p));
show image;
// Both groups have order 20. The image has order 10 and the kernel order 2.
rvJacOrder(image);
