rvunram W[3,1]{1,0,1};
let z=rvUnramElem[W]{};
let o=rvUnramElem[W]{1};
let i=rvUnramElem[W]{0,1};
rvhyperelliptic C[W]{z,neg(o),z,z,z,o};
let e1=rvJacPoint[C](z,z);
let f1=rvJacPoint[C](o,z);
let e2=add(add(e1,f1),rvJacPoint[C](neg(o),z));
let f2=add(add(e1,f1),rvJacPoint[C](i,z));
rvjactorsionbasis B[C,2](neg(o)){e1,e2,f1,f2};

let w=rvMukaiIdeal[B,2];
rvmukaimoduli M(w);
rvOG10Lattice[M];

let sigma=rvOG10Exc[M];
let G=rvOG10Basis[M,23];
let g=rvOG10[M,0](rvMukaiModuliGlue[M]);
assert add(G,G)==add(g,sigma);
assert G==rvOG10[M,1](zero[RVMukai(B)]);
assert add(sigma,neg(sigma))==zero[RVOG10(M)];

rvOG10Pair(sigma,sigma);
rvOG10Pair(sigma,G);
rvOG10Div(sigma);
rvOG10Fujiki[0](sigma,sigma);
rvOG10Fujiki[1](sigma,G);
rvOG10RR(zero[RVOG10(M)]);
rvOG10RR(sigma);
rvOG10RR(add(sigma,sigma));
