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

let w=rvAbMukai[B,1,-1]{0,0,0,0,0,0};
rvog6 M(w);
let delta=rvOG6[M,0](rvAbMukai[B,1,1]{0,0,0,0,0,0});
let e=rvOG6ExcHalf[M];
let exceptional=rvOG6Exc[M];
let exchange=add(delta,e);
let small=rvOG6[M,0](rvAbMukai[B,0,0]{0,1,0,0,-1,0});

rvOG6Orbit(delta);
rvOG6Orbit(e);
rvOG6SameOrbit(delta,e);
rvOG6SameOrbit(e,exceptional);

assert rvOG6Reflect(exchange,delta)==neg(e);
assert rvOG6Reflect(exchange,e)==neg(delta);
rvOG6Wall(exchange);
rvOG6Wall(exceptional);
rvOG6Wall(small);

let isotropic=rvOG6[M,1](rvAbMukai[B,0,0]{0,1,0,0,1,0});
rvOG6Orbit(isotropic);
