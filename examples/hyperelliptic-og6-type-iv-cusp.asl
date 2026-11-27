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

rvog6 M(rvAbMukai[B,1,-1]{0,0,0,0,0,0});
let e=rvOG6[M,0](rvAbMukai[B,0,0]{1,0,0,0,0,0});
let f=rvOG6[M,0](rvAbMukai[B,0,0]{0,1,0,0,0,0});
let eDual=rvOG6[M,0](rvAbMukai[B,0,0]{0,0,0,0,0,-1});
let fDual=rvOG6[M,0](rvAbMukai[B,0,0]{0,0,0,0,1,0});
let hE=rvOG6[M,0](rvAbMukai[B,0,0]{0,0,1,0,0,0});
let hF=rvOG6[M,0](rvAbMukai[B,0,0]{0,0,0,-1,0,0});
let h=add(hE,hF);
let delta=rvOG6[M,0](rvAbMukai[B,1,1]{0,0,0,0,0,0});

rvog6cusp P(h,e,f);
// The quotient has an integral basis of three square-minus-two classes.
rvOG6CuspGram[P];
let u=rvOG6CuspBasis[P,0];
rvOG6Pair(u,h);
rvOG6Pair(u,u);

let q=rvOG6CuspTransvection[P,1](neg(delta),eDual);
let q2=rvOG6CuspTransvection[P,0](neg(delta),q);
let q3=rvOG6CuspTransvection[P,1](delta,q2);
let commutator=rvOG6CuspTransvection[P,0](delta,q3);
let central=rvOG6CuspCentral[P,-2](eDual);
assert commutator==central;
assert central==add(eDual,neg(add(f,f)));
