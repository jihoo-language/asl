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

rvmukaimoduli M(rvMukaiIdeal[B,2]);
let Hk3=add(rvK3FromKummerDiv(rvKummerHyperplane[B]),rvK3FromKummerDiv(rvKummerNodeDiv[B,0]));
let H=rvOG10[M,0](rvMukai[0,0](Hk3));
let sigma=rvOG10Exc[M];

// The numerical classes in Mongardi's corrected degree-two example.
let D=rvOG10[M,0](rvMukai[3,3](add(Hk3,Hk3)));
let E=rvOG10[M,-1](rvMukai[2,2](Hk3));
let F=add(add(E,add(E,E)),sigma);
assert add(E,E)==add(D,neg(sigma));
rvOG10Wall(D);
rvOG10Wall(E);
rvOG10Wall(F);
rvOG10Wall(add(sigma,sigma));

let h0=add(add(H,H),sigma);
let h1=add(add(H,H),neg(sigma));
rvOG10Walls(h0,h1);
let middle=rvOG10WallAt[2](h0,h1);
assert middle==neg(sigma);
assert rvOG10Reflect(middle,h0)==h1;
rvOG10Walls(H,h0);
