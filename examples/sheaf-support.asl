// The link is a product of two four-point circles. The vertex lies below the link.
space Link {
    t00<t02, t00<t03, t00<t20, t00<t30, t01<t02, t01<t03,
    t01<t21, t01<t31, t02<t22, t02<t32, t03<t23, t03<t33,
    t10<t12, t10<t13, t10<t20, t10<t30, t11<t12, t11<t13,
    t11<t21, t11<t31, t12<t22, t12<t32, t13<t23, t13<t33,
    t20<t22, t20<t23, t21<t22, t21<t23, t30<t32, t30<t33,
    t31<t32, t31<t33
};
space Cone {
    t00<t02, t00<t03, t00<t20, t00<t30, t01<t02, t01<t03,
    t01<t21, t01<t31, t02<t22, t02<t32, t03<t23, t03<t33,
    t10<t12, t10<t13, t10<t20, t10<t30, t11<t12, t11<t13,
    t11<t21, t11<t31, t12<t22, t12<t32, t13<t23, t13<t33,
    t20<t22, t20<t23, t21<t22, t21<t23, t30<t32, t30<t33,
    t31<t32, t31<t33, north<t00, north<t01, north<t10, north<t11
};
let j = tab Link -> Cone {
    t00=>Cone.t00, t01=>Cone.t01, t02=>Cone.t02, t03=>Cone.t03,
    t10=>Cone.t10, t11=>Cone.t11, t12=>Cone.t12, t13=>Cone.t13,
    t20=>Cone.t20, t21=>Cone.t21, t22=>Cone.t22, t23=>Cone.t23,
    t30=>Cone.t30, t31=>Cone.t31, t32=>Cone.t32, t33=>Cone.t33
};

sheaf F[Link,3] {
    t00:1, t01:1, t02:1, t03:1, t10:1, t11:1, t12:1, t13:1,
    t20:1, t21:1, t22:1, t23:1, t30:1, t31:1, t32:1, t33:1
} {
    t00<t02:[1], t00<t03:[1], t00<t20:[1], t00<t30:[1],
    t01<t02:[1], t01<t03:[1], t01<t21:[1], t01<t31:[1],
    t02<t22:[1], t02<t32:[1], t03<t23:[1], t03<t33:[1],
    t10<t12:[1], t10<t13:[1], t10<t20:[1], t10<t30:[1],
    t11<t12:[1], t11<t13:[1], t11<t21:[1], t11<t31:[1],
    t12<t22:[1], t12<t32:[1], t13<t23:[1], t13<t33:[1],
    t20<t22:[1], t20<t23:[1], t21<t22:[1], t21<t23:[1],
    t30<t32:[1], t30<t33:[1], t31<t32:[1], t31<t33:[1]
};

sheaf Whole[Cone,3] {
    t00:1, t01:1, t02:1, t03:1, t10:1, t11:1, t12:1, t13:1,
    t20:1, t21:1, t22:1, t23:1, t30:1, t31:1, t32:1, t33:1,
    north:1
} {
    t00<t02:[1], t00<t03:[1], t00<t20:[1], t00<t30:[1],
    t01<t02:[1], t01<t03:[1], t01<t21:[1], t01<t31:[1],
    t02<t22:[1], t02<t32:[1], t03<t23:[1], t03<t33:[1],
    t10<t12:[1], t10<t13:[1], t10<t20:[1], t10<t30:[1],
    t11<t12:[1], t11<t13:[1], t11<t21:[1], t11<t31:[1],
    t12<t22:[1], t12<t32:[1], t13<t23:[1], t13<t33:[1],
    t20<t22:[1], t20<t23:[1], t21<t22:[1], t21<t23:[1],
    t30<t32:[1], t30<t33:[1], t31<t32:[1], t31<t33:[1],
    north<t00:[1], north<t01:[1], north<t10:[1], north<t11:[1]
};

let i = tab One -> Cone {pt=>Cone.north};
shextra Costalk[i,Whole];
shdeligne Upper[j,F,1];
shextracomplex PerverseCostalk[i,Upper];
shlocaldim[Costalk,pt,relative,0];
shlocaldim[Costalk,pt,relative,1];
shlocaldim[Costalk,pt,relative,2];
shlocaldim[Costalk,pt,relative,3];
shlocaldim[PerverseCostalk,pt,relative,2];
shlocaldim[PerverseCostalk,pt,relative,3];
let a = shlocalgen[Costalk,pt,punctured,1,0];
let b = shlocalgen[Costalk,pt,punctured,1,1];
show shLocalConnect(a);
show shLocalConnect(b);
let u = shLocalRep(a);
assert shLocalClass[Costalk,pt,punctured](u) == a;
assert shLocalClass[Costalk,pt,relative](shD(u)) == shLocalConnect(a);
assert shLocalClass[Costalk,pt,punctured](add(u,shD(shTerm[Whole,[t00],0]))) == a;
assert shLocalClass[Costalk,pt,relative](shD(add(u,shTerm[Whole,[north,t00],0]))) == shLocalConnect(a);
assert shLocalForget(shLocalConnect(a)) == zero[ShLocal(Costalk,pt,absolute,2)];
let top = shlocalgen[PerverseCostalk,pt,punctured,2,0];
show shLocalConnect(top);
