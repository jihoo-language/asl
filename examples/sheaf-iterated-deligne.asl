// The link is a product of two four-point circles. Both cone tips are below it.
space Link {
    t00<t02, t00<t03, t00<t20, t00<t30, t01<t02, t01<t03,
    t01<t21, t01<t31, t02<t22, t02<t32, t03<t23, t03<t33,
    t10<t12, t10<t13, t10<t20, t10<t30, t11<t12, t11<t13,
    t11<t21, t11<t31, t12<t22, t12<t32, t13<t23, t13<t33,
    t20<t22, t20<t23, t21<t22, t21<t23, t30<t32, t30<t33,
    t31<t32, t31<t33
};
space Susp {
    t00<t02, t00<t03, t00<t20, t00<t30, t01<t02, t01<t03,
    t01<t21, t01<t31, t02<t22, t02<t32, t03<t23, t03<t33,
    t10<t12, t10<t13, t10<t20, t10<t30, t11<t12, t11<t13,
    t11<t21, t11<t31, t12<t22, t12<t32, t13<t23, t13<t33,
    t20<t22, t20<t23, t21<t22, t21<t23, t30<t32, t30<t33,
    t31<t32, t31<t33, north<t00, north<t01, north<t10, north<t11,
    south<t00, south<t01, south<t10, south<t11
};
let j = tab Link -> Susp {
    t00=>Susp.t00, t01=>Susp.t01, t02=>Susp.t02, t03=>Susp.t03,
    t10=>Susp.t10, t11=>Susp.t11, t12=>Susp.t12, t13=>Susp.t13,
    t20=>Susp.t20, t21=>Susp.t21, t22=>Susp.t22, t23=>Susp.t23,
    t30=>Susp.t30, t31=>Susp.t31, t32=>Susp.t32, t33=>Susp.t33
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


space DoubleSusp {
    t00<t02, t00<t03, t00<t20, t00<t30, t01<t02, t01<t03,
    t01<t21, t01<t31, t02<t22, t02<t32, t03<t23, t03<t33,
    t10<t12, t10<t13, t10<t20, t10<t30, t11<t12, t11<t13,
    t11<t21, t11<t31, t12<t22, t12<t32, t13<t23, t13<t33,
    t20<t22, t20<t23, t21<t22, t21<t23, t30<t32, t30<t33,
    t31<t32, t31<t33, north<t00, north<t01, north<t10, north<t11,
    south<t00, south<t01, south<t10, south<t11, first<north, first<south,
    second<north, second<south
};
let k = tab Susp -> DoubleSusp {
    t00=>DoubleSusp.t00, t01=>DoubleSusp.t01, t02=>DoubleSusp.t02, t03=>DoubleSusp.t03,
    t10=>DoubleSusp.t10, t11=>DoubleSusp.t11, t12=>DoubleSusp.t12, t13=>DoubleSusp.t13,
    t20=>DoubleSusp.t20, t21=>DoubleSusp.t21, t22=>DoubleSusp.t22, t23=>DoubleSusp.t23,
    t30=>DoubleSusp.t30, t31=>DoubleSusp.t31, t32=>DoubleSusp.t32, t33=>DoubleSusp.t33,
    north=>DoubleSusp.north, south=>DoubleSusp.south
};

shdeligne FirstLower[j,F,0];
shdeligne FirstUpper[j,F,1];
shdelignecomplex Lower[k,FirstLower,1];
shdelignecomplex Upper[k,FirstUpper,1];

shhyperdim[Lower,0];
shhyperdim[Upper,4];
shhyperdim[Lower,1];
shhyperdim[Upper,3];
shhyperdim[Lower,2];
shhyperdim[Upper,2];
shhyperdim[Lower,3];
shhyperdim[Upper,1];
shhyperdim[Lower,4];
shhyperdim[Upper,0];
