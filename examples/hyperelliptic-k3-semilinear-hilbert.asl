// The resolved Kummer surface of y²=x⁵−x with full 2-torsion over F9.
rvunram F[3,1]{0,1};
rvunram W[3,1]{1,0,1};
let z=rvUnramElem[F]{};
let o=rvUnramElem[F]{1};
rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
rvjacbase L[C,W];
let zw=rvUnramElem[W]{};
let ow=rvUnramElem[W]{1};
let i=rvUnramElem[W]{0,1};
let e1=rvJacPoint[L](zw,zw);
let f1=rvJacPoint[L](ow,zw);
let e2=add(add(e1,f1),rvJacPoint[L](neg(ow),zw));
let f2=add(add(e1,f1),rvJacPoint[L](i,zw));
rvjactorsionbasis B[L,2](neg(ow)){e1,e2,f1,f2};

// The Enriques switch after Frobenius fixes 12 surface points and
// 232 points on the Hilbert square.
rvk3autom Enriques[B]{switch(2)};
rvK3FrobHilbFixed[Enriques,C,1,1];
rvK3FrobHilbFixed[Enriques,C,2,1];
rvK3FrobHilbTrace[Enriques,C,2,1,0];
rvK3FrobHilbTrace[Enriques,C,2,1,2];
rvK3FrobHilbTrace[Enriques,C,2,1,4];
rvK3FrobHilbTrace[Enriques,C,2,1,6];
rvK3FrobHilbTrace[Enriques,C,2,1,8];

// Frobenius moves the translating node of Shift. Hyperbolic has
// positive entropy and uses the same Hilbert-scheme calculation.
rvk3autom Shift[B]{translate(8)};
rvK3FrobHilbFixed[Shift,C,2,2];
rvk3autom Hyperbolic[B]{project(0),project(1),project(2)};
rvK3FrobHilbFixed[Hyperbolic,C,2,1];
