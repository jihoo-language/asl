// Push a twist-only quartic point through a Richelot correspondence to a branch chart.
rvunram F[7,1]{0,1};
let z=rvUnramElem[F]{};
let o=rvUnramElem[F]{1};
let t=rvUnramElem[F]{2};
let a=rvUnramElem[F]{3};
let h=rvUnramElem[F]{4};
let e=rvUnramElem[F]{5};
let s=rvUnramElem[F]{6};
rvhyperelliptic C[F]{z,h,t,e,t,o};
rvrichelot D[C]{[4],[2],[1]}{[],[1]}{[1],[],[1]};
rvjacchart N[D](t);

let point=rvKummer[C](z,o,a,t);
assert rvKummerLiftable(point)==S.lo;
let image=rvKummerRichelot[N](point);
assert rvKummerLiftable(image)==S.lo;

// Here k=3. The twist uses X=x/3 and the target chart uses t'=t/6.
rvhyperelliptic T[F]{z,o,e,s,a,o};
rvrichelot H[T]{[2],[3],[1]}{[],[1]}{[4],[],[1]};
rvjacchart M[H](a);
let p=rvJacPoint[T](o,h);
let expected=rvJacKummerPoint(rvRichelotPush[M](p));
let back=\u:RVKummer(M)=>rvKummer[N](rvKummerCoord[0](u),neg(rvKummerCoord[1](u)),rvKummerCoord[2](u),neg(rvKummerCoord[3](u)));
assert image==back(expected);
show image;
