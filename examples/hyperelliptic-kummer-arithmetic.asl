// Arithmetic on quartic points that lift only to the quadratic twist.
rvunram F[5,1]{0,1};
let z=rvUnramElem[F]{};
let o=rvUnramElem[F]{1};
let t=rvUnramElem[F]{2};
let three=rvUnramElem[F]{3};
let four=rvUnramElem[F]{4};
rvhyperelliptic C[F]{o,o,z,three,o,o};
rvhyperelliptic T[F]{three,o,z,t,three,o};

// The twist has X=x/2. Its Kummer coordinates return to C with weights 1,2,4,3.
let back=\k:RVKummer(T)=>rvKummer[C](rvKummerCoord[0](k),rvUnramMul(t,rvKummerCoord[1](k)),rvUnramMul(four,rvKummerCoord[2](k)),rvUnramMul(three,rvKummerCoord[3](k)));
let p=rvJacPoint[T](three,four);
let q=rvJacPoint[T](t,four);
let a=rvKummer[C](z,o,o,o);
let b=rvKummer[C](z,o,four,o);
assert rvKummerLiftable(a)==S.lo;
assert rvKummerLiftable(b)==S.lo;
assert a==back(rvJacKummerPoint(p));
assert b==back(rvJacKummerPoint(q));
let sum=back(rvJacKummerPoint(add(p,q)));
let diff=back(rvJacKummerPoint(add(p,neg(q))));
assert rvKummerDiffAdd(a,b,diff)==sum;
assert rvKummerDiffAdd(a,b,sum)==diff;
assert rvKummerScale[-3](a)==back(rvJacKummerPoint(rvJacScale[3](p)));
show sum;
show diff;
