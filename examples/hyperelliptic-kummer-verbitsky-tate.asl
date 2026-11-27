rvverbitskytate T{3};
rvVerbitskyTatePeriod[T,0,ordinary,2,0,0];
rvVerbitskyTateCheck[T,0,ordinary,0];
rvVerbitskyTateCheck[T,0,stable,0];
rvVerbitskyTatePeriod[T,0,stable,-2000000,1,0];
rvVerbitskyTateComplexMap[T,0,stable,-2000001,cochain,0];

let x=rvVerbitskyTateTerm[T,ordinary]{1}{0}{0};
let e=rvVerbitskyTateTerm[T,ordinary]{0}{1}{0};
let t=rvVerbitskyTateTerm[T,ordinary]{0}{0}{1};
assert rvVerbitskyTateBracket(e,x)==x;
assert rvVerbitskyTateBracket(e,t)==rvVerbitskyTateScale[-3,1](t);
let hi=rvVerbitskyTateTerm[T,stable]{0}{1}{-1000000};
show rvVerbitskyTateBV(hi);

rvhodgewick B[21]{0,1}{[2]};
let a=rvVerbitskyTateTerm[T,ordinary]{0}{1}{3};
let bar=rvVerbitskyTateBar[B,1](a);
assert rvVerbitskyHHBV(bar)==rvVerbitskyTateBar[B,1](rvVerbitskyTateBV(a));
assert rvVerbitskyHodgeProject[4](bar)==bar;

rvverbitskytate U{2,2};
rvhodgewick C[20]{0,1}{[0,1],[1,0]};
let e1=rvVerbitskyTateTerm[U,ordinary]{0,0}{1,0}{0,0};
let e2=rvVerbitskyTateTerm[U,ordinary]{0,0}{0,1}{0,0};
let both=rvVerbitskyTateCup(e1,e2);
assert rvVerbitskyTateCup(e2,e1)==neg(both);
assert rvVerbitskyHHBV(rvVerbitskyTateBar[C,1](both))==rvVerbitskyTateBar[C,1](rvVerbitskyTateBV(both));
show rvVerbitskyTateBV(both);
