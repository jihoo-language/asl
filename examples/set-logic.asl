setomega O[S];
setfamily T[S] {lo:1, hi:1} {lo<hi:[0]};
setfamily Future[S] {lo:0, hi:1} {lo<hi:[]};
setnat include[Future,T] {lo:[], hi:[0]};
let U = setClassify[O](include);
let bottom = setBottom[T,O];
let top = setTop[T,O];

setsieve[O,lo,1];
assert setNot(U) == bottom;
assert setOr(U,setNot(U)) == U;
assert setNot(setNot(U)) == top;
show U;
show top;

setsub Members(U);
assert setClassify[O](setInclusion[Members]) == U;
let exist = setExists(include,setTop[Future,O]);
assert exist == U;
assert setForall(include,setBottom[Future,O]) == bottom;
