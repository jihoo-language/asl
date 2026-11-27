space Diamond {a<b,a<c,b<d,c<d};
// The four simple sheaves have one nonzero stalk each.
sheaf A[Diamond,3] {a:1,b:0,c:0,d:0} {a<b:[],a<c:[],b<d:[],c<d:[]};
sheaf B[Diamond,3] {a:0,b:1,c:0,d:0} {a<b:[],a<c:[],b<d:[],c<d:[]};
sheaf C[Diamond,3] {a:0,b:0,c:1,d:0} {a<b:[],a<c:[],b<d:[],c<d:[]};
sheaf D[Diamond,3] {a:0,b:0,c:0,d:1} {a<b:[],a<c:[],b<d:[],c<d:[]};
shExtDim[A,B,1];
shExtDim[A,C,1];
shExtDim[A,D,1];
shExtDim[A,D,2];
let ab = shClass(shHomTerm[A,B,[a,b],0,0]);
let ac = shClass(shHomTerm[A,C,[a,c],0,0]);
let bd = shClass(shHomTerm[B,D,[b,d],0,0]);
let cd = shClass(shHomTerm[C,D,[c,d],0,0]);
let left = shExtCompose(bd,ab);
let right = shExtCompose(cd,ac);
assert add(left,right) == zero[ShExt(A,D,2)];
assert right == shExtGen[A,D,2,0];
show left;
show right;
show shRep(left);
// A declared sheaf map gives an actual natural 0-cochain.
shmap identityA[A,A] {a:[1],b:[],c:[],d:[]};
assert shMapC[identityA] == shHomId[A];
assert shExtCompose(left,shClass(shMapC[identityA])) == left;
assert shExtCompose(shClass(shHomId[D]),left) == left;
card \x : ShExt(A,D,2) => S.hi;
