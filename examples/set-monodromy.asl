// Synchronized GL(3,F2) monodromy on the points and lines of the Fano plane.
space Base { x, y, a, b, c, x<a, x<b, x<c, y<a, y<b, y<c };
setfamily Points[Base] {x:7, y:7, a:7, b:7, c:7} {
    x<a:[0,1,2,3,4,5,6], x<b:[0,1,2,3,4,5,6], x<c:[0,1,2,3,4,5,6],
    y<a:[0,1,2,3,4,5,6], y<b:[1,3,5,2,0,6,4], y<c:[0,2,1,3,4,6,5]
};
setfamily Lines[Base] {x:7, y:7, a:7, b:7, c:7} {
    x<a:[0,1,2,3,4,5,6], x<b:[0,1,2,3,4,5,6], x<c:[0,1,2,3,4,5,6],
    y<a:[0,1,2,3,4,5,6], y<b:[2,3,6,0,1,4,5], y<c:[2,1,0,3,6,5,4]
};
elements PointCover[Points];
elements LineCover[Lines];
betti PointCover;
betti LineCover;

setmonodromy M[Points, Lines, x];
setmonodromy Transposed[Lines, Points, x];
setmonodromy PP[Points, Points, x];
setmonodromy LL[Lines, Lines, x];
setmonosize[M];
setmonochars[M];
setmonostabilizer[M,0];
setmonomaps[M];
setmonoorbitals[M];

shlinearize P5[Points,5];
shlinearize L5[Lines,5];
shorbital N[M,P5,L5,1];
shorbital Nt[Transposed,L5,P5,1];
shorbitals Back[Transposed,L5,P5] {0:-1,1:2};
shorbitals Jp[PP,P5,P5] {0:1,1:1};
shorbitals Jl[LL,L5,L5] {0:1,1:1};
shmaprank[N,x];
shExtDim[P5,L5,0];

let n = shMapC[N];
let nt = shMapC[Nt];
let jp = shMapC[Jp];
let jl = shMapC[Jl];
let ip = shHomId[P5];
let il = shHomId[L5];
assert shExtCompose(nt,n) == add(shScale[2](ip),jp);
assert shExtCompose(n,nt) == add(shScale[2](il),jl);
let inverse = shExtCompose(add(shScale[3](ip),shScale[3](jp)),nt);
assert shMapC[Back] == inverse;
assert shExtCompose(inverse,n) == ip;
assert shExtCompose(n,inverse) == il;
assert shExtCompose(shClass(n),shClass(inverse)) == shClass(il);

shlinearize P2[Points,2];
shlinearize L2[Lines,2];
shorbitals Z[M,P2,L2] {};
shorbitals J[M,P2,L2] {0:1,1:-1};
shorbital Inc[M,P2,L2,1];
shorbital Off[M,P2,L2,0];
shExtDim[P2,L2,0];
shmaprank[Z,x];
shmaprank[J,x];
shmaprank[Off,x];
shmaprank[Inc,x];
