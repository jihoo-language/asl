// The base has two points in each of three levels. The fibres are A<C>B.
// The sheaf vanishes on B. The map from A to C is the stalk differential.
space Y {b0<b2, b0<b3, b1<b2, b1<b3, b2<b4, b2<b5, b3<b4, b3<b5};
space X {x0a<x2a, x0b<x2b, x0c<x2c, x0a<x3a, x0b<x3b, x0c<x3c, x1a<x2a, x1b<x2b, x1c<x2c, x1a<x3a, x1b<x3b, x1c<x3c, x2a<x4a, x2b<x4b, x2c<x4c, x2a<x5a, x2b<x5b, x2c<x5c, x3a<x4a, x3b<x4b, x3c<x4c, x3a<x5a, x3b<x5b, x3c<x5c, x0a<x0c, x0b<x0c, x1a<x1c, x1b<x1c, x2a<x2c, x2b<x2c, x3a<x3c, x3b<x3c, x4a<x4c, x4b<x4c, x5a<x5c, x5b<x5c};
let pi = tab X -> Y {x0a => Y.b0, x0b => Y.b0, x0c => Y.b0, x1a => Y.b1, x1b => Y.b1, x1c => Y.b1, x2a => Y.b2, x2b => Y.b2, x2c => Y.b2, x3a => Y.b3, x3b => Y.b3, x3c => Y.b3, x4a => Y.b4, x4b => Y.b4, x4c => Y.b4, x5a => Y.b5, x5b => Y.b5, x5c => Y.b5};
sheaf F[X,3] {x0a:1, x0b:0, x0c:1, x1a:1, x1b:0, x1c:1, x2a:2, x2b:0, x2c:2, x3a:2, x3b:0, x3c:2, x4a:2, x4b:0, x4c:2, x5a:2, x5b:0, x5c:2} {
    x0a<x2a:[1,0],
    x0b<x2b:[],
    x0c<x2c:[1,0],
    x0a<x3a:[1,0],
    x0b<x3b:[],
    x0c<x3c:[1,0],
    x1a<x2a:[1,0],
    x1b<x2b:[],
    x1c<x2c:[0,1],
    x1a<x3a:[1,0],
    x1b<x3b:[],
    x1c<x3c:[0,1],
    x2a<x4a:[1,1,0,1],
    x2b<x4b:[],
    x2c<x4c:[1,0,0,1],
    x2a<x5a:[1,0,0,1],
    x2b<x5b:[],
    x2c<x5c:[1,0,0,1],
    x3a<x4a:[1,0,0,1],
    x3b<x4b:[],
    x3c<x4c:[1,0,0,1],
    x3a<x5a:[1,0,0,1],
    x3b<x5b:[],
    x3c<x5c:[1,0,0,1],
    x0a<x0c:[0],
    x0b<x0c:[],
    x1a<x1c:[0],
    x1b<x1c:[],
    x2a<x2c:[0,1,0,2],
    x2b<x2c:[],
    x3a<x3c:[0,1,0,2],
    x3b<x3c:[],
    x4a<x4c:[0,1,0,2],
    x4b<x4c:[],
    x5a<x5c:[0,1,0,2],
    x5b<x5c:[]
};

shleray K[pi,F];
shpush R0[pi,F,0];
shpush R1[pi,F,1];
shdim[R0,0]; shdim[R0,1]; shdim[R0,2];
shdim[R1,0]; shdim[R1,1]; shdim[R1,2];
shedim[K,2,0,1]; shedim[K,2,2,0];
let a = shegen[K,2,0,1,0];
let hit = shED(a);
show hit;
assert shENext(hit) == zero[ShE(K,3,2,0)];
shedim[K,3,0,1]; shedim[K,3,2,0];

let u = shTerm[F,[x0a],0];
let v = shTerm[F,[x0a,x2a,x4c],0];
assert shD(shLerayIncl[K](u)) == shLerayIncl[K](shD(u));
assert shD(shLerayIncl[K](v)) == shLerayIncl[K](shD(v));
shdim[F,0]; shdim[F,1]; shdim[F,2]; shdim[F,3];
shhyperdim[K,0]; shhyperdim[K,1]; shhyperdim[K,2]; shhyperdim[K,3];
