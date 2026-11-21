space Circle {a<x,a<z,b<x,b<y,c<y,c<z};

// Constant coefficients over F2.
sheaf Constant[Circle,2] {a:1,b:1,c:1,x:1,y:1,z:1} {
    a<x:[1],a<z:[1],b<x:[1],b<y:[1],c<y:[1],c<z:[1]
};
shdim[Constant,0];
shdim[Constant,1];

// Over F3 the one nonidentity restriction gives unipotent rank-two holonomy.
sheaf Unipotent[Circle,3] {a:2,b:2,c:2,x:2,y:2,z:2} {
    a<x:[1,1,0,1],a<z:[1,0,0,1],b<x:[1,0,0,1],
    b<y:[1,0,0,1],c<y:[1,0,0,1],c<z:[1,0,0,1]
};
shdim[Unipotent,0];
shdim[Unipotent,1];
let section = shgen[Unipotent,0,0];
assert shD(shRep(section)) == zero[ShC(Unipotent,1)];
assert shClass(shRep(section)) == section;

// Over F5 holonomy 2 has neither fixed sections nor first cohomology.
sheaf Twisted[Circle,5] {a:1,b:1,c:1,x:1,y:1,z:1} {
    a<x:[2],a<z:[1],b<x:[1],b<y:[1],c<y:[1],c<z:[1]
};
shdim[Twisted,0];
shdim[Twisted,1];

// These coefficient maps have shape 1 by 0.
sheaf Edges[Circle,5] {a:0,b:0,c:0,x:1,y:1,z:1} {
    a<x:[],a<z:[],b<x:[],b<y:[],c<y:[],c<z:[]
};
shdim[Edges,0];
shdim[Edges,1];
show shTerm[Edges,[a,x],0,-1];
