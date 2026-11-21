// The lower target open sees both circles. The upper opens see only the second.
space TwoCircles {
    a<x,a<z,b<x,b<y,c<y,c<z,
    d<u,d<w,e<u,e<v,f<v,f<w
};
space Y {low<mid<top};
let project = tab TwoCircles -> Y {
    a=>Y.low,b=>Y.low,c=>Y.low,x=>Y.low,y=>Y.low,z=>Y.low,
    d=>Y.top,e=>Y.top,f=>Y.top,u=>Y.top,v=>Y.top,w=>Y.top
};
sheaf F[TwoCircles,5] {
    a:1,b:1,c:1,x:1,y:1,z:1,d:1,e:1,f:1,u:1,v:1,w:1
} {
    a<x:[1],a<z:[1],b<x:[1],b<y:[1],c<y:[1],c<z:[1],
    d<u:[1],d<w:[1],e<u:[1],e<v:[1],f<v:[1],f<w:[1]
};
shpush R0[project,F,0];
shpush R1[project,F,1];
shpush R2[project,F,2];
shdim[R0,0];
shdim[R1,0];
shdim[R2,0];

// The two local H1 generators restrict as zero and the surviving top circle.
assert shD(shTerm[R1,[low],0]) == zero[ShC(R1,1)];
assert shD(shTerm[R1,[low],1]) ==
    add(shTerm[R1,[low,mid],0,-1],shTerm[R1,[low,top],0,-1]);
assert shD(shTerm[R1,[mid],0]) ==
    add(shTerm[R1,[low,mid],0],shTerm[R1,[mid,top],0,-1]);
show shD(shTerm[R1,[low],1]);

shpull Pulled[project,R1];
shdim[Pulled,1];
