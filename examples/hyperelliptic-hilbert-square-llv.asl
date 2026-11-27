rvunram W[3,1]{1,0,1};
let z=rvUnramElem[W]{};
let o=rvUnramElem[W]{1};
let i=rvUnramElem[W]{0,1};
rvhyperelliptic C[W]{z,neg(o),z,z,z,o};
let e1=rvJacPoint[C](z,z);
let f1=rvJacPoint[C](o,z);
let e2=add(add(e1,f1),rvJacPoint[C](neg(o),z));
let f2=add(add(e1,f1),rvJacPoint[C](i,z));
rvjactorsionbasis B[C,2](neg(o)){e1,e2,f1,f2};

let a=rvHilbSquareDiv[0](rvK3Lattice[B]{[0,2,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]});
let b=rvHilbSquareDiv[0](rvK3Lattice[B]{[0,0,0,0,2,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]});
let delta=rvHilbSquareBasis[B,2,22];
let point=rvHilbSquareBasis[B,8,0];

// The two isotropic classes have q(a,b)=2.
assert rvHilbSquareLower(a,a)==zero[RVHilbSquare(B,0)];
assert rvHilbSquareRotate(a,b,a)==rvHilbSquareScale[2,1](a);

// The isotropic square is primitive relative to delta.
let a2=rvHilbSquareCup(a,a);
let delta2=rvHilbSquareCup(delta,delta);
assert rvHilbSquarePrimitive(delta,add(a2,delta2))==a2;
assert rvHilbSquareCup(delta,a2)==zero[RVHilbSquare(B,6)];

rvHilbSquareIntegral(rvHilbSquareCup(delta,rvHilbSquareLower(delta,point)));
rvHilbSquareIntegral(rvHilbSquareCup(delta,rvHilbSquareLambda(delta,point)));
