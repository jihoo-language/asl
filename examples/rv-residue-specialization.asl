// These row spaces have the same tropical Plücker data, with every minor a unit.
let trop2 = rvValMat[7,4]{[1,0,1,1],[0,1,1,2]};
let trop3 = rvValMat[7,4]{[1,0,1,1],[0,1,1,3]};
assert trop2 == trop3;
let x = rvResSpecialize[7,4,[0,0,0,0]]{[1,0,1,1],[0,1,1,2]};
let y = rvResSpecialize[7,4,[0,0,0,0]]{[1,0,1,1],[0,1,1,3]};
show x;
show y;
rvResSupport(x);
rvResSupport(y);
rvResLeq(x,y);
rvResContains[[0,1,1,2]](x);

// Filtration steps remain nested. Negative weights give the dual specialization.
rvphodge H[7] {
    {(1,1),(7,1),(49,1)}, {},
    {(2,[1,7,1]),(1,[0,1,0]),(0,[0,0,1])}
};
let step = rvResFil[H,2,[1,2,0]];
let plane = rvResFil[H,1,[1,2,0]];
rvResLeq(step,plane);
assert rvResIntersection(step,plane) == step;
assert rvResSum(step,plane) == plane;
assert rvResOrthogonal(plane) == rvResSpecialize[7,3,[-1,-2,0]]{[-1,0,1]};

// Cancellation at the chosen split prime occurs in the coefficient beta^8.
rvqpadic K[5,-1,2];
let z = zero[RVQpQuad(K)];
let one = rvQpQuad[K,1,0,1];
let beta = rvQpQuad[K,-1,-2,1];
let top = add(one,rvQpQuadPow[8](beta));
let wall = rvResSpecializeQuad[K,4,[0,0,8,8]]{[one,z,one,one],[z,one,one,top]};
show wall;
rvResSupport(wall);
let beyond = rvResSpecializeQuad[K,4,[0,0,9,9]]{[one,z,one,one],[z,one,one,top]};
show beyond;
rvResSupport(beyond);
