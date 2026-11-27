// Every possible first pivot in the first column has an unrepresentable inverse.
rvqpadic K[7,2,3];
let x = rvQpQuadPow[46](rvQpQuad[K,3,-1,1]);
let z = zero[RVQpQuad(K)];
let o = rvQpQuad[K,1,0,1];
let xx = rvQpQuadMul(rvQpQuad[K,2,0,1],x);
let xxx = rvQpQuadMul(rvQpQuad[K,3,0,1],x);
let M = rvValMatQuad[K,4]{[x,o,o,z],[xx,o,z,o],[xxx,z,o,o]};
assert M == rvValMatQuad[K,4]{[xxx,z,o,o],[x,o,o,z],[xx,o,z,o]};
rvValPlucker[[0,1,2]](M);
rvValPlucker[[0,1,3]](M);
rvValPlucker[[0,2,3]](M);
rvValPlucker[[1,2,3]](M);

rvqpadic Other[7,2,4];
let z2 = zero[RVQpQuad(Other)];
let o2 = rvQpQuad[Other,1,0,1];
let y = rvQpQuadPow[46](rvQpQuad[Other,3,-1,1]);
let yy = rvQpQuadMul(rvQpQuad[Other,2,0,1],y);
let yyy = rvQpQuadMul(rvQpQuad[Other,3,0,1],y);
let N = rvValMatQuad[Other,4]{[y,o2,o2,z2],[yy,o2,z2,o2],[yyy,z2,o2,o2]};
rvValPlucker[[0,1,2]](N);
rvValPlucker[[0,1,3]](N);
rvValPlucker[[0,2,3]](N);
rvValPlucker[[1,2,3]](N);

// All integer rows have large entries. The signed minors are 1,3,-2,1.
let R = rvValMat[2,4]{
    [29999999999999999999999998,30000000000000000000000001,1,90000000000000000000000003],
    [99999999999999999999999993,100000000000000000000000003,4,300000000000000000000000011],
    [69999999999999999999999995,70000000000000000000000002,4,210000000000000000000000011]
};
rvValPlucker[[0,1,2]](R);
rvValPlucker[[0,1,3]](R);
rvValPlucker[[0,2,3]](R);
rvValPlucker[[1,2,3]](R);
