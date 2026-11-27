// An exceptional point over the identity, then a tangent at an affine branch pair.
rvunram F[3,1]{0,1};
let z=rvUnramElem[F]{};
let o=rvUnramElem[F]{1};
let t=rvUnramElem[F]{2};
rvhyperelliptic C[F]{z,t,z,z,z,o};
let id=zero[RVKummer(C)];
let atIdentity=rvKummerExceptional(id,o,o,o,z);
assert atIdentity==rvKummerExceptional(id,t,t,t,o);
assert rvKummerBlowdown(atIdentity)==id;
assert rvKummerIsExceptional(atIdentity)==S.hi;
show atIdentity;

let divisor=add(rvJacPoint[C](z,z),rvJacPoint[C](o,z));
let node=rvJacKummerPoint(divisor);
let atNode=rvKummerExceptional(node,z,o,z,o);
assert rvKummerDirection[0](atNode)==z;
assert rvKummerDirection[3](atNode)==o;
assert rvKummerBlowdown(atNode)==node;
show atNode;
