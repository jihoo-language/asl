rvunram F[3,1]{0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
let two = rvUnramElem[F]{2};
rvhyperelliptic C[F]{z,o,z,z,z,o};
let center = zero[RVKummer(C)];

// The line through the center and (1:0:1:1) has residual point (1:0:1:2).
let x = rvKummerSmooth(rvKummer[C](o,z,o,o));
let y = rvKummerSmooth(rvKummer[C](o,z,o,two));
assert rvKummerNodeInvolution(center,x) == y;
assert rvKummerNodeInvolution(center,y) == x;
show y;

// The center's exceptional directions can be fixed or meet a smooth point.
let fixed = rvKummerExceptional(center,o,z,z,z);
assert rvKummerNodeInvolution(center,fixed) == fixed;
let direction = rvKummerExceptional(center,o,o,o,z);
let shifted = rvKummerExceptional(center,o,o,o,two);
let residual = rvKummerSmooth(rvKummer[C](o,o,o,z));
assert rvKummerNodeInvolution(center,direction) == residual;
assert rvKummerNodeInvolution(center,shifted) == residual;
assert rvKummerNodeInvolution(center,residual) == direction;
show residual;

// The other node stays fixed downstairs. Its exceptional conic carries
// the differential of the line involution.
let node = rvJacKummerPoint(rvJacPoint[C](z,z));
let h = rvKummerExceptional(node,o,z,o,o);
let shiftedH = rvKummerExceptional(node,o,two,o,o);
let image = rvKummerExceptional(node,o,z,o,two);
assert rvKummerNodeInvolution(center,h) == image;
assert rvKummerNodeInvolution(center,shiftedH) == image;
assert rvKummerNodeInvolution(center,image) == h;
let fixedH = rvKummerExceptional(node,o,z,z,z);
assert rvKummerNodeInvolution(center,fixedH) == fixedH;
show image;
