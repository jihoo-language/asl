// Trace the divisor of the point (2+ω, 2+2ω) over F_9 to F_3.
rvunram F[3,1]{0,1};
rvunram L[3,1]{1,0,1};
let z=rvUnramElem[F]{};
let o=rvUnramElem[F]{1};
let k=rvUnramElem[F]{2};
rvhyperelliptic C[F]{o,o,z,o,z,o};
let traced=rvJacQuadTrace[C](k){[1],[1]}{[2]}{[2]}{[2]};
rvjacbase D[C,L];
let p=rvJacPoint[D](rvUnramElem[L]{2,1},rvUnramElem[L]{2,2});
let norm=add(p,rvJacFrob[C,1](p));
assert rvJacExtend[D](traced)==norm;
assert traced==rvJacDescend[C](norm);
show traced;
