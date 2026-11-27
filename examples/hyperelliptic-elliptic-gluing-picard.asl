rvunram K[5,1]{0,1};
let z = rvUnramElem[K]{};
let o = rvUnramElem[K]{1};
let two = rvUnramElem[K]{2};
let three = rvUnramElem[K]{3};
let four = rvUnramElem[K]{4};
rvhyperelliptic E[K]{z,two,two,o};
rvhyperelliptic F[K]{z,three,four,o};

// The matching x ↦ (x+4)/(x+2) sends 0,1,2 to 2,0,4.
rvellglue C[E,F](o,four,o,two);
rvjacchart S[C](o);
rvJacCount[E];
rvJacCount[F];
rvJacCount[S];

let glue = \pq : RVJac(E) * RVJac(F) => rvEllGlueMap[S](fst(pq),snd(pq));
let dual = \d : RVJac(S) => rvEllGlueDual[S](d);
assert (\pq : RVJac(E) * RVJac(F) => dual(glue(pq)))
    == (\pq : RVJac(E) * RVJac(F) => (rvJacScale[2](fst(pq)),rvJacScale[2](snd(pq))));
assert (\d : RVJac(S) => glue(dual(d))) == (\d : RVJac(S) => rvJacScale[2](d));

// The kernel is the graph of the matching on two-torsion.
let identity = zero[RVJac(S)];
assert rvEllGlueMap[S](zero[RVJac(E)],zero[RVJac(F)]) == identity;
assert rvEllGlueMap[S](rvJacPoint[E](z,z),rvJacPoint[F](two,z)) == identity;
assert rvEllGlueMap[S](rvJacPoint[E](o,z),rvJacPoint[F](z,z)) == identity;
assert rvEllGlueMap[S](rvJacPoint[E](two,z),rvJacPoint[F](four,z)) == identity;

// The raw point (0,1) becomes (1,4) on S and is a pole of the right map.
let pole = rvJacPoint[S](o,four);
let leftNorm = add(rvEllGlueLeft[C](z,o),rvJacPoint[E](o,z));
assert rvEllGlueDual[S](pole) == (leftNorm,rvJacPoint[F](z,z));
