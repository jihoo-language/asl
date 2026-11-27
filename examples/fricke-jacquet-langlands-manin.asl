// Work on the actual rational modular symbols at the two levels.
rvmanin M23[23];
rvManinData[M23];
rvManinSpectrum[M23,3,2];
rvManinSpectrum[M23,3,3];

let edge=rvManinPath[M23]{[1,0],[0,1]};
assert edge==rvManinSymbol[M23,1,0];
assert rvManinSymbol[M23,0,1]==neg(edge);
assert rvManinSymbol[M23,1,1]==rvManinSymbol[M23,3,3];
assert add(add(rvManinSymbol[M23,1,1],rvManinSymbol[M23,1,-2]),rvManinSymbol[M23,-2,1])==zero[RVManin(M23)];
rvManinBoundary(edge);
let path=rvManinPath[M23]{[1,2],[3,23]};
assert path==add(rvManinPath[M23]{[1,2],[1,0]},rvManinPath[M23]{[1,0],[3,23]});
rvManinBoundary(path);

rvunram F23[23,1]{-5,0,1};
rvfrickess B23[F23];
rvfrickejl A[M23,B23];
rvFrickeJLData[A];
let s=rvManinBasis[M23,3,0];
let d=rvFrickeJLToBrandt[A](s);
assert rvFrickeJLToManin[A](d)==s;
assert rvManinJ(s)==s;
assert rvManinAtkin[2](s)==s;
assert rvFrickeJLToBrandt[A](rvManinHecke[3](s))==rvFrickeSSHecke[3](d);
assert rvFrickeJLToBrandt[A](rvManinHecke[5](s))==rvFrickeSSHecke[5](d);
assert rvFrickeJLToBrandt[A](neg(rvManinAtkin[23](s)))==rvFrickeSSFrob(d);
assert rvManinHecke[3](rvManinHecke[3](s))==add(rvManinHecke[9](s),rvManinScale[3](s));

rvmanin M41[41];
rvManinData[M41];
rvManinSpectrum[M41,3,2];
rvunram F41[41,1]{-3,0,1};
rvfrickess B41[F41];
rvfrickejl B[M41,B41];
rvFrickeJLData[B];
let difference=rvFrickeSSDiv[B41]{0,1,-1,0,0};
let anti=rvFrickeJLToManin[B](difference);
assert rvFrickeJLToBrandt[B](anti)==rvFrickeJLToBrandt[B](rvManinAtkin[41](anti));
assert rvFrickeJLToBrandt[B](anti)==neg(rvFrickeSSFrob(rvFrickeJLToBrandt[B](anti)));
assert rvFrickeOrderIntegral(rvFrickeJLToBrandt[B](anti))==difference;

// A congruent lift gives the same path, and the class supports rational scale.
assert edge==rvManinPath[M23]{[1,46],[0,1]};
let half=rvManinScale[1,2](edge);
assert add(half,half)==edge;
