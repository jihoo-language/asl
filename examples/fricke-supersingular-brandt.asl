// Supersingular 2-isogeny edges modulo the Fricke involution.
rvunram F7[7,1]{-3,0,1};
rvfrickess B7[F7];
rvFrickeSSData[B7];
rvFrickeSSMatrix[B7,3];
rvFrickeSSMatrix[B7,5];
rvFrickeSSMatrix[B7,9];
rvFrickeSSRelations[B7];

let h0=rvFrickeSSVertex[B7,0];
let h1=rvFrickeSSVertex[B7,1];
assert h0==rvUnramElem[F7]{1};
assert h1==rvUnramElem[F7]{5};
assert rvFrickeSSEndpoint[B7,0,0]==rvUnramElem[F7]{6};
let d0=rvFrickeSSBasis[B7,0];
let d1=rvFrickeSSBasis[B7,1];
assert rvFrickeSSHecke[3](d0)==rvFrickeSSDiv[B7]{0,4};
assert rvFrickeSSHecke[3](d1)==rvFrickeSSDiv[B7]{2,2};
assert rvFrickeSSHecke[3](add(d0,d1))==add(rvFrickeSSHecke[3](d0),rvFrickeSSHecke[3](d1));
assert rvFrickeSSHecke[3](rvFrickeSSHecke[3](d0))==add(rvFrickeSSHecke[9](d0),rvFrickeSSDiv[B7]{4,0});
assert add(d0,neg(d0))==zero[RVFrickeSSDiv(B7)];

// A characteristic-zero horizontal CM loop specializes to H=1.
rvcmclasspoly D4[-4];
rvfrickecmedge L4[D4,D4];
assert rvFrickeCMPolyEval[L4](h0)==rvUnramElem[F7]{0};

rvunram F23[23,1]{-5,0,1};
rvfrickess B23[F23];
rvFrickeSSData[B23];
rvFrickeSSMatrix[B23,3];
rvFrickeSSRelations[B23];

rvunram F41[41,1]{-3,0,1};
rvfrickess B41[F41];
rvFrickeSSData[B41];
rvFrickeSSMatrix[B41,3];
rvFrickeSSRelations[B41];
assert rvFrickeSSVertex[B41,1]==rvUnramElem[F41]{20,7};
assert rvFrickeSSVertex[B41,2]==rvUnramElem[F41]{20,34};
let anti=rvFrickeSSDiv[B41]{0,1,-1,0,0};
assert rvFrickeSSFrob(anti)==neg(anti);
assert rvFrickeSSHecke[3](anti)==rvFrickeSSDiv[B41]{0,-2,2,0,0};
assert rvFrickeSSFrob(rvFrickeSSBasis[B41,1])==rvFrickeSSBasis[B41,2];
assert rvFrickeSSFrob(rvFrickeSSHecke[5](rvFrickeSSBasis[B41,1]))==rvFrickeSSHecke[5](rvFrickeSSBasis[B41,2]);
