mixedhigher X[5] {x0:199999999,s0:0} {s0 => [(100000000,x0,-1)]};
mixedhigher Y[5] {y:34,t:-199999965} {t => [(100000000,y,1)]};
mixhom F[X,Y];
let f = mixU[17](mixunit[F,s0,y]);
let fc = mixHClass(f);
assert fc == mixHgen[F,0,0];
assert mixHRep(fc) == f;
assert mixU[99999982](fc) == mixHgen[F,-199999964,0];
assert mixU[99999983](fc) == zero[MixedH(F,-199999966)];
let fNull = mixHPrimitive(mixU[99999983](f));
assert mixD(fNull) == mixU[99999983](f);

let hx = mixHgen[X,199999999,0];
assert mixInduce(f,hx) == zero[MixedH(Y,199999999)];
mixcone K(f);
mixcone Split(zero[MixedU(F,0)]);

let xRaw = mixU[0](mixConeLift[K](mixchain[X,199999999][(x0,1)]));
let x = mixHClass(xRaw);
assert x == mixHgen[K,200000000,0];
let y = mixHgen[Y,34,0];
let yImage = mixConeI[K](y);
let xAt34 = mixU[99999983](x);
let short = add(yImage,xAt34);
assert short == mixHgen[K,34,0];
assert xAt34 == mixHgen[K,34,1];
assert mixU[16](short) == mixHgen[K,2,0];
assert mixU[17](short) == zero[MixedH(K,0)];
assert mixU[199999982](x) == mixHgen[K,-199999964,0];
assert mixU[199999983](x) == zero[MixedH(K,-199999966)];
let s = mixU[0](mixConeLift[K](mixchain[X,0][(s0,1)]));
let t = mixU[0](mixConeI[K](mixchain[Y,-199999965][(t,1)]));
let longNull = add(mixU[99999983](s),pscale[4](t));
assert mixD(longNull) == mixU[199999983](xRaw);
assert mixD(mixHPrimitive(mixU[199999983](xRaw))) == mixU[199999983](xRaw);

assert mixConeP[K](yImage) == zero[MixedH(X,33)];
assert mixConeP[K](short) == mixHgen[X,33,0];
assert mixConeP[K](x) == hx;
assert mixConeConnect[K](hx) == mixInduce(f,hx);
assert mixConeConnect[K](mixConeP[K](short)) == zero[MixedH(Y,33)];
assert mixConeI[K](mixU[17](y)) == pscale[4](mixU[100000000](x));
assert mixConeP[K](mixU[100000000](x)) == zero[MixedH(X,-1)];

show fc;
mixHdim[K,34];
mixHdim[K,0];
mixbars[K];
mixbars[Split];
