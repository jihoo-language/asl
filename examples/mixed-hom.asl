mixedcomplex V[5] { a:0, b:1, c:-1 }
    { b => [(a,1)] }
    { c => [(a,1)] };
mixhom End[V,V];
let a = mixchain[V,0][(a,1)];
let b = mixchain[V,1][(b,1)];
let c = mixchain[V,-1][(c,1)];
let f = mixunit[End,a,b];
let g = mixunit[End,a,a];
let e = mixId[End];
assert mixb(e) == zero[MixedChain(End,-1)];
assert mixB(e) == zero[MixedChain(End,1)];
assert mixb(f) == add(g,mixunit[End,b,b]);
assert mixB(f) == mixunit[End,c,b];
assert mixb(g) == pscale[4](mixunit[End,b,a]);
assert mixB(g) == pscale[4](mixunit[End,c,a]);
assert mixEval(f,a) == b;
assert mixEval(e,c) == c;
assert mixb(mixEval(f,a)) ==
    add(mixEval(mixb(f),a),pscale[4](mixEval(f,mixb(a))));
let shiftedF = mixU[2](f);
let shiftedB = mixU[1](b);
assert mixD(mixEval(shiftedF,shiftedB)) ==
    add(mixEval(mixD(shiftedF),shiftedB),pscale[4](mixEval(shiftedF,mixD(shiftedB))));
let bc = mixunit[End,b,c];
let ca = mixunit[End,c,a];
assert mixComp[End](ca,mixComp[End](bc,f)) == g;
assert mixComp[End](mixComp[End](ca,bc),f) == g;
assert mixComp[End](e,f) == f;
assert mixComp[End](f,e) == f;
assert mixEval(mixb(f),b) == mixEval(f,mixb(b));
assert mixEval(mixB(f),c) == mixEval(f,mixB(c));
assert mixComp[End](f,mixComp[End](g,e)) == mixComp[End](mixComp[End](f,g),e);
assert mixComp[End](f,g) == f;
assert mixb(mixComp[End](f,g)) ==
    add(mixComp[End](mixb(f),g),pscale[4](mixComp[End](f,mixb(g))));
assert mixB(mixComp[End](f,g)) ==
    add(mixComp[End](mixB(f),g),pscale[4](mixComp[End](f,mixB(g))));

mixedcomplex X[5] { x:0 } {} {};
mixedcomplex Y[5] { y:4, z:0 } {} {};
mixedcomplex Z[5] { q:8 } {} {};
mixhom XY[X,Y];
mixhom YZ[Y,Z];
mixhom XZ[X,Z];
mixhom XX[X,X];
let input = mixU[0](mixchain[X,0][(x,1)]);
let higher = mixU[2](mixunit[XY,x,y]);
let map = add(mixU[0](mixunit[XY,x,z]),higher);
let post = mixU[2](mixunit[YZ,y,q]);
let out = mixComp[XZ](post,map);
assert mixEval(map,input) ==
    add(mixU[0](mixchain[Y,0][(z,1)]),mixU[2](mixchain[Y,4][(y,1)]));
assert mixEval(out,input) == mixU[4](mixchain[Z,8][(q,1)]);
assert mixEval(out,input) == mixEval(post,mixEval(map,input));
assert mixComp[XY](map,mixU[0](mixId[XX])) == map;
assert mixD(out) == zero[MixedU(XZ,-1)];
assert mixD(mixU[0](mixId[XX])) == zero[MixedU(XX,-1)];
show mixEval(out,input);
