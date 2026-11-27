// Frobenius on the scalar and the first two nonzero Bockstein arrows.
let t = rhoClass[1,0](rmcword[0,1]{});
let t2 = rhoFrob(t);
let t4 = rhoFrob(t2);
assert t2 == rhoClass[2,0](rmcword[0,2]{});
assert t4 == rhoClass[4,0](rmcword[0,4]{});
assert rhoD(t2) == rhoFrob(rhoD(t));
assert rhoD(t4) == rhoFrob(rhoD(t2));
assert rhoD(t4) == rhobockgen[4,1,0,-4,4,0];
show rhoD(t4);

// The square of the second arrow uses the first arrow's right-unit transport.
let arrow = rmcword[0,0]{([1],[]),([],[1])};
let square = add(rmcword[0,1]{([2],[]),([1],[])},
    add(rmcword[1,0]{([2],[1]),([1],[])},
    add(rmcword[1,0]{([2],[]),([1],[1])},rmcword[1,0]{([2],[]),([],[0,1])})));
assert rhoFrob(arrow) == square;

let c = rhoReduce[1](rmcword[0,0]{([],[1]),([],[1])});
let fc = rhoFrob(c);
let want = add(rmcword[0,2]{([1],[]),([1],[])},
    add(rmcword[1,1]{([1],[]),([1],[1])},
    add(rmcword[1,1]{([1],[]),([],[0,1])},rmcword[1,1]{([],[0,1]),([1],[])})));
assert fc == rhoReduce[2](want);
let u = rhoFrobConnectNull(c);
let change = add(rhoFrob(rhoConnectRep(c)),rhoConnectRep(fc));
assert cobarD(u) == change;
assert change == rhoConnectRep(fc);
show extClass(change);

// Canonical cohomology representatives add a second, exact correction.
let x = extClass(c);
let fx = rhoFrob(x);
assert cobarD(rhoFrobConnectNull(x)) == add(rhoFrob(rhoConnectRep(x)),rhoConnectRep(fx));
assert rhoConnect(fx) == rhoFrob(rhoConnect(x));

let a = extClass(rhoReduce[32](rmcword[0,32]{}));
let h = rhoConnect(rhoFrob(a));
assert h == rhoFrob(rhoConnect(a));
assert rhoNull[64](h) == rmcword[0,64]{};
rhoorder(h);
