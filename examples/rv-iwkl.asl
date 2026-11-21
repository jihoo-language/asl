rvIwKL[[0,1,2,3],[2,3,0,1]];
rvIwKL[[0,1,2],[-4,3,4]];
rvIwKL[[1,0,2],[0,2,1]];

qmodule X[0] { x0:0,x1:0,x2:0,x3:0 } {};
let i = rvU[0](rvId[X,[0],3,[1]]);
let s0 = rvU[0](add(add(rvgen[X,X,[0],3,[1],1],rvgen[X,X,[0],3,[1],4]),
                     add(rvgen[X,X,[0],3,[1],10],rvgen[X,X,[0],3,[1],15])));
let s1 = rvU[0](add(add(rvgen[X,X,[0],3,[1],0],rvgen[X,X,[0],3,[1],6]),
                     add(rvgen[X,X,[0],3,[1],9],rvgen[X,X,[0],3,[1],15])));
let s2 = rvU[0](add(add(rvgen[X,X,[0],3,[1],0],rvgen[X,X,[0],3,[1],5]),
                     add(rvgen[X,X,[0],3,[1],11],rvgen[X,X,[0],3,[1],14])));
let first = ecomp(s1,ecomp(s0,ecomp(s2,s1)));
let second = ecomp(s1,ecomp(s2,ecomp(s0,s1)));
assert first == second;
rvIw(first);
rvIwLength(second);
rvIwSchubertKL[[2,3,0,1]](i);
rvIwSchubertKL[[2,3,0,1]](s0);
rvIwSchubertKL[[2,3,0,1]](first);
