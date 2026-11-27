space P { x0 < x1 < x2 < x3 < x4 < x5, x0 < z < x5 };
sheaf F[P,5] {x0:1,x1:1,x2:1,x3:1,x4:1,x5:1,z:0}
    {x0<x1:[0],x1<x2:[0],x2<x3:[0],x3<x4:[0],x4<x5:[0],x0<z:[],z<x5:[]};

let a1 = add(shHomTerm[F,F,[x0,x1],0,0],
    add(shHomTerm[F,F,[x1,x2],0,0],
    add(shHomTerm[F,F,[x2,x3],0,0],
    add(shHomTerm[F,F,[x3,x4],0,0],shHomTerm[F,F,[x4,x5],0,0]))));
let j0 = shDefStart[F];
assert shDefTerm[0](j0) == zero[ShHomC(F,F,1)];
let j1 = shDefAppend(j0,a1);
assert shDefTerm[1](j1) == a1;
let j2 = shDefLift(j1);
let j3 = shDefLift(j2);
let j4 = shDefLift(j3);
assert shD(shDefTerm[2](j2)) == shScale[-1](shDefCurvature(j1));
assert shD(shDefTerm[3](j3)) == shScale[-1](shDefCurvature(j2));
assert shD(shDefTerm[4](j4)) == shScale[-1](shDefCurvature(j3));
assert shDefObstruction(j3) == zero[ShExt(F,F,2)];

let c = add(shHomTerm[F,F,[x0,x1,x5],0,0],
    add(shHomTerm[F,F,[x0,x2,x5],0,0],
    add(shHomTerm[F,F,[x0,x3,x5],0,0],shHomTerm[F,F,[x0,x4,x5],0,0])));
assert shDefCurvature(j4) == c;
assert shD(c) == zero[ShHomC(F,F,3)];
assert shDefObstruction(j4) == shClass(c);
show shDefObstruction(j4);
