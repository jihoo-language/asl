mixedcomplex Left[5] { x0:1, x1:3, y1:2, y2:4 }
    { x1 => [(y1,-1)] }
    { x0 => [(y1,1)], x1 => [(y2,1)] };
mixedcomplex Right[5] { x0:1, x1:3, x2:5, y1:2, y2:4, y3:6 }
    { x1 => [(y1,-1)], x2 => [(y2,-1)] }
    { x0 => [(y1,1)], x1 => [(y2,1)], x2 => [(y3,1)] };
mixtensor T[Left,Right];
mixtensor S[Right,Left];

let A = mixBarRep[Left,0];
let B = mixBarRep[Right,0];
let X = mixBarNull[Left,0];
let Y = mixBarNull[Right,0];
let primary = mixTensorCycle[T,0,0];
let secondary = mixTensorTor[T,0,0];
assert primary == mixTensor[T](A,B);
assert secondary == add(mixTensor[T](A,Y),pscale[4](mixU[1](mixTensor[T](X,B))));
assert mixTensorTorNull[T,0,0] == mixTensor[T](X,Y);
assert mixD(primary) == zero[MixedU(T,9)];
assert mixD(secondary) == zero[MixedU(T,4)];
assert mixD(mixTensorNull[T,0,0]) == mixU[2](primary);
assert mixD(mixTensorTorNull[T,0,0]) == mixU[2](secondary);
assert mixPageClass[3,1](mixU[1](primary)) == mixgen[T,3,8,1,0];
assert mixPageClass[3,1](mixU[1](secondary)) == mixgen[T,3,3,1,0];
assert mixPageClass[3,2](mixU[2](secondary)) == zero[MixedE(T,3,1,2)];
assert mixSwap[S](primary) == mixTensorCycle[S,0,0];
assert mixSwap[S](secondary) == pscale[4](mixTensorTor[S,0,0]);
assert mixSwap[T](mixSwap[S](secondary)) == secondary;
assert mixD(mixSwap[S](mixTensor[T](X,Y))) == mixSwap[S](mixD(mixTensor[T](X,Y)));

let x = mixchain[Left,1][(x0,1)];
let y = mixchain[Right,1][(x0,1)];
let raw = mixTensor[T](x,y);
assert mixB(raw) == add(mixTensor[T](mixB(x),y),pscale[4](mixTensor[T](x,mixB(y))));
assert mixSwap[S](raw) == pscale[4](mixTensor[S](y,x));
assert mixSwap[S](mixB(raw)) == mixB(mixSwap[S](raw));
assert mixSwap[S](mixb(raw)) == mixb(mixSwap[S](raw));
assert mixTensor[T](mixU[2](x),mixU[3](y)) == mixU[5](raw);
mixbarrank[T,10,8];
mixbarrank[T,10,6];
mixbarrank[T,5,3];
mixbarrank[T,5,1];
mixbars[T];
