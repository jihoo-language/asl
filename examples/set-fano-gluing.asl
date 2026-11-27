// Synchronized GL(3,F2) monodromy on the points and lines of the Fano plane.
space Base { x, y, a, b, c, x<a, x<b, x<c, y<a, y<b, y<c };
setfamily Points[Base] {x:7, y:7, a:7, b:7, c:7} {
    x<a:[0,1,2,3,4,5,6], x<b:[0,1,2,3,4,5,6], x<c:[0,1,2,3,4,5,6],
    y<a:[0,1,2,3,4,5,6], y<b:[1,3,5,2,0,6,4], y<c:[0,2,1,3,4,6,5]
};
setfamily Lines[Base] {x:7, y:7, a:7, b:7, c:7} {
    x<a:[0,1,2,3,4,5,6], x<b:[0,1,2,3,4,5,6], x<c:[0,1,2,3,4,5,6],
    y<a:[0,1,2,3,4,5,6], y<b:[2,3,6,0,1,4,5], y<c:[2,1,0,3,6,5,4]
};
space Closed {cy,cx};
space Open {oc,ob,oa};
let i = tab Closed -> Base {cy=>Base.y,cx=>Base.x};
let j = tab Open -> Base {oc=>Base.c,ob=>Base.b,oa=>Base.a};
setrecollement R(i,j);
setrestrict CP[R,Points,closed];
setrestrict CL[R,Lines,closed];
setrestrict OP[R,Points,open];
setrestrict OL[R,Lines,open];
setgluelimit TP[R,OP];
setgluelimit TL[R,OL];
let alphaP = setAttach[CP,TP];
let alphaL = setAttach[CL,TL];
show alphaP;
show alphaL;
setfamily Probe[Closed] {cy:0,cx:1} {};
card (\u:SetNat(Probe,TP) => S.hi);
card (\u:SetNat(Probe,TL) => S.hi);
setglue GP[R](alphaP);
setglue GL[R](alphaL);
assert setComp(setGlueForget[GP,Points],setGlueRecover[GP,Points]) == setId[Points];
assert setComp(setGlueRecover[GL,Lines],setGlueForget[GL,Lines]) == setId[GL];

setnat a[CP,CL] {cx:[0,1,2,3,4,5,6],cy:[0,1,2,3,4,5,6]};
setnat b[OP,OL] {oa:[0,1,2,3,4,5,6],ob:[0,1,2,3,4,5,6],oc:[0,1,2,3,4,5,6]};
setnat inverse[OL,OP] {oa:[0,1,2,3,4,5,6],ob:[0,1,2,3,4,5,6],oc:[0,1,2,3,4,5,6]};
assert setComp(setGlueLimitMap[TL,TP](inverse),setGlueLimitMap[TP,TL](b)) == setId[TP];

setran JP[OP](j);
setran JL[OL](j);
setrestrict OPJ[R,JP,open];
setrestrict OLJ[R,JL,open];
let compare = setRanTranspose[JL](setComp(b,setRanCounit[JP,OPJ]));
let back = setRanTranspose[JP](setComp(inverse,setRanCounit[JL,OLJ]));
assert setComp(back,compare) == setId[JP];
assert setComp(compare,back) == setId[JL];

setmonodromy M[GP,GL,x];
setmonostabilizer[M,0];
setmonomaps[M];
card (\m:SetNat(GP,GL) => S.hi);
let fromStrata = \m:SetNat(GP,GP) =>
    setGlueMap[GP,GP](setGlueClosed[GP,GP](m),setGlueOpen[GP,GP](m));
let identity = \m:SetNat(GP,GP) => m;
assert fromStrata == identity;
card (\m:SetNat(GP,GP) => S.hi);
