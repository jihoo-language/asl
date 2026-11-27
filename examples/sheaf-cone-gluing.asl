// The link is an ordinal sum of four two-point antichains, a finite S^3.
space Link {a0, b0, a1, b1, a2, b2, a3, b3, a0<a1, a0<b1, b0<a1, b0<b1, a1<a2, a1<b2, b1<a2, b1<b2, a2<a3, a2<b3, b2<a3, b2<b3};
space Cone {u3, u2, tip, u4, u7, u5, u0, u6, u1, tip<u7, tip<u3, u7<u0, u7<u6, u3<u0, u3<u6, u0<u4, u0<u2, u6<u4, u6<u2, u4<u1, u4<u5, u2<u1, u2<u5};
let j = tab Link -> Cone {a0=>Cone.u7, b0=>Cone.u3, a1=>Cone.u0, b1=>Cone.u6, a2=>Cone.u4, b2=>Cone.u2, a3=>Cone.u1, b3=>Cone.u5};
let i = tab One -> Cone {pt=>Cone.tip};
sheaf F[Link,5] {a0:1, b0:1, a1:1, b1:1, a2:1, b2:1, a3:1, b3:1} {
    a0<a1:[1], a0<b1:[1], b0<a1:[1], b0<b1:[1], a1<a2:[1], a1<b2:[1], b1<a2:[1], b1<b2:[1], a2<a3:[1], a2<b3:[1], b2<a3:[1], b2<b3:[1]
};

let fundamental = shTerm[F,[a0,a1,a2,a3],0];
let h = shTerm[F,[a0,a1,a2],0,2];
shconeglue Split[j,tip,F,3] {zero[ShC(F,3)]};
shconeglue Attached[j,tip,F,3] {fundamental};
shconeglue Adjusted[j,tip,F,3] {add(fundamental,shD(h))};
shconehomotopy H[Attached,Adjusted] {h};

let a = shGlueApex[Attached,0];
let a2 = shConeTransfer[H](a);
assert shConeUntransfer[H](a2) == a;
assert shGlueTail[Adjusted](a2) == h;
assert shD(a2) == shConeTransfer[H](shD(a));
assert shGlueTail[Attached](shD(a)) == fundamental;
assert shD(shGlueLift[Attached](h)) == shGlueLift[Attached](shScale[-1](shD(h)));
let class = shhypergen[Attached__fiber,1,0];
assert shConeUntransfer[H](shConeTransfer[H](class)) == class;

shhyperdim[Split__fiber,1];
shhyperdim[Split__fiber,3];
shhyperdim[Split__fiber,4];
shhyperdim[Attached__fiber,1];
shhyperdim[Attached__fiber,3];
shhyperdim[Attached__fiber,4];

shextracomplex SplitAtTip[i,Split];
shextracomplex AttachedAtTip[i,Attached];
shlocaldim[SplitAtTip,pt,relative,1];
shlocaldim[SplitAtTip,pt,relative,3];
shlocaldim[SplitAtTip,pt,relative,4];
shlocaldim[AttachedAtTip,pt,relative,1];
shlocaldim[AttachedAtTip,pt,relative,3];
shlocaldim[AttachedAtTip,pt,relative,4];

// The stalk cohomology in degree three is the same skyscraper.
shcohomology SplitThree[Split,3];
shcohomology AttachedThree[Attached,3];
shmap Same[SplitThree,AttachedThree] {
    tip:[1],u0:[],u1:[],u2:[],u3:[],u4:[],u5:[],u6:[],u7:[]
};
shmap SameBack[AttachedThree,SplitThree] {
    tip:[1],u0:[],u1:[],u2:[],u3:[],u4:[],u5:[],u6:[],u7:[]
};
assert shExtCompose(shMapC[SameBack],shMapC[Same]) == shHomId[SplitThree];
assert shExtCompose(shMapC[Same],shMapC[SameBack]) == shHomId[AttachedThree];
