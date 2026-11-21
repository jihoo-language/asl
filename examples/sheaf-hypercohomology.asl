// The suspension of the four-point circle has two points at each of three levels.
space Oct {a0 < b0, a0 < b1, a1 < b0, a1 < b1, b0 < c0, b0 < c1, b1 < c0, b1 < c1};
// The nonsplit extension has a degree-two transgression.
sheaf QTwisted[Oct,3] {a0:1, a1:1, b0:3, b1:3, c0:5, c1:5} {
    a0<b0:[1,0,0],
    a0<b1:[1,0,0],
    a1<b0:[1,0,0],
    a1<b1:[1,0,0],
    b0<c0:[1,1,0,2,0,0,0,0,1,1,0,2,0,0,0],
    b0<c1:[1,0,0,0,1,0,0,0,1,0,2,2,0,0,0],
    b1<c0:[1,0,0,2,1,0,0,0,1,1,0,0,0,2,2],
    b1<c1:[1,0,0,0,1,0,0,0,1,0,0,0,0,2,2]
};
sheaf PTwisted[Oct,3] {a0:1, a1:1, b0:3, b1:3, c0:5, c1:5} {
    a0<b0:[1,0,0],
    a0<b1:[1,0,0],
    a1<b0:[0,1,0],
    a1<b1:[0,1,0],
    b0<c0:[1,0,0,0,1,0,0,0,1,0,0,0,0,0,0],
    b0<c1:[1,0,0,0,1,0,0,0,1,0,0,0,0,0,0],
    b1<c0:[1,0,0,0,1,0,0,0,0,0,0,1,0,0,0],
    b1<c1:[1,0,0,0,1,0,0,0,0,0,0,1,0,0,0]
};
shmap deltaTwisted[QTwisted,PTwisted] {
    a0:[0],
    a1:[0],
    b0:[0,2,0,0,0,2,0,1,1],
    b1:[0,2,0,0,0,2,0,1,1],
    c0:[2,2,0,0,0,0,0,2,0,0,1,0,0,2,0,0,0,0,0,2,0,1,1,1,1],
    c1:[0,2,0,0,0,0,0,2,0,0,0,0,0,2,0,0,0,0,0,2,0,1,1,1,1]
};
shcomplex KTwisted {QTwisted,deltaTwisted,PTwisted};
shhyperdim[KTwisted,0];
shhyperdim[KTwisted,1];
shhyperdim[KTwisted,2];
shhyperdim[KTwisted,3];

// The split extension retains all four classes.
sheaf QSplit[Oct,3] {a0:1, a1:1, b0:3, b1:3, c0:5, c1:5} {
    a0<b0:[1,0,0],
    a0<b1:[1,0,0],
    a1<b0:[1,0,0],
    a1<b1:[1,0,0],
    b0<c0:[1,0,0,0,1,0,0,0,1,0,2,2,0,0,0],
    b0<c1:[1,0,0,0,1,0,0,0,1,0,2,2,0,0,0],
    b1<c0:[1,0,0,0,1,0,0,0,1,0,0,0,0,2,2],
    b1<c1:[1,0,0,0,1,0,0,0,1,0,0,0,0,2,2]
};
sheaf PSplit[Oct,3] {a0:1, a1:1, b0:3, b1:3, c0:5, c1:5} {
    a0<b0:[1,0,0],
    a0<b1:[1,0,0],
    a1<b0:[0,1,0],
    a1<b1:[0,1,0],
    b0<c0:[1,0,0,0,1,0,0,0,1,0,0,0,0,0,0],
    b0<c1:[1,0,0,0,1,0,0,0,1,0,0,0,0,0,0],
    b1<c0:[1,0,0,0,1,0,0,0,0,0,0,1,0,0,0],
    b1<c1:[1,0,0,0,1,0,0,0,0,0,0,1,0,0,0]
};
shmap deltaSplit[QSplit,PSplit] {
    a0:[0],
    a1:[0],
    b0:[0,2,0,0,0,2,0,1,1],
    b1:[0,2,0,0,0,2,0,1,1],
    c0:[0,2,0,0,0,0,0,2,0,0,0,0,0,2,0,0,0,0,0,2,0,1,1,1,1],
    c1:[0,2,0,0,0,0,0,2,0,0,0,0,0,2,0,0,0,0,0,2,0,1,1,1,1]
};
shcomplex KSplit {QSplit,deltaSplit,PSplit};
shhyperdim[KSplit,0];
shhyperdim[KSplit,1];
shhyperdim[KSplit,2];
shhyperdim[KSplit,3];

let top = shhypergen[KTwisted,3,0];
assert shClass(shRep(top)) == top;
assert shD(shRep(top)) == zero[ShTot(KTwisted,4)];
