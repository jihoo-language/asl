// A signed direct sum in a homogeneously changed basis.
mixedcomplex M[3] {
    x0:-8, x1:-6, x2:-4, x3:-2, x4:-7, x5:-5,
    x6:-3, x7:-1, x8:-2, x9:0, x10:-1, x11:1,
    x12:-3, x13:-1, x14:1, x15:-2, x16:0, x17:2,
    x18:2, x19:1, x20:0, x21:-1, x22:-2, x23:3,
    x24:2, x25:1, x26:0, x27:-1, x28:4, x29:3,
    x30:2, x31:1, x32:0, x33:5, x34:4, x35:3,
    x36:2, x37:1
} {
    x1 => [(x4,-1)],
    x2 => [(x5,-1)],
    x3 => [(x6,-1)],
    x9 => [(x10,-1)],
    x13 => [(x15,-1)],
    x14 => [(x16,-1)],
    x22 => [(x6,-1)],
    x23 => [(x18,1)],
    x24 => [(x19,1)],
    x25 => [(x20,1)],
    x26 => [(x21,1), (x7,-1), (x10,-1)],
    x27 => [(x22,1), (x3,-1)],
    x31 => [(x16,-1)],
    x33 => [(x28,1)],
    x34 => [(x29,1)],
    x35 => [(x30,1), (x17,-1)],
    x36 => [(x31,1), (x14,-1)],
    x37 => [(x32,1), (x16,-1)]
} {
    x0 => [(x4,1)],
    x1 => [(x5,1)],
    x2 => [(x6,1)],
    x3 => [(x7,1)],
    x8 => [(x10,1)],
    x12 => [(x15,1)],
    x13 => [(x16,1)],
    x14 => [(x17,1)],
    x19 => [(x18,1)],
    x21 => [(x20,1)],
    x22 => [(x7,1)],
    x24 => [(x23,-1)],
    x26 => [(x25,-1), (x11,1)],
    x29 => [(x28,1)],
    x31 => [(x30,1)],
    x34 => [(x33,-1)],
    x36 => [(x35,-1)],
    x37 => [(x18,1)]
};

assert mixD(mixBarRep[M,0]) == zero[MixedU(M,-3)];
assert mixD(mixBarRep[M,1]) == zero[MixedU(M,0)];
assert mixD(mixBarNull[M,2]) == mixU[4](mixBarRep[M,2]);
assert mixD(mixBarNull[M,3]) == mixU[3](mixBarRep[M,3]);
assert mixD(mixBarRep[M,2]) == zero[MixedU(M,-2)];
assert mixD(mixBarRep[M,3]) == zero[MixedU(M,1)];

let source = mixPageClass[4,0](mixBarNull[M,2]);
let target = mixPageClass[4,4](mixU[4](mixBarRep[M,2]));
assert mixD(source) == target;
assert mixNext(target) == zero[MixedE(M,5,-9,4)];
mixbardim[M,-2];
mixbarrank[M,2,-2];
mixbarrank[M,2,-4];
mixbarperiodicdim[M,0];
mixbars[M];
