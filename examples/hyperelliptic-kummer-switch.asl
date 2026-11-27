rvunram F[3,1]{0,1};
let z = rvUnramElem[F]{};
let o = rvUnramElem[F]{1};
let two = rvUnramElem[F]{2};

// For y²=x+x³+x⁴+x⁵ the gradient at (1:0:2:1) is (0:2:1:1).
rvhyperelliptic C[F]{z,o,z,o,o,o};
let smooth = rvKummerSmooth(rvKummer[C](o,z,two,o));
let other = rvKummerSmooth(rvKummer[C](o,two,two,z));
assert rvKummerSwitch(smooth) == other;
assert rvKummerSwitch(other) == smooth;
show rvKummerSwitch(smooth);

// A point on the theta plane over F9 lifts only on the twist.
rvunram W[3,1]{1,0,1};
let zw = rvUnramElem[W]{};
let ow = rvUnramElem[W]{1};
let i = rvUnramElem[W]{0,1};
rvhyperelliptic D[W]{zw,neg(ow),zw,zw,zw,ow};
let x = add(ow,i);
let square = rvUnramMul(x,x);
let point = rvKummer[D](zw,ow,x,square);
let exceptional = rvKummerExceptional(zero[RVKummer(D)],ow,neg(x),square,zw);
assert rvKummerLiftable(point) == S.lo;
assert rvKummerSwitch(rvKummerSmooth(point)) == exceptional;
assert rvKummerSwitch(exceptional) == rvKummerSmooth(point);
let intersection = rvKummerExceptional(zero[RVKummer(D)],zw,zw,ow,zw);
assert rvKummerSwitch(intersection) == intersection;
show exceptional;
