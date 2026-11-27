// A degree-two point over F5, lifted through the unramified ring modulo 5^5.
rvunram W[5,5]{2,0,1};
let x = rvUnramElem[W]{0,1};
let z = rvLegendreParameter(x);
assert z == rvUnramElem[W]{0,2096};
assert rvUnramFrob[1](z) == rvUnramPow[5](z);
assert rvUnramPow[25](z) == z;
let a = rvLegendreFactor[0](x);
let b = rvLegendreFactor[1](x);
assert rvUnramFrob[1](a) == b;
assert rvUnramFrob[1](b) == a;
assert rvUnramMul(a,b) == rvLegendreUnit(x);
show a;
show b;
show rvLegendreUnit(x);
rvLegendreTrace(x);
rvLegendreHensel(x);

// Changing the polynomial lift changes coefficients in its named basis.
rvunram Other[5,5]{7,0,1};
let y = rvUnramElem[Other]{0,1};
let imageX = rvUnramElem[Other]{0,2756};
assert rvUnramMul(imageX,imageX) == rvUnramElem[Other]{-2};
assert add(rvUnramElem[Other]{3039},rvUnramMul(rvUnramElem[Other]{49},imageX))
    == rvLegendreFactor[0](y);
assert rvLegendreUnit(y) == rvUnramElem[Other]{2823};

// On a raw polynomial generator the lifted Frobenius differs from p-powering.
rvunram M[5,2]{2,-1,1};
let t = rvUnramElem[M]{0,1};
assert rvUnramFrob[1](t) == rvUnramElem[M]{1,-1};
assert rvUnramPow[5](t) == rvUnramElem[M]{6,-1};
assert rvUnramFrob[2](t) == t;
