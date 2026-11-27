rvhodgewick B[21]{0,1}{[2]};
rvVerbitskyHHDims[B,1,1,0];
rvVerbitskyHHDims[B,1,2,-3];
rvVerbitskyHHDims[B,1,4,-6];

let unit=rvVerbitskyHHBasis[B,1,0,0,0];
let x=rvVerbitskyHHBasis[B,1,0,1,0];
let deriv=rvVerbitskyHHClass(add(
    rvVerbitskyHHTerm[B,1,1,0,1]{1},
    rvVerbitskyHHTerm[B,1,1,0,2,2]{2}
));
assert rvVerbitskyHHBracket(deriv,x)==x;
assert rvVerbitskyHHBV(deriv)==rvVerbitskyHHScale[2,1](unit);
assert rvVerbitskyHHBV(rvVerbitskyHHCup(x,deriv))==x;

let carry=rvVerbitskyHHClass(add(
    add(rvVerbitskyHHTerm[B,1,2,-3,0]{1,2},rvVerbitskyHHTerm[B,1,2,-3,0]{2,1}),
    rvVerbitskyHHTerm[B,1,2,-3,1]{2,2}
));
assert rvVerbitskyHHBracket(deriv,carry)==rvVerbitskyHHScale[-3,1](carry);
let param=rvVerbitskyScale[-21,2](rvVerbitskyBasis[B,1,0,0]);
let tangent=rvVerbitskyTangent[B,1,-3]{(1,param)};
assert rvVerbitskyHHFromTangent(tangent)==carry;
rvVerbitskyHHCoeff[0]{1,2}(carry);
show rvVerbitskyHHCup(deriv,carry);
show rvVerbitskyHHCup(carry,carry);
