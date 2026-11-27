// At SU(3) level two the WZW root has order 60.
let vacuum = rvWzwBasis[3,2,[0,0,0]];
let adjoint = rvWzwBasis[3,2,[0,1,2]];
let fundamental = rvWzwBasis[3,2,[0,0,1]];
assert rvWzwGalois[7](vacuum) == adjoint;
assert rvWzwGalois[7](adjoint) == neg(vacuum);
assert rvWzwGalois[49](vacuum) == neg(vacuum);
show rvWzwGalois[7](vacuum);
show rvWzwGalois[7](adjoint);

// Compare an exact Galois-conjugate S entry with the signed vacuum row.
let entry = rvWzwSEntry[3,2,[0,1,2],[0,0,1]];
assert rvCycloGalois[7](entry) == neg(rvWzwSEntry[3,2,[0,0,0],[0,0,1]]);

// Coefficientwise transport and the signed label action are separate.
let v = add(rvWzwScale(rvCyclo[60,2],vacuum),
    add(rvWzwScale(rvCyclo[60,3],adjoint),fundamental));
let twisted = rvWzwGalois[7](v);
let semilinear = rvWzwGalois[7](rvWzwCoeffGalois[7](v));
assert rvWzwCoeffGalois[7](rvWzwS(v)) == rvWzwS(semilinear);
assert rvWzwCoord[[0,1,2]](twisted) == rvCyclo[60,2];
assert rvWzwCoord[[0,1,2]](semilinear) == rvCyclo[60,14];
assert rvWzwGalois[7](twisted) == neg(v);
assert rvWzwGalois[-1](v) == neg(rvWzwCharge(v));

// The twist changes by the square of the power and a vacuum phase.
let phase = rvWzwTEntry[3,2,[0,1,2]];
assert rvWzwT[1](twisted) == rvWzwScale(phase,rvWzwGalois[7](rvWzwT[49](v)));
show rvWzwCoord[[0,1,2]](semilinear);
