// On s=1 this cover is y^2+y=t^3, with one pole of order three.
rvascover Supersingular[[0,1]] {rvmon[[0,1]]{0,3},rvmon[[0,1]]{3,0}};
rvASGenus[Supersingular];
rvASPrank[Supersingular];
rvASLPolynomial[Supersingular];
rvASNewton[Supersingular];
rvASPoints[Supersingular,4];

// y^2+y=t+1/t has two simple rational poles.
rvascover Ordinary[[0,1]] {
    add(rvmon[[0,1]]{2,0},rvmon[[0,1]]{0,2}), rvmon[[0,1]]{1,1}
};
rvASGenus[Ordinary];
rvASPrank[Ordinary];
rvASLPolynomial[Ordinary];
rvASNewton[Ordinary];
rvASPoints[Ordinary,5];

rvascover GenusTwo[[0,1]] {rvmon[[0,1]]{0,5},rvmon[[0,1]]{5,0}};
rvASGenus[GenusTwo];
rvASLPolynomial[GenusTwo];
rvASNewton[GenusTwo];
rvASPoints[GenusTwo,4];

// The roots of t^2+t+1 are two poles over the algebraic closure.
rvascover QuadraticPoles[[0,1]] {
    rvmon[[0,1]]{2,0},
    add(rvmon[[0,1]]{2,0},add(rvmon[[0,1]]{1,1},rvmon[[0,1]]{0,2}))
};
rvASGenus[QuadraticPoles];
rvASPrank[QuadraticPoles];
rvASLPolynomial[QuadraticPoles];
rvASNewton[QuadraticPoles];
rvASPoints[QuadraticPoles,2];
