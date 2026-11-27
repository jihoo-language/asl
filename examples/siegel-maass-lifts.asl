// The weight-10 additive lift has seed Delta times phi_{-2,1}.
rvsiegelmaass Lift10[10]{minus(0,0,1,1)};
rvSiegelMaassJacobi[Lift10,1,0];
rvSiegelMaassCoeff[Lift10,4,4,0];
rvIgusaPhi10[4,4,0];
rvIgusaChi12[1,1,0];

// Both monomials have weight 16 and index 1.
rvsiegelmaass Lift16[16]{zero(1,0,1,3),minus(0,1,1,-2)};
rvSiegelMaassJacobi[Lift16,2,0];
rvSiegelMaassCoeff[Lift16,1,2,-1];
