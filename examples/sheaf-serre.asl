space A3 {a<b<c};
let atA = tab One -> A3 {pt=>A3.a};
let atB = tab One -> A3 {pt=>A3.b};
let atC = tab One -> A3 {pt=>A3.c};
sheaf F0[A3,3] {a:1,b:0,c:0} {a<b:[],b<c:[]};

// The simple at a moves to b, then c. Each move contributes one shift.
shserre K1[F0];
shserreshift[K1];
shcohomology F1[K1,1];
shpull F1a[atA,F1]; shdim[F1a,0];
shpull F1b[atB,F1]; shdim[F1b,0];
shpull F1c[atC,F1]; shdim[F1c,0];

shserre K2[F1];
shserreshift[K2];
shcohomology F2[K2,1];
shpull F2a[atA,F2]; shdim[F2a,0];
shpull F2b[atB,F2]; shdim[F2b,0];
shpull F2c[atC,F2]; shdim[F2c,0];

// The projective at c goes to the injective at c in unshifted degree.
shserre K3[F2];
shserreshift[K3];
shcohomology F3[K3,2];
shpull F3a[atA,F3]; shdim[F3a,0];
shpull F3b[atB,F3]; shdim[F3b,0];
shpull F3c[atC,F3]; shdim[F3c,0];

// The next step returns to the original simple. The total shift is two.
shserre K4[F3];
shserreshift[K4];
shcohomology F4[K4,2];
shpull F4a[atA,F4]; shdim[F4a,0];
shpull F4b[atB,F4]; shdim[F4b,0];
shpull F4c[atC,F4]; shdim[F4c,0];
