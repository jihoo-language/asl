// The real tau_0 arrow becomes the differential of the localized vertex.
let tau0 = rhoLocalize(rmcword[0,0]{([],[1])});
let primitive = rhoHomotopy(tau0);
assert primitive == rholocword[1,1]{};
assert cobarD(primitive) == tau0;
show primitive;

// The real xi_1 arrow maps to classical xi_1 with half its original degree.
let h = rhoLocalize(rmcword[0,0]{([1],[])});
assert rhoProject(h) == cword{[1]};
assert h == rhoSection[2](cword{[1]});
show extClass(rhoProject(h));

// Add a genuine localized boundary to a classical Ext representative.
let g = extgen[1,2,0];
let section = rhoSection[-7](extRep(g));
let cycle = add(section,cobarD(rholocword[-7,2]{}));
assert cobarD(cycle) == zero[RhoLoc(2,-7,-9)];
assert extClass(rhoProject(cycle)) == g;
assert add(rhoSection[-7](rhoProject(cycle)),cobarD(rhoHomotopy(cycle))) == cycle;
assert rhoHomotopy(section) == zero[RhoLoc(0,-7,-9)];
show extClass(rhoProject(cycle));
