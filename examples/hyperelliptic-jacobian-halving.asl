rvunram F[3,1]{0,1};
let zeroF = rvUnramElem[F]{};
let oneF = rvUnramElem[F]{1};
rvhyperelliptic C[F]{zeroF,oneF,oneF,oneF,zeroF,oneF};

let point = rvJacPoint[C](oneF,oneF);
let half = rvJacHalfPoint[C](oneF,oneF){[2],[2],[],[1]};
let other = rvJacHalfPoint[C](oneF,oneF){[2],[2],[1],[2],[1]};
assert half == rvjac[C]{[2],[],[1]}{[2]};
assert other == rvjac[C]{[2],[2],[1]}{[1],[1]};
assert add(half,half) == point;
assert add(other,other) == point;

let branch = rvJacPoint[C](zeroF,zeroF);
assert add(half,branch) == other;
let quarter = rvJacHalfPoint[C](zeroF,zeroF){[],[1],[2],[],[1]};
assert quarter == rvjac[C]{[2],[],[1]}{[],[2]};
assert add(quarter,quarter) == branch;
show half;
show other;
show quarter;
