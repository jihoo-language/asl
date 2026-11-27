rvunram F[7,1]{0,1};
let zeroF = rvUnramElem[F]{};
let oneF = rvUnramElem[F]{1};
let twoF = rvUnramElem[F]{2};
let fourF = rvUnramElem[F]{4};
let fiveF = rvUnramElem[F]{5};
rvhyperelliptic C[F]{zeroF,fourF,twoF,fiveF,twoF,oneF};
rvrichelot D[C]{[],[1]}{[1],[],[1]}{[-3],[2],[1]};
rvRichelotData[D];
rvHWWPoints[C,1];
rvHWWPoints[D,1];
rvHWWWeil[C];
rvHWWWeil[D];

rvrichelotdual E[D];
rvRichelotData[E];
rvHWWWeil[E];
