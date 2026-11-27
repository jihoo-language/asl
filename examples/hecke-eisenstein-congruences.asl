rvheckepacket Delta(rvHeckeT[12,2]);
rvheckeresidual R[Delta,691,-24];
rvHeckeEisensteinSturm[R];
assert rvHeckeResidualCoeff[R,7] == rvHeckeEisensteinCoeff[12,691,7];

rvheckepacket Weight16(rvHeckeT[16,2]);
rvheckeresidual S16[Weight16,3617,216];
rvHeckeEisensteinSturm[S16];

rvheckepacket W(rvHeckeT[24,2]);
rvheckeresidual Good[W,103,83];
rvheckeresidual Other[W,103,70];
rvheckeresidual Large[W,2294797,1504218];
rvHeckeEisensteinSturm[Good];
rvHeckeEisensteinSturm[Other];
rvHeckeEisensteinSturm[Large];
assert rvHeckeResidualReduce[Good](rvHeckeEigenGenerator[W]) == rvHeckeResidualCoeff[Good,2];
show rvHeckeResidualCoeff[Good,3];
