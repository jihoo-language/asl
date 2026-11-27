rvunram F[5,1]{0,1};
let z=rvUnramElem[F]{};
let o=rvUnramElem[F]{1};
rvhyperelliptic S[F]{z,o,z,z,z,o};
rvKummerK3ArtinTate[S,2];
rvKummerK3ArtinTate[S,10];
rvKummerK3Iwasawa[S,2,5,0];
rvKummerK3Iwasawa[S,2,5,1];
rvKummerK3Iwasawa[S,2,5,2];
rvKummerK3Iwasawa[S,2,5,1000000];

rvunram K[23,1]{0,1};
let kz=rvUnramElem[K]{};
let ko=rvUnramElem[K]{1};
let ktwo=rvUnramElem[K]{2};
let kfour=rvUnramElem[K]{4};
rvhyperelliptic C[K]{kz,kfour,neg(kfour),ko,neg(ktwo),ko};
rvKummerK3Iwasawa[C,2,23,1000000];

rvunram H[3,1]{0,1};
let hz=rvUnramElem[H]{};
let ho=rvUnramElem[H]{1};
rvhyperelliptic D[H]{hz,ho,ho,ho,ho,ho};
rvKummerK3Height[D];
rvKummerK3Iwasawa[D,1,3,1000000];
