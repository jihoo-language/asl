space Fork {root<u,root<v};
setfamily F[Fork] {root:0,u:2,v:3} {root<u:[],root<v:[]};
setj J[F];
setjcheck[F];
setjcheck[J];
setsections[F];
setsections[J];

setfamily G[Fork] {root:6,u:2,v:3}
  {root<u:[1,1,1,0,0,0],root<v:[2,1,0,2,1,0]};
setjcheck[G];
setnat local[F,G] {root:[],u:[1,0],v:[2,0,1]};
let extension = setJExtend[J](local);
assert setComp(extension,setJUnit[J]) == local;
assert setJExtend[J](setComp(extension,setJUnit[J])) == extension;
show extension;

setj JJ[J];
assert setComp(setJInverse[JJ],setJUnit[JJ]) == setId[J];
assert setComp(setJUnit[JJ],setJInverse[JJ]) == setId[JJ];

setomega O[Fork];
setsieve[O,root,3];
setdense[O,root,3];
setdense[O,root,1];
