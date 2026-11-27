// Two signs over F5 produce a rank-two system with nonabelian monodromy.
shdisk Sign[5,1,1] [1] [2];
shdiskseq Seed[Sign,Sign];
shdiskrigidity[Seed];
shconv MC[Seed,-1];
shconvrank[MC];
shdiskpoint Left[MC,0];
shdiskpoint Right[MC,1];
shdiskmat[Left,T];
shdiskmat[Right,T];
shdiskinfinity[MC];
shdiskgroup[MC];
shdiskrigidity[MC];

let raw=shConvRaw[MC]{[1],[-1]};
let v=shConvProject(raw);
assert shConvProject(shConvSection(v)) == v;
assert shConvRawT[-1](shConvRawT[1](raw)) == raw;
assert shConvRawT[2](shConvRawT[-2](raw)) == raw;
assert shConvProject(shConvRawT[1](raw)) == shConvT[1](v);
assert shConvProject(shConvRawT[-2](raw)) == shConvT[-2](v);
shconvcoords(shConvRawT[1](raw));
shconvcoords(shConvT[2](v));

shconv Back[MC,-1];
shconvrank[Back];
shdiskmonosame[Seed,Back];

// A different scalar and its inverse also change the nearby rank.
shconv Twist[Seed,2];
shconvrank[Twist];
shconv Untwist[Twist,3];
shconvrank[Untwist];
shdiskmonosame[Seed,Untwist];
