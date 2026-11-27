// The final types of the genus-twelve F3 curves with equal zeta factors.
let first = eoType{0,0,0,0,1,2,3,4,4,5,5,5};
let second = eoType{0,0,0,1,2,2,3,4,4,4,5,5};
eoDim(first);
eoDim(second);
show eoClosure(first,second);
show eoClosure(second,first);
eoClosureWitness(first,second);

// A codimension-one specialization in genus twelve beyond ordinary Bruhat.
let lower = eoType{1,2,3,4,5,6,7,7,7,8,9,9};
let upper = eoType{1,2,3,4,5,6,7,8,8,8,8,9};
eoDim(lower);
eoDim(upper);
show eoBruhat(lower,upper);
show eoClosure(lower,upper);
eoClosureWitness(lower,upper);
