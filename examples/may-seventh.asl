// This E_7 generator has leading word [xi_2 | xi_2 | xi_2 | xi_2].
let source = maygen[7, 4, 12, 0, 0];
let target = mayClass[7, 7](cword{[1], [1], [1], [1], [8]});
assert target == maygen[7, 5, 12, 7, 0];
assert mayD(source) == target;
show mayD(source);
maydim[8, 4, 12, 0];
maydim[8, 5, 12, 7];
