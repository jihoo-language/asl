let xi = mayClass[1, 3](pword[3]{([1], [])});
let tau0 = mayClass[1, 0](pword[3]{([], [0])});
let tau1 = mayClass[1, 2](pword[3]{([], [1])});
assert mayD(tau1) == mayMul(xi, tau0);
show mayD(tau1);

// A corrected lift of [tau_1 | tau_1 | tau_1] has nonzero d_5.
let source = pmaygen[3, 5, 3, 15, 6, 0];
let target = mayClass[5, 11](pword[3]{([], [0]), ([], [0]), ([], [0]), ([3], [])});
assert mayD(source) == target;
show mayD(source);
pmaydim[3, 6, 3, 15, 6];
pmaydim[3, 6, 4, 15, 11];
