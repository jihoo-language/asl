// One geometric theta graph, two arithmetic Frobenius actions.
let split = rvHyperArithmetic[7,1]{0,49,1,344,2,16809};
let nonsplit = rvHyperArithmetic[7,3]{0,49,1,344,2,16809};
assert rvHyperGeometric(split) == rvHyperGeometric(nonsplit);
rvHyperMonodromy(rvHyperGeometric(split));
rvHyperComponent(rvHyperGeometric(split));
rvHyperFrobenius(split);
rvHyperFrobenius(nonsplit);
rvHyperFixed(split);
rvHyperFixed(nonsplit);

// The two even tails have opposite Frobenius signs.
let mixed = rvHyperArithmetic[7,1]{0,1,2,345,3,16810};
let twist = rvHyperArithmetic[7,3]{0,1,2,345,3,16810};
assert rvHyperGeometric(mixed) == rvHyperGeometric(twist);
rvHyperMonodromy(rvHyperGeometric(mixed));
rvHyperComponent(rvHyperGeometric(mixed));
rvHyperFrobenius(mixed);
rvHyperFrobenius(twist);
rvHyperFixed(mixed);
rvHyperFixed(twist);
