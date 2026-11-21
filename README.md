# ASL · Alexandrov Section Lambda

ASL is a pure functional language in which every type is a finite T₀ topological space and every function is continuous. Its two-point space `S` has exactly three endomorphisms. You can return the low point, preserve the input, or return the high point. Try to exchange the two points and the interpreter reports the order relation your function breaks.

## Run a program

From the repository root, run the triangle example with Nix.

```sh
nix run path:. -- examples/descent.asl
```

```text
Circle.bc
betti(Circle) = [1, 1]
betti(Circle * Circle) = [1, 2, 1]
homotopy = no
```

The program rotates the boundary of a triangle by gluing three local maps. Its third iterate is the identity. The Betti numbers detect the one-dimensional hole in the circle and the two independent loops and two-dimensional class in its product with itself. The final line comes from asking whether the rotation can deform to the identity *in the finite function space*. [Triangle program](examples/descent.asl)

The flake also supplies a Rust development shell. The interpreter uses the Rust standard library and reads either a file or standard input.

```sh
nix develop path:.
cargo test --offline
cargo run --offline -- examples/observables.asl
printf '%s\n' 'show (\x : S => x)(S.hi);' | cargo run --quiet --offline -- -
```

## Types, terms and opens

A `space` declaration gives generators for a partial order. For example, `space Fork { root < left, root < right };` has three points. An open set containing `root` must also contain `left` and `right`. This is the **upper Alexandrov topology**. A function is continuous exactly when `x ≤ y` implies `f(x) ≤ f(y)`.

| Type | Points and order |
| --- | --- |
| `One` | The single point `One.pt`. |
| `S` | The Sierpiński space with `S.lo < S.hi`. |
| `X * Y` | Pairs, ordered in both coordinates. |
| `X -> Y` | All monotone maps, ordered pointwise. |
| `Star(X, p)` | The principal open consisting of points `q` with `p ≤ q`, with the inherited order. |
| `H(X, k)` | The discrete space of degree-`k` cohomology classes of the order complex of `X`, with coefficients in 𝔽₂. |
| `A(d)` | The discrete space of homogeneous degree-`d` elements of the mod-two Steenrod algebra. |
| `Susp(X)` | The suspension poset, with two incomparable poles above every point of `X`. |
| `Stable(X, s, d)` | The discrete space of reduced degree-`d` classes of the `s`-fold suspension of `X₊`. Adjoin a disjoint basepoint to `X` to form `X₊`. |

Products associate to the left. Arrows associate to the right and bind more loosely than products. Named spaces have distinct types, even when their orders are isomorphic. A star also has its own type. For `q` in the star, its local point is written `Star(X, p).q`. A declaration can introduce isolated points, chains such as `p < q < r`, or the empty space `space Void {};`. A cycle of strict comparisons is an error.

`X -> S` is the space of open subsets of `X`. The points mapped to `S.hi` form the open, and the pointwise order on its characteristic maps is inclusion. This makes evaluation at a point a higher-order program.

```text
space Fork { root < left, root < right };
let leftOpen = tab Fork -> S { root => S.lo, left => S.hi, right => S.lo };
let eta = \x : Fork => \u : Fork -> S => u(x);
show eta(Fork.left)(leftOpen);   // S.hi
show eta(Fork.root)(leftOpen);   // S.lo
```

`let` binds an immutable value. A lambda captures its lexical environment, and `f(x)(y)` applies a curried function twice. Construct a product with `(x, y)` and project with `fst(pair)` or `snd(pair)`. A `tab` gives one result expression for each named point of its domain. The interpreter checks types, totality and monotonicity. Two functions of the same type are equal when they agree on every input, including when their inputs or outputs are themselves functions.

## Building a function from local data

```text
space Joint { west < seam, east < seam };
let w = tab Star(Joint, west) -> S { west => S.lo, seam => S.hi };
let e = \x : Star(Joint, east) => S.hi;
let lamp = sew Joint -> S { west => w, east => e };
assert res[west](lamp) == w;
show lamp(Joint.seam);   // S.hi
```

`res[p](f)` restricts a map to the points above `p`. In `sew X -> Y`, the chart at `p` must have type `Star(X, p) -> Y`. The chosen stars must cover `X`. For every point in two charts, both maps must give the same value. These conditions determine a unique global continuous map. `tab`, `res` and `sew` use point labels from a declared space, a built-in space or a star.

The double-cover example has a local section on each vertex star of the triangle. Two of the overlaps force the sheet choices. On the third overlap one chart yields `Cover.ca0` and the other yields `Cover.ca1`. Switching the initial sheet exchanges both answers. Run `cargo run -- examples/rejected/double-cover.asl` to see that obstruction reported at the `sew`. [Double-cover program](examples/rejected/double-cover.asl)

## Cup squares of cohomology classes

The real projective plane and the wedge of a circle and a 2-sphere both have mod-two Betti numbers `[1, 1, 1]`. Each has exactly one nonzero class in degree one. An ASL function can take that class as input and square it.

```text
let squareRP = \u : H(RP2, 1) => cup(u, u);
let squareWedge = \u : H(Wedge, 1) => cup(u, u);
show squareRP(gen[RP2, 1, 0]);       // H(RP2, 2){g0}
show squareWedge(gen[Wedge, 1, 0]);  // H(Wedge, 2){0}
```

Run `nix run path:. -- examples/projective-plane.asl` for the space declarations and the calculation. On the projective plane the square generates degree two. On the wedge the degree-one class is pulled back from the circle, where degree two vanishes. This distinguishes the two cohomology rings. [ASL program](examples/projective-plane.asl) · [Hatcher, *Algebraic Topology*, §3.2](https://pi.math.cornell.edu/~hatcher/AT/ATch3.pdf)

| Expression | Value |
| --- | --- |
| `zero[X, k]` | The additive zero in `H(X, k)`. Even a group of rank zero contains this point. |
| `gen[X, k, i]` | Generator `i` of a computed basis, indexed from zero. Basis coordinates are deterministic for a fixed presentation. |
| `cocycle[X, k] { a < b, ... }` | A class specified by its representative. List the strict `k`-simplices with coefficient one using point labels of `X`. Repeated entries cancel in pairs. |
| `add(u, v)` | Sum of classes on the same space in the same degree. |
| `cup(u, v)` | Cup product of classes on the same space. Their degrees add. |
| `pull(f, u)` | For `f` of type `X -> Y` and `u` of type `H(Y, k)`, the induced class in `H(X, k)`. |

Degrees and generator indices are nonnegative decimal literals. Any ASL type can serve as the base of `H`. A `cocycle` literal uses point labels, so it accepts declared spaces, built-in spaces and stars. On a product space, its projections let you pull classes back from either factor.

A `k`-simplex is a chain of `k + 1` points. A cocycle has even total coefficient on the boundary of every `(k + 1)`-simplex. A failed check identifies such a simplex. Closed representatives differing by a coboundary give equal ASL values. For example, in the hexagonal order complex of `Circle`, marking either of the two edges incident to `ab` gives the same degree-one class. Adding the coboundary of the 0-cochain supported at `ab` exchanges the two markings. [Torus program](examples/torus-cup.asl)

`cup` uses the Alexander–Whitney formula. On a chain `x₀ < … < xₚ₊ᵣ`, it multiplies the first cochain's value on `x₀ < … < xₚ` by the second's value on `xₚ < … < xₚ₊ᵣ`. `pull` evaluates a cocycle on the image chain, with coefficient zero when vertices collapse. Both operations return classes modulo coboundaries. A printed `g0 + g1` denotes the sum of two basis generators, and `0` denotes the zero class. [Hatcher, *Algebraic Topology*, §3.2](https://pi.math.cornell.edu/~hatcher/AT/ATch3.pdf)

## The mod-two Steenrod algebra

A function can take a cohomology operation as an argument and return another operation. `sq[i]` denotes the Steenrod square of degree `i`. In `mul(a, b)`, the right factor acts first. `sq[0]` is the multiplicative identity. `add(a, b)` adds operations of equal degree, and `zero[A(d)]` is the additive zero in that degree.

```text
let afterTwo = \theta : A(3) => mul(sq[2], theta);
show afterTwo(sq[3]);   // A(5){Sq^5 + Sq^4 Sq^1}
assert mul(sq[1], sq[1]) == zero[A(2)];
```

Equality in `A(d)` is equality of operations in the algebra. The interpreter reduces words to a sum of admissible monomials. In each surviving word, successive indices satisfy $i_j \ge 2i_{j+1}$. An inadmissible pair rewrites by the Adem relation, with binomial coefficients taken modulo two. [Hatcher, *Algebraic Topology*, §4.L](https://pi.math.cornell.edu/~Hatcher/AT/ATch4.4.pdf)

$$
\mathrm{Sq}^a\mathrm{Sq}^b
= \sum_{j=0}^{\lfloor a/2\rfloor}
  \binom{b-j-1}{a-2j}\mathrm{Sq}^{a+b-j}\mathrm{Sq}^j
\qquad (0<a<2b).
$$

Run `nix run path:. -- examples/adem.asl` to see the seven terms in the normal form of $\mathrm{Sq}^{128}\mathrm{Sq}^{128}$. The interpreter reduces this product directly, so expressions in `A(256)` remain usable even though that graded piece exceeds the 4,096-point enumeration limit. [Adem example](examples/adem.asl)

`act(theta, u)` evaluates an operation on a class. For `u` of degree $n$, $\mathrm{Sq}^0 u = u$, $\mathrm{Sq}^n u = u \smile u$, and $\mathrm{Sq}^i u = 0$ for $i>n$. For $0<i<n$ the interpreter takes a cocycle representative $c$, computes $c \smile_{n-i} c$ on ordered simplices and reduces the answer modulo coboundaries. A formally nonzero operation can therefore evaluate to zero on a given input. [Medina-Mardones, *New formulas for cup-i products and fast computation of Steenrod squares*](https://arxiv.org/abs/2105.08025)

## Suspension and stable degrees

`north[X]` and `south[X]` are the poles of `Susp(X)`. The inclusion `equator[X]` has type `X -> Susp(X)`. The order complex consists of two cones on the equator. `susp(u)` raises the degree of a cocycle by one. On a chain ending at the north pole it reads the value of the old cocycle on the remaining face, and it assigns zero to all other chains.

The suspensions of the real projective plane and of $S^1 \vee S^2$ both have Betti numbers `[1, 0, 1, 1]`. Every product of two positive-degree classes vanishes in both. Applying $\mathrm{Sq}^1$ to the degree-two generator gives the nonzero degree-three class for the suspended projective plane, and zero for the suspended wedge. The program computes the operation on both order complexes. [Projective-plane example](examples/projective-plane.asl)

For large shifts, `Stable(X, s, d)` uses the suspension isomorphism
$\widetilde H^d(\Sigma^s X_+;\mathbb F_2) \cong H^{d-s}(X;\mathbb F_2)$ for $d\ge s$, where $X_+$ adjoins a disjoint basepoint. Below degree $s$ the group is zero. `shift[s](u)` increases the shift and total degree together. `act` uses the underlying class to compute the action, since Steenrod squares commute with suspension. This representation keeps the cost of a shift independent of the number of cones it denotes. [Hatcher, *Algebraic Topology*, §4.L](https://pi.math.cornell.edu/~Hatcher/AT/ATch4.4.pdf)

```text
let alpha = gen[RP2, 1, 0];
let distant = \v : Stable(RP2, 1000000, 1000001) => act(sq[1], v);
show distant(shift[1000000](alpha));   // Stable(RP2, 1000000, 1000002){g0}
```

The basepoint convention matters in degree zero. For a nonempty connected space, `susp` sends the constant unit to zero. Shifting the unit by one gives a nonzero class in `Stable(X, 1, 1)`. Stable classes support addition, zero, further shifts and the Steenrod action.

## Observing a computation

Observations are top-level statements. Every statement ends with `;` and `//` starts a line comment.

| Statement | Result |
| --- | --- |
| `show e;` | Print a point, cohomology class, Steenrod operation, pair or function table. |
| `assert e == f;` | Check extensional equality at a common type. Success is silent. |
| `homotopy f ~ g;` | Find a shortest zigzag of pointwise-comparable maps of the same type. Print its length, or `no` if the maps lie in different components. |
| `betti T;` | Compute ordinary Betti numbers over 𝔽₂ of the order complex of `T`, starting in degree zero. Trailing zeros are omitted. The empty space yields `[]`. |

A homotopy step uses a single direction across the whole domain, so it requires `f(x) ≤ g(x)` for every `x`, or the reverse inequality for every `x`. In the triangle example, any endomorphism above the identity must fix all three maximal edges. Each vertex is the unique common lower bound of its two incident edges, so it too is fixed. The dual argument handles endomorphisms below the identity. This isolates the identity in `Circle -> Circle` and explains the `no` printed for the rotation. [Stong, *Finite Topological Spaces*, §3–4](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/stong2.pdf)

For `betti`, a strict chain of `k + 1` points contributes a `k`-simplex. Dropping each point in turn gives its boundary over 𝔽₂ in positive degrees. The ranks of successive boundary matrices determine the homology. [Barmak and Minian, *Simple Homotopy Types and Finite Spaces*, §2](https://math.uchicago.edu/~may/FINITE/ArgentinianPapers/MinianSimple.pdf)

The evaluator tabulates a lambda on its domain. For a discrete `D` with `n` points, an individual map `D -> S` needs `n` entries, while enumerating the whole function space yields `2^n` maps. A total order of `n` points has `2^n − 1` nonempty chains. The interpreter admits up to 512 points per named space, 4,096 elements per enumerated type and 50,000 simplices per order complex. A cohomology group of rank `b` has `2^b` elements. Individual classes can be constructed and multiplied even when enumerating that space would exceed the limit. Work and storage budgets can stop a calculation earlier with a diagnostic.
