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

## Observing a computation

Observations are top-level statements. Every statement ends with `;` and `//` starts a line comment.

| Statement | Result |
| --- | --- |
| `show e;` | Print a point, pair or function table. |
| `assert e == f;` | Check extensional equality at a common type. Success is silent. |
| `homotopy f ~ g;` | Find a shortest zigzag of pointwise-comparable maps of the same type. Print its length, or `no` if the maps lie in different components. |
| `betti T;` | Compute ordinary Betti numbers over 𝔽₂ of the order complex of `T`, starting in degree zero. Trailing zeros are omitted. The empty space yields `[]`. |

A homotopy step uses a single direction across the whole domain, so it requires `f(x) ≤ g(x)` for every `x`, or the reverse inequality for every `x`. In the triangle example, any endomorphism above the identity must fix all three maximal edges. Each vertex is the unique common lower bound of its two incident edges, so it too is fixed. The dual argument handles endomorphisms below the identity. This isolates the identity in `Circle -> Circle` and explains the `no` printed for the rotation. [Stong, *Finite Topological Spaces*, §3–4](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/stong2.pdf)

For `betti`, a strict chain of `k + 1` points contributes a `k`-simplex. Dropping each point in turn gives its boundary over 𝔽₂ in positive degrees. The ranks of successive boundary matrices determine the homology. The projective-plane example presents six vertices, fifteen edges and ten faces, and evaluates to `[1, 1, 1]`. [Projective-plane program](examples/projective-plane.asl) · [Barmak and Minian, *Simple Homotopy Types and Finite Spaces*, §2](https://math.uchicago.edu/~may/FINITE/ArgentinianPapers/MinianSimple.pdf)

The evaluator tabulates a lambda on its domain. For a discrete `D` with `n` points, an individual map `D -> S` needs `n` entries, while enumerating the whole function space yields `2^n` maps. A total order of `n` points has `2^n − 1` nonempty chains. The interpreter admits up to 512 points per named space, 4,096 elements per enumerated type and 50,000 simplices per order complex. Work and storage budgets can stop a calculation earlier with a diagnostic.
