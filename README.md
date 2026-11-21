# ASL · Alexandrov Section Lambda

ASL is a pure functional language in which every type is a finite T₀ topological space and every function is continuous. Its two-point space `S` has exactly three endomorphisms. You can return the low point, preserve the input, or return the high point. Try to exchange the two points and the interpreter reports the order relation your function breaks.

The classifying spaces of the cyclic groups of orders 4 and 4,096 have isomorphic mod-two cohomology rings. ASL can distinguish them by how many copies of the degree-one generator it takes to produce a nonzero Massey product. [Run the higher-product experiment](#how-late-can-a-higher-product-appear)

Even four zero inputs can have 344 different answers. ASL computes the entire set of outcomes as a function that tests membership. [Run the complete Massey product](#the-complete-massey-product)

ASL computes a nonzero $d_7$ in internal degree 12 of the algebraic May spectral sequence and the transferred product $m_4(h_0,h_1,h_2,h_0)=h_2^2$ on the Adams $E_2$ page. [May differential](examples/may-seventh.asl) · [Transferred product and cobar lift](examples/ext-transfer.asl)

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
| `Cobar(s, t)` | The discrete space of normalized cobar cochains of length `s` and internal degree `t` for the dual mod-two Steenrod algebra. |
| `Ext(s, t)` | The discrete space of mod-two Adams $E_2$ classes in filtration `s` and internal degree `t`. |
| `XExt(X, s, t)` | The discrete space of mod-two Adams $E_2$ classes for the suspension spectrum of `X` with a disjoint basepoint. |
| `MapExt(X, Y, s, t)` | The discrete space of mod-two Adams $E_2$ classes for stable maps from `X` to `Y`, with disjoint basepoints. |
| `May(r, s, t, f)` | The discrete space of classes on the algebraic May page `r` in cobar length `s`, internal degree `t`, and filtration `f`. |
| `Susp(X)` | The suspension poset, with two incomparable poles above every point of `X`. |
| `Stable(X, s, d)` | The discrete space of reduced degree-`d` classes of the `s`-fold suspension of `X₊`. Adjoin a disjoint basepoint to `X` to form `X₊`. |
| `Bar(m)` | The discrete space of normalized mod-two 1-cochains of the cyclic group of order `m`. Its points are the functions on residues that vanish at zero. |
| `H2(m)` | The discrete space of classes in $H^2(BC_m;\mathbb F_2)$. It has two points for even `m` and one for odd `m`. |

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

## The normalized cobar algebra

The dual mod-two Steenrod algebra is $\mathbb F_2[\xi_1,\xi_2,\ldots]$ with $|\xi_i|=2^i-1$. In a literal such as `cword{[2], [0, 1]}`, each bracket lists the exponents of one monomial in order from $\xi_1$ onward. This word is $[\xi_1^2\mid\xi_2]$ in `Cobar(2, 5)`. Each factor has positive degree, and trailing zero exponents may be omitted. The empty word `cword{}` is the tensor unit in `Cobar(0, 0)`.

`add(a, b)` adds cochains of equal bidegree over $\mathbb F_2$, and `zero[Cobar(s, t)]` gives the zero cochain. `concat(a, b)` concatenates words bilinearly, adding their lengths and internal degrees. `cobarD(a)` replaces each factor in turn by its reduced coproduct and sums the results. It increases the length by one and preserves internal degree. The coproduct is multiplicative and uses $\Delta(\xi_n)=\sum_{i=0}^n \xi_{n-i}^{2^i}\otimes\xi_i$ with $\xi_0=1$. The reduced coproduct omits the two terms with a unit factor.

```text
let xi2 = cword{[0, 1]};
let d = \x : Cobar(1, 3) => cobarD(x);
assert d(xi2) == cword{[2], [1]};
assert cobarD(d(xi2)) == zero[Cobar(3, 3)];
```

Sparse operations work directly on words. A lambda over `Cobar(s, t)` tabulates all cochains of that bidegree and is subject to the finite enumeration limit. Run `nix run path:. -- examples/cobar.asl` for coproduct, unit, and derivation checks. [Cobar example](examples/cobar.asl)

## The Adams E₂ page

`Ext(s, t)` is $H^s(\mathrm{Cobar}(*,t),d)$, the classical mod-two Adams $E_2$ page of the sphere. The corresponding stem is $t-s$. ASL uses sparse elimination over $\mathbb F_2$ on cobar differentials. `extdim[s, t];` prints the rank, and `extgen[s, t, i]` selects the $i$th generator in the computed basis, starting from zero.

`extClass(c)` takes a closed cobar cochain in bidegree `(s, t)` to its class in `Ext(s, t)`. `extRep(x)` returns a deterministic closed representative of `x`. For a cobar boundary of positive length, `extPrimitive(c)` returns a deterministic `p` with `cobarD(p) == c`. `yoneda(x, y)` multiplies classes by concatenating their representatives and taking the result's class. Classes of the same bidegree support `add` and `zero[Ext(s, t)]`.

Run `nix run path:. -- examples/adams-e2.asl` to compute $h_0h_1=0$ and $h_1^3=h_0^2h_2\ne0$. [Adams E₂ example](examples/adams-e2.asl)

## The algebraic May spectral sequence

Give a generator power $\xi_i^{2^j}$ May weight $2i-1$. For an arbitrary monomial, count the set bits of each exponent and multiply that count by the corresponding generator weight. The weight $M$ of a cobar word is the sum over its factors. In internal degree $t$, the decreasing filtration $F^f$ is spanned by words with $t-M\ge f$. The cobar differential preserves this filtration, and concatenation adds filtration indices. The associated graded complex is page $E_0$. Its differential preserves $f$, and $d_r$ increases $s$ by one and $f$ by $r$, preserving $t$. [May, *The cohomology of restricted Lie algebras and of Hopf algebras*](https://doi.org/10.1016/0021-8693(66)90009-3)

| Expression or statement | Result |
| --- | --- |
| `maydim[r, s, t, f];` | Prints the rank of the page component. |
| `maygen[r, s, t, f, i]` | Returns generator `i` of a deterministic basis, indexed from zero. |
| `mayClass[r, f](c)` | Projects a cobar cochain in $F^f$ whose differential lies in $F^{f+r}$. |
| `mayRep(x)` | Returns a deterministic cobar lift of the page class. |
| `mayD(x)` | Evaluates $d_r$ on the page class. |
| `mayMul(x, y)` | Multiplies two classes on the same page, adding their cobar lengths, internal degrees, and filtration indices. |

`zero[May(r, s, t, f)]` and `add(x, y)` use the vector space structure in a fixed page component. For $r\ge1$, projection identifies lifts that differ by an $(r-1)$-cycle in $F^{f+1}$ or a cobar boundary in $F^f$ with a primitive in $F^{\max(0,f-r+1)}$. The computation stabilizes separately in each internal degree. The stabilized groups are the associated graded groups of `Ext(s, t)` for the induced filtration.

Run `nix run path:. -- examples/may.asl` to compute $d_1(h_{2,0})=h_1h_0$ and the nonzero $d_3$ of a corrected $b_{2,0}$ lift in internal degree 6. The page ranks record the effect on the source and target. A second program computes a nonzero $d_7$ in internal degree 12 and checks that both page components vanish on $E_8$. [First and third differentials](examples/may.asl) · [Seventh differential](examples/may-seventh.asl)

## Transferred homotopy operations on Ext

`ainfty(a1, ..., an)` gives the operation $m_n$ of the transferred minimal $A_\infty$ algebra on `Ext`. `alift(a1, ..., an)` gives the component $F_n$ of its $A_\infty$ quasi-isomorphism to the normalized cobar algebra. Both take $n\ge 2$ Ext classes. With input bidegrees $(s_i,t_i)$, their result types are `Ext(sum s + 2 - n, sum t)` and `Cobar(sum s + 1 - n, sum t)`, respectively. A negative result degree gives an error. The binary operation is `yoneda`, and the generator of `Ext(0, 0)` is a strict unit.

At fixed internal degree, let $P_s$ send each boundary of cobar length $s$ to its deterministic preimage by sparse elimination on the ordered words of length $s-1$. Write $i$ for `extRep`. For a cochain $c$ of length $s$, form the cycle $q=c+P_{s+1}(dc)$, let $p(c)$ be `extClass(q)`, and put $h(c)=P_s(q+i p(c))$. These maps give the contraction $d h+h d=1+i p$.

Starting with $F_1=i$, ASL computes $U_n=\sum_{j=1}^{n-1} F_j(a_1,\ldots,a_j)F_{n-j}(a_{j+1},\ldots,a_n)$ using cobar concatenation. It returns $m_n=p(U_n)$ and $F_n=h(U_n)$. The correction $P_{s+1}(dc)$ applies even when $U_n$ has nonzero differential. The calculation reuses lifts for each consecutive interval and caches the sparse splittings by bidegree. [Merkulov, *Strongly homotopy algebras of a Kähler manifold*](https://arxiv.org/abs/math/9809172)

The higher operations can depend on the chosen contraction. The Massey product collects classes from defining systems. In this transfer, $m_3(h_0,h_1,h_0)=h_1^2$ and $m_4(h_0,h_1,h_2,h_0)=h_2^2$. The adjacent product $h_2h_0$ obstructs a Massey defining system for the four inputs. The fourth lift is $[\xi_1\xi_3]+[\xi_1^2\xi_2^2]+[\xi_1^5\xi_2]$. Run `nix run path:. -- examples/ext-transfer.asl` to compute both operations and the lift. [Transfer example](examples/ext-transfer.asl)

## Milnor duality and the homology coaction

`milnor{[r1, r2, ...]}` is the element of the Steenrod algebra dual to the monomial $\xi_1^{r_1}\xi_2^{r_2}\cdots$. Its type is `A(d)`, where `d` is the monomial degree. Trailing zeros are optional, and `milnor{[0]}` is `sq[0]`. The interpreter converts Milnor basis operations to the existing admissible basis. `milnorPair(theta, c)` evaluates the pairing with a cochain in `Cobar(1, d)` and gives `S.hi` for one and `S.lo` for zero. In degree zero it also pairs against `Cobar(0, 0)`. [Milnor, *The Steenrod Algebra and Its Dual*](https://doi.org/10.2307/1969932)

The pairing of `mul(a, b)` against a monomial uses `a` on the left coproduct factor and `b` on the right. Thus `milnorPair(mul(sq[2], sq[1]), cword{[0, 1]})` is `S.hi`, and reversing the squares gives `S.lo`. Pairing a single square with a monomial reads the coefficient of the corresponding power of $\xi_1$. [Milnor duality example](examples/milnor-duality.asl)

`coact[n, i](u)` takes `u : H(X, q)` and pairs the coaction on the homology generator dual to `gen[X, n, i]` against `u`. The coefficient of $\xi^R$ is the `i`th coordinate of `act(milnor{[r1, r2, ...]}, u)`. For `n > q`, the result has type `Cobar(1, n - q)`. For `n = q`, it gives the `i`th coordinate of `u` as a multiple of the unit in `Cobar(0, 0)`.

## Adams E₂ with finite-space coefficients

`XExt(X, s, t)` computes $\operatorname{Ext}_{\mathcal A}^{s,t}(H^*(X;\mathbb F_2),\mathbb F_2)$, the mod-two Adams $E_2$ page of $\Sigma^\infty X_+$. It uses ordinary homology, so $H_0$ contributes one sphere summand per component. `XExt(One, s, t)` therefore has the rank of `Ext(s, t)`.

A cochain basis term $[a_1\mid\cdots\mid a_s]\otimes z$ has internal degree $t=|z|+\sum_i |a_i|$, where $z$ is dual to a canonical cohomology generator of `X`. The differential splits each word factor by the reduced coproduct. It also appends a factor from the reduced Milnor coaction on $z$, using the actual Steenrod action on $H^*(X;\mathbb F_2)$. For the projective plane the top generator has differential $[\xi_1]\otimes z_1$. The wedge in the projective-plane example has the same Betti numbers and a zero differential on its top generator. [Projective plane and wedge](examples/projective-plane.asl)

| Form | Meaning |
| --- | --- |
| `xextdim[X, s, t];` | Print the rank of `XExt(X, s, t)`. |
| `xextgen[X, s, t, i]` | Select basis generator `i`, starting at zero. |
| `xaction(a, x)` | Act by an `Ext(u, v)` class `a` on `x : XExt(X, s, t)` by prepending its cobar representative. The result is in `XExt(X, s + u, t + v)`. |
| `xmap(f, x)` | Apply the map induced by `f : X -> Y` to `x : XExt(X, s, t)`, using the dual of cohomology pullback. The result is in `XExt(Y, s, t)`. |

`zero[XExt(X, s, t)]`, `add(x, y)`, and equality use the vector space structure in a fixed bidegree. The example computes the two sphere summands of a discrete two-point space and checks that the collapse commutes with scalar action. [Finite-space Adams E₂ example](examples/space-adams.asl)

## Adams E₂ for stable maps

`MapExt(X, Y, s, t)` computes $\operatorname{Ext}_{\mathcal A}^{s,t}(H^*(Y;\mathbb F_2),H^*(X;\mathbb F_2))$ for maps from $\Sigma^\infty X_+$ to $\Sigma^\infty Y_+$. The Adams filtration is `s`, the internal degree `t` accepts signed integers, and the stem is $t-s$. A length-`s` cochain term sends a homology generator of `X` in degree `q` to a cobar word of internal degree `w` tensored with a generator of `Y` in degree `p`, where $p+w=q+t$. The differential combines the cobar differential and the Steenrod coactions at both endpoints. `MapExt(One, Y, s, t)` for nonnegative `t` recovers the groups computed by `XExt(Y, s, t)`.

| Expression or statement | Result |
| --- | --- |
| `mapextdim[X, Y, s, t];` | Print the rank of `MapExt(X, Y, s, t)`. |
| `mapextgen[X, Y, s, t, i]` | Select generator `i` in the computed basis, starting from zero. |
| `emap(f)` | Take the class of the map induced by `f : X -> Y` on homology, in `MapExt(X, Y, 0, 0)`. |
| `eid[X]` | Take the identity class in `MapExt(X, X, 0, 0)`. |
| `ecomp(g, f)` | Compose classes from `X` to `Y` and from `Y` to `Z`, applying `f` first and adding both bidegrees. |

`zero[MapExt(X, Y, s, t)]`, `add`, and equality act in a fixed type. The example composes maps of a two-point space and uses degree $-1$ to send the circle's degree-one homology generator to `One`. [Stable-map Adams E₂ example](examples/map-ext.asl)

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

## How late can a higher product appear?

For every cyclic group $C_m$ of order $m=2^r\ge 4$, the mod-two cohomology ring of $BC_m$ is $\mathbb F_2[x,y]/(x^2)$, with generators in degrees one and two. Take $x(g)=g\bmod 2$ and repeat it as the input to successively longer Massey products. Every result from arity 2 through $m-1$ is zero. At arity $m$ the result is $y$. The input stays in degree one and the answer stays in degree two throughout. [Totaro, Theorem 1.2](https://assets.cambridge.org/97811070/15777/excerpt/9781107015777_excerpt.pdf) · [Benson and Greenlees, §2](https://arxiv.org/pdf/2006.07160)

```sh
nix run path:. -- examples/massey-tower.asl
awk -v m=64 -f examples/massey-sweep.awk | nix run path:. -- -
awk -v m=4096 -v first=4096 -f examples/massey-sweep.awk | nix run path:. -- -
```

The first program supplies witnesses for $C_8$. The second emits ASL source for $C_{64}$ and runs every arity from 2 through 64. The first 62 output lines read `H2(64){0}`. Line 63 reads `H2(64){carry}`. The third command checks one system of arity 4,096, supplying 4,095 distinct cochains for a triangle with 8,390,655 proper intervals. [Worked program](examples/massey-tower.asl) · [Experiment generator](examples/massey-sweep.awk)

| Group | Arity requested in a single run | Observed result |
| --- | --- | --- |
| $C_{4096}$ | 4,095 | `H2(4096){0}` |
| $C_{4096}$ | 4,096 | `H2(4096){carry}` |
| $C_{4096}$ | 4,097 | Interval `[0,4096]` fails at residues `(1,4095)`. The coboundary is zero and the products sum to one. |
| $C_{8192}$ | 8,192 | Constructing the binomial witnesses reaches the interpreter's finite-value storage limit. |

To reproduce one row, set both `first` and `last` to the requested arity. The default range is 2 through `m`.

**Supplying a defining system.** Write `massey(a, b, c){0 < 2 => u, 1 < 3 => v}` for a triple product. The inputs and witnesses all have type `Bar(m)`, with $m\ge 2$. Write any normalized cochain by listing its support, as in `bar[4]{1, 3}`. Repeated residues cancel. `binom[m, k]` gives the cochain $g\mapsto\binom{g}{k}\bmod 2$ for the representative $0\le g<m$ and $k\ge1$. Both new types support `add` and `zero[T]`. The finite group supplies an index for these types. Its classifying space can have cohomology in infinitely many degrees.

For $n$ inputs, put input $i$ in entry $a_{i,i+1}$, starting at $i=0$. Supply one witness for each longer interval except $[0,n]$. At every residue pair the interpreter checks

$$
a_{ij}(g)+a_{ij}(h)+a_{ij}((g+h)\bmod m)
=\sum_{i<k<j} a_{ik}(g)a_{kj}(h) \pmod 2.
$$

For adjacent entries the right side is zero, so each input must be a cocycle. At the omitted corner the right side gives a 2-cocycle $\omega$. The answer is its class for the supplied witnesses. For even $m$ its coordinate is $\sum_g\omega(g,1)\bmod 2$. This sum telescopes to zero for every coboundary. The carry cocycle, which is one when $g+h\ge m$, pairs to one because only $g=m-1$ contributes. Its class is `carry[m]`. For odd $m$ every 2-cocycle is exact.

Witness choices can change the answer. Over $C_2$ take $y(1)=1$ and use inputs $(y,0,y)$. Both inner products are zero, so either $0$ or $y$ is a valid witness on the right. With the left witness fixed at zero, the answer is respectively zero or the nonzero square of $y$.

```text
let y = bar[2]{1};
let z = zero[Bar(2)];
let triple = \w : Bar(2) => massey(y, z, y) { 0 < 2 => z, 1 < 3 => w };
show triple(z);   // H2(2){0}
show triple(y);   // H2(2){carry}
```

Changing one choice can change thousands of witnesses while all proper interval equations continue to hold. This experiment emits two ordinary `massey` calls with identical lists of 127 inputs. Each call supplies 8,000 inner witnesses. Of those, 1,805 change between the two calls. The output is `H2(64){carry}` followed by `H2(64){0}`. [Defining-system generator](examples/massey-choices.awk)

```sh
awk -f examples/massey-choices.awk | nix run path:. -- -
```

The generator uses a directed graph with 64 steps from its initial vertex to its terminal vertex. Every internal layer has two vertices. Between consecutive internal layers its adjacency matrix is $\left(\begin{smallmatrix}1&1\\1&0\end{smallmatrix}\right)$. Removing the edge $61\to64$ preserves every adjacent input. If $N$ is the adjacency matrix of the whole graph, the proper entries of $(I+N)^g$ give the cochains at residue $g$. The top corner of $(I+N)^{64}=I+N^{64}$ counts full paths modulo two. The first graph has 17,167,680,177,565 such paths and the second has 12,422,650,078,084. The interpreter receives just the cochains and checks both defining systems.

**Checking thousands of witnesses.** In `powerdef(b1, ..., bL)` every interval of length $j$ uses the same cochain $b_j$. The omitted corner has arity $L+1$. Set $P_h(t)=1+\sum_{j=1}^{L} b_j(h)t^j$ in $\mathbb F_2[t]/(t^{L+1})$. The proper interval equations are exactly the coefficients of $P_{g+h}=P_gP_h$. The checker tests $P_{h+1}=P_1P_h$ for every residue, including the wrap to zero. Since $P_0=1$ and $1$ generates $C_m$, these checks imply every pair equation. This works for any cochains the program supplies. Explicit `massey` checks each interval across residue pairs directly.

**Why the first obstruction can be arbitrarily late.** The binomial witnesses in the experiment are the entries of $(I+J)^g$, where $J$ has ones on the first superdiagonal and zeros elsewhere. More generally, each defining system gives a homomorphism to the group of upper unitriangular matrices modulo its top-right entry. At arity $n$, a lift of the generator has size $n+1$ and form $U=I+N$ with ones on the first superdiagonal. Since $m$ is a power of two, $U^m=I+N^m$. For $n<m$ there is no path of $m$ strict steps, so $N^m=0$ and the lift respects the group relation. For $n=m$ exactly one such path connects $0$ to $m$. Its product is one for every choice of the higher entries. No choice of witnesses removes that obstruction. [Dwyer's lifting criterion, Theorem 2.4 in Merkurjev and Scavia](https://www.math.ucla.edu/~merkurev/papers/massey4.pdf)

## The complete Massey product

`bracket(u1, ..., un)` takes at least two positive-degree cohomology classes of the same finite space $X$. It computes the ordinary mod-two Massey product on the order complex. For input degrees $d_i$, the result has type `H(X, D) -> S` with $D=2+\sum_i(d_i-1)$. Applying this function to a class returns `S.hi` exactly when some defining system gives that class. Two inputs give the singleton containing their usual cup product. An empty set means there is no defining system. A singleton containing zero means there is a defining system and every choice gives zero. [May, *Matric Massey Products*, §1](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/May-matric.pdf)

```sh
nix run path:. -- examples/massey-locus.asl
awk -f examples/bracket-relator.awk | nix run path:. -- -
```

**Four zeros and 344 answers.** Let $G$ be the graph with two vertices and four parallel edges. Its first cohomology is $V=\mathbb F_2^3$. On $G\times G$, take $z$ to be zero in degree one. In `bracket(z, z, z, z)` the three gap-two witnesses may be arbitrary cocycles, and the gap-three equations vanish. The only surviving top term is the cup product of the first and last gap-two witnesses. Thus the complete answer is the image of a bilinear map.

Identify $H^2(G\times G)$ with the $3\times3$ matrices over $\mathbb F_2$. Multiplying degree-one classes $(u,v)$ and $(w,t)$ gives $ut^\mathsf T+wv^\mathsf T$. These are exactly the matrices of rank at most two. Of the 512 matrices, $(8-1)(8-2)(8-4)=168$ are invertible. The interpreter prints `card = 344 / 512`. It accepts a diagonal matrix of rank two and the matrix supported on its missing diagonal entry. It rejects their sum, the identity. Even the sum of two attainable answers can be unattainable. [Matrix program](examples/massey-locus.asl)

**Two triples can vanish separately and still obstruct a fourfold product.** The second command builds the presentation complex $K$ of the group with generators $a,b,c,s$ and relator $[a,s][[a,b],c]$, using $[x,y]=xyx^{-1}y^{-1}$. Its two degree-two classes are $0$ and $y$. The program gives

| Inputs, in order | Full result |
| --- | --- |
| $a,b,c$ | $\{0,y\}$ |
| $b,c,a$ | $\{0,y\}$ |
| $a,b,c,a$ | $\varnothing$ |

Both triples need a witness for the middle pair $b,c$. In the upper unitriangular matrix description, write $e$ for the $(1,3)$ entry of the generator $s$, with rows indexed from zero. Evaluating the relator gives the two triple obstruction coordinates $1+e$ and $e$. Either can vanish separately. A fourfold defining system requires the same $e$ to annul both. The interpreter obtains this result from simplicial cochains. The generator attaches a cone to the source of a finite mapping cylinder and emits the corresponding poset and cocycles. [Presentation generator](examples/bracket-relator.awk) · [Barmak, §2](https://mate.dm.uba.ar/~jbarmak/thmafinal.pdf) · [Isaksen, *When is a fourfold Massey product defined?*](https://arxiv.org/abs/1309.4637)

Vary the attaching word and the ordered inputs with `awk -v word=aaaa -v queries='aa aaa aaaa aaaaa' -f examples/bracket-relator.awk | nix run path:. -- -`. This attaches a 2-cell by a degree-four map of circles, giving the Moore space $M(\mathbb Z/4,1)$. The four outputs are $\{0\},\{0\},\{y\},\varnothing$. For `word=abAbaBAB` and `queries='abb bab'`, the relator $[[a,b],b]$ gives $\{y\}$ and $\{0\}$. Uppercase letters denote inverses. [Fenn and Sjerve, *Massey Products and Lower Central Series of Free Groups*](https://doi.org/10.4153/CJM-1987-015-5)

**How the search is complete.** Each non-adjacent proper interval $[i,j]$ needs a cochain $a_{ij}$ of degree $1+\sum_{i\le k<j}(d_k-1)$ solving $d a_{ij}=\sum_{i<k<j}a_{ik}\smile a_{kj}$. The interpreter finds a primitive by sparse elimination and varies it by all linear combinations of cohomology representatives. A change by a coboundary $d\lambda$ is an upper unitriangular change of basis $I+\lambda E_{ij}$ on the defining system. It fixes shorter intervals and changes the terminal obstruction only by a coboundary after clearing the top entry. For four inputs ASL can solve the two proper triple conditions as a joint linear system. General searches that reach a resource limit report a diagnostic. For example, five zero degree-one inputs on the wedge of three circles and a 2-sphere reach the current per-statement work budget. A completed `bracket` contains the full answer.

## Observing a computation

Observations are top-level statements. Every statement ends with `;` and `//` starts a line comment.

| Statement | Result |
| --- | --- |
| `show e;` | Print a point, cohomology class, Steenrod operation, cobar cochain, Ext class, pair or function table. |
| `assert e == f;` | Check extensional equality at a common type. Success is silent. |
| `homotopy f ~ g;` | Find a shortest zigzag of pointwise-comparable maps of the same type. Print its length, or `no` if the maps lie in different components. |
| `betti T;` | Compute ordinary Betti numbers over 𝔽₂ of the order complex of `T`, starting in degree zero. Trailing zeros are omitted. The empty space yields `[]`. |
| `extdim[s, t];` | Print the dimension of `Ext(s, t)` over 𝔽₂. |
| `xextdim[X, s, t];` | Print the dimension of `XExt(X, s, t)` over 𝔽₂. |
| `mapextdim[X, Y, s, t];` | Print the dimension of `MapExt(X, Y, s, t)` over 𝔽₂. |
| `card p;` | For a predicate `p : X -> S`, print the number of points sent to `S.hi` and the size of `X`. On the empty domain this gives `card = 0 / 0`. |

A homotopy step uses a single direction across the whole domain, so it requires `f(x) ≤ g(x)` for every `x`, or the reverse inequality for every `x`. In the triangle example, any endomorphism above the identity must fix all three maximal edges. Each vertex is the unique common lower bound of its two incident edges, so it too is fixed. The dual argument handles endomorphisms below the identity. This isolates the identity in `Circle -> Circle` and explains the `no` printed for the rotation. [Stong, *Finite Topological Spaces*, §3–4](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/stong2.pdf)

For `betti`, a strict chain of `k + 1` points contributes a `k`-simplex. Dropping each point in turn gives its boundary over 𝔽₂ in positive degrees. The ranks of successive boundary matrices determine the homology. [Barmak and Minian, *Simple Homotopy Types and Finite Spaces*, §2](https://math.uchicago.edu/~may/FINITE/ArgentinianPapers/MinianSimple.pdf)

The evaluator tabulates a lambda on its domain. For a discrete `D` with `n` points, an individual map `D -> S` needs `n` entries, while enumerating the whole function space yields `2^n` maps. A total order of `n` points has `2^n − 1` nonempty chains. The interpreter admits up to 512 points per named space, 4,096 elements per enumerated type and 50,000 simplices per order complex. A cohomology group of rank `b` has `2^b` elements. Individual classes can be constructed and multiplied even when enumerating that space would exceed the limit. Work and storage budgets can stop a calculation earlier with a diagnostic.
