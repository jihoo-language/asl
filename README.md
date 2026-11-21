# ASL · Alexandrov Section Lambda

ASL is a pure functional language that starts from finite T₀ spaces and continuous maps. It also computes with graded cochains and sparse polynomial rings. Its two-point space `S` has exactly three endomorphisms. You can return the low point, preserve the input, or return the high point. Try to exchange the two points and the interpreter reports the order relation your function breaks.

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
| `PCobar(p, s, t)` | The corresponding discrete space of normalized cobar cochains over an odd prime field $\mathbb F_p$. |
| `PCobarAn(p, n, s, t)` | The odd-prime cobar cochains for the finite Steenrod subalgebra $\mathcal A_p(n)$. |
| `Ext(s, t)` | The discrete space of mod-two Adams $E_2$ classes in filtration `s` and internal degree `t`. |
| `PExt(p, s, t)` | The discrete space of odd-prime Adams $E_2$ classes over $\mathbb F_p$. |
| `PExtAn(p, n, s, t)` | The discrete space of Ext classes over $\mathcal A_p(n)$. |
| `PMay(p, r, s, t, f)` | The discrete space of odd-prime algebraic May classes on page `r` in filtration `f`. |
| `PMayAn(p, n, r, s, t, f)` | The algebraic May page of the finite odd-primary profile $\mathcal A_p(n)$. |
| `PCEAn(p, n, r, s, t, f)` | The Cartan–Eilenberg page for the central extension by $E(Q_n)$ in $\mathcal A_p(n)$. |
| `CobarAn(n, s, t)` | The discrete space of normalized cobar cochains for the dual of $\mathcal A(n)$. |
| `ExtAn(n, s, t)` | The discrete space of Ext classes over $\mathcal A(n)$ in filtration `s` and internal degree `t`. |
| `XExt(X, s, t)` | The discrete space of mod-two Adams $E_2$ classes for the suspension spectrum of `X` with a disjoint basepoint. |
| `MapExt(X, Y, s, t)` | The discrete space of mod-two Adams $E_2$ classes for stable maps from `X` to `Y`, with disjoint basepoints. |
| `MapMay(X, Y, r, s, t, f)` | The discrete space of algebraic May page classes of bivariant cobar morphisms. |
| `MapCobar(X, Y, s, t)` | The discrete space of length-`s` bivariant cobar cochains of signed internal degree `t`. |
| `QKoszul(X, Y, I, s, t)` | The discrete space of bivariant Koszul cochains for a finite set of Milnor primitives. |
| `QExt(X, Y, I, s, t)` | The corresponding cohomology classes for the exterior subalgebra generated by those primitives. |
| `VKoszul(X, Y, I, s, t)` | Bivariant Koszul cochains with the highest variable inverted. |
| `VExt(X, Y, I, s, t)` | Cohomology classes of the localized complex. |
| `RVHom(X, Y, I, modulus, point)` | Linear maps between the Margolis homologies of two endpoints at an affine field point. |
| `RVLaurent(X, Y, I, modulus, point)` | Laurent-polynomial matrices between the Margolis homology spaces. |
| `May(r, s, t, f)` | The discrete space of classes on the algebraic May page `r` in cobar length `s`, internal degree `t`, and filtration `f`. |
| `MayAn(n, r, s, t, f)` | The corresponding May page for the finite Steenrod subalgebra `A(n)`. |
| `CEAn(n, r, s, t, f)` | The Cartan–Eilenberg page for the central exterior generator of `A(n)`. |
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

## Sheaf cohomology on a finite poset

`sheaf F[X,p] {a:2,b:1} {a<b:[1,-1]};` declares a sheaf of vector spaces over the prime field $\mathbb F_p$ on the named poset `X`. Give a dimension for every point and one restriction $F(a)\to F(b)$ for each Hasse cover $a<b$. List matrix entries by rows, with the target dimension giving the number of rows. Zero-dimensional sources or targets require an empty matrix list. Signed decimal entries are reduced modulo $p$. The declaration composes along the Hasse diagram in topological order and checks that all paths with the same endpoints give the same matrix.

`ShC(F,k)` is the normalized Roos cochain space in degree $k$. The factor for a strict chain $x_0<\cdots<x_k$ is the terminal stalk $F(x_k)$. Use `shTerm[F,[x0,...,xk],i]` for its basis coordinate $i$, numbered from zero, or add a signed coefficient as the fourth argument. `shD(c)` takes the alternating sum over faces. When a face drops the terminal point, the contribution uses the restriction from its old terminal stalk to the new one. In particular, on a two-point chain the differential of a section is $c(b)-r_{ab}c(a)$.

`ShH(F,k)` is the resulting cohomology space. `shClass(c)` projects a closed cochain to its class and `shRep(h)` returns a deterministic closed representative. Use `shdim[F,k];` for the dimension and `shgen[F,k,i]` for cohomology basis element $i$. The basis calculation orders chains by their point indices in the space declaration, then by terminal stalk coordinate. Both types support `zero`, `add`, and signed scalar multiplication `shScale[c](v)`. Sparse cochains can be formed and differentiated directly. Cohomology uses bounded prime-field reduction of the adjacent degrees.

The circle example compares constant coefficients, rank-two unipotent holonomy over $\mathbb F_3$, scalar holonomy $2$ over $\mathbb F_5$, and a sheaf supported on the three maximal points. [Finite-poset sheaf example](examples/sheaf-cohomology.asl)

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

## The complex-motivic cobar algebra

`MCobar(s, t, w)` is the length, topological degree, and signed weight component of the normalized dual Steenrod cobar over $\mathbb C$. The coefficient ring is $\mathbb F_2[\tau]$, with $\tau$ of degree zero and weight minus one. In `mcword[k]{([xi-powers], [tau-powers]), ...}`, each pair gives the powers of the $\xi_i$ starting at one and the $\tau_i$ starting at zero. The entire word has coefficient $\tau^k$. A factor must have positive topological degree. Powers of $\tau_i$ normalize by $\tau_i^2=\tau\xi_{i+1}$, contributing to the same global coefficient.

Use `cobarD`, `concat`, `add`, and `zero` on motivic cochains. The differential preserves topological degree and weight. `tauShift[k](x)` multiplies by $\tau^k$ and lowers weight by `k`. `realize(x)` sends the coefficient $\tau$ to one, $\xi_i$ to the square of the classical Milnor generator with index $i$, and $\tau_i$ to the one with index $i+1$. It commutes with products and the differential.

A normalized monomial $\prod_{i\ge1}\xi_i^{b_i}\tau_{i-1}^{\epsilon_{i-1}}$ realizes to the classical exponents $2b_i+\epsilon_{i-1}$. If a word has $o$ odd exponents across its slots, its intrinsic weight is $W=(t-o)/2$, and its global coefficient exponent is $W-w$. This gives a finite basis at every tridegree. In a coproduct the number of odd exponents can drop when two $\tau_i$ factors meet. For example the reduced coproduct of $\tau_0\tau_1$ contains $\tau[\xi_1|\xi_1]$. [Complex-motivic cobar example](examples/motivic-cobar.asl)

### The complex-motivic Adams E2

`MExt(s, t, w)` is cohomology of the trigraded motivic cobar. Here `s` is Adams filtration and `t-s` is the stem. `mextdim[s, t, w];` computes the dimension and `mextgen[s, t, w, i]` selects generator `i`. Use `extClass(x)` on a closed motivic cochain, `extRep` to read a cocycle, and `extPrimitive` to find a cochain bounding a given motivic boundary. `yoneda` multiplies classes by concatenating representatives. `tauShift[k]` and `realize` act on classes by applying their cochain maps and taking the resulting class.

The elements $h_0=[\tau_0]$, $h_1=[\xi_1]$, and $h_2=[\xi_1^2]$ satisfy $h_0 h_1=0$ and $\tau h_1^3=h_0^2 h_2\ne0$. The class $h_1^4$ is nonzero in weight four. Multiplying by $\tau$ gives a boundary in weight three, where the explicit five-term primitive can be checked with `cobarD`. In even topological degree $t$ at weight $t/2$, every exponent is even under realization. Halving those exponents identifies that slice with the classical cobar in half the topological degree. These computations concern the motivic Adams $E_2$ page. [Complex-motivic Ext example](examples/motivic-ext.asl)

## The real-motivic cobar algebra

`RMCobar(s, t, w)` is the normalized real-motivic dual Steenrod cobar over $\mathbb F_2[\tau,\rho]$. Length `s` is nonnegative, and topological degree `t` and weight `w` are signed. The degrees of $\tau$ and $\rho$ are $(0,-1)$ and $(-1,-1)$. Use `rmcword[a,b]{([xi-powers], [tau-powers]), ...}` for a word with global coefficient $\rho^a\tau^b$. Its factors specify $\xi_i$ powers starting at one and $\tau_i$ powers starting at zero.

The right unit sends $\tau$ to $\tau+\rho\tau_0$ and fixes $\rho$. The source coefficient in arrow slot $j$ is therefore $\tau+\rho\sum_{k<j}\tau_{0,k}$. At that slot the square relation is $\tau_i^2=(\tau+\rho\sum_{k\le j}\tau_{0,k})\xi_{i+1}+\rho\tau_{i+1}$, with the unindexed generators interpreted in slot $j$. The interpreter expands these squares into normal monomials. Applying `cobarD` transports the coefficient in the first coface and splits the arrow factors in the internal cofaces. Unit arrows are discarded after normalizing each coface. `concat` transports the second word's coefficient across the first word's arrows before multiplying. The operations `add`, `zero`, and `rhoShift[k]` act on cochains, with the last multiplying by $\rho^k$.

If a normal word has topological degree $T$ before its global coefficient and has $o$ odd exponents in its encoding by classical words, set $W=(T-o)/2$. A monomial of tridegree $(s,t,w)$ has coefficient $\rho^{T-t}\tau^{W-(T-t)-w}$. Its possible degrees $T$ lie from $\max(0,t)$ to $2(t-w)$. This bounds the basis even at negative weights. The chain map `realComplex(x)` sets $\rho=0$ and has type `MCobar(s, t, w)` at the same tridegree. Complex `MCobar` and `MExt` have zero groups in negative topological degree. The classical `realize` operation requires a topological degree in the unsigned classical range.

For example $d\tau=\rho[\tau_0]$ and $d(\tau^2)=\rho^2[\tau_0^2]$, with the latter expanding to three normal terms. [Real-motivic cobar example](examples/real-motivic-cobar.asl)

## The odd-prime cobar algebra

`PCobar(p, s, t)` uses coefficients in $\mathbb F_p$ for an odd prime `p` that fits in an unsigned 32-bit integer. The dual Steenrod algebra is polynomial on $\xi_i$ of degree $2(p^i-1)$ for $i\ge1$ and exterior on $\tau_i$ of degree $2p^i-1$ for $i\ge0$. The literal `pword[3]{([1], [0]), ([], [1])}` denotes $[\xi_1\tau_0\mid\tau_1]$. Each factor gives the polynomial exponents in generator order and the strictly increasing exterior indices. Every factor has positive degree. The empty word `pword[p]{}` gives the tensor unit.

The coproducts are $\Delta(\xi_n)=\sum_{i=0}^n\xi_{n-i}^{p^i}\otimes\xi_i$ and $\Delta(\tau_n)=\tau_n\otimes1+\sum_{i=0}^n\xi_{n-i}^{p^i}\otimes\tau_i$, where $\xi_0=1$. Multiplying tensor factors uses the sign $(-1)^{|b||c|}$ in $(a\otimes b)(c\otimes d)$. Exterior factors use increasing order.

For the reduced coproduct $\overline\Delta u=\sum u'\otimes u''$, the cobar differential is $d[u]=\sum(-1)^{|u'|}[u'\mid u'']$. Its graded Leibniz rule uses the stem $t-s$. For example, at $p=3$, $d[\xi_1^2]=2[\xi_1\mid\xi_1]$ and $d[\tau_1]=[\xi_1\mid\tau_0]$.

`add`, `concat`, `cobarD`, and `zero` operate on these cochains. Use `pscale[k](c)` to multiply by the residue of the nonnegative integer `k` in the cochain's field. Sparse operations work directly on their supports. A lambda over `PCobar(p, s, t)` tabulates all $p^r$ cochains for word basis size $r$ subject to the finite enumeration limit. [Odd-prime cobar example](examples/odd-cobar.asl)

## The Adams E₂ page

`Ext(s, t)` is $H^s(\mathrm{Cobar}(*,t),d)$, the classical mod-two Adams $E_2$ page of the sphere. The corresponding stem is $t-s$. ASL uses sparse elimination over $\mathbb F_2$ on cobar differentials. `extdim[s, t];` prints the rank, and `extgen[s, t, i]` selects the $i$th generator in the computed basis, starting from zero.

`extClass(c)` takes a closed cobar cochain in bidegree `(s, t)` to its class in `Ext(s, t)`. `extRep(x)` returns a deterministic closed representative of `x`. For a cobar boundary of positive length, `extPrimitive(c)` returns a deterministic `p` with `cobarD(p) == c`. `yoneda(x, y)` multiplies classes by concatenating their representatives and taking the result's class. Classes of the same bidegree support `add` and `zero[Ext(s, t)]`.

Run `nix run path:. -- examples/adams-e2.asl` to compute $h_0h_1=0$ and $h_1^3=h_0^2h_2\ne0$. [Adams E₂ example](examples/adams-e2.asl)

### Restriction to finite Steenrod subalgebras

For $n\geq 0$, the dual of $\mathcal A(n)$ is $\mathbb F_2[\xi_1,\ldots,\xi_{n+1}]/(\xi_i^{2^{n+2-i}})$. `CobarAn(n, s, t)` uses words in the surviving positive-degree monomials. `ExtAn(n, s, t)` is the cohomology in length `s` and internal degree `t`. `extandim[n, s, t];` prints the rank, and `extangen[n, s, t, i]` selects generator `i`, starting from zero.

`resA[n](c)` takes a full cobar cochain or a cochain over $\mathcal A(m)$ with $m\geq n$ to `CobarAn(n, s, t)`. It kills any word with a factor outside the quotient monomial basis. On an Ext class, it restricts a closed representative and takes its class in `ExtAn(n, s, t)`. Use `resA[n](cword{...})` to construct cochains. The operations `cobarD`, `concat`, `extClass`, `extRep`, `extPrimitive`, and `yoneda` preserve the input algebra profile, and `concat` and `yoneda` require both inputs to use the same profile.

The corrected full-algebra May lift `mayRep(maygen[7, 4, 12, 0, 0])` becomes a cocycle over $\mathcal A(1)$ because each term of its differential is killed by restriction. Its class generates `ExtAn(1, 4, 12)` and represents $v_1^4$. The full group `Ext(4, 12)` and the restricted group `ExtAn(0, 4, 12)` both have rank zero. [Chromatic restriction example](examples/chromatic-restriction.asl)

## The odd-prime Adams E₂ page

`PExt(p, s, t)` computes the cohomology of the signed `PCobar(p, s, t)` complex. The statement `pextdim[p, s, t];` prints the dimension over $\mathbb F_p$, and `pextgen[p, s, t, i]` selects basis generator `i` starting at zero. A component of rank $r$ has $p^r$ points, including zero.

The operators `extClass(c)`, `extRep(x)`, and `extPrimitive(c)` work on odd-prime cochains and classes. They project closed cochains, choose deterministic lifts, and solve the cobar boundary equation over $\mathbb F_p$. The operators `pscale[k](x)` and `add(x, y)` act on field coefficients. The Yoneda product `yoneda(x, y)` concatenates lifts at a common prime. Its graded commutation sign uses the stem $t-s$.

### Restriction to finite subalgebras at odd primes

The dual of $\mathcal A_p(n)$ is $\mathbb F_p[\xi_1,\ldots,\xi_n]/(\xi_i^{p^{n+1-i}})\otimes\Lambda(\tau_0,\ldots,\tau_n)$. `PCobarAn(p, n, s, t)` uses words in the quotient monomials, and `PExtAn(p, n, s, t)` computes their signed cobar cohomology. The case $n=0$ has only the exterior generator $\tau_0$. Use `pextandim[p, n, s, t];` to find the rank and `pextangen[p, n, s, t, i]` to select generator `i`.

`resA[n]` preserves the odd prime and projects full cochains or cochains from a larger finite profile by keeping the words whose factors all survive. On classes it performs the same operation on a chosen representative. The cochain, Ext, and field operations preserve the profile. Binary products require the same prime and profile on both inputs.

At $p=3$, the full `PMay(3, 5, 3, 15, 6)` lift supports a nonzero $d_5$. Over $\mathcal A_3(1)$ its projection becomes a closed representative of $v_1^3$ in `PExtAn(3, 1, 3, 15)`. The full Ext group and its further restriction to $\mathcal A_3(0)$ both have rank zero. [Odd-primary restriction example](examples/odd-profile-restriction.asl)

The example computes powers of $a_0=[\tau_0]$ and the divided binomial class $b_{1,0}$ at $p=3$ and $p=5$. Its cobar coefficient on $[\xi_1^i\mid\xi_1^{p-i}]$ is the integer $\binom{p}{i}/p$ reduced modulo $p$. It also checks the nonzero anticommuting products of $[\xi_1]$ and $[\xi_1^9]$ at $p=3$. [Odd-prime Adams E₂ example](examples/odd-adams.asl)

### Odd-prime Massey defining systems

`pMassey(x1, ..., xn){i<j => w, ...}` reads its adjacent inputs from left to right. Each input is a `PExt` class or a closed `PCobar` cochain at the same odd prime. It uses the canonical cobar lifts of the classes. In the defining triangle, positions run from $0$ to $n$, and an omitted proper witness is zero. On the interval $[i,j]$, the witness is a `PCobar` cochain of cobar length $1+\sum s_k-(j-i)$ and internal degree $\sum t_k$ for the inputs in that interval.

Write $\overline A=(-1)^{t-s+1} A$ for a cochain of bidegree $(s,t)$. Each proper interval must satisfy $d A_{ij}=\sum_{i<k<j}\overline A_{ik} A_{kj}$ with cobar concatenation. `pCurvature` takes the same inputs and returns the closed upper-right sum in `PCobar`. `pMassey` takes its class in `PExt`. Both evaluate the supplied defining system.

`pPowerDef(b1, ..., b_(n-1))` uses `b1` for each adjacent input and `b_k` on every proper interval of length $k$. For $n=p$ and $b_k=[\xi_1^k]/k!$, the binomial coproduct gives the required proper equations. The upper-right sum represents the negative of the nonzero $b_{1,0}$ class because $(p-1)!=-1$ in $\mathbb F_p$. [Odd-prime power examples at 3 and 5](examples/odd-massey-power.asl)

## The odd-prime May spectral sequence

`PMay(p, r, s, t, f)` computes the algebraic May page of the signed `PCobar` complex over $\mathbb F_p$. The weight of $\xi_i^e$ is $(2i-1)$ times the sum of the base-$p$ digits of $e$, and the weight of $\tau_j$ is $2j+1$. Sum the weights over the factors of a word to get $M$. In internal degree $t$, filtration $F^f$ contains the words with $t-M\ge f$.

`pmaydim[p, r, s, t, f];` prints a rank and `pmaygen[p, r, s, t, f, i]` selects basis generator `i` starting at zero. `mayClass[r, f](c)` projects a raw `PCobar` lift. Page zero is the associated graded complex. For $r>0$, a lift in $F^f$ must have differential in $F^{f+r}$. The projection takes the leading piece of a filtered cocycle modulo the leading pieces of boundaries from $F^{\max(0,f-r+1)}$.

`mayRep(x)` returns a chosen raw lift and `mayD(x)` raises $s$ by one and $f$ by $r$. `mayNext(x)` takes a $d_r$-cycle to its class on page $r+1$, and `mayPrev(x)` lifts a class to a cycle on the preceding page. `mayMul(x, y)` concatenates lifts at a common prime and page. The differential satisfies the Leibniz rule with the sign of the stem $t-s$, and `pscale[k](x)` and `add(x, y)` act on field coefficients.

The example computes $d_1([\tau_1])=[\xi_1\mid\tau_0]$ at $p=3$. It also computes a nonzero $d_5$ from bidegree $(3,15)$ in filtration $6$ to $(4,15)$ in filtration $11$. Both components vanish on $E_6$. [Odd-prime May example](examples/odd-may.asl)

### Odd-primary May pages over finite Steenrod subalgebras

`PMayAn(p, n, r, s, t, f)` applies the same weight and filtered quotient to `PCobarAn(p, n, s, t)`. The statement `pmayandim[p, n, r, s, t, f];` prints a rank, and `pmayangen[p, n, r, s, t, f, i]` selects generator `i`. The May operations and `pscale` preserve the prime and profile.

`resA[n]` restricts a full odd-primary May class or a class over a containing finite profile at the same prime. Each projected word keeps its weight, so restriction commutes with page products, differentials, and passage to the next page on cycles. At $p=3$ the full $d_5$ on the corrected $[\tau_1]^3$ lift has nonzero target. Over $\mathcal A_3(1)$ that target vanishes and the restricted lift is a cocycle representing the nonzero class in `PExtAn(3, 1, 3, 15)`. [Odd finite May example](examples/odd-profile-may.asl)

## Algebraic Steenrod squares on Adams E₂

`cCup[i](a, b)` is the normalized cobar cup-$i$ product. The inputs have a common Steenrod profile, either full `Cobar` or the same finite `CobarAn`. For input bidegrees $(s,t)$ and $(u,v)$ the result has bidegree $(s+u-i,t+v)$. An index greater than $\min(s,u)$ gives zero when the result length is nonnegative. Index zero is concatenation.

For result length $N$, the calculation takes increasing subsets $U\subseteq\{0,\ldots,N\}$ of size $N-i$. It assigns the $j$th vertex $u_j$ to the first input when $u_j+j$ is even and to the second otherwise, with $j$ starting at one. It keeps assignments with $N-s$ and $N-u$ cofaces, respectively. The cofaces lift both words to length $N$ using the full coproduct for an inner coface and a unit for an end coface. ASL multiplies the lifted words factor by factor in the dual algebra and projects away the words containing units. Finite profiles use their quotient multiplication.

`extSq[k](x)` applies cup-$(s-k)$ to a canonical representative of $x$ and itself, then takes the class in bidegree $(s+k,2t)$. It works on `Ext` and `ExtAn`. The result is independent of the representative. Indices above $s$ give zero. The top square is the Yoneda square. Here $\mathrm{Sq}^0$ is Frobenius on the Milnor exponents, so $\mathrm{Sq}^0(h_i)=h_{i+1}$ in the full algebra. The example computes the nonzero intermediate square $\mathrm{Sq}^1(h_0h_2)=h_0^2h_3$. [Cobar squares example](examples/cobar-squares.asl)

### Corrected cochain squares and May lifts

On a length-`s` binary cobar cochain `a`, write $i=s-k$. The expression `cNakamura[k](a)` computes $N^k(a)=a\smile_i a+a\smile_{i+1}da$ and has bidegree $(s+k,2t)$. It satisfies $dN^k(a)=N^k(da)$ on arbitrary cochains. At $k=0$ it squares every factor, at $k=s+1$ it gives `concat(a, cobarD(a))`, and higher indices give zero. On a closed input its class is `extSq[k](extClass(a))`. It also works over a finite binary Steenrod profile.

`mayLiftTo[r, f](a)` takes `a` in $F^f$ and solves by deterministic linear reduction for a correction in $F^{f+1}$. The corrected cochain has differential in $F^{f+r}$, and its filtration-`f` terms agree with the input. The operation returns the original cochain if it already meets the condition and reports an obstruction if no correction exists. It works on `Cobar`, `CobarAn`, and `MapCobar`.

Starting from $u=[\xi_2]$, the example computes $b=N^1(u)$ and $c=N^2(b)$. The differential of $b$ gives $d_3$. The raw $c$ satisfies the $E_6$ condition. Correcting it to satisfy the $E_7$ condition exhibits the nonzero $d_7$. Its restriction to $A(1)$ lifts to a closed representative of $v_1^4$. [Corrected square and May lift example](examples/nakamura-may.asl)

## The algebraic May spectral sequence

Give a generator power $\xi_i^{2^j}$ May weight $2i-1$. For an arbitrary monomial, count the set bits of each exponent and multiply that count by the corresponding generator weight. The weight $M$ of a cobar word is the sum over its factors. In internal degree $t$, the decreasing filtration $F^f$ is spanned by words with $t-M\ge f$. The cobar differential preserves this filtration, and concatenation adds filtration indices. The associated graded complex is page $E_0$. Its differential preserves $f$, and $d_r$ increases $s$ by one and $f$ by $r$, preserving $t$. [May, *The cohomology of restricted Lie algebras and of Hopf algebras*](https://doi.org/10.1016/0021-8693(66)90009-3)

| Expression or statement | Result |
| --- | --- |
| `maydim[r, s, t, f];` | Prints the rank of the page component. |
| `maygen[r, s, t, f, i]` | Returns generator `i` of a deterministic basis, indexed from zero. |
| `mayClass[r, f](c)` | Projects a cobar cochain in $F^f$ whose differential lies in $F^{f+r}$. |
| `mayRep(x)` | Returns a deterministic cobar lift of the page class. |
| `mayD(x)` | Evaluates $d_r$ on the page class. |
| `mayNext(x)` | Projects a $d_r$-cycle to the next page. |
| `mayPrev(y)` | Chooses a cycle on the previous page representing `y`. |
| `mayMul(x, y)` | Multiplies two classes on the same page, adding their cobar lengths, internal degrees, and filtration indices. |

`zero[May(r, s, t, f)]` and `add(x, y)` use the vector space structure in a fixed page component. For $r\ge1$, projection identifies lifts that differ by an $(r-1)$-cycle in $F^{f+1}$ or a cobar boundary in $F^f$ with a primitive in $F^{\max(0,f-r+1)}$. The computation stabilizes separately in each internal degree. The stabilized groups are the associated graded groups of `Ext(s, t)` for the induced filtration.

Run `nix run path:. -- examples/may.asl` to compute $d_1(h_{2,0})=h_1h_0$ and the nonzero $d_3$ of a corrected $b_{2,0}$ lift in internal degree 6. The page ranks record the effect on the source and target. A second program computes a nonzero $d_7$ in internal degree 12 and checks that both page components vanish on $E_8$. [First and third differentials](examples/may.asl) · [Seventh differential](examples/may-seventh.asl)

`mayNext(x)` uses the quotient $E_{r+1}=H(E_r,d_r)$. It sends an incoming $d_r$ boundary to zero and reports an error when $d_r(x)$ is nonzero. A class can survive even when its current cobar representative needs a higher-filtration correction. `mayPrev(y)` uses the canonical cobar lift of a class on page $r>0$ to choose an earlier cycle, determined modulo incoming $d_{r-1}$ boundaries. Thus `mayNext(mayPrev(y)) == y`. Both operations preserve the bidegree, filtration, and the specified profile or map endpoints. [Survival example](examples/may-survival.asl)

### May pages over finite Steenrod subalgebras

`MayAn(n, r, s, t, f)` applies the same filtered construction to `CobarAn(n, s, t)` and stabilizes to the associated graded of `ExtAn(n, s, t)`. The statement `mayandim[n, r, s, t, f];` prints the rank, and `mayangen[n, r, s, t, f, i]` gives basis generator `i` starting at zero. Use `mayClass[r, f]` and `mayRep` to pass between classes and finite-profile cobar lifts. The operations `mayD` and `mayMul` preserve the profile. A product takes two classes on the same page over the same subalgebra.

On a scalar May class, `resA[n]` induces restriction to `A(n)` in the same page and degrees. It restricts the cobar lift along the dual Hopf quotient and projects to the target page. Surviving words keep their May indices, so restriction commutes with page differentials and products. Restricting through an intermediate subalgebra gives the direct restriction.

The example computes the nonzero full $d_7$ on $b_{2,0}^2$. Its restriction in `A(1)` has zero $d_7$ and survives to the nonzero class $v_1^4$ in `ExtAn(1, 4, 12)`. Restriction to `A(0)` kills the class. [Finite-profile May example](examples/profile-may.asl)

## The Cartan–Eilenberg spectral sequence of A(n)

The central extension $E(Q_n)\to\mathcal A(n)\to\mathcal A(n)//E(Q_n)$ filters the finite-profile cobar complex. In a word of length `s`, let `q` count factors containing $\xi_{n+1}$. This top generator is exterior in $\mathcal A(n)_*$, so its exponent is zero or one in each factor. The base filtration is `f = s - q`. A cobar split either keeps the top generator or replaces it with lower generators, raising `f` by one or two. Concatenation adds filtrations.

`CEAn(n, r, s, t, f)` uses total cobar length `s` and internal degree `t` on page `r`. Page zero is the associated graded complex. Its differential is zero, so page one has the same groups and its differential records the filtration-one splits. Page two is the Ext of the lower-generator truncated Hopf coalgebra tensored with the polynomial algebra on the fiber class $v_n=[\xi_{n+1}]$. Its $d_2(v_n)$ is represented by the reduced coproduct of $\xi_{n+1}$ in the lower generators.

`ceandim[n, r, s, t, f];` prints the rank and `ceangen[n, r, s, t, f, i]` selects generator `i`. `ceClass[r, f]` takes a `CobarAn` cochain satisfying the page cycle condition and `ceRep` returns a corrected lift. `ceD` increases `s` by one and `f` by `r`. `ceMul` multiplies classes on the same page over the same profile. `ceNext` takes a $d_r$-cycle to its next-page class, and `cePrev` uses the chosen lift to represent a class on the previous page.

For $\mathcal A(1)$ the example computes $d_2(v_1)=h_0h_1$ and $d_3(v_1^2)=h_1^3$. The class of $v_1^4$ survives and its closed lift generates `ExtAn(1, 4, 12)`. [Cartan–Eilenberg example](examples/cartan-eilenberg.asl)

### Odd-primary central extensions

`PCEAn(p, n, r, s, t, f)` uses the signed cobar complex of $\mathcal A_p(n)$ with the filtration for $E(Q_n)\to\mathcal A_p(n)\to\mathcal A_p(n)//E(Q_n)$. If $q$ counts the factors containing $\tau_n$ in a word of length $s$, its filtration is $f=s-q$. Splitting a factor increases $f$ by one when its $\tau_n$ count stays the same and by two when it loses $\tau_n$. The page-zero differential is therefore zero.

`pceandim[p, n, r, s, t, f];` prints a rank and `pceangen[p, n, r, s, t, f, i]` selects a generator. The `ceClass[r, f]` and `ceRep` operations project and lift cochains in `PCobarAn`. The operations `ceD`, `ceMul`, `ceNext`, `cePrev`, and `pscale` retain coefficients in $\mathbb F_p$ and preserve the profile.

Over $\mathcal A_p(1)$ the class of $[\tau_1]^{p-1}\mid[\xi_1]$ survives to $E_p$ at $p=3$ and $p=5$. Its differential is $-a_0^{p-1}b_{1,0}$, with the $\tau_0$ factors first. Both components vanish on $E_{p+1}$. The class with leading term $[\tau_1]^p$ survives and gives a nonzero class in bidegree $(p,p(2p-1))$ of `PExtAn` at these primes. [Odd-primary Cartan–Eilenberg example](examples/odd-ce.asl)

## Brown–Peterson cobar coordinates

`BP(p, n, a)` is the integral polynomial ring in `v1` through `vn` and the coordinates `t1_1` through `ta_n` of `a` composable arrows. Here `p` is a prime that fits in 32 bits. Height `n = 0` gives the integers. Use `bpv[p, n, a, i]` for $v_i$, `bpt[p, n, a, j, i]` for the $i$th coordinate of arrow $j$, and `bpint[p, n, a, k]` for a signed integer. Arrow and generator indices start at one. `add`, `bpNeg`, `bpMul`, and `bpPow[k]` perform integral arithmetic within one ring. `zero[BP(p, n, a)]` gives its zero polynomial.

`bpCoface[j](f)` goes from `BP(p, n, a)` to `BP(p, n, a + 1)`. Coface zero transports the base through the first arrow and shifts the old arrow slots. Cofaces `1` through `a` merge adjacent arrows over their transported base. Coface `a + 1` forgets the last arrow. `bpD(f)` sums all cofaces of the unnormalized complex with signs $(-1)^j$. In height two at $p=2$, for example, $d(t_{1,2})=v_1t_{1,1}t_{2,1}-t_{1,1}t_{2,1}^2$.

The coordinate formulas use the Hazewinkel recursion $pL_k=\sum_{i=0}^{k-1}L_i v_{k-i}^{p^i}$ with $L_0=1$ and the right unit $\eta_R(L_k)=\sum_{i=0}^k L_i t_{k-i}^{p^i}$ with $t_0=1$. Intermediate log coefficients use exact powers of $p$ as denominators. Each public coordinate formula must divide integrally. ASL bounds heights by 8 and arrow counts by 16. Coefficients and intermediate numerators and denominators use checked signed 128-bit arithmetic. Large expansions can reach the term storage or computation limit. [Brown–Peterson cobar example](examples/bp-cobar.asl)

### Invariant ideals and connecting cocycles

`BPI(p, n, a, h)` is the quotient of `BP(p, n, a)` by the invariant ideal $I_h=(p,v_1,\ldots,v_{h-1})$, for $1\le h\le n+1$. `bpReduce[h](f)` takes an integral cochain or a cochain from a shallower quotient to this quotient. It fixes elements already at depth `h`. Addition, `bpNeg`, `bpMul`, `bpPow`, `bpCoface`, and `bpD` act in the quotient. `bpLift(f)` uses integer coefficients in $0,\ldots,p-1$ and the monomials that survive the quotient. It selects a representative in each cochain degree.

For a closed `BPI` cochain `x`, `bpDelta(x)` takes the differential of its canonical lift in the preceding quotient and divides exactly by $v_{h-1}$ for $h\ge2$, or by $p$ in `BP` at $h=1$. The result is a closed cochain one arrow longer. It represents the connecting map for multiplication by this invariant element. A noncycle gives an error. On integral `BP` cochains `bpBockP[e](f)` computes $d(f)/p^e$ if every coefficient is divisible by $p^e$.

At $p=3$ the first connecting cocycle of $v_2$ from $I_2$ is $2v_1^2t+t^3$ modulo $p$. Its next connecting representative is $\beta=4v_1tu+5t^2u-tu^2$, with $t$ and $u$ the two arrows' first coordinates. The differential of the integral lift $2v_1^2t+t^3$ is $3\beta$. In weight 12 the normalized one-cochains have basis $v_1^2t$, $v_1t^2$, and $t^3$. Their differentials have $t^2u$ coefficients $9$, $0$, and $-3$, respectively, while $\beta$ has coefficient $5$. Thus $\beta$ has order three in the integral cobar cohomology of this weight.

At $p=2$ the analogous canonical iterated lifts of $v_2$ and $v_3$ vanish. At $p=3$ the three-step lift of $v_3$ is a closed 27-term cochain. At $p=2$ and height one, `bpBockP[3]` applied to $v_1^4$ gives $\alpha_{4/3}$ with $8\alpha_{4/3}=d(v_1^4)$. Here $v_1^4$ spans the zero-cochains of weight eight, so the nonzero cocycle detects order eight at that height. [Greek-letter cocycle example](examples/bp-greek.asl)

### Cobar cup products and Massey cocycles

`bpCup(f, g)` takes cochains with the same prime, height, and ideal depth. If `f` and `g` have `a` and `b` arrows, the result has `a + b` arrows. The Alexander–Whitney formula evaluates `f` on the first `a` arrows and `g` on the last `b` arrows over the base transported through the first `a`. It then multiplies the two values. The product is associative and satisfies $d(f\smile g)=df\smile g+(-1)^a f\smile dg$. For example, $t\smile v_1=(v_1+pt)t$ and $v_1\smile t=v_1t$ when $v_1$ is a zero-cochain and $t$ is the first coordinate of one arrow.

`bpDegeneracy[j](f)` inserts an identity arrow at slot `j + 1`. It sets that slot's coordinates to zero, shifts the later slots, and lowers the number of arrows by one. `bpNormalize(f)` projects onto the normalized cochains by applying $x\leftarrow x-d^{j+1}s^j(x)$ in increasing order of `j`. Its result lies in every degeneracy kernel, and it commutes with `bpD`. The cup product preserves normalized cochains.

`bpMassey(A, B, C, U, V)` checks that `A`, `B`, and `C` are closed and that $dU=A\smile B$ and $dV=B\smile C$. It returns the cocycle $A\smile V+(-1)^{a+1}U\smile C$, where `a` is the arrow count of `A`. The result uses this chosen defining system. Different witnesses can give different cohomology classes.

At $p=3$, take $a=t$ and $w=v_1t+t^2$. Then $dw=a\smile a$ and the chosen triple value is $M=2v_1tu+4t^2u+tu^2$. The equality $M+\beta=d(v_1^2t)$ identifies its class with $-\beta$, of order three. [Brown–Peterson Massey cocycle example](examples/bp-massey.asl)

### Integral normalized BP Ext

`bpext[p, n, a, t];` computes the integral cohomology of the normalized finite-height BP cobar complex in arity `a` and internal topological degree `t`. It reports the free rank and the positive torsion invariant factors in divisibility order. Every monomial of a normalized cochain uses at least one coordinate from each arrow slot. The generators $v_i$ and $t_{j,i}$ both have degree $2(p^i-1)$. The differential is `bpD`.

`bporder[t](f);` takes a homogeneous closed normalized cochain and reports its exact additive class order, including one for a boundary and `infinite` for a class of infinite order. In positive arity, `bpNull[t](f)` returns a normalized cochain whose differential is the minimal positive multiple of `f` that is a boundary. After inverting $p$ the Hazewinkel coordinates identify the complex with that of a pair groupoid. Inserting the origin contracts it in positive internal degree. This certifies the rational ranks and makes the integral torsion $p$-primary. Smith reduction uses enough modular precision to find every pivot of the certified rank. The null cochain uses the origin contraction and an exact integral correction, and ASL verifies its differential with `bpD`. The calculation describes a fixed finite-height integral cobar complex.

At $p=2$ the arity-one group in degree eight is cyclic of order eight at height one and of order sixteen at height two. The height-two class $d(v_1^4+8v_1v_2)/16$ generates the latter group. At $p=3$ and height two in degree twelve, the arity-one and arity-two groups have orders nine and three. [Integral BP Ext example](examples/bp-ext.asl)

### First-class BP Ext classes and products

`BPExt(p, n, a, t)` is the discrete additive group of that cohomology. `bpextgen[p, n, a, t, j]` selects a generator in invariant-factor order, starting at zero. `bpExtClass[t](f)` sends a homogeneous closed normalized integral cochain to its exact class, and `bpExtRep(x)` gives a deterministic normalized representative. Projection of the representative returns `x`. The modular coordinate frame uses a modulus equal to the square of the largest invariant factor and checks the modulus bound. Integer multiplication `bpScale[k](x)`, addition, and `bpNeg` reduce torsion coordinates modulo their invariant factors.

`bpCup(x, y)` multiplies classes at a common prime and height by the Alexander–Whitney product on their representatives and projects to the sum of their arities and internal degrees. `BPExt(p, n, 0, 0)` is the free group $\mathbb Z$, with the generator one as the cup unit and checked signed integer coefficients. Finite groups can be enumerated within the usual limit. The free group gives an enumeration error when used as a domain for an operation that enumerates all its elements.

The example checks that $[t_1]^2$ is nonzero of order two, enumerates all four classes in `BPExt(2, 2, 2, 8)`, and identifies the signed Massey representative at $p=3$ with $-\beta$ in integral Ext. [BP Ext ring example](examples/bp-ext-ring.asl)

### The relative BP cobar complex

For heights `n < m`, `BPRel(p, n, m, a)` is degree `a` of the normalized height-`m` integral cobar modulo the embedded height-`n` cobar. The operation `bpRel[n](f)` takes a normalized integral upper cochain to the quotient. It retains the monomials containing a coordinate of index above `n`. The operation `bpRelLift(q)` reads that representative in the upper complex. The relative `bpD` differentiates the lift and projects the result. Addition, zero, and integer scaling act on the quotient.

`bprelext[p, n, m, a, t];` reports the free rank and torsion factors in internal degree `t`. Over $\mathbb Z[1/p]$ insertion of the origin preserves both heights and contracts the positive-degree quotient. Its degree-zero quotient vanishes. The same exact modular Smith reduction computes the integral relative groups. The operations `bporder[t]` and `bpNull[t]` give orders and minimal-multiple witnesses for closed relative cochains.

For $p=2$ and $n=1$, $m=2$, the relative class of $t_{1,2}$ has order two in degree six. The relative differential of $v_2$ is twice that cochain. [Relative BP cobar example](examples/bp-relative.asl)

### The BP height long exact sequence

`BPRelExt(p, n, m, a, t)` is the discrete additive group of relative classes. Use `bprelgen[p, n, m, a, t, j]` for the generator of invariant factor `j`, starting at zero, and `bpRelClass[t](q)` and `bpRelRep(x)` to project cycles and read their canonical representatives. The operations `add`, `zero`, `bpNeg`, and `bpScale` act on relative classes.

On absolute `BPExt` classes at height `m`, `bpRel[n]` gives the map to relative cohomology. The connecting map `bpHeightDelta(x)` differentiates the canonical upper lift, reindexes the resulting lower-height cochain, and projects to absolute `BPExt` at height `n` with arity increased by one. It also acts on closed raw relative cochains and returns the lower raw cochain. Together with `bpRaise[m]` these maps give the long exact sequence of the height inclusion.

At $p=2$, the connecting image of the relative class of $t_{1,2}$ in internal degree six is the nonzero height-one class $[v_1tu-tu^2]$. Its height-two image vanishes. In degree eight the order-sixteen class $\sigma$ projects to the relative order-two generator. The kernel of that projection is generated by $2\sigma$, the image of the order-eight height-one class. [BP height exact sequence example](examples/bp-height-les.asl)

### The integral BP height spectral sequence

In ambient height $N$, use the increasing filtration $F_h C_N=C_h$ of the normalized integral cobar. It is zero below height zero and constant from height $N$ onward. For $r\ge 1$, put $Z_r(h,a)=\{x\in F_h C_N^a : dx\in F_{h-r} C_N^{a+1}\}$ and $Z_0(h,a)=F_h C_N^a$. The page is the quotient of $Z_r(h,a)$ by $Z_{r-1}(h-1,a)+dZ_{r-1}(h+r-1,a-1)$. Page zero is the associated graded. Page one is the cohomology of consecutive height quotients, and a page-$r$ differential lowers height by $r$ and increases arity by one.

`bpHeightPage[p, N, r, h, a, t];` gives the free rank and torsion factors. The cycle lattice is saturated. To compute its quotient, ASL first computes the rational constraint ranks in the relative complexes, then uses sufficiently precise modular kernels of the constraints and Smith factors of the denominator. In positive internal degree the pages after page zero are finite $p$-primary groups. `bpHeightImage[p, N, r, h, a, t];` gives the order of the differential image, or the free rank when the image is infinite on page zero. Pages with $r>N$ have stabilized.

In degree six at $p=2$, the height-two relative class gives a $d_1$ of order two into height one. In degree eight, the two stable height pieces have orders eight and two, and their extension in the height-two integral group is cyclic of order sixteen. [Integral BP height page example](examples/bp-height-pages.asl)

The page type is `BPHeight(p, N, r, h, a, t)`. The generator `bpHeightgen[p, N, r, h, a, t, i]` follows the invariant-factor order for finite groups and the monomial order on page zero. These classes support `zero`, `add`, `bpNeg`, and `bpScale[k]`. Use `bpHeightClass[r, h, t](x)` on a homogeneous normalized cochain at ambient height `N`. It tests the height bound and the cycle condition. `bpHeightRep` gives an integral upper-height representative. The Smith generator residues are lifted to actual cycles by a relative cobar primitive when necessary.

`bpHeightD` gives the differential, with a zero group at negative target height. `bpHeightNext` requires a cycle for that differential and corrects the chosen representative by a lower-height cochain before projecting to the next page. `bpHeightPrev` projects a chosen representative to the preceding page. Applying `bpHeightNext` to that result recovers the original class. In degree eight at height two, the page-one generator has a nonzero raw differential in height one. Subtracting a lower-height cochain gives a closed representative with integral order sixteen. [Typed BP height pages example](examples/bp-height-classes.asl)

### Bockstein and algebraic Novikov pages

`bpBockdim[p, n, r, a, t, f];` and `bpNovdim[p, n, r, a, t, f];` give dimensions over $\mathbb F_p$ on page `r`, in normalized cobar degree `a`, internal degree `t`, and filtration `f`. The Bockstein calculation filters the integral normalized cobar by powers of $p$. The Novikov calculation uses powers of the invariant ideal $(p,v_1,\ldots,v_n)$. For a monomial of Novikov weight $w$, the sum of the source $v_i$ exponents, membership in $F^f$ means its coefficient is divisible by $p^{\max(0,f-w)}$. The Bockstein calculation uses $w=0$ for every monomial.

Page zero is $F^f/F^{f+1}$. If $Z_r^{a,f}$ consists of $x\in F^f C^a_t$ with $dx\in F^{f+r} C^{a+1}_t$, page `r` for `r ≥ 1` is the quotient of $Z_r^{a,f}$ by $Z_{r-1}^{a,f+1}+dZ_{r-1}^{a-1,f-r+1}$. The calculator imposes the congruences on an integral basis and reduces the denominator in cycle coordinates modulo $p$. Exact arithmetic uses checked signed 128-bit coefficients after removing common filtration powers when possible. These pages describe the stated finite-height complex.

In internal degree eight at $p=2$, the Bockstein starting in cobar degree zero and filtration zero survives to page four at height two and to page three at height one. In internal degree 16 the height-two Novikov $d_2$ from $(a,f)=(1,1)$ to $(2,3)$ has rank one. [Filtered Brown–Peterson page example](examples/bp-pages.asl)

The page types are `BPBock(p, n, r, a, t, f)` and `BPNov(p, n, r, a, t, f)`. Use `bpBockgen[p, n, r, a, t, f, i]` or `bpNovgen[p, n, r, a, t, f, i]` for generator `i`. They support `zero`, `add`, and `pscale` over $\mathbb F_p$. Use `bpBockClass[r, t, f](x)` or `bpNovClass[r, t, f](x)` to project an integral homogeneous normalized filtered cycle. The `Rep` operations give a chosen integral cycle representative, including its filtration factor. The `D` operations apply the normalized cobar differential and project in degree $(a+1,t,f+r)$ on the same page.

`bpBockNext` and `bpNovNext` solve for a corrected lift of a cycle to its next-page class. `bpBockPrev` and `bpNovPrev` project the chosen representative onto the preceding page. In the height-two Bockstein example, $v_1^4$ gives a page-three cycle. The corrected lift $v_1^4+8v_1v_2$ represents its page-four class, whose $d_4$ is nonzero. The implementation can remove common powers of $p$ for page differentials and transitions. A request for an absolute integral representative remains subject to the signed 128-bit coefficient bound.

`bpCup` multiplies page classes in a common family, page, prime, and height. It applies the Alexander–Whitney product to their corrected integral lifts and projects in degree $(a+a',t+t',f+f')$. On either family, $d_r$ satisfies the signed Leibniz rule, and `Next` preserves products of cycles. The computation keeps common powers of $p$ factored out. For example, the filtration-one class of $p$ shifts the nonzero height-two $d_4$ in weight eight from filtration $0\to4$ to filtration $1\to5$. At $p=3$, on the Bockstein page zero, $d_0(t\smile t^2)=-t\smile d_0(t^2)$ is nonzero. [Multiplicative Brown–Peterson page example](examples/bp-page-products.asl)

`bpRaise[m](x)` includes height `n` into height `m ≥ n` on integral and invariant-ideal cochains, Bockstein and Novikov page classes, and integral `BPExt` classes. It preserves each $v_i$ and each arrow coordinate $t_{j,i}$ while expanding the slots. It preserves the prime, ideal depth, internal degree, and filtration. The inclusions commute with cofaces, differentials, and cup products.

At $p=2$ the height-one class of $v_1^4$ supports $d_3$. Its height-two image survives $d_3$, and its height-one $d_3$ target maps to zero. The corrected height-two lift $g=v_1^4+8v_1v_2$ supports $d_4$. Integrally, the inclusion maps $\alpha=d(v_1^4)/8$ to $2\sigma$ in cohomology, where $\sigma=dg/16$. The representative identity is $2\sigma-\operatorname{bpRaise}(\alpha)=d(v_1v_2)$. [Brown–Peterson height inclusion example](examples/bp-height.asl)

### Full integral triple Massey cosets

`bpTriple(A, B, C)` takes three positive-cobar-degree homogeneous `BPExt` classes at the same prime and height. It is defined exactly when `bpCup(A, B)` and `bpCup(B, C)` vanish. `bpTripleDefined(A, B, C);` tests this condition. ASL lifts the adjacent products to exact integral primitives and applies the signed defining-system formula. It then adds the indeterminacy $A\,H^{b+c-1}+H^{a+b-1}\,C$ in the two corresponding internal degrees.

The result has type `BPTriple(p, n, a, t)`, where `a` is the target cobar degree $a_A+a_B+a_C-1$ and `t` is the sum of the internal degrees. `bpTripleBase(bracket)` gives its canonical base class. `bpTripleContains(bracket, x);` tests membership, and `bpTripleSize(bracket);` prints the number of classes as `1` or an exact power of `p`. Coset equality uses the integral Hermite lattice generated by the cyclic relations and the indeterminacy directions, together with the reduced base class. Membership uses this lattice without enumerating the ambient group. `bpMassey` gives the raw cocycle of the defining system supplied by the caller.

At $p=3$, the triple of the first generator in bidegree $(1,4)$ is the singleton $-\beta$. At $p=2$ and height two, take the outer classes $\eta$ in $(1,2)$ and a zero middle class in $(1,4)$. The target group in $(2,8)$ has two order-two generators. The bracket contains zero and the first generator. Neither the second generator nor its sum with the first belongs. With the nonzero middle class $\alpha_2=[v_1t+t^2]$ in $(1,4)$, the same outer classes give the affine coset $g_1+\langle g_0\rangle$. [Integral BP triple coset example](examples/bp-triple.asl)

### Integral higher BP defining systems

`bpHigherSystem[n](A1, ..., An, witnesses...)` checks one defining system on `n ≥ 3` positive-degree homogeneous `BPExt` classes at a common prime and height. Number vertices from zero to `n`. Supply integral normalized BP witnesses in increasing interval length and then increasing starting vertex. For `n=4` their order is $X_{02},X_{13},X_{24},X_{03},X_{14}$. A witness on $[i,j)$ has cobar degree $1+\sum_{i\le l<j}(a_l-1)$ and internal degree $\sum_{i\le l<j}t_l$.

The required differential of a proper witness is the sum of $X_{ik}\smile X_{kj}$ over $i<k<j$, with sign $(-1)^{\sum_{i\le l<k-1}(a_l-1)}$. The terminal sum is projected to `BPExt`. At arity three this gives the signed convention of `bpMassey`. Changing a proper witness can change the returned class. `bpChosenHigher(A1, ..., An)` builds the witnesses in the same order using deterministic integral primitives. It reports an obstruction to this particular choice when an intermediate sum represents a nonzero class.

For the integral class $\alpha=[t_1]$, the chosen repeated systems at $p=3,5,7,11$ have zero values before arity $p$ and class $-\beta$ at arity $p$, where $\beta$ is the two-step connecting class of $v_2$. The height-two example checks a fivefold result using explicitly supplied integral witnesses at $p=5$. [Higher BP defining system example](examples/bp-higher.asl)

### Divided Greek-letter connecting maps

`BPThick(p, n, a, h, j)` is the quotient of `BPI(p, n, a, h−1)` by $v_{h-1}^j$, for $2\le h\le n+1$ and $j\ge1$. Apply `bpThick[j](x)` to a cochain in that `BPI` ring. Its canonical polynomials retain powers of $v_{h-1}$ below `j` and use coefficients in $0,\ldots,p-1$. The generator $v_{h-1}$ is strictly invariant in the underlying `BPI`, so arithmetic, cup, cofaces, codegeneracies and differential descend. `bpThickLift` returns the chosen monomial lift in the underlying `BPI`. Reduction to a deeper invariant-ideal quotient is available once that quotient kills $v_{h-1}$.

`bpThickDelta(x)` differentiates this lift and divides exactly by $v_{h-1}^j$. A failed division detects a noncycle. The output is closed in `BPI(p, n, a+1, h−1)`. For `j=1` the operation agrees with `bpDelta` after reduction to `BPI` depth `h`. Changing a lift by $v_{h-1}^j y$ changes the quotient of its differential by the boundary $dy$.

At height two, apply `bpDelta(bpThickDelta(...))` to the class of $v_2^s$ in the thick quotient to obtain an integral representative of $\beta_{s/j}$. Its internal degree is $2s(p^2-1)-2j(p-1)$. The exact $v_1$ divisibility of $d(v_2^s)$ modulo $p$ is $p^{\nu_p(s)}$ for $s>0$. The example checks order two for $\beta_{2/2}$ at $p=2$ in weight eight and order three for $\beta_{3/3}$ at $p=3$ in weight 36. [Divided Greek-letter example](examples/bp-divided-greek.asl)

### Cohomology of the BP quotients

`BPIExt(p, n, a, h, t)` and `BPThickExt(p, n, a, h, j, t)` are the actual cohomology groups of the normalized `BPI` and `BPThick` cobar complexes in internal degree `t`. The calculator uses their homogeneous monomial bases and reduces cycles modulo boundaries over $\mathbb F_p$. `bpIExtdim[p, n, a, h, t];` and `bpThickExtdim[p, n, a, h, j, t];` report the dimensions. The corresponding `gen` forms append a generator index. `bpIExtClass[t]` and `bpThickExtClass[t]` project homogeneous normalized cocycles. Their `Rep` operations return chosen cocycle representatives. Both groups have `zero`, `add`, `pscale`, and the signed Alexander–Whitney `bpCup` product.

`bpThick[j]` and `bpReduce` also act on quotient cohomology. The `bpThickDelta` class operation goes from the thick quotient with parameter `h` to the underlying `BPIExt` depth `h−1`, raises `a` by one, and subtracts $2j(p^{h-1}-1)$ from `t`. On `BPIExt` at depth `h>1`, `bpDelta` lowers the depth by one and subtracts $2(p^{h-1}-1)$ from the internal degree. At depth one it maps to integral `BPExt` by division by $p$, which preserves the internal degree. These maps apply the connecting calculations to cocycle representatives and then project in the target. [First-class divided Greek-letter example](examples/bp-divided-classes.asl)

### The chromatic Bockstein for an invariant generator

`BPVPage(p,n,r,a,h,t,f)` is page `r` of the $v_h$-adic filtration of the normalized `BPI(p,n,a,h)` cobar, for `1 <= h <= n`. Filtration `f` consists of the cochains divisible by $v_h^f$. The associated graded in that filtration is the cobar of `BPI` at depth `h+1`, shifted in internal degree by $2f(p^h-1)$. In particular `E1` in that grade has the dimension of the corresponding `BPIExt` group. In any fixed internal degree the filtration is finite and the stable page is the associated graded of `BPIExt` at depth `h`.

Use `bpVdim[p,n,r,a,h,t,f];` for the dimension over $\mathbb F_p$ and `bpVgen[p,n,r,a,h,t,f,k]` for a basis class. `bpVClass[r,t,f](cochain)` requires a normalized homogeneous `BPI` cochain in filtration `f` whose differential lies in `f+r`. `bpVRep(x)` returns a chosen correcting lift, and `bpVD(x)` raises arity by one and filtration by `r`. `bpVNext(x)` projects a $d_r$-cycle to the next page, and `bpVPrev(x)` chooses a lift to the previous page. The pages support `zero`, `add`, `pscale`, and the signed `bpCup` product. The computation reduces sparse differential columns over $\mathbb F_p$ in the relevant thick quotients, subject to the usual work and storage bounds.

At height two and depth one, $d(v_2)=v_1t_1^p-v_1^pt_1$ modulo $p$. The example detects nonzero $d_2(v_2^2)$ at $p=2$ and $d_9(v_2^9)$ at $p=3$. It compares the $d_2$ target with the divided connecting cochain multiplied back by $v_1^2$. [Chromatic Bockstein example](examples/bp-vadic.asl)

## Transferred homotopy operations on Ext

`ainfty(a1, ..., ak)` gives the operation $m_k$ of the transferred minimal $A_\infty$ algebra on `Ext`. `alift(a1, ..., ak)` gives the component $F_k$ of its $A_\infty$ quasi-isomorphism to the normalized cobar algebra. Both take $k\ge 2$ Ext classes in a common algebra profile. With input bidegrees $(s_i,t_i)$, the results have bidegrees $(\sum s_i+2-k,\sum t_i)$ and $(\sum s_i+1-k,\sum t_i)$, respectively. They are Ext classes and cobar cochains over the input algebra. A negative result degree gives an error. The binary operation is `yoneda`, and the unit in bidegree $(0,0)$ is strict.

At fixed internal degree, let $P_s$ send each boundary of cobar length $s$ to its deterministic preimage by sparse elimination on the ordered words of length $s-1$. Write $i$ for `extRep`. For a cochain $c$ of length $s$, form the cycle $q=c+P_{s+1}(dc)$, let $p(c)$ be `extClass(q)`, and put $h(c)=P_s(q+i p(c))$. These maps give the contraction $d h+h d=1+i p$.

Starting with $F_1=i$, ASL computes $U_n=\sum_{j=1}^{n-1} F_j(a_1,\ldots,a_j)F_{n-j}(a_{j+1},\ldots,a_n)$ using cobar concatenation. It returns $m_n=p(U_n)$ and $F_n=h(U_n)$. The correction $P_{s+1}(dc)$ applies even when $U_n$ has nonzero differential. The calculation reuses lifts for each consecutive interval and caches the sparse splittings by bidegree. [Merkulov, *Strongly homotopy algebras of a Kähler manifold*](https://arxiv.org/abs/math/9809172)

The higher operations can depend on the chosen contraction. The Massey product collects classes from defining systems. In this transfer, $m_3(h_0,h_1,h_0)=h_1^2$ and $m_4(h_0,h_1,h_2,h_0)=h_2^2$. The adjacent product $h_2h_0$ obstructs a Massey defining system for the four inputs. The fourth lift is $[\xi_1\xi_3]+[\xi_1^2\xi_2^2]+[\xi_1^5\xi_2]$. Run `nix run path:. -- examples/ext-transfer.asl` to compute both operations and the lift. [Transfer example](examples/ext-transfer.asl)

### Signed transfer over odd primes

The same `ainfty` and `alift` forms accept `PExt` or `PExtAn` classes at one common prime and Steenrod algebra profile. They return an Ext class or cobar cochain of that profile in the bidegrees above. In stem grading $|a|=t-s$, let $\sigma$ raise degree by one. The signed binary map on suspended cochains is $b_2(\sigma a,\sigma b)=(-1)^{|a|}\sigma(ab)$. The cobar contraction closes a cochain by $q=c-P(dc)$, projects by $p(c)=[q]$, and uses $h(c)=P(q-i p(c))$. Its identity is $dh+hd=1-ip$.

On suspended inputs the recursion sums $b_2$ over split pairs of previously computed lifts. It applies $p$ to obtain the higher product and $-h$ to obtain the lift. The public operations undo suspension with the sign $(-1)^{\sum_{j=1}^n(n-j)|a_j|}$ for both outputs. Thus the binary product is `yoneda`, and the binary lift satisfies $dF_2(a,b)=i(m_2(a,b))-i(a)i(b)$.

For $x=[\xi_1]$ the first nonzero repeated higher product is $m_p(x,\ldots,x)=(-1)^{p(p-1)/2} b_{1,0}$. At $p=5$ the transfer gives $b_{1,0}$ and the explicit signed defining system in the example gives $-b_{1,0}$. The example also computes a nonzero mixed ternary product with $[\tau_0]$ at $p=3$. [Odd-prime transfer example](examples/odd-transfer.asl)

### Coherent restriction

`resHigher[n](a1, ..., ak)` computes the component $R_k$ of the $A_\infty$ morphism induced by restriction to $\mathcal A(n)$. Its inputs are ordered as for `ainfty` and all use the full algebra or one common $\mathcal A(m)$ with $m\geq n$. The result lies in `ExtAn(n, sum s + 1 - k, sum t)`, or `PExtAn(p, n, sum s + 1 - k, sum t)` for odd-primary inputs of a common prime `p`. The unary component is `resA[n]`, which preserves `yoneda`. For $k>1$, a unit input or an unchanged profile gives zero. ASL restricts canonical lifts of consecutive blocks and applies the target cobar contraction to the resulting tensors.

The example computes $R_3(h_2,h_0^2,h_1^3)=v_1^4$ over $\mathcal A(1)$. The individual restrictions of $h_2$ and $h_1^3$ both vanish. It computes additional nonzero components in arities two and four. [Coherent restriction example](examples/ainfty-restriction.asl)

Over an odd prime the target projection on suspended tensors uses the tensor contraction $H_0=\sum_j(ip)^{\otimes(j-1)}\otimes h\otimes 1^{\otimes(k-j)}$ with Koszul signs. The projection satisfies $P_1=p$ and $P_k=-P_{k-1}b_2H_0$, where $b_2$ acts on adjacent pairs with the suspended prefix sign. The public components undo suspension with the same sign as `ainfty`.

At $p=3$, write $a_0=[\tau_0]$, $h_{1,1}=[\xi_1^3]$, and $b_{1,0}$ for the degree-$(2,12)$ relation class. Both $R_2(a_0^3,h_{1,1})$ and $R_3(a_0^2,b_{1,0},a_0)$ give $v_1^3$ over $\mathcal A_3(1)$. At $p=5$, $R_2(a_0^5,[\xi_1^5])$ gives $v_1^5$ over $\mathcal A_5(1)$. [Signed coherent restriction example](examples/odd-coherent-restriction.asl)

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

### Bivariant cobar cochains

`mword[X, q, i, Y, p, j](u)` sends the homology basis element dual to `gen[X, q, i]` to `u` tensored with the element dual to `gen[Y, p, j]`. It sends every other basis element to zero. For `u : Cobar(s, d)` its type is `MapCobar(X, Y, s, p + d - q)`. These cochains support `add`, `zero`, and equality.

`eD(c)` applies the bivariant differential, raising `s` by one and preserving `t`. `ecomp(g, f)` also composes two `MapCobar` values. It pairs matching intermediate homology coordinates and concatenates the word from `f` followed by the word from `g`. The differential satisfies $d^2=0$ and the Leibniz identity over $\mathbb F_2$.

`eClass(c)` takes a closed cochain to its class in `MapExt` of the same endpoints and bidegree. `eRep(x)` gives the deterministic closed representative of a `MapExt` class, so `eClass(eRep(x)) == x`. The example checks the differential of the projective-plane identity term by term. [Bivariant cobar example](examples/map-cobar.asl)

## Algebraic Massey products of Adams maps

`eBracket(fn, ..., f1)` takes $n\ge3$ composable `MapExt` classes $f_i:X_{i-1}\to X_i$, with the outermost map first. It returns a membership predicate on `MapExt(X0, Xn, sum s + 2 - n, sum t)`. A negative result filtration gives an error. `eInBracket(fn, ..., f1, x)` checks membership for one class `x` without enumerating the target group.

Three inputs use a sparse affine calculation. For the chosen representatives $F,G,H$, the interpreter finds primitives $U,V$ with $dU=GF$ and $dV=HG$. If either has no primitive, the predicate is empty. Otherwise the answer is the affine coset of $[HU+VF]$ by the sum of the spaces of composites $h\circ u$ and $v\circ f$, where $u,v$ range over the `MapExt` groups in the endpoints and bidegrees of $U,V$.

With more inputs, the interpreter fills entries $A_{ij}:X_i\to X_j$ in order of increasing span $j-i$. The adjacent entries are canonical cobar representatives. For each other proper interval it solves $dA_{ij}=\sum_{i<k<j}A_{kj}A_{ik}$ in filtration $1-(j-i)+\sum_{i<k\le j}s_k$ and internal degree $\sum_{i<k\le j}t_k$. A negative filtration allows only the zero cochain. For an exact obstruction the interpreter tries one primitive plus each canonical representative of the corresponding cohomology group. At the omitted top corner it takes the class of $\sum_{0<k<n}A_{kn}A_{0k}$.

If a proper entry changes by a boundary $dB$, the unitriangular change of basis $1+B E_{ij}$ leaves the shorter entries fixed and adjusts the longer ones. The curvature of the defining system is supported in the top corner, where unitriangular conjugation acts trivially. Discarding the top entry that conjugation introduces changes the curvature by a boundary. Thus searching one representative of each cohomology class gives every top class. The search reports a resource-limit diagnostic if it cannot finish the requested calculation.

On `One` the algebraic triple Massey product $\langle h_0,h_1,h_0\rangle$ consists of the nonzero class $h_1^2$. A four-input example begins and ends at a two-point space and uses zero maps. The possible outcomes are the matrices of rank at most one, 10 out of 16, since the top term factors through a one-dimensional homology group. The construction uses Adams $E_2$ classes and cobar primitives. A topological Toda bracket uses stable maps and their nullhomotopies. [Triple example](examples/map-massey.asl) · [Four-input example](examples/map-massey-four.asl)

### Explicit defining systems

`eCurvature(fn, ..., f1){i < j => w, ...}` checks a chosen defining system and returns its top cochain. `eMassey(fn, ..., f1){i < j => w, ...}` returns the class of that same cochain. Both forms take at least two adjacent arrows, with the outermost first. Vertices are numbered from the source $X_0$ to the target $X_n$.

Adjacent entries may be `MapExt` classes or chosen closed `MapCobar` cochains. Classes use their `eRep` representatives. A proper nonadjacent entry $A_{ij}$ is a raw cochain from $X_i$ to $X_j$ of bidegree $(\sum_{i<k\le j}s_k+1-(j-i),\sum_{i<k\le j}t_k)$, where $(s_k,t_k)$ is the degree of the adjacent input $f_k$. An omitted entry is zero, including when that entry would have negative cobar length.

For every interval except $[0,n]$, ASL checks $dA_{ij}+\sum_{i<k<j}A_{kj}\circ A_{ik}=0$ and identifies the interval if the equation fails. The top entry is omitted from the defining system. Its curvature is the closed cochain $\sum_{0<k<n}A_{kn}\circ A_{0k}$ in bidegree $(\sum_k s_k+2-n,\sum_k t_k)$. The example gives explicit primitives for the two proper intervals of $\langle h_0,h_1,h_0\rangle$ and obtains the class $h_1^2$. [Explicit defining system](examples/map-system.asl)

## May pages for bivariant cobar morphisms

`MapMay(X, Y, r, s, t, f)` applies the algebraic May construction to `MapCobar(X, Y, s, t)`. A term that sends a degree-`q` source generator to a word tensored with a degree-`p` target generator has word degree $d=q+t-p$. Its filtration is $f=d-M$, where $M$ is the May weight of the word. Splitting a word and appending or prepending the endpoint coactions preserve the decreasing filtration. At each signed internal degree `t`, the filtration is finite and the stabilized page is the associated graded of `MapExt(X, Y, s, t)`.

`mapmaydim[X, Y, r, s, t, f];` prints a rank, and `mapmaygen[X, Y, r, s, t, f, i]` gives basis generator `i` starting at zero. Apply `mayClass[r, f](c)` to a `MapCobar` cochain satisfying the page cycle condition. Then `mayRep(x)` gives its canonical `MapCobar` lift, and `mayD(x)` raises `s` by one and `f` by `r`. The usual `zero`, `add`, and equality operate in each page component.

On a fixed page `r`, `ecomp(g, f)` composes classes with matching endpoints and adds their lengths, signed internal degrees, and filtrations. The page differentials satisfy the Leibniz identity. A finite-space map `u` induces a filtration-zero class on each page via `mayClass[r, 0](eRep(emap(u)))`.

The example detects the projective plane’s Steenrod coaction at both the source and the target on $E_0$. It checks functorial composition and shows two stabilized slices whose ranks sum to a signed `MapExt` rank. With both endpoints `One`, the pages recover scalar May. [Bivariant May example](examples/map-may.asl)

## Transferred operations on stable maps

For composable `MapExt` classes $f_i:X_{i-1}\to X_i$ of bidegrees $(s_i,t_i)$, pass the maps in composition order with the outermost map first. For $n\ge2$, `eHigher(f_n, ..., f_1)` gives the transferred $m_n$ in `MapExt` from $X_0$ to $X_n$ of bidegree $(\sum_i s_i+2-n,\sum_i t_i)$. `eLift(f_n, ..., f_1)` gives the corresponding $F_n$ in `MapCobar` of bidegree $(\sum_i s_i+1-n,\sum_i t_i)$. An intermediate lift of negative cobar length contributes zero, and a negative output length gives an error.

The transfer uses the ordered sparse cobar reduction. Write $P$ for its selected preimage of a boundary and $r$ for the chosen representative of a class. For a cochain $c$, $q=c+P(dc)$ is a cycle. The projection and homotopy are $p(c)=[q]$ and $h(c)=P(q+r(p(c)))$. Starting with $F_1=r$, the computation sums the compositions of the lifts over consecutive splits of the inputs, then applies $p$ for $m_n$ and $h$ for $F_n$. The higher products depend on this contraction.

The binary product is `ecomp` and `eid[X]` is a strict unit. On `One`, the operations recover the sphere transfer with the arguments reversed. The example composes a path of matrices through four finite spaces and obtains the square of $h_2$ at arity four. [Stable-map transfer example](examples/map-transfer.asl)

## Koszul complexes for finite exterior Steenrod subalgebras

A list such as `[0, 1]` selects the Milnor primitives $Q_0,Q_1$. Lists may be empty and have strictly increasing indices. For a list $I$, the bivariant `QKoszul(X, Y, I, s, t)` complex computes Ext over the sub-Hopf algebra $E(I)=\Lambda(Q_i\mid i\in I)$. A basis term maps a degree-$q$ homology generator of `X` to a degree-$p$ generator of `Y` times a polynomial monomial $\prod_i v_i^{a_i}$. Its length is $s=\sum_i a_i$, and its signed internal degree is $t=p-q+\sum_i a_i(2^{i+1}-1)$.

The differential `qD(F)` sums $v_i(Q_i^YF+FQ_i^X)$. It uses the duals of the cohomology actions obtained from $Q_0=\mathrm{Sq}^1$ and $Q_i=\mathrm{Sq}^{2^i}Q_{i-1}+Q_{i-1}\mathrm{Sq}^{2^i}$. The operation `qClass(c)` projects a closed cochain to `QExt(X, Y, I, s, t)` and `qRep(x)` returns a deterministic representative. Use `qextdim[X, Y, I, s, t];` to print the rank and `qextgen[X, Y, I, s, t, k]` for generator `k` starting at zero. The operation `ecomp(g, f)` composes matching endpoints at a common list $I$ and adds polynomial exponents, lengths, and signed degrees.

The restriction `qRes[I](z)` acts on mod-two scalar and bivariant cobar cochains and Ext classes. A scalar input uses `One` at both endpoints. A word whose letters are all single generators $\xi_{i+1}$ for $i\in I$ becomes the product of the corresponding $v_i$. All other words map to zero. This map commutes with the differential and composition. From a finite profile `A(n)`, the selected indices must lie between zero and `n`. Restriction from a Koszul type to a subset sends each omitted variable to zero.

The example computes the polynomial ring for `One`, checks composition of maps, and detects the $A(1)$ class $v_1^4$ after restriction to $Q_1$. [Koszul restriction example](examples/q-koszul.asl)


### Localizing bivariant Koszul Ext

For a nonempty list $I$ with $n=\max I$, `VKoszul(X, Y, I, s, t)` and `VExt(X, Y, I, s, t)` invert $v_n$. Both $s$ and $t$ are signed. The basis is still finite in each bidegree. If a term maps a degree-$q$ generator to a degree-$p$ generator and $D_i=2^{i+1}-1$, the lower exponents solve

$$\sum_{i\in I,\ i<n}(D_n-D_i)a_i=D_ns-(t+q-p),\qquad a_i\ge0.$$

Its highest exponent is $s-\sum_{i<n}a_i$ and may be negative. The differential and composition use the same endpoint actions as the polynomial Koszul complex. The operations `qD`, `qClass`, `qRep`, `add`, `zero`, and `ecomp` apply in the localized types. `vLocal(x)` embeds a polynomial cochain or sends a polynomial class to its localized class. `vShift[k](x)` multiplies by $v_n^k$ and shifts bidegree by $(k,kD_n)$. For one selected primitive, the localized dimensions are the corresponding bivariant Margolis dimensions with Laurent shifts.

Use `vextdim[X, Y, I, s, t];` and `vextgen[X, Y, I, s, t, k]` for ranks and basis generators. Over $I=[0,1]$ and with target $\mathbb{RP}^4$, the example verifies $[z_1]=(v_0/v_1)[z_3]\ne0$ and $(v_0/v_1)^2[z_3]=0$. Localizing the one-primitive complex for $I=[1]$ kills $[z_1]$. [Localized Koszul example](examples/v-koszul.asl)

## Declared graded Koszul modules

A `qmodule` declaration specifies a finite graded homology module over the exterior algebra on a selected list of Milnor primitives. The first block lists basis generators and their nonnegative homological degrees. The second block lists the images under selected primitives.

```text
qmodule Flash[0, 1] {
    a1:2, b1:1, a2:4, b2:3
} {
    0(a1)=>[b1], 0(a2)=>[b2], 1(a2)=>[b1]
};
```

Each image is the sum of its listed generators over $\mathbb F_2$. Repeated summands cancel and omitted rules give zero. ASL checks that every $Q_i$ drops degree by $2^{i+1}-1$, squares to zero, and commutes with the other declared primitives. Module names can appear as either endpoint of the ordinary or localized Koszul types. The list in a Koszul type must be a subset of the declared primitives at each module endpoint. A morphism can have a finite-space endpoint at one end and a declared module at the other.

The literal `qterm[X,q,i,Y,p,j,I]{a0,...}` gives a `QKoszul` cochain sending generator `i` of degree `q` to generator `j` of degree `p` times the monomial with the specified exponents in the variables of `I`. Coordinates start at zero within each degree. Module coordinates follow declaration order within that degree. The class `qid[X,I]` is the identity at endpoint `X` for that primitive selection.

The example computes a lightning-flash module with four pairs of generators. The differential has an explicit witness for the fourth power of $v_0$ annihilating the last class, and the third power is nonzero. [Declared module example](examples/q-module.asl)

### Bockstein pages of the localized Koszul category

`VBock(X, Y, [I], r, s, t, f)` is the page-$r$ group of the localized Koszul complex in signed bidegree $(s,t)$. The nonempty list `I` selects the Milnor primitives, with $n=\max(I)$. The filtration of a term $\prod_i v_i^{a_i}$ is $f=\sum_{i<n}a_i$. Only $v_n$ is inverted. In each signed bidegree, the positive differences between the internal degree of $v_n$ and those of the lower variables bound the possible lower exponents.

The differential $v_n Q_n+\sum_{i<n}v_iQ_i$ acts on coefficient maps. The first summand preserves filtration, so $E_1$ takes Margolis $Q_n$ homology of the coefficients at each monomial. The other summands raise filtration by one. Higher differentials use corrected lifts.

| Expression or statement | Result |
| --- | --- |
| `vbockdim[X, Y, [I], r, s, t, f];` | Print the rank of a page component. |
| `vbockgen[X, Y, [I], r, s, t, f, i]` | Select basis generator `i`, starting from zero. |
| `vbClass[r, f](a)` | Take a `VKoszul` cochain in $F^f$ with differential in $F^{f+r}$ to the page. |
| `vbRep(x)` | Return the chosen `VKoszul` lift. |
| `vbD(x)` | Evaluate $d_r$, raising `s` by one and `f` by `r`. |
| `vbNext(x)` | Advance a $d_r$-cycle to the next page. Incoming boundaries give zero. |
| `vbPrev(x)` | Choose a cycle on the preceding page using the canonical raw lift. |
| `ecomp(g, f)` | Compose classes on the same page with the same primitive selection, adding both degrees and filtration. |
| `vShift[k](x)` | Multiply by $v_n^k$, preserving filtration. |

The example over real projective four-space uses $I=[0,1]$. Here $Dz_2=v_0z_1$ and $Dz_4=v_1z_1+v_0z_3$. The lift $z_2+(v_0/v_1)z_4$ gives the nonzero differential $d_2(z_2)=(v_0^2/v_1)z_3$. Both components vanish on $E_3$. The same calculation after multiplying by $v_1^{-1}$ starts in negative cohomological degree. [Localized Bockstein example](examples/vbock-rp4.asl)

For a differential on a prescribed page $m\ge 1$, take a declared module with generators $a_j$ of degree $2j$ and $b_j$ of degree $2j-1$ for $1\le j\le m$. Its nonzero primitive actions are $Q_0a_j=b_j$ and $Q_1a_j=b_{j-1}$ for $j>1$. The lift of $a_1$ is $\sum_{j=1}^m(v_0/v_1)^{j-1}a_j$. Its differential cancels in every filtration below $m$ and leaves $v_0^m v_1^{1-m}b_m$. The program verifies that this term represents a nonzero $d_m$ from $(s,t,f)=(0,2,0)$ to $(1,2,m)$ and becomes zero on $E_{m+1}$. [Flash module generator](examples/flash-bockstein.awk)

```sh
awk -v m=97 -f examples/flash-bockstein.awk | nix run path:. -- -
```

The output is `VBock(One, Flash, [0, 1], 97, 1, 2, 97){g0}`. The generator checks the interpreter's source-size limit before emitting the program.

## Rank varieties over binary extension fields

`rvdim[X, Y, I, modulus, point];` takes a list of Milnor primitives $I$ and a coefficient vector $\lambda$ in the same order. After forgetting the internal grading, it computes Margolis homology for $Q(\lambda)=\sum_{i\in I}\lambda_iQ_i$ over the selected field. The output is the dimension of the space of linear maps from the source homology to the target homology. Both endpoints can be finite spaces or declared Koszul modules. The affine zero vector gives the product of the original homology dimensions.

The bits of `modulus` encode the coefficients of an irreducible polynomial over $\mathbb F_2$, from the constant term upward. For example, `3` gives $\mathbb F_2$, `7` gives $\mathbb F_4$ using $x^2+x+1$, and `11` gives $\mathbb F_8$ using $x^3+x+1$. For a modulus of degree $m$, each coefficient of `point` is an integer from zero to $2^m-1$ encoding a polynomial of degree below $m$. The interpreter reduces the action matrices over that field and uses $\dim H=\dim C-2\operatorname{rank}Q(\lambda)$.

`rvspectrum[X, Y, I, modulus];` reports each Hom dimension with its multiplicity among projective points. The first nonzero coordinate is normalized to one, so every projective point is counted once. A spectrum can enumerate up to 4,096 projective points. An empty selection gives an empty spectrum.

For the example module, the determinant of the three-by-three action block is $\lambda_1^3+\lambda_0^2\lambda_2$. In the chart $\lambda_2=1$, the equation parametrizes as $(\lambda_0,\lambda_1)=(r^3,r^2)$. At a nonzero point the inverse parameter is $r=\lambda_0/\lambda_1$. The cusp has Margolis dimension four, and the remaining $q-1$ affine curve points and the one point at infinity have dimension two over $\mathbb F_q$. The other $q^2$ projective points have dimension zero. [Rank variety example](examples/rank-variety.asl)

### Frobenius orbits and zeta coefficients

`rvorbit[modulus, point];` requires a nonzero projective vector. It normalizes the first nonzero coordinate to one and reports the orbit under coordinatewise squaring. The list starts at its lexicographically least vector, and its length is the degree of the closed point over $\mathbb F_2$.

`rvclosed[X, Y, I, m];` reports exact counts $C_1,\ldots,C_m$ of closed points, grouped by Margolis Hom dimension. The interpreter chooses an irreducible modulus for each $\mathbb F_{2^r}$ and evaluates the rank at the least representative of every Frobenius orbit of exact length $r$. Both endpoints accept finite spaces or declared Koszul modules. The report lists the strata observed through degree $m$, so $m=0$ or an empty selection gives an empty report. The enumeration uses the interpreter's work and storage budgets.

`rvzeta[X, Y, I, m];` reports the coefficients $a_0,\ldots,a_m$ of the Hasse–Weil zeta series for each observed stratum. The counts over the extension fields are $N_r=\sum_{d\mid r}d C_d$. The coefficients use $a_0=1$ and $n a_n=\sum_{r=1}^n N_r a_{n-r}$, with exact integer arithmetic.

The TwistedCusp action block has determinant $y^6+x^2y^3z+x^4z^2$. Its dimension-two stratum first appears over $\mathbb F_4$, where $[1:1:\omega]$ and $[1:1:\omega^2]$ form a degree-two orbit. Through degree six the stratum has closed-point counts `[0, 3, 0, 6, 0, 20]` and zeta coefficients `[1, 0, 3, 0, 12, 0, 48]`, agreeing with $(1-T^2)/(1-4T^2)$. [Rank arithmetic example](examples/rank-arithmetic.asl)

### Specialization to Margolis homology

`RVHom(X, Y, I, modulus, point)` is the discrete space of linear maps from the Margolis homology of `X` to that of `Y`. It uses the same field and coefficient encoding as `rvdim` and retains the entire affine coefficient vector. The maps support `zero`, `add`, equality, and bounded finite enumeration over the full extension field. Use `rvgen[X, Y, I, modulus, point, k]` for the $k$th basis map, indexed from zero, and `rvId[X, I, modulus, point]` for the identity. The operation `rvScale[c](f)` multiplies a map by the encoded field element `c`. Composition `ecomp(g, f)` requires matching endpoints at the same affine field point.

`rvSpec[modulus, point](z)` evaluates each Koszul variable $v_i$ at $\lambda_i$ and takes the resulting map on Margolis homology. It accepts `QKoszul` cochains and `QExt` classes. For a raw cochain, the evaluated map must commute with $Q(\lambda)$ at the given point. Koszul boundaries specialize to zero, and specialization preserves identities and composition. The same operation accepts `VKoszul` cochains and `VExt` classes when the coefficient for the inverted variable is nonzero. Multiplication by $v_n^k$ then specializes to multiplication by $\lambda_n^k$, including negative $k$.

Over $\mathbb F_4$ at $(1,3,1)$, the Cusp generators give $[b_2]=3[b_1]$ and $[b_3]=2[b_1]$, with $[b_1]$ nonzero. All three classes vanish at $(2,3,1)$. [Specialization example](examples/rv-specialization.asl)

### Frobenius, field extension, and trace on maps

`rvFrob[k](f)` sends each affine coordinate and each coefficient of an `RVHom` map through $x\mapsto x^{2^k}$. It preserves identities and composition and agrees with applying Frobenius to the point before specialization.

`rvExtend[modulus, root](f)` embeds the field of `f` in the field given by `modulus`. The value `root` is encoded in the larger field and must satisfy the source field's polynomial. The smaller degree must divide the larger degree. `rvDescend[modulus, root](f)` uses the indicated embedding of the smaller field into the field of `f`. It recovers the affine point and the map coefficients in the smaller field. Each must belong to the image of the embedding. The interpreter recovers preimages by binary elimination on the powers of `root`.

`rvTrace[modulus, root](f)` uses the same embedding as descent and requires the affine point to lie in the smaller field. It applies the relative field trace $a\mapsto\sum_{j=0}^{[L:K]-1}a^{|K|^j}$ to each map coefficient before expressing it in the smaller field. The basis coordinates agree under field extension because the ordered elimination uses the same leading coordinates and normalization commutes with the embedding. Thus trace is additive, satisfies the projection formula with maps extended from the smaller field, and composes down a tower of fields.

For modulus 7 embedded in modulus 19 using root 6, the cusp point $(1,3,1)$ becomes $(1,7,1)$. The trace of the encoded scalar 2 from modulus 19 to modulus 7 is 1, and the trace of the identity is zero because the extension has degree two. [Field descent example](examples/rv-descent.asl)

### Constructive Hilbert 90 and twisted descent

`rvInverse(f)` inverts an isomorphism of Margolis homology spaces and exchanges its endpoints. The unique endomorphism of the zero-dimensional space is its identity and its own inverse.

For an endomorphism $A$ defined over $L$ at an affine point in $K$, set $\sigma(x)=x^{|K|}$ and $d=[L:K]$. The operation `rvHilbert90[modulus, root](A)` specifies the inclusion of $K$ into $L$ as in `rvDescend`. It checks the ordered cocycle norm $A\sigma(A)\cdots\sigma^{d-1}(A)=I$ and returns an invertible frame $B$ with $A\sigma(B)=B$.

The interpreter applies $T=A\sigma$ to standard basis vectors to verify $T^d=I$. It tests the fixed vectors $\sum_{i=0}^{d-1}T^i(\alpha^k e_j)$ for the ambient polynomial generator $\alpha$ and $0\le k<d$. Independence of the field automorphisms makes the matrix $(\sigma^i(\alpha^k))_{i,k}$ invertible, so these vectors span over $L$. Keeping the first independent vectors gives a deterministic frame by elimination.

`rvTwistedDescend[modulus, root](A_Y, F, A_X)` takes a map between the two descent data. Both data use the same affine field point in $L$. The operation tests $A_Y\sigma(F)=F A_X$ by expressing $F$ as $B_Y^{-1}FB_X$ in the fixed frames and checking its coordinates lie in $K$. It returns the map with those coordinates over $K$. The inverse operation `rvTwistedExtend[modulus, root](A_Y, f, A_X)` takes a map at the corresponding point over $K$ and returns $B_Y\iota(f)B_X^{-1}$. The formulas allow rectangular maps and preserve identities and composition.

The example computes a frame for the unipotent two-by-two datum over $\mathbb F_4/\mathbb F_2$, which requires a power of the extension generator to obtain the second independent column. It also descends rectangular maps and verifies a degree-three norm whose first two factors have different products in the two orders. [Hilbert 90 example](examples/rv-hilbert90.asl)

### Shintani descent of Frobenius norms

For an invertible `RVHom` endomorphism $A$ over $L$ at an affine point in $K$, `rvFrobNorm[modulus, root](A)` returns the ordered norm $N=A\sigma(A)\cdots\sigma^{d-1}(A)$ over $L$, with $\sigma(x)=x^{|K|}$ and $d=[L:K]$. The modulus and root specify the inclusion of $K$ into $L$ as in `rvDescend`. Applying $A\sigma$ a total of $d$ times to each standard basis vector computes the norm columns.

The relation $A\sigma(N)A^{-1}=N$ makes the monic invariant factors of $N$ fixed by $\sigma$. The interpreter computes these factors by polynomial Smith reduction of $TI-N$ over $L$ and expresses their coefficients over $K$. `rvInvariantFactors[modulus, root](A);` prints the nonunit factors in divisibility order, with coefficients in ascending powers. `rvShintani[modulus, root](A)` returns their block companion matrix $C$ over $K$ in the same order. Each block has ones on the subdiagonal and the lower coefficients of its invariant factor in its final column in characteristic two.

`rvShintaniFrame[modulus, root](A)` constructs an invertible matrix $B$ over $L$ with $NB=B\iota(C)$. The reduction tracks the inverse Smith row operations on the cyclic generators by evaluating the polynomial operations at $N$. Their Krylov bases give the columns of $B$. Under $A\mapsto S^{-1}A\sigma(S)$ the norm becomes $S^{-1}NS$, so the canonical descended matrix stays the same.

Over $\mathbb F_4/\mathbb F_2$ the example descends a norm with distinct eigenvalues $\omega$ and $\omega^2$ to the companion matrix of $T^2+T+1$. The identity norm and a unipotent norm both have characteristic polynomial $T^2+1$, while their invariant factors are $(T+1,T+1)$ and $(T^2+1)$, respectively. [Shintani example](examples/rv-shintani.asl)

### Equal-characteristic Laurent isocrystals

`RVLaurent(X, Y, I, modulus, point)` consists of finite Laurent-polynomial matrices between the two Margolis homology spaces. The affine point and field encoding are the same as for `RVHom`. Use `rvU[k](f)` to form $u^k f$ from an `RVHom` map. Laurent maps support equality, `zero`, `add`, `rvScale`, and composition `ecomp` at a common field point. Exponents are checked signed 64-bit integers. The operator `rvFrob[k]` raises the field coefficients and point coordinates to $2^k$ and fixes $u$.

An invertible Laurent endomorphism over $L((u))$ gives a semilinear isocrystal over a specified subfield $K$. Its affine point must lie in $K$. The parameters `modulus, root` describe the inclusion $K\hookrightarrow L$. With $d=[L:K]$ and $\sigma=\operatorname{Frob}^{[K:\mathbb F_2]}$, `rvUNorm[modulus, root](A)` computes the ordered product $N=A\sigma(A)\cdots\sigma^{d-1}(A)$ over $L$. The determinant of $N$ must be nonzero, including when that determinant has more than one Laurent term.

`rvNewtonCoeff[modulus, root, j](A)` gives the coefficient $c_j(u)$ of $T^j$ in $\det(TI-N)$. It returns a scalar `RVLaurent` map from `One` to `One` at the descended point. The calculation uses division-free Berkowitz reduction over the full Laurent coefficient ring. The coefficient polynomials descend since $N$ and $\sigma(N)$ are conjugate.

`rvNewton[modulus, root](A);` takes the lower convex hull of $(j,\nu_u(c_j))$ for the nonzero characteristic coefficients. It reports the negative edge slopes divided by $d$ in increasing order, each with its horizontal width as multiplicity. The example computes slopes $0,1$ for an ordinary diagonal matrix, $1/2,1/2$ for a quadratic Frobenius cycle, and $1/3,1/3,1/3$ for a three-cycle. It checks cancellation in the determinant of a matrix with all entries of valuation zero. Composing with $u^{-2}I$ shifts the three-cycle slopes to $-5/3$. [Newton polygon example](examples/rv-newton.asl)

### Hodge exponents and framed affine Deligne–Lusztig loci

An invertible `RVLaurent` endomorphism $X$ represents the lattice $X L[[u]]^n$ in $L((u))^n$. `rvHodge(A);` reports the Hodge exponents of $A$ in increasing order. The sum of its first $k$ exponents is the least $u$-valuation of a nonzero $k$-by-$k$ minor. The interpreter uses exact Laurent arithmetic and valuation-pivoted Bareiss reduction to obtain those determinantal valuations.

Use `rvPosition(X, Y);` for the Hodge exponents of $X^{-1}Y$ when both frames have the same Margolis endpoint and exact affine point. The calculation applies the adjugate of $X$ by Cayley–Hamilton and subtracts $\nu_u(\det X)$ from each exponent. This allows all invertible frames over $L((u))$, including frames whose inverse has infinite power-series entries.

For $\sigma$ acting on coefficients and fixing $u$, `rvMuAt[modulus, root](A, X);` reports the position of $X L[[u]]^n$ relative to its image under $A\sigma$. The affine point must lie in the specified subfield. `rvInADLV[modulus, root, [mu0,...]](A, X)` evaluates to `S.hi` when this position equals the supplied weakly increasing coweight. `rvInADLVClosure` takes the same arguments and tests whether the determinant totals agree and each ascending prefix sum of the actual position is at least the corresponding sum for the supplied coweight.

For the two-cycle $A=\begin{pmatrix}0&u\\1&0\end{pmatrix}$, the Newton slopes are both $1/2$. The standard lattice has position $(0,1)$, and the lattice framed by $\operatorname{diag}(u^2,1)$ has position $(-1,2)$. The coweight $(0,1)$ accepts the standard lattice for exact membership. The closure with coweight $(-1,2)$ contains both. The example also changes the latter frame by a basis of determinant $1+u$ and computes a Hodge exponent arising from determinant cancellation. [Affine Deligne–Lusztig example](examples/rv-adlv.asl)

### Rank-two Bruhat–Tits shells

For an invertible rank-two `RVLaurent` endomorphism $A$, `rvTreeShell[modulus, root, r](A);` reports how many vertices of the Bruhat–Tits tree at distance $r$ from the standard lattice have each position under $A\sigma$. The affine point must lie in the specified subfield, with $\sigma$ the corresponding relative Frobenius on the field $L$ of $A$. Vertex positions are the ascending Hodge exponents of $X^{-1}A\sigma(X)$ for their lattice frames $X$. The difference of those exponents is the tree displacement.

Every vertex at radius $r>0$ has a unique normalized lattice with elementary divisors $(0,r)$. Its image modulo $u^r$ is a primitive line over $L[u]/(u^r)$. The interpreter uses the disjoint charts $(1,a)$, with $a$ any residue, and $(b,1)$, with $b$ divisible by $u$. Thus for $Q=|L|$ the sphere has $Q^r+Q^{r-1}$ vertices. Radius zero has just the standard lattice. If $Y=\operatorname{adj}(X)A\sigma(X)$, the first exponent is the minimum valuation of an entry of $Y$ minus $r$. Their sum is $\nu_u(\det A)$.

`rvTreeBall[modulus, root, R](A);` collects the same histogram for radii up to $R$. `rvTreeCount[modulus, root, R, [mu0,mu1]](A);` counts vertices in that ball with exact position $\mu$. `rvTreeClosureCount` uses the ascending-prefix closure condition from `rvInADLVClosure`. These counts run over the finite $L$-residue vertices in the chosen shells or ball. Each query can visit up to 4,096 vertices.

For $L=\mathbb F_4$ over $\mathbb F_2$ at radius three, relative Frobenius fixes 12 vertices. The other displacements $2$, $4$, and $6$ occur at 12, 24, and 32 vertices, respectively. The supersingular two-cycle exchanges the endpoints of the standard edge and gives positions $(-2,3)$ and $(-3,4)$ at 16 and 64 vertices. The example also checks invariance under a change of basis fixing the standard lattice and the effect of a negative Tate twist. [Bruhat–Tits tree example](examples/rv-tree.asl)

### Affine flags and Iwahori relative position

For rank $n$, the standard descending lattice chain has its first $i$ basis vectors multiplied by $u$ in $\Lambda_i$, for $0\le i<n$. It continues by $\Lambda_{i+n}=u\Lambda_i$. A lattice frame acts on the entire chain. The stabilizer is lower triangular modulo $u$.

`rvIw(A);` reports the Iwahori position of the standard chain and its image by an invertible `RVLaurent` endomorphism $A$. It uses the affine permutation window $[w(0),\ldots,w(n-1)]$ with distinct residues modulo $n$. If $A e_j=u^k e_i$ in the monomial representative, then $w(j)=i+nk$. In the Schur complement the pivot minimizes $n\nu_u(a_{ij})+i-j$ in the original row and column indices. This choice permits row and column elimination in the two lower-triangular Iwahori stabilizers. Fraction-free Laurent elimination obtains the pivot valuations exactly. `rvIwLength(A);` reports $\sum_{i<j}|\lfloor(w(j)-w(i))/n\rfloor|$.

`rvIwPosition(X, Y);` takes the window of $X^{-1}Y$, and `rvIwAt[modulus, root](A, X);` takes the window of $X^{-1}A\sigma(X)$. Both frames use the same Margolis endpoint and affine point. For the semilinear position this point must lie in the specified subfield. The adjugate calculation accounts for the scalar determinant of $X$ by translating the entire window.

`rvInIwADLV[modulus, root, window](A, X)` returns `S.hi` when the computed semilinear position equals the supplied window. `rvInIwADLVClosure` checks extended affine Bruhat order. The components $\sum_j(w(j)-j)/n$ must agree and $R_w(i,j)=\#\{a<j:w(a)\ge i\}$ must be at most the same rank for the supplied window at every $i,j$. For fixed residues of $i$ and $j$, these ranks are sums of truncated affine functions. The calculation compares their finitely many breakpoints.

The identity and the finite simple reflection both preserve the standard lattice. Their affine windows are $(0,1)$ and $(1,0)$. The affine simple reflection has window $(-1,2)$, and both reflections have length one. A coefficient Frobenius over $\mathbb F_4/\mathbb F_2$ changes the flag of an upper-triangular unit frame. [Affine flag example](examples/rv-iwahori.asl)

### Affine Kazhdan–Lusztig polynomials

`rvIwKL[[x0,...],[w0,...]];` computes the affine type-$A$ Kazhdan–Lusztig polynomial $P_{x,w}(q)$ for the supplied windows. The output lists coefficients in ascending powers of $q$, with `[0]` for zero. Equal windows give one. The Bruhat order and determinant component decide whether distinct windows are comparable. A common length-zero translation is removed before the calculation.

The interpreter constructs the normalized Hecke expansion $B_w=\sum_x h_{x,w}(v)H_x$. For $w=ys$ with a right ascent from $y$, it multiplies $B_y$ by $B_s=H_s+v$ and subtracts the lower canonical-basis terms $[v]h_{z,y}\,B_z$ for which $zs<z$. The desired coefficient is read from $h_{x,w}(v)=v^{\ell(w)-\ell(x)}P_{x,w}(v^{-2})$. The calculation uses exact signed integers and checks its storage and arithmetic bounds.

`rvIwSchubertKL[[w0,...]](A);` first takes the affine flag position $x$ of the supplied Laurent endomorphism $A$, then reports the affine Schubert intersection-cohomology polynomial $P_{x,w}$. The finite permutation $3412$ gives $P_{e,w}=1+q$. The affine rank-three window $(-4,3,4)$ also gives $1+q$ and has a reduced expression using all three affine simple reflections. [Affine Schubert example](examples/rv-iwkl.asl)

### The generic extended affine Iwahori–Hecke algebra

`RVIwHecke(n)` consists of finite sums $\sum_w c_w(q)T_w$ with signed integer polynomial coefficients and affine permutation windows of rank $n$. Use `rvHeckeT[[w0,...]]` for $T_w$, `rvHeckeId[n]` for $T_e$, and `rvHeckeQ[n]` for $qT_e$. The operations `zero`, `add`, `neg`, and `rvHeckeMul` act on these algebra elements.

The right multiplication rule is $T_xT_s=T_{xs}$ on an ascent and $T_xT_s=qT_{xs}+(q-1)T_x$ on a descent. For a right length-zero rotation $\omega^c$, the new window is $x(j+c)$ using the period rule $x(j+n)=x(j)+n$. The interpreter applies this rotation and then the right simple steps of a reduced expression.

`rvHeckeKL[[w0,...]]` constructs $B_w=\sum_{x\le w}P_{x,w}(q)T_x$ from one cached canonical expansion. This is $q^{\ell(w)/2}$ times the normalized $C'_w$. The statement `rvHeckeCoeff[[x0,...]](h);` prints the coefficient of $T_x$ in ascending powers of $q$. The statement `rvHeckeAt[q,[x0,...]](h);` evaluates the coefficient at the supplied integer.

For finite $\mathrm{GL}_3$ flags over $\mathbb F_q$, these evaluated structure coefficients count intermediate complete flags with the prescribed relative positions. The example computes $T_s^2=qT_e+(q-1)T_s$ and verifies the braid relation and $B_sB_tB_s=B_{sts}+qB_s$ for the two adjacent finite generators. It also checks the noncentral rotation and displays a cubic coefficient in $T_s^4$. [Affine Hecke example](examples/rv-hecke.asl)

### The spherical Hecke algebra and Hall counts

`RVSphHecke(n)` is the spherical Hecke algebra of $\mathrm{GL}_n(L((u)))$ with polynomial coefficients in $q=|L|$ and the maximal compact subgroup $K=\mathrm{GL}_n(L[[u]])$ assigned volume one. Use `rvSphBasis[[lambda0,...]]` for $S_\lambda$, with a weakly increasing signed coweight. The identity `rvSphId[n]` is $S_0$, and `rvSphQ[n]` is $qS_0$. The elements support `zero`, `add`, `neg`, and `rvSphMul(a,b)`.

The map `rvSphIwahori` sends $S_\lambda$ to the sum of the distinct $T_w$ in the double coset $S_n t^\lambda S_n$. Its windows are obtained by assigning a distinct permutation of the exponents $\lambda$ and a permutation of the row residues. If $i$ is this map, it satisfies $i(a)i(b)=P_n(q)i(ab)$ with $P_n(q)=\prod_{j=1}^n(1+q+\cdots+q^{j-1})$. The interpreter checks that each orbit in the raw product is complete with constant coefficient, then divides that coefficient exactly by $P_n(q)$.

`rvSphCoeff[[nu0,...]](h);` prints the coefficient of $S_\nu$ in ascending powers of $q$. `rvSphAt[q,[nu0,...]](h);` evaluates it at the supplied integer. For nonnegative coweights these values count intermediate lattices with prescribed quotient types, equivalently submodules with prescribed submodule and quotient partitions in a finite nilpotent Jordan module. In rank two the example computes $S_{(0,1)}^2=S_{(0,2)}+(q+1)S_{(1,1)}$. In rank three it gives $1+q+q^2$ at $(1,1,1)$ in $S_{(0,1,1)}S_{(0,0,1)}$, and $q^2+q$ at $(0,2,2)$ in $S_{(0,0,2)}^2$. [Spherical Hecke example](examples/rv-spherical-hecke.asl)

### Normalized Satake transform

`RVSatake(n)` is the ring of symmetric Laurent polynomials in $X_1,\ldots,X_n$ with coefficients in $\mathbb Z[v^{\pm1}]$. Use `rvSatV[n,k]` for the scalar $v^k$. The elements support `zero`, `add`, `neg`, equality, and `rvSatakeMul`.

`rvSatake(h)` takes a spherical Hecke element to this ring by sending $q$ to $v^2$. For a basis element $S_\lambda$, write the coweight in descending order, take $c=\min\lambda$, and put $\nu=\lambda-c(1,\ldots,1)$. The image is
$$
(X_1\cdots X_n)^c\,v^{\sum_i(n+1-2i)\nu_i}P_\nu(X_1,\ldots,X_n;v^{-2}).
$$
The interpreter expands the Hall–Littlewood $P$ polynomial by interlacing partitions. At each step its branching multiplier is $\prod_{a>0,\ m_a(\mu)>m_a(\nu)}(1-t^{m_a(\mu)})$ with $t=v^{-2}$. Separating the determinant power gives an exact Laurent image for negative coweights.

`rvSatakeCoeff[[a0,...]](f);` reports the coefficient of $X_1^{a_0}\cdots X_n^{a_{n-1}}$ as a map from signed powers of $v$ to integer coefficients. `rvSatakeAt[v,[x0,...]](f);` substitutes integers for $v$ and the $X_i$ and prints a reduced exact rational. A base of zero requires its Laurent powers to be nonnegative.

In rank two the images of $S_{(0,1)}$, $S_{(0,2)}$, and $S_{(1,1)}$ are $v(X_1+X_2)$, $v^2(X_1^2+X_2^2)+(v^2-1)X_1X_2$, and $X_1X_2$. The example checks multiplication in ranks two and three and evaluates a negatively translated coweight. [Satake example](examples/rv-satake.asl)

### Symbolic rank equations

`RVPoly(I, d)` is the discrete space of homogeneous polynomials of degree `d` over $\mathbb F_2$, with one variable $x_i$ for each selected primitive. The literal `rvmon[I]{a0,...}` has one exponent per primitive and infers its degree from their sum. Use `zero`, `add`, and `rvMul` for the polynomial operations. `rvHasse[a0,...](p)` extracts the indicated coefficient of $p(x+u)$ in the formal variables $u$. `rvVanish[modulus, point](p)` evaluates at a field point encoded the same way as `rvdim`, returning `S.hi` when the value is zero.

`rvMinor[X, I, rows, columns]` computes a minor of the universal operator $\sum_i x_iQ_i$ on $H_*(X)$, with the grading forgotten. Each coordinate is a `(degree, index)` pair, with index starting at zero in that degree. The two lists must have the same size, with distinct coordinates within each list. The empty minor is one. `rvminors[X, I, k];` prints the reduced linear basis of the span of the size-`k` minors. This basis generates the determinantal ideal for rank less than `k`, and `rvInMinors[X, I, k](p)` tests membership in its degree-`k` component.

For the Cusp module the size-three minors span $x_1^3+x_0^2x_2$. Its first Hasse derivatives are $0$, $x_1^2$, and $x_0^2$, giving the unique projective singular point $[0:0:1]$. The second Hasse derivative in $x_0$ is $x_2$. The same equation detects rank loss at points over both $\mathbb F_4$ and $\mathbb F_8$. [Symbolic rank equations example](examples/rank-equations.asl)

### Homogeneous rank ideals

`rvideal J[I] { f0, f1, ... };` declares an ideal of the polynomial ring over $\mathbb F_2$ with variables indexed by `I`. Its homogeneous generators may have different degrees. `rvfitting J[X, I, k];` declares the ideal generated by the canonical linear basis of size-`k` minors, in the order reported by `rvminors`.

`rvgroebner[J];` prints the reduced monic Gröbner basis in increasing graded lexicographic order, with earlier variables dominant. `rvNormal[J](f)` reduces the homogeneous polynomial `f` in this order. `rvInIdeal[J](f)` returns `S.hi` exactly when the residue is zero. For a member `f`, `rvWitness[J, j](f)` gives the coefficient of original generator `j` in a deterministic expression $f=\sum_j c_j f_j$. Indices start at zero. The coefficient has degree $\deg(f)-\deg(f_j)$, and the operation requires this degree to be nonnegative.

`rvhilbert[J, d];` counts the degree-`d` monomials outside the initial monomial ideal. This is the dimension of the corresponding graded quotient, computed by inclusion–exclusion over the distinct least common multiples of the initial generators. For $(xy+z^2, x^2+yz)$ the Gröbner basis includes the cubic $xz^2+y^2z$ with witness coefficients $x$ and $y$. The degree-$d$ quotient has dimension four for $d\ge2$. For the Cusp size-three ideal, the quotient has dimension $3d$ for $d\ge1$, so the example computes 3,000 at degree 1,000. [Gröbner and fitting ideal example](examples/rank-ideals.asl)

#### Intersection and saturation

| Declaration | Ideal assigned to `K` |
| --- | --- |
| `rvintersect K[J, H];` | The intersection $J\cap H$, for ideals in the same ring. |
| `rvcolon K[J](f);` | The colon $J:f=\{g\mid fg\in J\}$. |
| `rvsaturate K[J](f);` | The saturation $J:f^\infty=\{g\mid f^n g\in J\text{ for some }n\}$. |
| `rvprojective K[J];` | The saturation by the irrelevant ideal, computed as $\bigcap_i(J:x_i^\infty)$. |

The polynomial `f` must be homogeneous in the variables of `J`. Zero gives the unit ideal for both colon and saturation, and one gives `J`. In the ring with no variables, projective saturation gives the unit ideal. The original generators for witness indices of each derived ideal are its reduced Gröbner basis in the order printed by `rvgroebner`.

Projective saturation sends $(x^2,xy,xz,y^2)$ to $(x,y^2)$, a doubled point with Hilbert function two in positive degrees. The singular scheme of the cuspidal cubic $y^3+x^2z$ has ideal $(x^2,y^2)$ and eventual Hilbert function four. Both nilpotent coordinates survive as classes at the single projective point $[0:0:1]$. [Intersection and saturation example](examples/rank-schemes.asl)

#### Rees algebras and normal cones

`rvrees Blow[J, [u0,...]] { h0, ... };` presents the Rees algebra $R[HT]\subset R[T]$ for $R=\mathbb F_2[x_i]/J$ and $H=(h_j)\subset R$. Each new index names a coordinate $U_j$ for $h_jT$. Use one new index per listed generator, in strictly increasing order above the indices of `J`. The generators are `RVPoly` values in the original variables of one common positive degree $d$. An empty generator list leaves $R$, and a zero generator of that degree gives the relation $U_j=0$.

The defining relations come from eliminating $T$ from $(J,U_j+Th_j)$. They are bihomogeneous for degrees $(1,0)$ on $x_i$ and $(d,1)$ on $U_j$. ASL uses their ordinary homogeneous grading with every coordinate of degree one. `rvnormalcone Cone[J, [u0,...]] { h0, ... };` additionally kills the images of the $h_j$ in degree zero of the Rees grading, giving $\operatorname{gr}_H(R)$. Both declarations give named ideals for operations such as `rvgroebner`, `rvWitness`, and `rvbetti`.

For the cuspidal cubic and $H=(x,y)$, the normal cone in the chart $z\ne0$ has equations $(x,y,U^2)$, with $U$ nonzero. For $(x^2,xy,y^2)$ in the plane, the Rees equations include $AC+B^2$. This equation remains in the exceptional fiber after the two linear syzygies vanish at the origin. [Rees algebra and normal cone example](examples/rank-blowups.asl)

#### The Gauss graph and projective dual

`rvgauss G[[u0,...]](f);` defines the closed graph of the Gauss map of a hypersurface. The polynomial `f` must be homogeneous, nonzero, and of positive degree. Supply one strictly increasing list of fresh dual indices, one per original variable in order. The graph uses the union of the two variable lists. The construction eliminates $T$ from $(f,U_i+T\partial_i f)$ using the first Hasse partials, then saturates separately by the ideals of all base coordinates and all dual coordinates.

`rvdual D[[u0,...]](f);` eliminates the base coordinates from the graph and saturates in the dual coordinates. It gives the scheme-theoretic image over $\mathbb F_2$. A linear polynomial with a nonempty projective hypersurface gives the point of its constant gradient. When every first partial vanishes, the smooth Gauss graph and dual image are empty, and the resulting ideals are the unit ideals.

`rvGenerator[J, j]` is the polynomial at index `j`, starting at zero, in the reduced Gröbner basis printed by `rvgroebner[J]`. It can be used directly in another declaration. For the cusp $y^3+x^2z=0$, the Gauss map on its normalization is $[s:t]\mapsto[0:t^2:s^2]$. Its dual is the line $U=0$ and its bidual is the point $[1:0:0]$. For the smooth conic $xz+y^2=0$, the dual is the pencil of lines through the nucleus $[0:1:0]$. The smooth Fermat cubic has Frobenius Gauss map $[x:y:z]\mapsto[x^2:y^2:z^2]$ and dual equation $U^3+V^3+W^3=0$. [Gauss graph and duality example](examples/rank-duality.asl)

#### Hasse orders and osculating curves

`rvcurve C[[s_idx,t_idx]] { f0, ... };` declares a map from the projective line by binary forms of a common positive degree $d$. The forms must have no common projective basepoint over the algebraic closure. ASL checks the gcd in the affine chart $s=1$ and the values at $[0:1]$.

`rvorders[C];` prints the least increasing orders $\epsilon_j$ for which the Hasse derivative rows of $f_i(1,T)$ are independent over $\mathbb F_2(T)$. The length of the sequence is the coordinate rank. Lucas parity restricts the search to submasks of the exponents present in the forms. `rvOscMinor[C, k, [i0,...,ik]]` takes the determinant in the chosen coordinate positions for the rows $\epsilon_0,\ldots,\epsilon_k$. It returns the binary form of degree $(k+1)d-2\sum_{j=0}^k\epsilon_j$, including its ramification factors. Coordinate positions start at zero.

`rvosculating Image[C,k,[u0,...]];` gives the homogeneous ideal of the projective image of those Plücker coordinates. Supply one fresh index per increasing $(k+1)$-tuple of coordinate positions, in lexicographic order of the tuples. Exact substitution and elimination compute the relations, followed by projective saturation. For $q=2^e$ and forms $[s^{q+1},s^qt,st^q,t^{q+1}]$, the order sequence is $[0,1,q,q+1]$. The tangent coordinates are $[s^{2q},0,s^qt^q,s^qt^q,0,t^{2q}]$ and their image is a conic. The four plane coordinates give the original parametrization. [Hasse orders and osculating images example](examples/rank-osculation.asl)

#### Frobenius orders and local ramification

`rvFrobOrders[C, e];` prints the Stöhr–Voloch Frobenius orders for $q=2^e$ with $e>0$. Use the lexicographically first independent subset of the coordinate forms. For $r+1$ such forms, prepend their $q$-th powers as the first row and choose the least $r$ independent Hasse orders $\nu_0,\ldots,\nu_{r-1}$. `rvFrobWronskian[C, e]` returns the determinant as a binary form of degree $(q+r)d-2\sum_j\nu_j$, including its ramification factors.

`rvOrdAt[P,[s,t]](F);` prints the vanishing multiplicity of the binary form `F` at the projective point $[s:t]$ over the binary field with modulus `P`. The field coordinates use the same encoding as `rvVanish`. This includes the point at infinity. The zero form has order `infinite` at every valid point. The affine calculation uses exact local substitution and splits the support at powers of two, so a sparse multiple root can have large order without a search through each lower order.

The Ballico–Hefez parametrization $[s^5:s^4t+st^4:t^5]$ has ordinary Hasse orders $[0,1,4]$. For $q=4$ its Frobenius orders are $[0,4]$ and its Wronskian is $(st^4+s^4t)(st^{16}+s^{16}t)$. It vanishes twice at each point of $\mathbb P^1(\mathbb F_4)$ and once at each other point of $\mathbb P^1(\mathbb F_{16})$. For $q=2$ or $16$ the Frobenius orders are $[0,1]$. [Frobenius order and ramification example](examples/rank-frobenius.asl)

#### Artin–Schreier covers

`rvascover C[[s_idx,t_idx]] { F,G };` defines the smooth projective normalization of $y^2+y=F/G$. The binary forms must have the same positive degree $d$ and no common geometric zero. Each geometric pole of $F/G$ must have odd order. If $q(T)=G(1,T)$, the infinity order is $d-\deg q$. The quotient $q/\gcd(q,q')$ gives the finite odd-pole radical. Square-free layers of the derivative gcd check that this radical contains every finite pole.

For $n$ geometric poles, `rvASGenus[C];` prints $(d+n)/2-1$ and `rvASPrank[C];` prints $n-1$. `rvASPoints[C,m];` counts points over $\mathbb F_{2^m}$. A rational pole contributes one point. At every other rational parameter, the fiber has two points when the absolute trace of $F/G$ is zero and zero points when it is one.

`rvASLPolynomial[C];` prints the signed coefficient list of the zeta numerator. It computes the point counts through extension degree equal to the genus, then uses Newton identities and the functional equation. `rvASNewton[C];` gives the slopes and horizontal multiplicities of the lower hull of $(i,v_2(c_i))$ for its nonzero coefficients. The genus-zero numerator is `[1]`. Point queries may use higher extension degrees, subject to the field enumeration and computation budgets.

The covers $y^2+y=t^3$ and $y^2+y=t+1/t$ have genus one and $L$-polynomials $1+2T^2$ and $1+T+2T^2$, respectively. The genus-two cover $y^2+y=t^5$ has $L=1+4T^4$. The denominator $t^2+t+1$ of the fourth example has two geometric poles and gives $p$-rank one even though neither pole is rational over $\mathbb F_2$. [Artin–Schreier cover example](examples/rank-artin-schreier.asl)

#### Cartier on holomorphic differentials

For a declared cover with finite denominator $q(x)$ and pole radical $R(x)$, put $Q(x)=\sqrt{q(x)R(x)}$. The differentials $\omega_j=x^j\,dx/Q(x)$ for $0\le j<g$ form a basis of the holomorphic differentials. `ASDiff(C, P)` is this vector space over the binary field with modulus `P`. Use `rvASOmega[C,P,j]` for $\omega_j$ and `rvASDiffScale[c](v)` to scale coefficients using the field encoding of `c`.

`rvASCartier(v)` applies Cartier. On $\omega_j$ the coefficient of $\omega_i$ is the coefficient of $x^{2i+1-j}$ in $Q$. The operator acts on scalar coefficients by square root, the inverse of absolute Frobenius. `rvASCartierRank[C,k];` prints the rank of its $k$-th power, and `rvASaNumber[C];` prints $g$ minus the rank of Cartier. The rank for $k=0$ is $g$ and the eventual rank is the $p$-rank.

For $y^2+y=x^5$, Cartier sends $\omega_1$ to $\omega_0$ and has square zero. The $a$-number is one. For $y^2+y=x^9$ the ranks are $4,2,1,0$ and the $a$-number is two. For $y^2+y=x+1/x$, Cartier fixes $\omega_0$. [Cartier operator example](examples/rank-cartier.asl)

#### Cyclic quartic Artin–Schreier–Witt covers

`rvwittcover W[[x_idx]] { f,g };` defines the smooth projective curve of $y^2+y=f(x)$ and $z^2+z=g(x)+f(x)y$. Use unary `RVPoly` monomials or binary homogeneous forms containing `x_idx`. Setting the other variable to one gives an ordinary polynomial in $x$. The degree $m$ of $f$ must be positive and odd. The second polynomial is zero, constant, or of odd positive degree $n$, with $n=0$ for zero or constant $g$.

The deck generator $(y,z)\mapsto(y+1,z+y)$ has order four. The single pole at infinity is totally ramified. `rvWittBreaks[W];` prints the upper jumps $m,u$ with $u=\max(2m,n)$, and `rvWittGenus[W];` prints $(m+2u-3)/2$.

`rvWittPoints[W,k];` counts the one point at infinity and four points for each affine parameter whose length-two Witt trace is zero. The trace uses componentwise Frobenius and Witt addition $(a,b)+(c,d)=(a+c,b+d+ac)$. `rvWittLPolynomial[W];` uses the point counts to recover the integral zeta numerator, and `rvWittNewton[W];` prints its dyadic slopes and multiplicities.

`rvWittCharacterL[W];` gives the degree-$(u-1)$ Artin factor for $\chi(a,b)=i^{a+2b}$ on $W_2(\mathbb F_2)$. Each coefficient pair `(a, b)` represents $a+bi$. The quartic numerator is the product of this factor, its complex conjugate, and the $L$-polynomial of the quadratic subcover $y^2+y=f(x)$. For $f=x$ and $g=0$ the character factor is $1+(1+i)T$ and the quartic numerator is $1+2T+2T^2$. The point counts at higher extension degrees in the example give independent checks on these polynomials. [Artin–Schreier–Witt cover example](examples/rank-witt.asl)

#### Quartic Witt differentials and the Weierstrass semigroup

Let $m,u$ be the upper jumps of a declared quartic Witt cover and $b_2=2u-m$ its second lower jump. The quadratic subfield has $\operatorname{ord}_\infty(x)=-2$ and $\operatorname{ord}_\infty(y)=-m$. Rewrite the second equation as $w^2+w=A(x)+B(x)y$, starting with $w=z$. Whenever the leading term is an even pole $2a$ from $A$, add $x^{a/2}$ to $w$ if $a$ is even, or $x^{(a-m)/2}y$ if $a$ is odd. The final leading pole in that subfield is odd and has order $b_2$.

The pole orders of $x,y,w$ on the quartic curve are $4,2m,b_2$. The four monomials $1,y,w,yw$ give a free basis over $\mathbb F_2[x]$ with different pole orders modulo four, so the Weierstrass semigroup is $\langle4,2m,b_2\rangle$. Use `rvWittGaps[W];` to list its gaps. Since $\operatorname{div}(dx)=(2g-2)\infty$, `rvWittBasis[W];` prints the triples $(a,b,c)$ for holomorphic $x^a y^b w^c\,dx$ in increasing order of $4a+2mb+b_2c$, with $b,c\in\{0,1\}$.

`WittDiff(W,P)` is the space of holomorphic differentials over the binary field with modulus `P`. `rvWittOmega[W,P,j]` makes basis element $j$, numbered from zero, and `rvWittDiffScale[c](v)` scales its coordinates. `rvWittCartier(v)` takes the square root of each coefficient and applies Cartier to the basis elements. For $H\,dx$, it computes $D_xH$ in the normal form and recovers its square root by removing successive leading poles. `rvWittCartierRank[W,k];` prints the rank of the $k$-th iterate and `rvWittaNumber[W];` prints $g$ minus its rank for $k=1$. A single ramified pole gives $p$-rank zero.

For $f=x^3$ and $g(x)=0$, the genus is six and the Cartier ranks are $6,4,3,2,1,0$. The basis contains a Cartier chain from $w\,dx$ through $xy\,dx$, $y\,dx$, $x\,dx$, and $dx$. [Witt holomorphic differential example](examples/rank-witt-cartier.asl)

#### One-pole Witt towers and cyclotomic characters

`rvwitttower W[[x_idx]] { f0,...,f_{s-1} };` defines a cyclic degree-$2^s$ cover of the projective line for $1\le s\le8$. Its equation is $F(\mathbf y)-\mathbf y=(f_0,\ldots,f_{s-1})$ in the length-$s$ Witt ring. Each coordinate uses the univariate or binary-homogenized syntax of `rvwittcover`. The first has positive odd degree. Later coordinates may be zero, constant, or of positive odd degree. Infinity is its only pole and is totally ramified.

`rvWittTowerBreaks[W];` prints the upper jumps $u_1=\deg f_0$ and $u_j=\max(2u_{j-1},\deg f_{j-1})$, with zero or constant degree taken as zero. `rvWittTowerGenus[W];` gives the genus from $2g-2=-2\cdot 2^s+\sum_{j=1}^s 2^{j-1}(u_j+1)$. For extension degree $e$, `rvWittTowerPoints[W,e];` counts the one point at infinity and $2^s$ points at every finite parameter with zero Witt trace. The calculation lifts finite-field elements to Teichmüller representatives in the unramified ring modulo $2^s$ and sums their Frobenius conjugates.

`rvWittTowerCharacterL[W,j];` prints the primitive Artin $L$-factor of conductor $2^j$ and degree $u_j-1$. Each signed array lists one coefficient in the power basis $1,\zeta,\ldots,\zeta^{2^{j-1}-1}$ for $\zeta=\zeta_{2^j}$. The character sends a Witt trace $a$ in $\mathbb Z/2^j$ to $\zeta^a$. `rvWittTowerLPolynomial[W];` multiplies the integral cyclotomic norms of these factors. `rvWittTowerNewton[W];` prints the dyadic Newton slopes and their multiplicities. Computing the numerator uses counts of character values only through extension degree $u_s-1$, subject to the field, work, storage, and signed coefficient bounds. The leading coefficient is $2^g$, so an integral numerator requires genus below 127 in the signed 128-bit coefficient type.

The length-four tower $(x,0,0,0)$ has upper jumps $[1,2,4,8]$ and genus $35$. Its degree-$70$ numerator uses extension degrees at most seven. The slopes are $1/8,1/4,3/8,1/2,5/8,3/4,7/8$ with multiplicities $8,12,8,14,8,12,8$. [Witt tower example](examples/rank-witt-tower.asl)


### Koszul syzygies and their products

For a declared rank ideal $J\subset S=\mathbb F_2[x_i:i\in I]$, `RVKoszul(J, s, t)` is the discrete space of chains in $(S/J)_{t-s}\otimes\Lambda^s(e_i)$. The constructor `rvkterm[J, [i,j,...]](p)` reduces the homogeneous coefficient $p$ modulo $J$ and wedges the listed generators in increasing primitive order. The internal degree is the polynomial degree plus the wedge length. The empty wedge is the unit of the exterior algebra.

`rvkD(a)` replaces each wedge factor $e_i$ by $x_i$ and sums the resulting terms. It requires positive homological degree and preserves internal degree. `rvkProd(a,b)` multiplies the coefficients modulo $J$ and takes the wedge product. A repeated exterior generator gives zero. Both chains and classes support `zero` and `add`.

`RVTor(J, s, t)` is the corresponding homology group $\operatorname{Tor}^S_s(S/J,\mathbb F_2)_t$. Use `rvbetti[J,s,t];` for its rank and `rvtorgen[J,s,t,j]` for its zero-based basis generators. `rvTorClass(a)` projects a cycle to its class and `rvTorRep(x)` chooses a cycle representative. `rvTorProd(x,y)` multiplies represented cycles and takes their class. The calculation uses standard monomials of the reduced Gröbner basis, so individual sparse chains can be multiplied without enumerating the ambient finite space.

For the complete intersection $(xy+z^2,x^2+yz)$, the two classes in $(1,2)$ multiply to the generator at $(2,4)$. The example also multiplies a cycle by a nonzero chain boundary and verifies that the product is a boundary. For $(x^2,xy,y^2)$ the ranks in $(1,2)$ and $(2,3)$ are three and two, and all products of positive homological degree vanish. [Syzygies example](examples/rank-syzygies.asl)

#### Multigraded Koszul slices

For a monomial ideal, `RVMKoszul(J,s,[mu])` fixes the full multidegree $\mu$, the sum of the coefficient exponents and the exterior indicators in primitive order. Its homology is `RVMTor(J,s,[mu])`. The interpreter enumerates only wedge subsets of the support of $\mu$ and takes the standard coefficient $x^{\mu-1_W}$ for each surviving wedge $W$.

`rvmkterm[J,[wedge]](rvmon[I]{powers})` uses a literal monomial so its multidegree can be inferred. `rvmbetti[J,s,[mu]];` reports the local rank, and `rvmtorgen[J,s,[mu],i]` selects a local basis class by zero-based index. The Koszul differential, primitive, products and class operations also work in these types. Products add multidegrees.

`rvMComponent[[mu]](x)` projects a total chain or class of internal degree $|\mu|$ to one multidegree. `rvMForget(x)` includes a multigraded chain or class in the total group. Including a class requires the total Tor slice to fit the computation limits. For example, over $\mathbb F_2[x]/(x^2)$ the cycle $x e_x$ generates `RVMTor(J,1,[2])`.

#### Hochster cochains

For a squarefree monomial ideal $J$, `rvhochster[J,[W],q];` gives the reduced cohomology dimension of the induced Stanley–Reisner complex on the selected vertices $W$. It computes `RVMTor` at homological degree $|W|-q-1$ and multidegree $1_W$. Degree $q=-1$ is the augmented degree.

`rvhface[J,[W],[F]]` gives the cochain of an actual face $F\subseteq W$ as $x_F e_{W\setminus F}$, including the empty face. The Koszul differential adds one allowed vertex to $F$, which gives the augmented simplicial coboundary. A repeated exponent gives an acyclic slice because toggling its exterior factor preserves the support of the coefficient. In the Flag example the cochain of edge $15$ evaluates to one on the cycle $1-4-2-5-1$ and represents both the local triple value and the corresponding transferred product. [Hochster example](examples/rank-hochster.asl)

### Triple Massey products in Tor

`rvTBracket(a,b,c)` returns the predicate for the triple Massey product of three `RVTor` classes over the same ideal. If either adjacent product is nonzero in Tor, the predicate is empty. Otherwise the interpreter chooses cycle representatives $A,B,C$ and solves $dU=AB$, $dV=BC$. The class of $AV+UC$ lies in homological degree $s_a+s_b+s_c+1$ and internal degree $t_a+t_b+t_c$. All possible values form its coset modulo $a\operatorname{Tor}_{s_b+s_c+1,t_b+t_c}+\operatorname{Tor}_{s_a+s_b+1,t_a+t_b}c$.

`rvTInBracket(a,b,c,x)` tests membership using this linear span directly. `rvkPrimitive(w)` returns a chosen chain whose differential is the supplied boundary $w$, raising homological degree by one. It reports an error if $w$ has no such primitive.

In the flag Stanley–Reisner example with ideal $(x_1x_2,x_3x_4,x_5x_6,x_1x_3,x_4x_5)$, the classes of $x_1e_2$, $x_3e_4$, and $x_5e_6$ have a triple product with exactly one value, the nonzero class of $x_1x_5e_2e_3e_4e_6$. The first adjacent product is zero at the chain level and the second is $d(x_5e_3e_4e_6)$. Both indeterminacy groups are in $(3,4)$ and vanish. [Tor Massey example](examples/rank-massey.asl)

#### Explicit defining systems

`rvTCurvature(a1,...,an){i < j => w, ...}` takes chosen witnesses and returns the top Koszul chain. `rvTMassey(a1,...,an){...}` takes its homology class. There are at least two adjacent inputs in order from left to right, with vertices numbered $0,\ldots,n$. An adjacent input is either an `RVTor` class with its canonical representative or a chosen closed `RVKoszul` chain. All chains and classes use the same ideal.

For each proper interval $[i,j]$ the witness $A_{ij}$ has homological degree $\sum_{i<k\le j}s_k+j-i-1$ and internal degree $\sum_{i<k\le j}t_k$. Omitted witnesses are zero. The interpreter checks $d A_{ij}=\sum_{i<k<j}A_{ik}A_{kj}$ on every proper interval. The top chain is $\sum_{0<k<n}A_{0k}A_{kn}$, of homological degree $\sum_k s_k+n-2$. These expressions evaluate the supplied defining system.

The fourfold example uses the classes of $x_1e_5$, $x_2e_6$, $x_3e_7$, and $x_4e_8$. It computes a nonzero class in $(6,8)$ and verifies that the raw curvature differs from $x_4x_5e_1e_2e_3e_6e_7e_8$ by an explicit boundary. [Fourfold defining system example](examples/rank-massey-four.asl)

#### Transferred higher products

`rvTHigher(a1,...,an)` evaluates the minimal $A_\infty$ operation $m_n$ on `RVTor` or `RVMTor` for $n\ge2$. It adds $n-2$ to the total homological degree and preserves the total internal degree. `rvTLift(a1,...,an)` evaluates the corresponding Koszul chain $F_n$, which adds $n-1$ to the homological degree. The unary lift $F_1$ is the chosen cycle representative, $m_1=0$, and $m_2$ is `rvTorProd`. On multigraded classes the operation and the lift add the full multidegrees and use the corresponding local slices.

Let $P$ be the linear boundary preimage given by the ordered incoming Koszul columns. For an arbitrary chain $w$, set $q=w+P(dw)$. Then $p(w)=[q]$ and $h(w)=P(q+i(p(w)))$, with $i$ the chosen cycle inclusion. Over $\mathbb F_2$ the signs disappear, so $dh+hd=1+ip$ and $h^2=hi=ph=0$. For each interval the interpreter forms $U=\sum F_{\rm left}F_{\rm right}$ and takes $m=p(U)$ and $F=h(U)$. The closure step handles $U$ even when $dU\ne0$. Each interval lift is computed once per call.

For a monomial ideal, the reduction preserves the exponent vector of the coefficient plus the exterior indicators. The staircase example imposes $x_jy_k=0$ for $k\le j$ except $(j,k)=(n,1)$. On the inputs $[x_i e_{y_i}]$, every proper contiguous product vanishes and the top transferred product is the nonzero class of $x_ny_1$ times the exterior product of all generators except $e_{x_n}$ and $e_{y_1}$. The calculation fixes a contraction and preserves the variable multigrading. Unrestricted Massey products range over all defining systems, including choices that mix variable multidegrees. [Staircase generator](examples/staircase-transfer.awk)

```sh
awk -v n=6 -f examples/staircase-transfer.awk | nix run path:. -- -
```

The generated program checks the proper products and prints `RVTor(Stair, 10, 12){g0}`. The checks run at arities three through six with the standard computation limits.

Use multigraded slices for the longer calculation.

```sh
awk -v local=1 -v n=20 -f examples/staircase-transfer.awk | nix run path:. -- -
```

The result is the generator of `RVMTor` in homological degree 38 and multidegree $(1,\ldots,1)$ on forty variables. This program checks all 189 proper contiguous products and the nonzero twentieth product within the same computation limits.

### The directed Koszul bar complex

`rvbarpath P[J] { [mu1], ..., [mun] };` takes nonzero multidegrees in the primitive order of a monomial ideal $J$. It defines an ordered directed $A_\infty$ category with vertices $0,\ldots,n$. From $i$ to $j$ with $i<j$, the morphisms are the multigraded Tor groups at $\mu_{i+1}+\cdots+\mu_j$. Composition and all higher operations use the transferred product on those groups. Each vertex has a formal identity. The reduced bar complex uses only the positive intervals.

`RVBar(P,i,j,b)` sums the ordered tensors over all cuts $i=i_0<\cdots<i_k=j$, using suspended homological degree $b=\sum_{v=1}^k(s_v+1)$. Build a tensor with `rvBarWord[P,[i0,...,ik]](c1,...,ck)`. Each class must have the multidegree of its interval, and the constructor expands multilinearly. Printed factors show the interval, unsuspended degree and Tor basis index.

`rvBarD(x)` replaces each consecutive block of $m\ge2$ factors by its $m_m$ product. Every replacement lowers $b$ by one and tensor length by $m-1$. The Stasheff identities give $d^2=0$. Use `rvBarLength[k](x)` for the length-$k$ part and `rvbardim[P,i,j,b];` for the basis dimension, with an optional final length to count just that part. The basis calculation counts feasible cuts using cached interval homology. For squarefree flag ideals, dominated vertices can be removed when finding homology dimensions. The actual morphisms and products use the original Koszul slices.

For the three-square example, two distinct binary contractions of the length-three word have the same nonzero differential and cancel. The staircase generator also supports `-v bar=1`, which checks that the full word maps to the singleton top class. [Bar differential example](examples/rank-bar.asl) · [Staircase generator](examples/staircase-transfer.awk)

#### Length spectral sequence

`RVBarE(P,i,j,r,b,l)` is the $(b,l)$ group on page $r$ of the increasing filtration $F_l$ by words of length at most $l$. The associated graded is $E_0$. A representative on page $r>0$ lies in $F_l C_b$ and has differential in $F_{l-r}C_{b-1}$. The calculation quotients by lower-filtration cycles and the incoming boundaries of the appropriate filtration.

`rvbaredim[P,i,j,r,b,l];` prints the dimension and `rvbaregen[P,i,j,r,b,l,k]` selects a basis class. The operation `rvBarClass[r,l](x)` checks the representative condition and `rvBarRep` returns a corrected raw lift. On page classes `rvBarD` evaluates $d_r$ in bidegree $(b-1,l-r)$. Filtration lengths outside the range of word lengths have zero page groups. Use `rvBarPrev` for the chosen preceding-page lift and `rvBarNext` for a $d_r$-cycle. A class hit by an incoming differential maps to zero on the next page.

The full staircase interval has two bar basis words. The adjacent input word lies in $(2n,n)$ and the singleton top lies in $(2n-1,1)$. The first nonzero differential is $d_{n-1}$ between them, and both vanish on $E_n$.

```sh
awk -v pages=1 -v n=20 -f examples/staircase-transfer.awk | nix run path:. -- -
```

This checks the preceding differentials, the nonzero $d_{19}$ and the two zero groups on $E_{20}$ using the standard computation limits.

#### Dual bar cochains and path Ext

`RVBarDual(P,i,j,b)` is the finite dual of `RVBar(P,i,j,b)`. The operation `rvBarDualize(x)` uses the ordered bar basis to convert a chain support to its coordinate functional. `rvBarPair(f,x)` evaluates that functional, returning `S.hi` for one. The transposed differential `rvBarDualD` raises $b$ by one and satisfies $\langle D^*f,x\rangle=\langle f,Dx\rangle$.

`RVPathExt(P,i,j,b)` is the cohomology of the dual bar complex of the directed path category. Use `rvPathClass` on a closed dual cochain and `rvPathRep` for a chosen representative. `rvpathextdim[P,i,j,b];` prints the dimension and `rvpathextgen[P,i,j,b,k]` selects a basis class. The vertex units are formal and the types describe positive intervals.

For consecutive intervals $i<k<j$, `rvPathCup(f,g)` concatenates dual words at $k$ and adds cohomological degrees. It is dual to splitting the raw word at that cut. Terms contracting across the cut have zero pairing, which gives the cochain Leibniz formula. The induced product on `RVPathExt` is associative. The Squares example computes the dual differential of a singleton as a sum of two different cuts. [Dual bar example](examples/rank-path-ext.asl)

```sh
awk -v dual=1 -v n=20 -f examples/staircase-transfer.awk | nix run path:. -- -
```

Every proper consecutive product of the adjacent classes is nonzero. The full raw cup is the dual of the long bar word, and the dual of the singleton top has exactly that coboundary. The full path product therefore prints `RVPathExt(StairBar, 0, 20, 40){0}`.

#### Bar–cobar rectification

`RVBarCobar(P,i,j,b)` is the positive-interval part of the dg bar–cobar category $\Omega B(A)$ of the directed Tor $A_\infty$ category. Its words have outer factors, each a nonempty bar word desuspended once. Thus $b$ is the sum of the suspended degrees of all letters minus the number of outer factors. `rvBCDesuspend(x)` puts each supported bar word in a single outer factor, and `rvBCMul(x,y)` concatenates the outer lists at their common endpoint.

The differential `rvBCD` sums the bar differential inside each outer factor and all splittings of one factor at a proper inner gap. Over $\mathbb F_2$ this is a strict differential satisfying the Leibniz rule for `rvBCMul`.

The inclusion $i$ puts a Tor class in a singleton bar and desuspends it. For the split differential, the elementary homotopy $h_0$ removes the outer separator at the first gap of the flattened word if that separator is present. The projection $p$ retains singleton words. With $\delta$ the internal bar differential, `rvBCProject` and `rvBCHomotopy` apply
$ p'=p\sum_{j\ge0}(\delta h_0)^j$ and $h'=h_0\sum_{j\ge0}(\delta h_0)^j$. Each application of $\delta h_0$ lowers the flattened letter count. These maps satisfy $p'i=1$ and $Dh'+h'D=1+ip'$. On singleton inclusions $p'(i(a)i(b))=m_2(a,b)$. [Bar–cobar example](examples/rank-barcobar.asl)

```sh
awk -v omega=1 -v n=20 -f examples/staircase-transfer.awk | nix run path:. -- -
```

The differential of the one-factor word is the singleton top plus all nineteen outer splits. The example checks the retract identity both on that word and on its first-gap split, which projects to the twentieth operation.

### The Hochschild–Connes mixed complex of Tor

For a monomial ideal $J$, `RVHHChain(J,b,[mu])` is the normalized Hochschild chain group of its augmented multigraded $A_\infty$ Tor algebra over $\mathbb F_2$. The word `rvHHWord(a0,...,an)` denotes $a_0[a_1|\cdots|a_n]$. Each input is an `RVMTor` class over $J$. The coefficient $a_0$ may be the unit at zero multidegree. Every other factor has positive multidegree and contributes its suspension. The weight is the sum of the exponent vectors, $b=\sum_i s_i+n$, and `rvHHLength[n]` selects normalized bar length $n$. Negative homological degrees have zero chain groups.

The differential `rvHHd` uses the transferred $m_k$ on each cyclically consecutive block, for $k\ge2$. A block containing the coefficient becomes the new coefficient after contraction, in its cyclic order. Blocks away from the coefficient contract in the bar. Products with the unit have the strict unital values. The Connes operator `rvConnesB` inserts the unit as the coefficient and sums the rotations of all factors as reduced inputs. A word with unit coefficient maps to zero. Repeated terms cancel over $\mathbb F_2$, and the two operations satisfy $d^2=B^2=dB+Bd=0$.

For $J=(x)$ there is one exterior Tor generator $e$. In weight $N>0$ the two normalized words have degrees $2N$ and $2N-1$. Their Hochschild differential vanishes and $B(e[e|\cdots|e])$ is $N$ times the word with unit coefficient and $N$ reduced factors. [Hochschild–Connes example](examples/rank-hochconnes.asl)

```sh
awk -v hh=1 -v n=20 -f examples/staircase-transfer.awk | nix run path:. -- -
```

The length-one part of the differential is the unit with the nonzero twentieth transferred product in its bar, in degree $39$ and weight $(1,\ldots,1)$ on forty variables. The two remaining terms come from multiplication by the unit at the cyclic ends. The generator checks the full differential and its Connes image at arity twenty. It also checks the second Hochschild differential through arity nine within the standard computation limits.

#### Hochschild homology and the induced Connes operator

`RVHH(J,b,[mu])` is the homology of the normalized Hochschild differential in degree $b$ and weight $\mu$. Use `rvHHdim[J,b,[mu]];` for its dimension and `rvHHgen[J,b,[mu],k]` for a generator of the computed basis. `rvHHClass(x)` checks that the chain $x$ is a cycle and reduces modulo the incoming boundaries. `rvHHRep` gives a chosen cycle representative. The Connes operator also acts on these classes, raising $b$ by one. The identity $dB+Bd=0$ sends cycles to cycles and boundaries to boundaries, so the result is independent of the representative.

For the exterior generator example, the induced operator in weight $N>0$ from degree $2N-1$ to $2N$ is nonzero precisely when $N$ is odd. In the Flag example the dimensions at degrees four, five and six are one, three and two. The $m_3$ term in the differential of $1[a|b|c]$ identifies two different representatives of the nonzero Connes image of the degree-four class. [Hochschild homology example](examples/rank-hochhom.asl)

#### Finite negative cyclic homology

`RVNCChain(J,N,b,[mu])` is the degree-$b$, weight-$\mu$ part of the negative cyclic complex modulo $u^N$, for $N\ge1$. The variable $u$ has homological degree $-2$ and leaves the weight unchanged. The differential `rvNCd` is $d+uB$. Its square vanishes by the mixed-complex identities. In degree $b$ the coefficient of $u^j$ comes from ordinary Hochschild degree $b+2j$.

Use `rvNCTerm[N,j](x)` to insert an ordinary Hochschild chain as the coefficient of $u^j$, and `rvNCPart[j](z)` to read it. The power $j$ is below $N$. The homology is `RVNC(J,N,b,[mu])`. Use `rvNCClass` on a cycle and `rvNCRep` for a chosen representative. The query `rvNCdim[J,N,b,[mu]];` prints the dimension, and `rvNCgen[J,N,b,[mu],k]` selects a basis class.

On chains and homology, `rvNCU[k]` multiplies by $u^k$ and lowers $b$ by $2k$. Components with powers at least $N$ vanish. The map `rvNCTruncate[M]` to $1\le M\le N$ discards all components of powers at least $M$. At $N=1$ the constant coefficient identifies the complex with ordinary Hochschild chains. These are individual finite truncations.

For the exterior generator, odd weight $n$ leaves the two classes represented by $1[e|\cdots|e]$ at $u^0$ and $e[e|\cdots|e]$ at $u^{N-1}$. In even positive weight every power of both words survives, for total dimension $2N$. [Finite negative cyclic example](examples/rank-negative.asl)

### Rank-two localized torsion bars

For $I=[i,n]$ let $w=v_i/v_n$ and $L=D_n-D_i$. Multiplication by $v_n^{-s}$ identifies the slice at filtration $s$ with filtration zero and gives normalized grade $g=t-sD_n$. Multiplication by $w$ lowers this grade by $L$. The filtration-zero complex is a free graded $\mathbb F_2[w]$ module on the coefficient matrices between the endpoints. Its differential $v_n^{-1}qD$ sums the two endpoint actions for $Q_n$ and $wQ_i$.

`vbars[X, Y, [i, n]];` computes the decomposition of its cohomology and reports the free generator grades and the torsion pairs $(d,e)$ for a generator of grade $d$ killed by $w^e$. The calculation uses graded Smith reduction on the coefficient matrix. Minimum-exponent pivots allow homogeneous row and column operations, so each entry stays a single monomial. The tracked target and source basis changes give the torsion representatives and their boundary witnesses. The cokernel projects onto the free image by the normalized differential of degree $-D_n$. This split quotient subtracts one free generator in the original grade of each pivot source column from the unpivoted target grades.

`vbardim[X, Y, [i, n], s, t];` recovers the rank at any signed bidegree from this decomposition. A bar of grade $d$ contributes when $d-g=kL$ for an integer $k\ge0$. A torsion bar also requires $k<e$.

`vbarRep[X, Y, [i, n], d, e, k]` gives the raw representative in bidegree $(0,d)$ of the $k$th torsion bar with pair $(d,e)$, with indexing from zero among equal pairs. `vbarNull` takes the same arguments and gives the witness in bidegree $(-1,d-eL)$. The differential of that witness equals `vRatio[e]` of the representative. The operation `vRatio[k](x)` multiplies a rank-two localized cochain or class by $w^k$ for a nonnegative integer $k$. The Flash example gives a bar of grade 7 and exact order 4. [Graded Smith example](examples/v-bars.asl)

### Margolis homology over truncated polynomial rings

`rvjet[X, Y, I, modulus, N, coefficients];` takes one list of coefficients for each selected primitive in `I`. The entries of each list give $\lambda_i(u)$ in ascending powers of $u$ over the same binary field $K$ as `rvdim`. The order `N` must be positive. Missing entries are zero and entries beyond degree $N-1$ must be zero.

Set $R=K[u]/(u^N)$. After forgetting the internal grading, the interpreter forms the $R$-linear differential $D(F)=\sum_i\lambda_i(u)(Q_i^Y F+FQ_i^X)$ on the module of coefficient matrices from `X` to `Y`. The result gives the dimension of $\ker D/\operatorname{im}D$ over $K$ and the multiplicity of each Jordan block length for multiplication by $u$. A block of length $j$ is a summand $R/(u^j)$. At `N = 1` the dimension agrees with `rvdim` at the constant coefficients.

The three arcs $(u,0,1)$, $(0,u,1)$, and $(u^3,u^2,1)$ through the cusp have the same residue. At order seven from `One` to `Cusp` they give block multiplicities `{1: 4}`, `{1: 2, 2: 2}`, and `{2: 2, 7: 2}`. The example also computes the blocks for endomorphisms of `Cusp`. [Margolis jet example](examples/rank-jets.asl)

### Named Margolis arcs and cochains

`rvpath Axis[[0,1,2], 7, 7, [[],[0,1],[1]]];` names the arc $(0,u,1)$ at order seven over $\mathbb F_4$. Use the name with `rvjet[X, Y, Axis];` or with `JetHom(X, Y, Axis, parity)`. The two parities record target degree minus source degree modulo two, and $u$ has even parity.

`jetterm[X,q,i,Y,p,j,Axis,k]` gives the elementary cochain that sends generator `i` in degree `q` of `X` to $u^k$ times generator `j` in degree `p` of `Y`. The exponent `k` must be below the arc order. Use `add`, `zero`, and `rvScale[c]` for the field operations. `jetD` applies the two endpoint actions and flips the parity. `ecomp(g,f)` composes cochains on one arc, adding their parities and multiplying the powers of $u$ modulo the truncation. The differential obeys the Leibniz rule.

On the $y$-axis in the Cusp example, the corrected lift $a_1+u a_3$ has differential $u^2 b_3$. The lift $a_1+u a_2+u^2 a_3$ is closed along the cusp parametrization. [Named arc example](examples/jet-hom.asl)

### The u-adic Bockstein pages

`JetPage(X, Y, Axis, r, parity, f)` is the filtration-`f` slice on page `r` of the two-periodic complex over the named arc. The decreasing filtration is $F^f=u^f C$, where $C$ is the raw jet complex. Page zero is the associated graded complex and its differential uses the constant coefficients of the arc. Page one is its Margolis homology. The differential $d_r$ flips parity and increases filtration by $r$. At arc order $N$, page $N$ has stabilized, and summing its ranks over parity and filtration gives the dimension from `rvjet`.

`jetdim[X, Y, Axis, r, parity, f];` prints a rank and `jetgen[X, Y, Axis, r, parity, f, i]` selects a generator. `jetClass[r, f](c)` projects a raw cochain in $F^f$ whose differential lies in $F^{f+r}$. The operation `jetRep` returns a canonical lift. Use `jetD` for the page differential, `jetNext` to advance a cycle to the next page, and `jetPrev` to take a chosen lift to the preceding page. Composition on a common page and arc adds filtrations and obeys the Leibniz rule. The operations `add`, `zero`, and `rvScale` act in each slice.

On $(0,u,1)$, the class of $a_2$ supports $d_1$ with target $u b_2$, and the corrected class of $a_1$ supports $d_2$ with target $u^2 b_3$. On $(u,0,1)$ both remaining differentials from filtration zero occur on page one. Along $(u^3,u^2,1)$ the corrected $a_1+u a_2+u^2 a_3$ gives a nonzero stabilized class. [Jet page example](examples/jet-pages.asl)

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
| `pextdim[p, s, t];` | Print the dimension of `PExt(p, s, t)` over 𝔽ₚ. |
| `pextandim[p, n, s, t];` | Print the dimension of `PExtAn(p, n, s, t)` over 𝔽ₚ. |
| `pmaydim[p, r, s, t, f];` | Print the dimension of `PMay(p, r, s, t, f)` over 𝔽ₚ. |
| `pmayandim[p, n, r, s, t, f];` | Print the dimension of `PMayAn(p, n, r, s, t, f)` over 𝔽ₚ. |
| `pceandim[p, n, r, s, t, f];` | Print the dimension of `PCEAn(p, n, r, s, t, f)` over 𝔽ₚ. |
| `extandim[n, s, t];` | Print the dimension of `ExtAn(n, s, t)` over 𝔽₂. |
| `mayandim[n, r, s, t, f];` | Print the dimension of a finite-profile May page component over 𝔽₂. |
| `ceandim[n, r, s, t, f];` | Print the dimension of a Cartan–Eilenberg page component over 𝔽₂. |
| `xextdim[X, s, t];` | Print the dimension of `XExt(X, s, t)` over 𝔽₂. |
| `mapextdim[X, Y, s, t];` | Print the dimension of `MapExt(X, Y, s, t)` over 𝔽₂. |
| `qextdim[X, Y, I, s, t];` | Print the dimension of `QExt(X, Y, I, s, t)` over 𝔽₂. |
| `vextdim[X, Y, I, s, t];` | Print the dimension of `VExt(X, Y, I, s, t)` over 𝔽₂. |
| `vbars[X, Y, [i, n]];` | Print free generator grades and torsion bars of rank-two localized Koszul cohomology. |
| `vbardim[X, Y, [i, n], s, t];` | Read a signed bidegree rank from the rank-two decomposition. |
| `mapmaydim[X, Y, r, s, t, f];` | Print the dimension of a bivariant May page component over 𝔽₂. |
| `card p;` | For a predicate `p : X -> S`, print the number of points sent to `S.hi` and the size of `X`. On the empty domain this gives `card = 0 / 0`. |

A homotopy step uses a single direction across the whole domain, so it requires `f(x) ≤ g(x)` for every `x`, or the reverse inequality for every `x`. In the triangle example, any endomorphism above the identity must fix all three maximal edges. Each vertex is the unique common lower bound of its two incident edges, so it too is fixed. The dual argument handles endomorphisms below the identity. This isolates the identity in `Circle -> Circle` and explains the `no` printed for the rotation. [Stong, *Finite Topological Spaces*, §3–4](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/stong2.pdf)

For `betti`, a strict chain of `k + 1` points contributes a `k`-simplex. Dropping each point in turn gives its boundary over 𝔽₂ in positive degrees. The ranks of successive boundary matrices determine the homology. [Barmak and Minian, *Simple Homotopy Types and Finite Spaces*, §2](https://math.uchicago.edu/~may/FINITE/ArgentinianPapers/MinianSimple.pdf)

The evaluator tabulates a lambda on its domain. For a discrete `D` with `n` points, an individual map `D -> S` needs `n` entries, while enumerating the whole function space yields `2^n` maps. A total order of `n` points has `2^n − 1` nonempty chains. The interpreter admits up to 512 points per named space, 4,096 elements per enumerated type and 50,000 simplices per order complex. A cohomology group of rank `b` has `2^b` elements. Individual classes can be constructed and multiplied even when enumerating that space would exceed the limit. Work and storage budgets can stop a calculation earlier with a diagnostic.
