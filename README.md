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

## Finite set families on a poset

A `setfamily F[X] {x:2, y:3} {x<y:[0,2]};` gives a covariant finite set
diagram on a named site with two comparable points. Supply the cardinality
of every fiber, including empty fibers, and the image list for every cover.
Indices start at zero. A list for any other strict comparable pair is
optional and must agree with every composite. Identity maps are automatic.

`elements E[F];` defines the usual finite total poset. Its points are
`E.x_0`, `E.x_1`, and so on. A point `E.x_i` lies below `E.y_j`
exactly when the transition from `x` to `y` takes `i` to `j`.
The total uses names `point_x_i` with numeric site indices if appending the
fiber indices would exceed the identifier length limit. `setProjection[E]`
is an ordinary monotone map from `E` to `X`.

`setnat u[F,G] {x:[...], y:[...]};` declares a natural transformation of
type `SetNat(F,G)` between diagrams on the same site. Every component must
commute with the transitions. `setId[F]` and `setComp(g,f)` give identities
and composition. If `E` and `D` are the declared total spaces for `F` and
`G`, `setTotal[E,D](u)` is an ordinary monotone map over the site.

`setsections[F];` counts global sections, and `setsect[E,i]` gives an
ordinary section of the projection in lexicographic order of the fiber
indices. See [examples/set-families.asl](examples/set-families.asl) for a
connected three-sheeted circle cover with no section and the trivial cover
with three sections.

### Pullback and dependent sums and products

For a monotone map `f : Y -> X` of named sites, `setpullback PH[H](f);`
pulls a family `H` on `X` back to `Y`. Given a natural map `u : H -> K`,
`setPullNat[PH,PK](u)` pulls its components back along the common map.

If `E` is the elements space of `F` and `G` is a family on `E`,
`setsum Sum[E,G];` and `setproduct Pi[E,G];` give families on the base of
`F`. The sum fiber at `x` is the disjoint union of `G(x,i)` over
`i` in `F(x)`. The product fiber at `x` consists of compatible sections
of `G` over all points of `E` lying above `x` in the base, including the
points lying over `x` itself. Product transitions restrict sections to
the smaller future. If the current fiber is empty, the future can still
admit multiple sections. An entirely empty future has one section.

Let `PH` and `PP` be the pullbacks of `H` and `Pi` along
`setProjection[E]`. Then `setEval[Pi,PP]` has type `SetNat(PP,G)`.
`setTranspose[Pi](u)` turns `u : SetNat(PH,G)` into a map
`SetNat(H,Pi)`, and `setUntranspose[PH](v)` goes back. For the sum,
`setSumLift[Sum](u)` turns `u : SetNat(G,PH)` into a map
`SetNat(Sum,H)`, and `setSumRestrict[PH](v)` goes back. See
[examples/dependent-sets.asl](examples/dependent-sets.asl) for both
adjunctions, the explicit evaluation map, and an empty current fiber.

### Pointwise Kan extensions of finite set families

For a monotone map `f : X -> Y` of named sites and a family `F` on `X`, declare `setlan L[F](f);` and `setran R[F](f);`. The fiber of `L` at `y` is the disjoint union of the fibers `F(x)` for which $f(x)\le y$, modulo the identifications from every source transition in that subposet. Its transition sends each class to the class of the same occurrence. The fiber of `R` consists of sections over the subposet $y\le f(x)$, and its transitions restrict sections. An empty indexing subposet gives an empty left fiber and a singleton right fiber. Left classes are indexed by their first occurrence in point and fiber order. Right sections are indexed by their lexicographically ordered coordinates.

Let `PH` be `setpullback PH[H](f)` for a family `H` on `Y`. Then `setLanLift[L](u)` takes `u : SetNat(F,PH)` to `SetNat(L,H)`, and `setLanRestrict[PH](v)` goes back. In the other adjunction, `setRanTranspose[R](u)` takes `u : SetNat(PH,F)` to `SetNat(H,R)`, and `setRanUntranspose[PH](v)` is its inverse. If `PL` and `PR` pull back `L` and `R` along `f`, the maps `setLanUnit[L,PL]` and `setRanCounit[R,PR]` insert an occurrence into its left fiber and evaluate a right section at that source point. For extensions `LH` and `RH` of `PH`, `setLanCounit[LH]` and `setRanUnit[RH]` give the other counit and unit.

For a commuting square with horizontal maps `gPrime : XPrime -> X` and `g : YPrime -> Y`, first pull back `F` along `gPrime` and extend it along `fPrime : XPrime -> YPrime` to obtain `LPrime` and `RPrime`. Pull back `L` and `R` along `g` to obtain `LPull` and `RPull`. Then `setLanBC[LPrime,LPull]` maps a left class to the class of its image occurrence. The right comparison `setRanBC[RPrime,RPull]` goes from `RPull` to `RPrime` by precomposing sections. ASL checks the four maps commute on every point of `XPrime`.

Even a strict pullback square can give a comparison that fails to be invertible. Over the chain $0<1$, take the two distinct inclusions of a point and put two elements on the source. When the source lands at $0$, the left comparison maps an empty fiber to two elements. When the source lands at $1$, the right comparison maps two elements to a singleton. The example computes both maps and checks the adjunctions on a fork where the left quotient identifies elements from different branches. [Kan extension example](examples/set-kan.asl)

### Artin gluing on an open and its closed complement

Use `setrecollement R(i,j);` for named sites and monotone maps `i : C -> X` and `j : U -> X`. The maps must be full order embeddings whose images partition `X`, with `j(U)` upward closed. Points can have different names and orderings in the three declarations. Empty strata are allowed. Use `setrestrict A[R,F,closed];` and `setrestrict B[R,F,open];` to pull the ambient family `F` back along the two inclusions.

For any open family `B`, `setgluelimit T[R,B];` constructs $i^*\operatorname{Ran}_j B$. Over a point $c$ in the closed stratum its fiber consists of compatible sections on $\{u\in U:i(c)\le j(u)\}$. This means that each open transition in the selected subposet carries the assigned source value to the assigned target value. There is one section of the empty subposet. Transitions between closed points restrict these section tuples, indexed in lexicographic order in the declared open site. The construction uses the same section solver and finite limits as right Kan extension.

A map `alpha : SetNat(A,T)` specifies the attachment. Use `setglue G[R](alpha);` to build the ambient family with those two strata. A crossing arrow $i(c)\to j(u)$ evaluates the chosen section $\alpha_c(a)$ at $u$. The compatibility of the open section handles composition with a later open arrow, and naturality of $\alpha$ handles composition with an earlier closed arrow.

When `A` and the open family of `T` are restrictions of the same ambient `F`, `setAttach[A,T]` constructs the attachment using every crossing transition of `F`. After gluing this attachment, `setGlueRecover[G,F]` and `setGlueForget[G,F]` are the canonical maps $F\to G$ and $G\to F$. They identify corresponding local indices and are inverse. Supplying another attachment with the same two strata is allowed, but the reconstruction maps require the canonical attachment of the given original family.

For two glued families `G,G2` on the same recollement, take maps `a : SetNat(A,A2)` and `b : SetNat(B,B2)`. The induced map `setGlueLimitMap[T,T2](b)` applies $b$ to each coordinate of each compatible section. The expression `setGlueMap[G,G2](a,b)` checks $\alpha_2a=(i^*\operatorname{Ran}_j b)\alpha$ at every closed occurrence and joins the two arrays of components. Conversely, `setGlueClosed[G,G2](m)` and `setGlueOpen[G,G2](m)` extract this compatible pair from an ambient map. These constructions are inverse, since every component is on exactly one stratum. The first example verifies both inverse laws with a noninjective open transition. It also splits the open stratum a second time and gives the maps that undo this finer gluing. [Artin gluing example](examples/set-artin-gluing.asl)

In the Fano example, the closed stratum is the two minima and the open stratum is the three maxima of $K_{2,3}$. The point and line families have identical restrictions, and each matching fiber has $7^3=343$ sections. The attachments agree at one closed point and differ at the other. The program exhibits an isomorphism and inverse between the two right Kan extensions of the open restrictions. Both the natural-map type and the source stabilizer detect the empty map set from points to lines. The natural-map calculation propagates equalities in both directions along invertible transitions. The unique endomorphism of the point family is recovered from its two restrictions. Consequently there is no compatible pair of stratum maps in that direction. [Fano gluing example](examples/set-fano-gluing.asl)

### Joint monodromy and permutation sheaves

`setmonodromy M[F,G,x];` compares two set families on the connected Hasse component of the chosen base point `x`. All cover maps in that component must be bijections. A spanning tree transports each fiber to the base point. Each remaining edge gives a pair of permutations of the base fibers, obtained by traversing the same loop in both families. ASL enumerates the subgroup generated by those pairs. The computation has explicit work and storage limits.

Use `setmonosize[M];` for the image order and `setmonoat[M,i];` for the two permutation arrays of the indexed image element. The identity has index zero. `setmonochars[M];` counts elements by their pair of fixed-point counts. The two permutation characters agree precisely when all the bins lie on the diagonal. `setmonostabilizer[M,i];` prints the order of the stabilizer of source point `i`, then its numbers of common fixed points in the source and target. Each source orbit can map equivariantly to any target point fixed by its stabilizer. `setmonomaps[M];` multiplies these counts to obtain the exact number of set maps on the component, with an overflow diagnostic if the answer exceeds an unsigned 128-bit integer.

`shlinearize LF[F,p];` replaces every finite fiber by its free vector space over $\mathbb F_p$ and extends transition maps linearly. The result is an ordinary sheaf, and the construction also allows nonbijective transitions. For a monodromy comparison, `setmonoorbitals[M];` lists the sizes of the orbits in the product of the target and source base fibers. Order the pairs by target index first and source index second, and index each orbit by its first encounter. Their indicator matrices form a basis of the space of linear maps commuting with both loop actions.

`shorbital N[M,LF,LG,i];` declares the sheaf map given by the $i$th orbital, transported to all stalks of the chosen component. Both sheaves must be the indicated linearizations over the same field. `shorbitals L[M,LF,LG] {i:a,j:b};` takes a signed linear combination, and the maps are zero off the chosen component. Use `shmaprank[N,x];` for the rank of a stalk matrix. `shMapC[N]` gives the actual degree-zero cochain for composition and passage to `ShExt`.

In the Fano example two connected seven-sheeted covers have Betti ranks $[1,8]$. The common monodromy image has order 168 and the permutation characters agree. A point stabilizer fixes one source point and no target line, so the covers admit no map over the base. The incidence orbital matrix $N$ satisfies $NN^{\mathsf T}=2I+J$ for the all-ones matrix $J$ and has absolute determinant 24. The example constructs its inverse over $\mathbb F_5$ and verifies both composites in native sheaf Hom. Over $\mathbb F_2$ the two orbitals span Ext in degree zero and the four maps have ranks $0,1,3,4$, so the permutation sheaves are nonisomorphic there. [Fano monodromy example](examples/set-monodromy.asl)

### Internal logic of finite set families

`setomega O[X];` defines the truth classifier on a named site. Its fiber
at `x` consists of upward subsets of the future of `x`, also called
cosieves. Transports restrict an upward subset to the new future. Use
`setsieve[O,x,i];` to inspect the subset at a particular fiber index.

A predicate on `F` is a natural transformation to `O`.
`setClassify[O](i)` takes a componentwise injection into `F` to the
predicate that records where each element's future images lie in the
subfamily. `setsub M(P);` forms the subfamily on which `P` is true now.
Its inclusion `setInclusion[M]` has classifying predicate `P`.

`setTop[F,O]`, `setBottom[F,O]`, `setAnd(P,Q)`, and `setOr(P,Q)`
give the full future, empty future subset, intersection, and union.
`setImplies(P,Q)` records the points where every later instance of
`P` implies `Q`. `setNot(P)` is implication into bottom. A proposition
that becomes true at the later point of a chain can be dense while
excluded middle still fails.

For `p : SetNat(A,B)`, predicate pullback is `setComp(Q,p)`.
`setExists(p,P)` and `setForall(p,P)` quantify a predicate on `A`
along `p`. At each stage, existential truth needs a witness at that
stage. Universal truth needs every fiber element over the transported
target at every later stage to satisfy the predicate. See
[examples/set-logic.asl](examples/set-logic.asl) for strict failure of
excluded middle and a round trip through comprehension.

### Double-negation sheafification

`setdense[O,x,i];` tests whether a cosieve is double-negation dense.
On a finite poset, density means that every future point has a later
point in the cosieve. Every maximal future must therefore belong to it,
and those points already suffice. The smallest dense cosieve at `x`
contains exactly the maximal points above `x`.

`setj J[F];` defines the double-negation sheafification. Its fiber at
`x` consists of tuples of values of `F` at the maximal futures of `x`.
Transport restricts tuples and `setJUnit[J]` sends each input to its
images at all those maximal points. Matching families on the smallest
dense cosieve are exactly these tuples. A family is a double-negation
sheaf when these unit components are all bijections, which
`setjcheck[F];` tests.

For a natural map `u : SetNat(F,G)` to a double-negation sheaf `G`,
`setJExtend[J](u)` is the unique extension from `J` to `G`. It applies
`u` at the maximal future points and glues those images using the
bijection for `G`. Precomposing the extension with `setJUnit[J]`
recovers `u`. If `F` already satisfies the sheaf condition,
`setJInverse[J]` gives the inverse of the unit. See
[examples/set-sheafification.asl](examples/set-sheafification.asl) for a
fork with six sections after reflection, explicit extension, and
idempotence.

### Lawvere–Tierney topologies and sheafification

Given a classifier `setomega O[X];`, declare `settopology J[O] {s1,...};` for a subset $S$ of the site. At $x$ its local cosieve operator is
$j_{S,x}(U)=\{y\ge x:\ S\cap\uparrow y\subseteq U\}$.
The typed endomorphism `setTopMap[J]` has type `SetNat(O,O)`. The statement `settopologyMap J(f);` accepts a natural endomorphism of the classifier. It checks truth, extensivity and idempotence at every fiber, recovers $S$ by testing $\uparrow x\setminus\{x\}$, and verifies the reconstructed formula on every cosieve. In that formula a selected successor lies in an intersection exactly when it lies in both cosieves, which certifies the meet law. Use `settopologysupport[J];` to print the recovered subset in site order.

Use `settopologymeet M[J,K];` and `settopologyjoin L[J,K];` for the meet and join of two nuclei. Their subsets are respectively the union and intersection of the two subsets. The statement `settopologyjoinsteps[J,K];` checks the actual alternating closure table and gives the number of passes of $j_J\circ j_K$ needed to reach the join on all local cosieves. On a chain with alternating retained points, the order of the two operators matters. The tests include alternating supports on chains of thirty and eighty-nine points, which need sixteen and forty-five passes.

The declaration `setreflect R[J,F];` builds the $J$-sheafification of $F$. Its fiber at $x$ is the compatible sections of $F$ over the full subposet $S\cap\uparrow x$. This is right Kan extension of the restriction of $F$ to $S$ along the inclusion, and its transport forgets the coordinates outside the new future. The unit `setReflectUnit[R]` takes each value to its transports to the retained points. The query `setreflectcheck[J,F];` tests whether every component of that unit is bijective. In that case `setReflectInverse[R]` gives the inverse.

For a natural map $f:F\to G$ into a $J$-sheaf, `setReflectExtend[R](f)` applies $f$ to each retained coordinate and reconstructs a value of $G$ by its bijective unit. This implements the unique extension with composite $f$ along the unit of $F$. The diamond example has two compatible choices at the excluded bottom point. It verifies restriction and extension on all four natural maps in the adjunction, and the unit and inverse on the reflected family. The full subset gives the identity reflection and the empty subset gives the terminal one. For the subset of maximal points, the result is naturally isomorphic to double-negation sheafification. [Topology and reflection example](examples/set-topologies.asl)

## Sheaf cohomology on a finite poset

`sheaf F[X,p] {a:2,b:1} {a<b:[1,-1]};` declares a sheaf of vector spaces over the prime field $\mathbb F_p$ on the named poset `X`. Give a dimension for every point and one restriction $F(a)\to F(b)$ for each Hasse cover $a<b$. List matrix entries by rows, with the target dimension giving the number of rows. Zero-dimensional sources or targets require an empty matrix list. Signed decimal entries are reduced modulo $p$. The declaration composes along the Hasse diagram in topological order and checks that all paths with the same endpoints give the same matrix. Sheaves on the same named poset share its order data. Strict restrictions are stored in runs of equal nonzero matrices.

`ShC(F,k)` is the normalized Roos cochain space in degree $k$. The factor for a strict chain $x_0<\cdots<x_k$ is the terminal stalk $F(x_k)$. Use `shTerm[F,[x0,...,xk],i]` for its basis coordinate $i$, numbered from zero, or add a signed coefficient as the fourth argument. `shD(c)` takes the alternating sum over faces. When a face drops the terminal point, the contribution uses the restriction from its old terminal stalk to the new one. In particular, on a two-point chain the differential of a section is $c(b)-r_{ab}c(a)$.

`ShH(F,k)` is the resulting cohomology space. `shClass(c)` projects a closed cochain to its class and `shRep(h)` returns a deterministic closed representative. Use `shdim[F,k];` for the dimension and `shgen[F,k,i]` for cohomology basis element $i$. The basis calculation orders chains by their point indices in the space declaration, then by terminal stalk coordinate. Both types support `zero`, `add`, and signed scalar multiplication `shScale[c](v)`. Sparse cochains can be formed and differentiated directly. Cohomology uses bounded prime-field reduction of the adjacent degrees.

The circle example compares constant coefficients, rank-two unipotent holonomy over $\mathbb F_3$, scalar holonomy $2$ over $\mathbb F_5$, and a sheaf supported on the three maximal points. [Finite-poset sheaf example](examples/sheaf-cohomology.asl)

### Sheaf pullback and higher direct images

`shpull G[f,F];` pulls a sheaf on $Y$ back along a named monotone function $f:X\to Y$. The stalk at $x$ is $F(f(x))$. A cover whose points have the same image acts by the identity. Other covers use the full restriction along the two images, including when they are separated by intermediate points.

`shpush G[f,F,q];` declares the higher direct image $R^qf_*F$ of a sheaf on $X$. Its stalk at $y$ is the Roos cohomology $H^q(X_y;F)$ for $X_y=f^{-1}\{z\mid y\le z\}$. If $y\le z$, then $X_z\subseteq X_y$. The restriction takes a cocycle representative, drops all chains outside $X_z$, and projects to its class there. Empty opens give zero stalks, and equal opens use the same basis and restriction identity.

For the map to a point, the stalk of $R^qf_*F$ is the global cohomology in degree $q$. The example sends two circles to different points in a chain. A restriction of $R^1f_*F$ forgets the circle lost when passing to the smaller open. [Sheaf direct-image example](examples/sheaf-direct-image.asl)

### Sheaf maps and hypercohomology

`shmap f[F,G] {a:[...],b:[...]};` declares a natural transformation. The sheaves must use the same named poset and prime field. Give a row-major matrix $f_a:F(a)\to G(a)$ at every point. The checker verifies $G(a\to b)f_a=f_bF(a\to b)$ at each cover. Use `shApply[f](v)` for the induced map on Roos cochains or their cohomology.

A declaration `shcomplex K {F,d,G,e,H};` puts the listed sheaves in successive stalk degrees, starting in degree zero. Each map must have the indicated source and target, and every adjacent composite must vanish at every stalk. The total space `ShTot(K,n)` sums the Roos cochains of bidegree $(a,q)$ with $a+q=n$. Its differential is $D=\delta_{\mathrm{Roos}}+(-1)^a d_{\mathrm{stalk}}$. Use `shTotEmbed[K,q](v)` to include a cochain of the degree-$q$ sheaf and `shTotPart[q](z)` to project a total cochain to that component. Displayed total terms carry a `q` prefix for the stalk degree.

`ShHyper(K,n)` is total hypercohomology. The same `shD`, `shClass`, `shRep`, and `shScale` operations apply to the total complex and its classes. `shhyperdim[K,n];` prints its rank and `shhypergen[K,n,i]` returns a canonical basis class. On the six-point model of the sphere, the example compares two complexes whose stalk cohomology sheaves are constant in degrees zero and one. The split example has hypercohomology ranks $1,1,1,1$ in total degrees zero through three. The nonsplit example has ranks $1,0,0,1$. [Sheaf hypercohomology example](examples/sheaf-hypercohomology.asl)

### The sheaf hypercohomology spectral sequence

`ShE(K,r,p,q)` is the page-$r$ group for the filtration of `ShTot(K,p+q)` by Roos degree. `shedim[K,r,p,q];` prints its rank and `shegen[K,r,p,q,i]` returns a basis class. Page zero carries the signed stalk differential. Page one is the Roos complex on stalk cohomology, so page two is $H^p(X;H^q(K))$. `shED(z)` has target bidegree $(p+r,q-r+1)$. A signed stalk degree $q$ outside the complex gives a zero group.

`shEClass[r,p](z)` projects a total cochain of filtration at least $p$ when $Dz$ has filtration at least $p+r$. Its stalk degree is inferred from the total degree. `shERep(x)` recovers a corrected total representative. `shELiftTo[r,p](z)` keeps the terms of degree $p$ and solves for terms in higher Roos degrees to meet the same cycle condition. It returns the input when the condition already holds and reports an obstruction when no correction exists. These computations use sparse finite filtration windows over the prime field.

`shENext(x)` sends a $d_r$-cycle to its next-page class. `shEPrev(y)` uses the canonical total representative to select an earlier-page cycle, so `shENext(shEPrev(y)) == y` when $r>0$ on the source of `shEPrev`. Incoming boundaries map to zero under `shENext`.

The transgression generator builds a length-$r$ complex on a $2r+2$-point model of the $r$-sphere. One entry in a restriction matrix distinguishes the split and twisted versions. In the twisted version the class in bidegree $(0,r-1)$ survives to page $r$, where $d_r$ maps it nontrivially to $(r,0)$. Both groups then vanish. Set `untwisted=1` to recover the split version, where they survive. The generator supports larger declarations within the interpreter's limits. Over $\mathbb F_3$ the twisted and split examples both run through page seven under the default limits. [Transgression generator](examples/sheaf-transgression.awk)

```sh
awk -v r=6 -v prime=3 -f examples/sheaf-transgression.awk | cargo run --quiet -- -
```

### Certified cellular resolutions

`shresolution R[X,p] {u:(0,a),v:(0,b),e:(1,b)} {e:[u:1,v:-1]} {u:1,v:1};` declares an augmented chain complex resolving the constant sheaf on `X` over $\mathbb F_p$. Each generator `(degree,point)` contributes the projective sheaf $P_x$, with a one-dimensional stalk at every $y\ge x$ and identity restrictions between those stalks. The second block lists sparse boundaries and the third lists the augmentation. A boundary term must lower the degree by one and be anchored at or below the point of its source. Coefficients are signed integers reduced modulo $p$, and omitted entries are zero. The checker verifies that the augmented differential squares to zero and tests exactness at every point $y$ on the generators anchored at $x\le y$.

`shcellular C[R,K];` uses the certified resolution for a complex `K` on the same named space and field. The horizontal degree is the generator's resolution degree. The summand at a generator anchored at $x$ in stalk degree $q$ is $K^q(x)$. Use `shCellTerm[C,q,e,i]` for the $i$-th coordinate there, optionally followed by a signed coefficient inside the brackets. The horizontal differential combines the boundary coefficients with the restriction maps of `K`. The vertical differential carries the sign $(-1)^a$ in horizontal degree $a$.

Use `ShTot`, `ShHyper`, `ShE` and their usual operations with the name `C`. `shTotEmbed` and `shTotPart` use the chain coordinates of ordinary Roos cochains. Comparing the two projective resolutions of the constant sheaf gives the same hypercohomology and identifies the spectral sequences from $E_2^{a,q}=H^a(X;H^q(K))$ onward. [Cellular cochain example](examples/sheaf-cellular.asl)

With `cellular=1`, the transgression generator certifies a resolution with two generators in each degree. Over $\mathbb F_3$ and $\mathbb F_5$, the twisted example has rank-one $d_{64}$ from $(0,63)$ to $(64,0)$, and both groups vanish on page 65. With `untwisted=1` the differential is zero and both groups survive. Both versions run under the default limits. [Transgression generator](examples/sheaf-transgression.awk)

```sh
awk -v r=64 -v prime=3 -v cellular=1 -f examples/sheaf-transgression.awk | cargo run --quiet -- -
```

### The Leray complex of a finite map

`shleray K[f,F];` constructs the derived direct image of a sheaf $F$ on the named source of $f:X\to Y$. Its stalk in degree $q$ at $y$ is the actual Roos cochain space $C^q(X_y;F)$ with the Roos differential. The restriction maps keep the chains in the smaller inverse-image open. The declared complex uses the component names `K__0`, `K__1`, and so on. Its spectral sequence has $E_2^{p,q}=H^p(Y;R^q f_*F)$.

For an ordinary cochain $a$, `shLerayIncl[K](a)` places its family of restrictions in horizontal degree zero. This is a cochain map. For each chain $(x_0,\ldots,x_q)$ in $X$, the horizontal coefficient complex is the constant-coefficient Roos complex of the downset $\{y\mid y\le f(x_0)\}$. The downset has a largest point, so only horizontal degree zero contributes. Consequently the inclusion induces the comparison $H^n(X;F)\cong\mathbb H^n(Y;K)$.

The example projects an eighteen-point product onto a six-point model of the sphere. Both $R^0f_*F$ and $R^1f_*F$ are rank-one constant sheaves, and $d_2$ takes the generator of $E_2^{0,1}$ nontrivially to $E_2^{2,0}$. Replacing the upper-right entry of the matrix at `x2a<x4a` by zero splits the example and makes this differential zero. [Leray spectral sequence example](examples/sheaf-leray.asl)

### Smart truncation and a Deligne extension

`shtruncate T[K,c];` forms the smart truncation $\tau_{\le c}K$. It retains the sheaves below stalk degree $c$ and, when the cutoff meets an outgoing differential, makes `T__c` from the stalkwise kernels of that differential. The restriction matrices preserve these kernels. The preceding differential factors through them because its composite with the outgoing differential vanishes.

`shdeligne T[j,F,c];` performs one extension step $\tau_{\le c}Rj_*F$ for an ordinary sheaf on the source of the named open embedding $j$. The embedding must be injective, reflect the order, and have upward-closed image. The construction forms Roos cochains through degree $c+1$ when that degree exists, then computes the terminal kernels. The result is a complex of sheaves with the usual total cochains, hypercohomology, and spectral sequence.

The example suspends the finite product of two circles. For the constant sheaf on the open torus, the two tip stalks see the entire torus. Cutoff zero gives hypercohomology ranks $1,0,2,1$, agreeing with ordinary cohomology. Cutoff one gives $1,2,0,1$. These ranks agree in complementary degrees on the three-dimensional suspension. [Deligne extension of a torus suspension](examples/sheaf-deligne.asl)

### Derived direct images of bounded complexes

`shleraycomplex P[f,K];` uses a bounded complex $K$ on the source of $f$ as coefficients. In degree $q$ the stalk at $y$ is the degree-$q$ total Roos space of $K$ on $X_y$. Its differential includes the signed stalk differential of $K$. The cover restrictions keep the same chain coordinates on the smaller open. `shLerayIncl[P](z)` accepts a `ShTot(K,q)` cochain and puts all its restrictions in horizontal degree zero. The downset argument for an ordinary sheaf applies to each chain and coefficient degree of $K$, giving a quasi-isomorphism on total complexes.

`shdelignecomplex T[j,K,c];` extends a bounded complex along a named open embedding and applies the same smart cutoff. This permits successive Deligne steps with different cutoffs. The double-suspension example adds two tips at each step to the finite torus. Cutoffs $(0,1)$ give ranks $1,0,0,2,1$ and cutoffs $(1,1)$ give $1,2,0,0,1$. The ranks agree in complementary degrees on the four-dimensional model. [Iterated Deligne example](examples/sheaf-iterated-deligne.asl)

### Derived right Beck–Chevalley comparison

`shbcsquare Q[f,u,v,fp];` checks named monotone maps $f:X\to Y$, $u:X'\to X$, $v:Y'\to Y$ and $f':X'\to Y'$ with $fu=vf'$. For a sheaf $F$ on $X$, `shbc B[Q,F];` constructs the right comparison $v^*Rf_*F\to Rf'_*u^*F$. Use `shbccomplex B[Q,K];` for a bounded complex. The named complexes `B__from` and `B__to` live on $Y'$. Their stalks at $y'$ are the total Roos cochains on $P=f^{-1}(\uparrow v(y'))$ and $P'=f'^{-1}(\uparrow y')$ with the indicated coefficients.

For a strict flag in $P'$, the comparison reads the cochain on its image in $P$. A repeated image makes the flag degenerate and gives zero. For a strict image flag, the terminal stalk coordinate is unchanged. This is a cochain map and commutes with restriction on $Y'$. `shBCApply[B](a)` applies it to `ShTot(B__from,n)` or to `ShHyper(B__from,n)` via class representatives. `shbcmap M[B,q];` names the ordinary map of the two degree-$q$ cochain sheaves. `shbchom M[B,q];` names the induced map on the degree-$q$ cohomology sheaves `M__source` and `M__target`. [Normalized comparison example](examples/sheaf-base-change.asl)

`shbcrank[B,y,q];` gives the source dimension, target dimension, rank, kernel dimension and cokernel dimension at stalk $y$. `shbckernel[B,y,q];` and `shbccokernel[B,y,q];` give a nonzero witness when the indicated obstruction exists and an empty list otherwise. Witnesses use the source or target cohomology basis of `shbchom` and list pairs of index and coefficient.

`shbcuniversal[Q,p];` tests the derived comparison for **every** coefficient sheaf over $\mathbb F_p$. For each $y'\in Y'$ and $x\in P$, it computes the augmented homology of the comma poset $\{z\in P'\mid u(z)\le x\}$. All groups must vanish. This homological initiality criterion follows by resolving coefficients by injective diagrams. On the injective diagram that is $\mathbb F_p$ at points below $x$, the two local complexes compute the constant cohomology of the downset ending at $x$ in $P$ and of the comma poset, respectively. The former has a terminal point. Thus a nonzero augmented group gives a concrete coefficient obstruction, including a degree-zero kernel if the comma poset is empty. The report gives the first failing pair, the lowest failed augmented degree for that pair, and a signed nonboundary cycle.

`shbccomma[Q,p,y,x];` lists all nonzero augmented ranks for that pair. `shbcwitness[Q,p,y,x,d,i];` gives an actual chain representative of the indexed homology class in degree $d$. The chain is printed as pairs of coefficients and ordered point lists. The empty flag represents the augmentation coordinate in degree $-1$.

The generator takes the cone on a finite sphere or a triangulated projective plane and selects its link by a strict pullback square. For the sphere of dimension six, ordinary right Kan base change is invertible on the constant set family, and the first derived failure occurs in degree six. For the projective plane the universal comparison is exact at odd primes and has augmented obstructions in degrees one and two over $\mathbb F_2$. [Derived base-change generator](examples/sheaf-base-change.awk)

```sh
awk -v kind=sphere -v dim=6 -v prime=5 -f examples/sheaf-base-change.awk | cargo run --quiet -- -
awk -v kind=rp2 -v prime=2 -f examples/sheaf-base-change.awk | cargo run --quiet -- -
```

### Derived gluing on a finite cone

`shconeglue G[j,v,F,r] {c0,...};` takes an ordinary sheaf `F` on the named open site and a named open embedding `j` onto the complement of a point `v`. The point must be strictly below the whole open image. For a bounded complex `B`, use `shconegluecomplex G[j,v,B,r] {c0,...};`. Supply closed degree-$r$ cochains of the open input, as `ShC` or `ShTot` respectively. Each column is the image of a separate apex basis vector. Their number is the dimension of the apex stalk, concentrated in degree $r$.

The ambient complex `G` uses total Roos cochains on $U_{\ge u}$ at each open point $j(u)$, and the specified vector space at $v$. Between open points the matrices restrict cochains. The map from $v$ to $j(u)$ sends each basis vector to its chosen cochain restricted to $U_{\ge u}$. The unit from the original open complex puts a stalk element in horizontal degree zero by transporting it to all points of that principal open. Since the open has initial point $u$, this is a stalkwise quasi-isomorphism. The named complex works with `shcohomology`, `shextracomplex`, and the other sheaf-complex operations.

`G__fiber` is a second actual complex on `One` representing the costalk. The ambient cone has initial point $v$, so its derived global sections identify with the apex complex. The restriction to the punctured cone identifies with the attachment map under the open quasi-isomorphism. Its fiber therefore computes cohomology supported at the apex. In degree $n$ its coordinates are $A^n\oplus C^{n-1}$, where $C$ is the open total complex and $A$ is the degree-$r$ apex vector space. The differential is $d(a,b)=(0,\alpha a-Db)$. The coordinates are ordered with apex vectors before the open cochain basis. Use `shGlueApex[G,i]` for an apex basis vector, `shGlueLift[G](b)` to put an open cochain in the tail, and `shGlueTail[G](x)` to read that tail in positive degree. The native total-cochain and hypercohomology operations compute classes and representatives on this fiber.

For two gluings on the same cone and open complex, `shconehomotopy H[G,G2]{h0,...};` checks each degree-$r-1$ equation $Dh_j=c'_j-c_j$. Then `shConeTransfer[H](x)` is the chain isomorphism $(a,b)\mapsto(a,b+ha)$ between their fibers, and `shConeUntransfer[H](y)` is its inverse. They also send fiber hypercohomology classes by applying the map to a representative. In degree zero the list is empty and the two attachments must be identical.

On the finite $S^3$ link, attaching by zero or by a fundamental cocycle gives the same degree-three cohomology sheaf. The costalk ranks in degrees $1,3,4$ are $1,1,1$ for the zero attachment and $1,0,0$ for the fundamental attachment. The example computes both the explicit fiber and the relative Roos complex at the apex, and verifies the signed chain isomorphism for an exact change of cocycle over $\mathbb F_5$. [Derived cone gluing example](examples/sheaf-cone-gluing.asl)

### Perverse sheaves on a marked disk

`shdisk P[p,n,m] [u entries] [v entries];` declares the punctured-disk quiver over $\mathbb F_p$. The nearby space $V$ has dimension $n$ and the vanishing space $W$ has dimension $m$. Give the matrices $u:V\to W$ and $v:W\to V$ in row-major order. The declaration requires invertible nearby monodromy $T=I-vu$. The vanishing monodromy is $S=I-uv$. Use `shdiskmat[P,u];`, `shdiskmat[P,v];`, `shdiskmat[P,T];`, and `shdiskmat[P,S];` to see the matrices.

`shdiskmap f[P,Q] [A entries] [B entries];` checks the two naturality equations $u_QA=Bu_P$ and $v_Q B=Av_P$. The resulting `ShDiskHom(P,Q)` value uses the basis of the actual kernel of these simultaneous linear equations. Use `shdiskhomdim[P,Q];` and `shDiskGen[P,Q,i]` for its dimension and indexed basis, and `shdiskmapmat(f);` for the pair of matrices. Use `shDiskId[P]` and `shDiskComp(g,f)` for identities and composition. Hom values support `zero`, `add`, `pscale`, equality and finite enumeration. A matrix operator can also be declared as an ordinary native sheaf map between matching vector sheaves on `One` with `shdiskoperator M[P,T,V,V];`.

The declarations `shdiskshriek L[P];`, `shdiskstar R[P];`, and `shdiskmiddle I[P];` retain the same nearby representation. Write $K=I_V-T$. The first has $(u,v)=(I_V,K)$ and the second $(K,I_V)$. For the middle extension the vanishing space is $\operatorname{im}K$, with the corestriction of $K$ followed by inclusion. The basis of this image consists of the first independent columns of $K$. The maps `shDiskCanonical[L,I]` and `shDiskCanonical[I,R]` give the canonical epimorphism and monomorphism. Their composite is `shDiskCanonical[L,R]`. Canonical maps also go from a shriek extension to any quiver with the same based nearby data, and from such a quiver to the star extension. Use `shDiskExtMap[E_P,E_Q](f)` on two declared extensions of the same kind.

`shdiskdual D[P];` gives the dual quiver, with $u_D=v^{\mathsf T}$, $v_D=-T^{-\mathsf T}u^{\mathsf T}$, and $T_D=T^{-\mathsf T}$. Given declared duals of the endpoints, `shDiskDualMap[D_Q,D_P](f)` transposes both matrices. For the declared double dual `DDP`, the maps `shDiskBidual[P,DDP]` and `shDiskBidualInv[P,DDP]` have nearby components $-T$ and $-T^{-1}$ and identity on the vanishing space.

The stalk is $[V\xrightarrow{u}W]$ in perverse degrees $-1,0$, and the costalk is $[W\xrightarrow{v}V]$ in degrees $0,1$. The query `shdiskranks[P];` gives their dimensions and cohomology ranks. `shdiskstalk A[P];` and `shdiskcostalk B[P];` declare these as native complexes on `One`. The storage degrees are $0,1$ for the stalk and $1,2$ for the costalk, so both use the same shift.

With `shdiskrealize G[P,j,apex,a,b,c,d];` the named open embedding `j` exhibits a cone with that apex over the four-point circle whose only strict arrows run from $a,b$ to $c,d$. The generated sheaf `G__nearby` takes $V$ at each open point and uses $T$ on $b\to d$, with identity elsewhere. Its degree-zero diagonal and the cochains supported on $b\to d$ give the cellular complex $[V\xrightarrow{I-T}V]$ inside the actual Roos cochains. The apex complex attaches by $(I_V,v)$ because $vu=I-T$. The ambient complex is strictified on principal opens as in finite-cone gluing. Its named `G__costalk` model uses $v$. The map `shDiskCostalkLift[G](x)` takes its raw cochains to the relative ambient cochains supported at the apex. It sends $W$ to the apex flag and the last $V$ to the negative of the $b\to d$ coordinate on the apex-to-$b$ flag. The sign gives the stated differential $v$.

The example has two quivers with $T=I$ and the same stalk and costalk ranks. One has a vanishing transvection and the other has $S=I$, which distinguishes the extensions. It also verifies the factorization through the middle image, duality on a non-endomorphism, the actual relative costalk differential, and the braid relation between two rank-one Picard–Lefschetz monodromies. [Perverse disk example](examples/sheaf-disk.asl)

### Hurwitz moves and parabolic cohomology

`shdiskseq S[P,Q,...];` declares an ordered tuple of disk quivers over the same prime field, with one common based nearby space $V$. The vanishing spaces may differ. `shdiskprod[S];` reports the ordered product of the nearby monodromies. `shhurwitz H[S] [1,-2,...];` acts successively from left to right by the indicated signed braid generators, numbered from one. A positive generator takes $(P,Q)$ to $(Q,P')$, with $u_{P'}=u_P T_Q$ and $v_{P'}=T_Q^{-1}v_P$. The product is preserved and the corresponding vanishing monodromies are permuted. `shdiskpoint R[H,i];` gives the quiver at zero-based strand $i$. `shdiskseqsame[S,H];` compares the based quivers, and `shdiskmonosame[S,H];` compares the ordered based nearby matrices.

`shdiskseqmap f[S,T]{f0,f1,...};` takes one native `ShDiskHom` morphism per strand. All must have the same nearby linear map. If `HS` and `HT` apply the same signed word to `S` and `T`, `shdiskseqmapbraid fh[HS,HT,f];` transports the entire map to the new quivers. `shdiskseqmapat[fh,i];` reports both of its matrices at the specified strand.

If the product is the identity, `ShPar(S)` is the parabolic cohomology of the punctured sphere. Define
$$
V\xrightarrow{\partial}\bigoplus_i\operatorname{im}(I-T_i)\xrightarrow{\Sigma} V,\qquad
\partial x=((I-T_i)x)_i,\qquad
\Sigma(a)=\sum_i T_1\cdots T_{i-1}a_i .
$$
The middle cohomology is $\ker\Sigma/\operatorname{im}\partial$. Write a nearby vector with `shParVec[S][entries]` or a raw tuple with `shParRaw[S]{[entries at strand 0],...}`. Each tuple entry is given in the common nearby basis and must lie in its local image. The types are `ShParVec(S)` and `ShParRaw(S)`. Use `shParBoundary` and `shParSum` for the two arrows, `shParClass` on a cocycle, and `shParRep` for the chosen representative. `shParPrimitive` gives an actual nearby vector whose boundary is the input tuple, or diagnoses failure. `shparrank[S];` gives the cohomology dimension, `shParGen[S,i]` a basis class, and `shparrows(x);` and `shparvector(w);` report ambient coordinates. All three types support ordinary finite functions, `zero`, `add`, and `pscale`.

`shParBraid[H](x)` takes vectors, raw tuples, or classes from the declared source of `H` to `H`. On an adjacent pair the positive move sends $(a,b)$ to $(b,T_Q^{-1}(a+(T_P-I)b))$. This preserves the transported sum and takes each coupled boundary to the boundary of the same nearby vector. `shParMap[f](x)` uses the common nearby matrix of a sequence map. `shParIdentify[T](x)` uses the identity on ambient coordinates when the ordered based monodromies are equal.

For four scalar monodromies $-1$ over $\mathbb F_5$, the example checks the explicit action on a two-dimensional parabolic group. The first two generators both have order five and satisfy the braid relation, their product cubed is $-I$, and they do not commute. The example constructs the vanishing-basis isomorphism that returns the first braided tuple to its source. [Parabolic Hurwitz example](examples/sheaf-parabolic-hurwitz.asl)

### Multiplicative middle convolution

`shconv C[S,lambda];` applies multiplicative middle convolution to the ordered affine monodromies $A_1,\ldots,A_r$ in `S` over $\mathbb F_p$. The signed integer `lambda` represents a field element distinct from zero and one. The finite product may be arbitrary. On $V^r$, the generator $B_i$ fixes all rows except row $i$, whose blocks are $\lambda(A_j-I)$ for $j<i$, $\lambda A_i$ for $j=i$, and $A_j-I$ for $j>i$.

The resulting nearby space is the exact quotient $V^r/(K+L)$, where $K$ is the direct sum of the local fixed spaces in their respective blocks. The common fixed space $L$ of the $B_i$ is formed from $\ker(\lambda A_1\cdots A_r-I)$ by sending $v$ to $(A_2\cdots A_r v,\ldots,A_r v,v)$. These two relation spaces intersect trivially. The implementation checks that the quotient rank is $\sum_i\operatorname{rank}(A_i-I)+\operatorname{rank}(\lambda A_1\cdots A_r-I)-\dim V$. It declares an ordinary reusable disk sequence `C` with nearby operators induced by the $B_i$ and vanishing spaces given by the images of one minus those operators.

Use `shconvrank[C];` for the new rank, `shConvRaw[C]{[entries of first block],...}` for a vector of type `ShConvRaw(C)`, and `shConvGen[C,i]` for a zero-based basis vector of `ShConvVec(C)`. `shConvProject(x)` maps an ambient vector to the quotient and `shConvSection(y)` chooses its canonical lift. The generators `shConvRawT[i](x)` and `shConvT[i](y)` act on their respective spaces. Their indices are signed and one-based, with a negative index giving the inverse. `shconvcoords(x);` reports either ambient blocks or quotient coordinates. Both vector types support `zero`, `add`, `pscale`, and finite functions.

If `f` is a declared map between the original sequences, `shconvmap g[C,D,f];` descends its common nearby map acting diagonally on the blocks to a new disk sequence map. Both output sequences must use the same scalar. `shConvMap[g](x)` acts on raw and quotient vectors, and `shdiskseqmapat[g,i];` reads the induced nearby and vanishing maps at any strand.

For any disk sequence `S`, `shdiskinfinity[S];` reports the inverse ordered product. `shdiskrigidity[S];` gives the rigidity index, the sum of the dimensions of the local matrix centralizers including infinity, minus $(r-1)(\dim V)^2$. `shdiskgroup[S];` computes the order of the generated finite matrix group within the usual work and storage limits. The example takes two scalar signs over $\mathbb F_5$ to two rank-two transvections generating $\mathrm{SL}_2(\mathbb F_5)$ of order 120. The rigidity index remains two. Convolution by the inverse scalar returns the seed in this example, for both $\lambda=-1$ and $\lambda=2$. [Middle convolution example](examples/sheaf-middle-convolution.asl)

### Extraordinary inverse image and local supports

`shextra K[i,F];` and `shextracomplex K[i,A];` construct $Ri^!$ of an ordinary sheaf `F` or a bounded complex `A` on the named target of $i:Z\to X$. The map must be injective, reflect the order, and have downward-closed image. At $z$ let $U$ be the principal open above $i(z)$ and let $U^\circ$ be its complement of the image of $Z$. The stalk model consists of Roos cochains on $U$ supported on chains meeting $Z$. Since $Z$ is closed, these chains begin in $Z$. The differential includes the signed stalk differential of the coefficients. Restrictions keep the same chains on smaller opens. The resulting complex on $Z$ has components `K__0`, `K__1`, and so on, and `shcohomology H[K,n];` takes its cohomology sheaf.

`ShLocal(K,z,mode,n)` is the cohomology in degree $n$ for `absolute` ($U$), `punctured` ($U^\circ$), or `relative` ($(U,U^\circ)$). Use `shlocaldim[K,z,mode,n];` for its dimension and `shlocalgen[K,z,mode,n,i]` for a basis class. `shLocalClass[K,z,mode](c)` restricts an ambient cochain to the indicated model and takes its class. A relative input must vanish on the punctured open. `shLocalRep(h)` gives a canonical ambient cochain by extending its local representative by zero. Closure is understood on the relevant local complex. Local classes support `zero`, `add` and `shScale`.

`shLocalForget` maps relative classes to absolute classes, and `shLocalRestrict` maps absolute classes to punctured classes. The underlying cochain restriction to $U^\circ$ is surjective in each degree with the relative complex as its kernel. The connecting map `shLocalConnect(h)` takes the absolute differential of a punctured representative extended by zero, then takes the relative class in degree $n+1$. Changing the representative or the lift gives the same class.

On the cone of a finite torus, the relative ranks at the vertex in degrees zero through three are $0,0,2,1$ for constant coefficients. The degree-one connecting map takes a punctured basis to a relative basis in degree two. Applying the cutoff-one Deligne extension makes the degree-two relative group vanish while the degree-three class survives. [Local supports on the cone of a torus](examples/sheaf-support.asl)

### Degree-one higher Yoneda products

Give composable degree-one `ShExt` classes to `shMassey(a,b,c){0<2=>u,1<3=>v}` in the order of the maps $F_0\to\cdots\to F_n$. A closed degree-one `ShHomC` cochain can also be used as an input. The braces specify the proper nonadjacent intervals and omitted witnesses are zero. An entry $A_{ij}$ lies in `ShHomC(Fi,Fj,1)`. The evaluator checks $dA_{ij}+\sum_{i<k<j} A_{kj}\circ A_{ik}=0$ on every proper interval, using the signed composition of `shExtCompose`. `shCurvature` returns the resulting closed degree-two sum on the full interval. `shMassey` takes its class in `ShExt(F0,Fn,2)`.

`shMasseySolve(a,b,c)` finds the first defining system in a fixed search order and returns its class. `shMasseyWitness[i,j](a,b,c)` returns a proper entry of that system. `shInMassey(a,b,c,x)` tests membership of the final class `x` in the full set of outcomes. `shMasseyFind[i,j](a,b,c,x)` returns the indicated entry of the first system realizing `x`. The search solves each interval equation by prime-field elimination. It explores the cycle choices in degree-one cohomology, since a degree-zero change of basis absorbs exact differences in that interval and adjusts only longer intervals. `S.lo` means the finite search exhausted all choices. Reaching a computation or storage bound raises a diagnostic.

Two disjoint arms of lengths $r$ and $s$ between a common bottom and top give a one-dimensional $\operatorname{Ext}^2$ between the two endpoint skyscraper sheaves. The cover classes along each arm give a unique nonzero higher product. The two products differ by $(-1)^{s-r+1}$. In a proper interval $[i,j]$, the explicit witness has the long-edge coefficient $(-1)^{j-i-1}$. The `explicit=1` option checks these witnesses and compares the resulting class with the solver. The 64- and 35-input example runs over $\mathbb F_5$ with this option under the default limits. [Higher Yoneda generator](examples/sheaf-higher-yoneda.awk)

```sh
awk -v left=3 -v right=5 -v prime=5 -v explicit=1 -f examples/sheaf-higher-yoneda.awk | cargo run --quiet -- -
```

### Stalk cohomology and the Serre functor

`shcohomology H[K,q];` takes stalk cohomology of a bounded complex of sheaves. At each point it computes the kernel of the outgoing differential modulo the image of the incoming differential. A restriction transports a representative and projects it into the target stalk cohomology.

`shserre K[F];` computes the incidence-algebra derived Serre functor on an ordinary sheaf. Write $h$ for the maximum number of strict arrows in a chain, printed by `shserreshift[K];`. For each chain $x_0<\cdots<x_\ell$, the complex has a summand $I_{x_\ell}\otimes F(x_0)$ in declared degree $h-\ell$. Here $I_x(y)=\mathbb F_p$ for $y\le x$ and is zero elsewhere. The bar boundary alternates over deleted points. The first face applies $F(x_0\to x_1)$, the interior faces retain the coefficient, and the last face projects $I_{x_\ell}$ to $I_{x_{\ell-1}}$. The generated sheaves are named `K__0`, `K__1`, and so on. Declared degree $q$ corresponds to cohomological degree $q-h$ in $S(F)$.

On the linearly ordered site $A_n$, the Serre image of the projective interval $[a,n]$ is the injective interval $[1,a]$. If $b<n$, the image of $[a,b]$ is $[a+1,b+1][1]$. The example follows the simple at the first point of $A_3$ through four Serre steps and recovers it with accumulated shift two. [Sheaf Serre example](examples/sheaf-serre.asl)

`shserrecomplex T[K];` applies the same construction to a bounded complex of sheaves. A chain of length $\ell$ with coefficients in $K^q(x_0)$ contributes in declared degree $q+h-\ell$. The intrinsic degree is $q-\ell$, and `shserreshift[T];` prints $h$. The bar differential has the same faces as for an ordinary sheaf. The coefficient differential acts at the initial point $x_0$ with sign $(-1)^\ell$. The result retains the maps of the input complex, including the data that distinguishes nontrivial extensions from their stalk cohomology.

On a four-point diamond, the complex $P_b\oplus P_c\to P_a$ and its split cohomology complex both have stalk cohomology $P_d$ in degree zero and the simple at $a$ in degree one. The Serre image of the first has one intrinsic degree, with one-dimensional stalks at $a,b,c$ and zero at $d$. The split image has one-dimensional stalks at $d$ in intrinsic degrees $-1$ and $0$. These can be read using `shcohomology` and `shpull` to the named point. [Serre complex on the diamond](examples/sheaf-serre-complex.asl)

On the outward-oriented $D_4$ star, three successive applications of `shserrecomplex` recover the four simples, the projective at the center, and the highest-root representation, each in intrinsic degree $-2$ as in $S^3(F)\cong F[2]$. The checks compute the three complexes before taking stalk cohomology. [Tests over three prime fields](tests/sheaf_serre_complex.rs)

### Derived sheaf Hom and Yoneda products

For two sheaves `F` and `G` on the same named poset over the same prime field, `ShHomC(F,G,n)` has one matrix factor $\operatorname{Hom}(F(x_0),G(x_n))$ for every strict chain $x_0<\cdots<x_n$. This is Hom out of the normalized bar projective resolution of $F$. The differential uses the source restriction on the first face, the target restriction on the last face, and the unchanged matrix on the other faces, with alternating signs. If `F` is constant of rank one, this is the Roos complex of `G`.

`shHomTerm[F,G,[x0,...,xn],column,row]` gives a single basis matrix coordinate. An optional final signed coefficient is reduced in the prime field. The source column and target row are zero-based. `shMapC[m]` includes a declared natural map in degree zero, and `shHomId[F]` gives the identity. Use `shD`, `shClass`, `shRep`, and `shScale` on the appropriate grades. `ShExt(F,G,n)` is the resulting Ext group, with `shExtDim[F,G,n];` and `shExtGen[F,G,n,i]` giving its dimension and basis classes.

`shExtCompose(psi,phi)` composes two raw cochains or two classes with matching middle sheaf. For $\phi$ of degree $n$ and $\psi$ of degree $m$, it applies $\phi$ on the initial chain, then $\psi$ on the remaining chain, and multiplies by $(-1)^{mn}$. The sign gives $d(\psi\circ\phi)=d\psi\circ\phi+(-1)^m\psi\circ d\phi$. On classes this is Yoneda composition. In the diamond over $\mathbb F_3$, the two products along the branches are nonzero and opposite in $\operatorname{Ext}^2$. [Sheaf Yoneda example](examples/sheaf-yoneda.asl)

### Formal deformations of sheaf restrictions

`ShDef(F,n)` is the type of restrictions on the stalks of $F$ over $\mathbb F_p[t]/(t^{n+1})$ that reduce to those of $F$. The coefficients $A_e$ for $e>0$ are in `ShHomC(F,F,1)` and include every strictly comparable pair, including transitive arrows. Identity arrows stay fixed. Use `shDefStart[F]` for order zero, `shDefAppend(J,A)` to check and append the next coefficient, and `shDefTerm[k](J)` to read a known coefficient. Coefficient zero is the original restriction cochain. Equality compares every known coefficient. Small jet types can be tabulated.

For a valid prefix through order $n$, `shDefCurvature(J)` gives the coefficient of the next functoriality defect with the unknown coefficient set to zero. In the convention for `shExtCompose` it is $C=-\sum_{i+j=n+1} A_i\circ A_j$ for $i,j>0$. The interpreter checks $dA+C=0$ when appending. Associativity makes $C$ a closed $2$-cochain. `shDefObstruction(J)` is its class in `ShExt(F,F,2)`. When it vanishes, `shDefLift(J)` chooses a coefficient with $dA=-C$ by tracked elimination of the incoming differential. The obstruction concerns the specified prefix, and the lift makes one deterministic choice among the possible next coefficients.

An elementary change of stalk bases $G=I+t^e Q$ with $Q$ in `ShHomC(F,F,0)` and $1\le e\le n$ acts by $R'_{xy}=G_y R_{xy}G_x^{-1}$. Use `shDefGauge[e](J,Q)` for the transformed jet. Products of these changes give all truncated changes of basis that reduce to the identity. The operation `shDefGaugeWitness[e](J,Q)` gives a $1$-cochain whose differential is $C(J)-C(J')$, certifying the invariance of the obstruction class.

In the parallel-arm example every original nonidentity restriction is zero. The long arm has five covers with first-order coefficient $1$ and the other arm passes through a zero stalk. The program lifts through order four, then computes a nonzero obstruction in order five. The first possible product of the long-arm cover series is $t^5$ with coefficient $1$, independent of their higher coefficients. The other product is zero. A change in the next direct restriction affects both paths by the same coefficient, so it cannot eliminate the difference. [Parallel-arm deformation example](examples/sheaf-deformation.asl)

The generator checks all chosen lower lifts and the final class for any arm length that fits the interpreter's computation limits. Over $\mathbb F_5$, this command computes an obstruction in order thirty. [Parallel-arm generator](examples/sheaf-deformation-arms.awk)

```sh
awk -v r=30 -v p=5 -f examples/sheaf-deformation-arms.awk | nix run path:. -- -
```

### Kuranishi equations for sheaf deformations

`shkuranishi K[F,n];` computes the normalized deformation chart through total degree $n$ in the basis of $\operatorname{Ext}^1(F,F)$. To choose an ordered linear slice, give independent degree-one classes in braces after the declaration. The command `shkuranishiDims[K];` reports the chosen tangent dimension, the dimension of $\operatorname{Ext}^2(F,F)$, and the truncation. The coordinates $x_0,x_1,\ldots$ follow the chosen order.

The chart splits the Roos Hom complex using linear preimages of boundaries. For a $2$-cochain $q$, it subtracts the chosen preimage of $dq$ before projecting the resulting cycle. This gives a projection $\pi$ and homotopy $h$ on arbitrary cochains with $dh+hd=1-i\pi$. In each degree it computes $A=i(x)-hQ$, where $Q=-A\circ A$ is the quadratic functoriality defect using the sign of `shExtCompose`. The two outputs are the actual curvature $R=dA+Q$ and the obstruction $\kappa=\pi Q$.

The calculation verifies $dQ=A\circ R-R\circ A$ and $R=i\kappa+h(dQ)$ coefficientwise. The term $h(dQ)$ is linear in coefficients of $R$ from lower orders, so the residual and the obstruction coordinates generate the same ideal through the computed truncation. The split also gives $\pi R=\kappa$. All products use sparse commutative monomials and ordinary finite-field coefficients.

Use `shKurA[K,[e0,...]]`, `shKurR[K,[e0,...]]`, and `shKurKappa[K,[e0,...]]` for coefficients in the raw lift, raw curvature, and $\operatorname{Ext}^2$, respectively. The statement `shkuranishiEq[K,j];` prints the polynomial in obstruction coordinate $j$. The command `shKurArc[K,[c0,...],[w0,...],r]` substitutes $x_i=c_i t^{w_i}$ with positive weights, checks the equations through order $r$, and verifies the resulting restriction jet with the deformation append operation. Use `shKurArcKappa` with the same arguments for the obstruction coefficient at order $r$.

The two arms of the example give $x_0x_1-x_2x_3x_4$. A common tail makes a cubic curvature coefficient with nonzero ordinary differential. A weighted arc equalizes the two path products at order three and gives their common tail composite in order four. The diagonal slice gives an actual square term, including in characteristic two. [Kuranishi example](examples/sheaf-kuranishi.asl)

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

### The rho-Bockstein spectral sequence

`RhoBock(r, s, t, w, f)` is the component at filtration `f` of page `r` of the real-motivic cobar filtered by powers of $\rho$. Use `rhobockdim[r, s, t, w, f];` for its dimension and `rhobockgen[r, s, t, w, f, i]` for generator `i`. Page zero is the associated graded. The page-one block is $\rho^f\operatorname{MExt}(s,t+f,w+f)$. Its real-motivic tridegree includes the degree $(-1,-1)$ of each factor of $\rho$. The filtration is finite at fixed `t,w`, since the exponent of $\rho$ is at most $t-2w$.

`rhoClass[r,f](c)` checks the requested filtration and cycle condition on an `RMCobar` cochain. `rhoRep` gives a corrected real cochain representative, and `rhoD` applies the actual cobar differential and projects at $(s+1,t,w,f+r)$. `rhoNext` takes a cycle for that page differential to its next-page class. `rhoPrev` projects a chosen lift to the preceding page. `rhoMul` concatenates the corrected lifts of classes on the same page and adds their filtrations. It satisfies the Leibniz rule on each page.

The leading word of $d_1(\tau)$ is $\rho[\tau_0]$. For $k\ge1$ the leading word of $d_{2^k}(\tau^{2^k})$ is $\rho^{2^k}\tau^{2^{k-1}}[\xi_1^{2^{k-1}}]$. The example follows these classes through $d_{16}$, checks that each target is nonzero on its page, and checks that the target maps to zero on the next page. It also corrects the lift of $\rho^2\tau^3[\xi_1]$ by adding $\rho^3\tau^2([\xi_1\tau_0]+[\tau_1])$. [Rho-Bockstein example](examples/rho-bockstein.asl)

### Inverting rho on the real-motivic cobar

`RhoLoc(s, t, w)` is the normalized cobar after inverting $\rho$. Write $T_0=\tau$ for the initial vertex and $U_j=T_j-T_{j-1}$ for the increment across arrow $j$. The literal `rholocword[t,a]{(u,[xi-powers]), ...}` encodes $T_0^a$ followed by factors $U_j^{u_j}\xi^{M_j}$ with classical $\xi$ generators. Each factor has positive degree. If $Q$ is the sum of their classical $\xi$ degrees, then $t-w=a+\sum_j u_j+Q$ and the implied coefficient is $\rho^{2Q-t}$. This leaves a finite basis at each signed tridegree. On localized cochains `rhoShift[k]` accepts signed exponents.

The first coface transports $T_0$ to $T_0+U_1$. Internal cofaces split $U_j$ primitively and use the classical Milnor coproduct on $\xi$. The normalized `cobarD` discards unit arrows. `concat` transports the right coefficient across the left arrows by $T_0\mapsto T_0+\sum U_j$. `rhoLocalize` sends the original real cochains into these coordinates using $\tau_{0,j}=U_j/\rho$ and $\tau_{i+1,j}=(\tau_{i,j}^2+T_j\xi_{i+1,j})/\rho$.

`rhoProject` sets $T_0$ and all $U_j$ to zero, mapping into `Cobar(s,t-w)`. The classical target requires an unsigned degree. `rhoSection[t]` embeds classical words as pure $\xi$ words at topological degree `t`. The maps commute with `cobarD` and `concat`. Writing them as $P$ and $i$, one has $Pi=1$.

On positive-length cochains `rhoHomotopy` gives a homotopy $H$ lowering length by one. A monomial contributes when its initial vertex exponent is zero and its first factor with positive $U_j$ exponent has a unit $\xi$ factor. The homotopy drops that slot and substitutes $T_0+\sum_{k<j}U_k$ for its $U_j$. It satisfies $dH+Hd=1+iP$ and $PH=Hi=H^2=0$ over $\mathbb F_2$. This gives the localized cohomology identification with classical Ext. For example the real arrow $[\tau_0]$ becomes the boundary of the vertex cochain $\rho^{-1}T_0$, and the real arrow $[\xi_1]$ maps to the classical generator in internal degree one. [Explicit rho-localized contraction](examples/rho-localization.asl)

### Real-motivic Ext and exact rho torsion

`RMExt(s, t, w)` is the cohomology of the normalized real-motivic cobar at signed tridegree $(s,t,w)$. `rmextdim[s,t,w];` computes its dimension and `rmextgen[s,t,w,i]` selects a generator. Use `extClass` on a real cocycle, `extRep` to read a representative, and `extPrimitive` on a real cobar boundary. The Yoneda product uses balanced concatenation of representatives. `rhoShift[j]` multiplies by $\rho^j$ for a nonnegative exponent. `realComplex` specializes a class to `MExt(s,t,w)`, and `rhoProject` gives its localized image in classical `Ext(s,t-w)`.

`rhoorder(x);` accepts a real cocycle or an `RMExt` class. It reports the least exponent $j$ with $\rho^j x=0$ in cohomology, zero for the zero class, or `infinite` if the localized image is nonzero. Use `rhoNull[j](x)` at positive cobar length to obtain a witness $y$ at degree $(s-1,t-j,w-j)$ with $dy=\rho^j x$, using the chosen representative when $x$ is a class. The calculator checks that the requested exponent is minimal.

At fixed $k=t-w$, a raw basis word of intrinsic degree $T$ has $\rho$ exponent $T-t$ and $\tau$ exponent $W-T+k$. Once $t\le 0$, the condition $T-t\ge0$ is automatic. Every further multiplication by $\rho$ is an isomorphism of these normalized cochain slices. Thus finite torsion is killed by exponent $\max(t,0)$. The exact tests use the incoming differential, with a direct scalar test in length one. In the example the real classes $h_0$ and $\rho h_1$ have the same tridegree and different behavior under the two comparison maps. The real class $h_0^2$ has exponent one, and the raw cochain $[\tau_0^{32}]$ has exponent 32 with witness $\tau^{32}$. [Real-motivic Ext example](examples/real-motivic-ext.asl)

### Graded rho-module bars

`rmbars[s,k];` decomposes real-motivic Ext as a graded $\mathbb F_2[\rho]$-module at fixed $k=t-w$. Each entry gives birth degree $B$ and length $r$, or `infinite` for a free interval. A finite interval contributes at $B-r<t\le B$, and a free interval contributes at every $t\le B$. Multiplication by $\rho$ lowers $t$ and $w$ by one. `rmbarRank[s,k,t,j];` counts intervals containing both $t$ and $t-j$, giving the rank of multiplication by $\rho^j$ in that component.

The calculation uses the normalized cobar at $t=0,w=-k$. Filtering its basis by intrinsic topological degree $T\ge t$ recovers the real slice at degree $t$. The two adjacent differential matrices are reduced in decreasing order of $T$. A reduced incoming column pairs its source degree $D$ with its pivot degree $B$, giving a finite interval of length $B-D$ when that difference is positive. An outgoing dependency with unpaired leading coordinate gives a free interval. The tracked column combinations give real cocycles and nulls, with their own topological degrees.

Use `rmbarGen[s,k,i]` for the bar's `RMExt` generator and `rmbarRep[s,k,i]` for its traced real cocycle. For a finite bar, `rmbarNull[s,k,i]` gives the paired cochain $u$ with $du=\rho^r z$, where $z$ is the traced representative. The null starts in filtration zero in topological degree $B-r$ and supports the differential $d_r$ to the filtration-$r$ class of $\rho^r z$. The example verifies this correspondence for length six in $(s,k)=(2,5)$ and length seven in $(3,6)$. [Real-motivic rho bars example](examples/real-motivic-bars.asl)

### Multiplicative rho extensions

`rmbarCoords(x)` expresses an `RMExt(s,t,w)` class in the basis of bars at fixed $k=t-w$ that live at $t$. Its `RMBar(s,t,w)` result has terms $\rho^{B-t}b_i$, where $i$ is the index in `rmbars[s,k]` and $B$ is that bar's birth. `rmbarTerm[s,t,w,i]` selects a living term. `rmbarClass` reconstructs the original class and `rmbarLift` gives the actual sum of shifted bar cocycles. `RMBar` supports `add`, `zero`, and finite enumeration.

`rhodepth(x);` prints the smallest exponent among the occupied bars, which is the largest $a$ such that $x$ is in the image of multiplication by $\rho^a$. It prints `infinite` for zero. `rhoDivide[a](x)` chooses a preimage in topological degree $t+a$ and weight $w+a$ using the same bar indices. This choice omits any bars that die under the shift, and the interpreter checks the result with `rhoShift[a]`.

`rhoEinf[f](x)` takes a class of depth at least $f$ to $F^fH/F^{f+1}H$ on the stable `RhoBock` page. For a nonempty slice it uses page $1+\max(0,t-2w)$. Use `rhoEinf[r,f](x)` to request a later or equal page when comparing grades with `rhoMul`. `rhohidden(a,b);` reports a `hidden` product when the nonzero Yoneda product has depth larger than the sum of the two input depths. It reports `zero` for a vanishing product and `visible` when the depths add. For a permissible division of the product, `rhoProductNull[j](a,b)` gives a cochain with differential equal to the sum of the product of `extRep` lifts and $\rho^j$ times the `rmbarLift` of the chosen quotient's coordinates.

The example gives $[\tau_0][\tau_0^2]$ exact depth one and order one. The two leading classes multiply to zero in filtration zero on the stable page, while their real-motivic product has a nonzero filtration-one class. It also checks depth two for $[\xi_1][\tau_0^4]$ and $[\tau_0^2][\tau_0^4]$. [Real-motivic multiplicative extensions](examples/real-rho-extensions.asl)

### Real-motivic Ext modulo powers of rho

`RMRhoCobar(e,s,t,w)` is the normalized real cobar modulo $\rho^e$ for $e\ge1$. `RMRhoExt(e,s,t,w)` is its cohomology, with dimension `rmrhoextdim[e,s,t,w];` and basis class `rmrhoextgen[e,s,t,w,i]`. `rhoReduce[e]` maps real cochains and classes to these types. It also reduces from a larger or equal truncation to a smaller cutoff. The differential and balanced concatenation in the quotient discard terms at exponent $e$ or above. Use `extClass`, `extRep`, `yoneda`, `add`, `zero`, and `rhoShift` with the quotient types. At $e=1$, `realComplex` identifies these types with `MCobar` and `MExt`. `rhoReduce[1]` on complex-motivic cochains and classes gives the inverse identification.

`rhoTruncLift` reads the normalized real cochain representing a truncated cochain or class. For a class $x$ modulo $\rho^e$, `rhoConnectRep(x)` computes $d\widetilde x/\rho^e$ from that lift. Its class `rhoConnect(x)` is the connecting map $\delta_e$ in `RMExt(s+1,t+e,w+e)`. In the long exact sequence, this map follows reduction and is followed by multiplication by $\rho^e$. For $a\le e$ one has $\delta_a(\operatorname{red}_a x)=\rho^{e-a}\delta_e(x)$. The internal Bockstein `rhoBeta(x)` is `rhoReduce[e](rhoConnect(x))`. It satisfies the Leibniz rule and squares to zero.

`rhoBetaRep(x)` returns the actual truncated cochain $d\widetilde x/\rho^e$ reduced modulo $\rho^e$. When this cochain is a boundary, `rhoBetaNull(x)` returns a truncated primitive with differential `rhoBetaRep(x)`. The primitive has type `RMRhoCobar(e,s,t+e,w+e)`. Use `extPrimitive` to find a primitive of any exact truncated cochain of positive cobar length. Both operations solve against the incoming differential and report a diagnostic for nonboundaries. A closed lift modulo $\rho^{2e}$ also gives a null by dividing its cochain difference from the canonical lift by $\rho^e$. The interpreter checks the difference and the resulting differential explicitly. [Explicit internal Bockstein null](examples/real-rho-bockstein-null.asl)

The example computes a three-term non-scalar cycle with nonzero $\beta_2$. The scalar $\tau^{32}$ is closed modulo $\rho^e$ for $1\le e\le32$. Its typed connecting class has representative $\rho^{32-e}[\tau_0^{32}]$ and exact real $\rho$-order $e$. The example checks the connecting class at $e=3$ and $e=32$, and computes the nonzero quotient differential of $\tau^{32}$ at $e=33$. [Real-motivic rho quotient example](examples/real-rho-quotient.asl)

### Algebraic Frobenius in the motivic cobar

`rhoFrob(x)` squares each term in the commutative tensor algebra of its cobar slots. It acts on real and complex motivic cochains and Ext classes. The cobar length stays fixed and the topological degree and weight double. In the real cobar it squares the scalar coefficients and normalizes every arrow square with the right-unit transport from the preceding arrows. The operation is additive, commutes with the normalized differential and balanced product, and sends a shift by $\rho^a$ to a shift by $\rho^{2a}$. It commutes with `realComplex` and induces `extSq[0]` after complex realization.

On `RMRhoCobar` and `RMRhoExt` Frobenius doubles the cutoff $e$. It commutes with reduction across doubled cutoffs, the internal `rhoBeta`, and `rhoConnect`. On `RhoBock(r,s,t,w,f)` it squares the lift and projects to `RhoBock(2r,s,2t,2w,2f)`. In these degrees `rhoD(rhoFrob(x))` equals `rhoFrob(rhoD(x))`, including page zero. The nonzero class $[\xi_1]^4$ maps to zero both in real motivic Ext and from the first to the second page.

`rhoConnectRep(x)` also accepts a closed truncated cochain with its termwise real lift. For a truncated cochain or class $x$, `rhoFrobConnectNull(x)` gives a real cochain $U$ with $dU=F(\delta_e^{\rm rep}x)+\delta_{2e}^{\rm rep}(F x)$. For cochains it divides the discarded part of the squared lift by $\rho^{2e}$. For a class, it also lifts a primitive of the boundary between the reduced square and the canonical target representative. If $B$ is that lift, the numerator is $F(\widetilde x)+\widetilde{F x}+dB$, which is divisible term by term by $\rho^{2e}$. The example checks nonzero page differentials from $\tau$, $\tau^2$, and $\tau^4$ and an explicit two-slot homotopy. [Motivic Frobenius example](examples/rho-frobenius.asl)

### Real-motivic triple Massey products

Use `rmTriple(a,b,c)` for the affine coset of the ordered triple Massey product. Each input may be an `RMExt` class or a closed `RMCobar` cochain. The adjacent balanced products $ab$ and $bc$ must be boundaries. `rmMasseyNull(a,b)` gives an actual primitive of $ab$ in positive total cobar length. `rmMassey(a,b,c)` gives the corresponding cocycle $uc+av$, where $du=ab$ and $dv=bc$. This raw operation needs only the defining system. For an adjacent product of length zero the product must already vanish and its negative-length primitive contributes zero.

The `RMTriple(s,t,w)` value has indeterminacy $aH^{s_b+s_c-1}+H^{s_a+s_b-1}c$ in the relevant signed tridegrees. `rmTripleBase(M)` gives a class in the coset and `rmTripleDir[i](M)` gives a direction in the indeterminacy basis. `rmTripleDim(M);` prints its affine dimension and `rmTripleContains(M,x);` checks membership. A candidate may be an `RMExt` class or a closed cochain in the output grade.

For example $\langle\rho,[\tau_0],[\xi_1]\rangle$ contains $[\tau_0^2]$ and $[\tau_0^2]+\rho^2[\xi_1^2]$. For $N$ a power of two the defining cocycle for $\langle[\tau_0^N],\rho^N,[\tau_0^N]\rangle$ is $\rho^N[\tau_0^{2N}]$, with both primitives equal to $\tau^N$. The scalar indeterminacy vanishes and the class has exact $\rho$ exponent $N$. The example computes the full coset at $N=4$ and checks the sparse defining cocycle and its minimal null at $N=32$. [Real-motivic Massey product example](examples/real-motivic-massey.asl)

### Higher operations on real-motivic Ext

On homogeneous `RMExt` classes, `ainfty(a1,...,an)` computes the transferred operation $m_n$, and `alift(a1,...,an)` computes the comparison map $I_n$ into the real cobar. Their lengths are $\sum s_i+2-n$ and $\sum s_i+1-n$, respectively. Both add the topological degrees and the weights. The real overload accepts one input, giving $m_1=0$ and $I_1=\texttt{extRep}$. The binary operation is Yoneda multiplication. The unit is strict, and scalar classes such as $\rho$ are ordinary inputs. Negative intermediate lengths contribute zero. A requested value in negative length is outside the supported types.

The transfer uses a trigraded contraction over $\mathbb F_2$. For the linear section $q$ of the differential on boundaries, choose the first linearly independent columns in the ordered real cobar basis. Projection first replaces $x$ by the cycle $x+q(dx)$ and then takes its `RMExt` coordinates. The homotopy applies $q$ to this cycle plus the representative of its class. The higher operations sum the ordered binary trees using the balanced cobar product. The contraction has $dh+hd=1+ip$ and $h^2=hi=ph=0$. A section linear over $\mathbb F_2[\rho]$ would send the nonzero $\rho$-torsion class $h_0$ into a torsion-free cochain group, so the transfer is taken over $\mathbb F_2$.

When the two adjacent products vanish, the chosen $m_3$ lies in the coset computed by `rmTriple`. The example checks both the coset with positive indeterminacy and the singleton $\langle[\tau_0^4],\rho^4,[\tau_0^4]\rangle$. With the stated contraction, $m_4(\rho,h_0,\rho,h_0^2)$ and $m_5(\rho,h_0,h_0,h_1,h_0^2)$ are nonzero. These values depend on the contraction. [Real-motivic transfer example](examples/real-motivic-transfer.asl)

### The rho-linear minimal complex

`RhoMin(s,t,w)` represents a sparse graded polynomial $\sum_f \rho^f x_f$ with $x_f$ in `MExt(s,t+f,w+f)`. Use `rhoMin[f](x)` to insert a coefficient and `rhoMinCoef[f](x)` to extract it. `zero`, `add` and `rhoShift` act on these polynomials. The differential `rhoMinD` increases $s$ by one and strictly increases filtration. Its homology is real-motivic Ext. The module uses the normal-word section from the complex cobar and extends the complex-motivic cohomology contraction coefficientwise over $\mathbb F_2[\rho]$.

Write $d_R=d_0+Q$, where $d_0$ preserves $\rho$ exponent and $Q$ strictly raises it. With the complex contraction $i,p,h$, the minimal differential is $\Delta=pQ(1+hQ)^{-1}i$. `rhoMinI` and `rhoMinP` are the corrected inclusion and projection, $(1+hQ)^{-1}i$ and $p(1+Qh)^{-1}$. `rhoMinH` on positive-length real cochains gives $H=h(1+Qh)^{-1}$. The sums terminate at fixed tridegree. They satisfy $\Delta^2=0$, $PI=1$, $I\Delta=d_RI$, $\Delta P=Pd_R$ and $d_RH+Hd_R=1+IP$, with $H^2=HI=PH=0$.

A coefficient of $\Delta$ gives a page differential only after the lower-filtration obstructions have been removed. The example compares the actual $d_6$ and $d_7$ with the images of corrected real lifts. In tridegree $(3,3,-4)$ a term at filtration seven cancels a differential component at filtration eight, leaving the length-seven differential from filtration three to filtration ten. The same example computes the sole $\rho^{32}$ coefficient of $\Delta(\tau^{32})$ and lifts it to the actual real boundary. [Rho-linear minimal complex example](examples/rho-minimal.asl)

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

## Dwork unit roots of Legendre elliptic curves

For an odd prime $p$ and $\lambda\in\mathbb F_p\setminus\{0,1\}$, `legendreTrace[p,lambda];` counts points of $y^2=x(x-1)(x-\lambda)$ over $\mathbb F_p$ using the quadratic character and includes the point at infinity. It reports the trace $a_p=p+1-\#E$ of Frobenius on $H^1(E)$. The parameter uses its representative in $0,\ldots,p-1$.

At precision $s\ge1$, put $z=\operatorname{Teich}(\lambda/16)\pmod{p^s}$ and
$B_r(t)=\sum_{0\le n<p^r}\binom{2n}{n}^2t^n$.
The query `legendreDwork[p,lambda,s];` evaluates both $B_s(z)$ and $B_{s-1}(z)$ modulo $p^s$. It gives the ratio and its signed Frobenius unit root $(-1)^{(p-1)/2}B_s(z)/B_{s-1}(z)$ for an ordinary curve. For a supersingular curve the output reports the two sums and the trace and marks the ratio undefined. Independently, `legendreHensel[p,lambda,s];` lifts the unit root of $T^2-a_pT+p$ obtained from the point count. The Dwork query checks equality with that lift at the requested precision.

`dworkPeriod[p,z,r,s];` evaluates the raw truncation $B_r(z)$ for a signed $128$-bit integer representative $z$ and keeps truncation order $r$ separate from precision $s$. Here $B_0=1$. The central binomial recurrence splits off powers of $p$ before dividing by $n+1$. It retains the unit part when a term vanishes and reuses it when the valuation later decreases. Moduli fit in $u64$, with widened products and the normal computation limits for enumerating terms.

`dworkCongruence[p,s];` checks every coefficient of
$B_{s+1}(t)B_{s-1}(t^p)\equiv B_s(t)B_s(t^p)\pmod{p^s}$.
The form `dworkCongruence[p,s,z];` computes and reports the actual two residues at an integer $z$, including values with $z^p\ne z$. The coefficientwise check uses the polynomial storage limit.

For $p=7$ and $\lambda=3$ the example gives trace $4$, raw ratio $10$ and signed unit root $39\pmod{49}$. The unit at precision six reduces to the same residue. The change $\lambda\mapsto1-\lambda$ gives the quadratic twist by $-1$ and changes the trace to $-4$. The example also computes the unit $169738\pmod{5^8}$ for $p=5,\lambda=2$. [Legendre Dwork example](examples/legendre-dwork.asl)

### Unramified coefficients and Legendre closed points

`rvunram W[p,s]{h0,...,1};` declares the truncated unramified ring $(\mathbb Z/p^s)[X]/(H)$, where $H$ is monic with the listed coefficients in increasing degree and has irreducible reduction. This presents $W_s(\mathbb F_{p^f})$ with $f=\deg H$. Use `rvUnramElem[W]{a0,...}` for elements of type `RVUnram(W)` and omit trailing zero coordinates. Values support `zero`, `add`, `neg`, equality, `rvUnramMul`, `rvUnramInv` on units, and `rvUnramPow[k]` for nonnegative $k$. The degree is at most 16, and $p^s$ must fit in 64 unsigned bits. Small finite rings can be enumerated within the usual limit.

`rvUnramFrob[k](x)` applies the $k$th power of the unramified Frobenius. It fixes integer scalars, and its image of $X$ is the Hensel lift of $X^p\bmod p$ among the roots of $H$. For $H=X^2-X+2$ at $p=5$, this image is $1-X$, whereas $X^5=6-X\pmod{25}$. `rvUnramTeich(x)` returns the unique lift of the input residue satisfying $x^{p^f}=x$. On those lifts Frobenius is $p$th powering. `rvUnramTrace` and `rvUnramNorm` give the sum and product of the Frobenius conjugates as constant elements of the same type. `rvUnramVal(x);` gives the minimum coefficient valuation, displaying `>=s` for zero.

The Legendre operations take a value of `RVUnram(W)` whose residue $\lambda$ is different from $0$ and $1$ at odd $p$. `rvLegendreParameter(lambda)` lifts $\lambda/16$ to a Teichmüller parameter $z$. With $B_s(z)=\sum_{j=0}^{p^s-1}\binom{2j}{j}^2z^j$, `rvLegendreFactor[i](lambda)` gives the signed semilinear factor $(-1)^{(p-1)/2}B_s(\sigma^i z)/B_{s-1}(\sigma^{i+1}z)$ for an ordinary curve. Frobenius permutes its $f$ factors cyclically. `rvLegendreUnit(lambda)` multiplies them to descend the unit root to the scalar ring. `rvLegendreData(lambda);` reports the two period truncations at $z$, their valuations at the given precision, and the ordinary or supersingular status. Singular residue parameters and supersingular unit-factor requests have distinct diagnostics.

`rvLegendreTrace(lambda);` counts the curve independently over $\mathbb F_{p^f}$ by a quadratic character sum. `rvLegendreHensel(lambda);` lifts the unit root of $T^2-a_qT+q$ from that trace when it exists, using $q=p^f$. These queries enumerate the residue field. The period calculation uses the exact central binomial recurrence with the existing statement work limits. For the degree-two point $\lambda=X$ at $p=5$, take $H=X^2+2$. Modulo $5^5$ the two semilinear factors are $3039+49X$ and $3039+3076X$, with scalar product $2823$. The example compares this value with the independent point-count unit and checks transport of a nonconstant factor under a change of polynomial lift. [Closed-point Legendre example](examples/legendre-closed-points.asl)

### Split quadratic Legendre Frobenius orders

`rvlegendrecm K(lambda);` counts the curve at an ordinary Legendre parameter in the residue field $\mathbb F_q$ of its `RVUnram` type. If $a$ is the trace, the declaration registers `RVQpQuad(K)` with radicand $D=a^2-4q<0$ and the split embedding $\sqrt D\equiv a\pmod p$. Use `rvCMPi[K]` for $\pi=(a+\sqrt D)/2$, `rvCMPi[K,m]` for $\pi^m$, and `rvQpQuadConj` for the conjugate. These elements work in `rvphodgequad`. `rvQpQuadSpecialize[W](x)` maps any element integral at the chosen embedding to a scalar in `RVUnram(W)` of the same residue characteristic, at the precision of `W`. This specialization agrees on $\pi$ with the independent Dwork norm `rvLegendreUnit(lambda)` in the source ring.

`RVCMOrder(K,m)` for $m\ge1$ uses the integral basis $(1,\pi^m)$ of $\mathbb Z[\pi^m]$. Use `rvCMTerm[K,m,b,c]` for $b+c\pi^m$, with `add`, `neg`, `zero`, equality, and `rvCMMul` for multiplication. `rvCMLift(x)` embeds the order element in `RVQpQuad(K)`. If $m$ divides $n$, `rvCMInclude[m](x)` takes an element of `RVCMOrder(K,n)` to its exact coordinates in `RVCMOrder(K,m)`. These inclusions compose along towers and preserve multiplication. The order is the subring generated by the specified Frobenius power. Computing the full endomorphism ring of the curve is a separate arithmetic calculation.

Write $u_0=0$, $u_1=1$, and $u_{m+1}=a u_m-q u_{m-1}$. In the original basis one has $\pi^m=u_m\pi-q u_{m-1}$. Thus `rvCMOrder[K,m];` reports the trace, index $|u_m|$ inside $\mathbb Z[\pi]$, and discriminant $D u_m^2$. `rvCMIndex[K,m,n];` reports $|u_n/u_m|$ when $m$ divides $n$. `rvCMTrace[K,m];` and `rvCMNorm[K,m];` query the trace and norm separately. The interpreter computes powers and inclusion indices by quadratic arithmetic without forming $q^m$. Integer coefficients and reports are bounded by the supported signed 128-bit range.

For $\lambda=X$ over $\mathbb F_{25}$ with $X^2=-2$, the chosen embedding gives $2823\bmod 5^5$. At powers $m=3,4,6$ the traces are $142,866,-11086$ and the absolute indices are $21,92,2982$. In the order at power three the generator from power six has coordinates $(-15625,142)$. The example checks the tower and its multiplication and also compares the two unit calculations for $\lambda=2$ over $\mathbb F_5$. [Legendre Frobenius order example](examples/legendre-cm-orders.asl)

### Higher Hasse–Witt matrices of hyperelliptic curves

`rvhyperelliptic C[p]{f0,...,1};` declares the chosen integral lift $f(x)$ for the smooth odd-degree model $y^2=f(x)$ over $\mathbb F_p$. Coefficients are signed integers in increasing degree. The interpreter checks that $p$ is odd, $f$ is monic of degree $2g+1$, and $f\bmod p$ is squarefree. The declaration retains the integral lift for the coefficient calculation. The point counts depend on its reduction.

`rvHW[C,r,s]` returns $A_r\bmod p^s$, where
$$
A_r(i,j)=[x^{p^rj-i}]f(x)^{(p^r-1)/2},\qquad 1\le i,j\le g,
$$
with $A_0=I$. The calculation uses exact modular polynomial powering, truncating to the necessary coefficient range. `rvHWStatus[C];` reports the genus and whether $A_1\bmod p$ is invertible. Raw matrices are available on smooth nonordinary curves. On ordinary curves `rvHWUnit[C,s]` computes the ordered matrix $M_s=A_sA_{s-1}^{-1}$ with **both** matrices at precision $s$. `rvHWConverges[C,s];` compares its reduction with the independently computed $M_{s-1}$ for $s\ge2$.

The matrix type is `RVHWMatrix(p,g,s)`, with row-major literal `rvHWMat[p,g,s]{entries}`. Use `rvHWCompose(A,B)` for $AB$ and `rvHWInv(A)` when the determinant is a unit modulo $p$. `rvHWEntry[i,j](A)` uses one-based indices. `rvHWChar[k](A)` gives the coefficient of $T^{g-k}$ in $\det(TI-A)$ by a division-free recurrence. Both return values in `RVHWScalar(p,s)`, with literal `rvHWScalar[p,s,c]`. Both types support zero, addition, negation, equality and precision reduction `rvHWReduce[s](x)`. Scalars also support `rvHWInv` on units and `rvHWVal(x);` for the $p$-adic valuation, with `>=s` for zero. The modulus $p^s$ fits in unsigned 64 bits, and the polynomial and matrix calculations use the interpreter's work and storage limits.

`rvHWPoints[C,r];` counts points over $\mathbb F_{p^r}$ by the quadratic character of $f(x)$ and includes the single point at infinity. It uses deterministically selected irreducible polynomials for the finite extensions. `rvHWWeil[C];` computes the descending coefficients of the Frobenius characteristic polynomial $P(T)$ from the first $g$ point counts using exact Newton identities and the reciprocal functional equation. `rvHWCertify[C,s];` verifies by monic division that $\det(TI-M_s)$ is a unit factor of $P(T)\bmod p^s$ and checks that the complementary factor reduces to $T^g$ modulo $p$. It prints the complementary factor and remainder. This certification uses a point-count calculation independent of the coefficient matrices. Finite-field enumeration and exact integer arithmetic remain subject to their resource limits.

For the genus-two curve in the example, $M_3=\begin{pmatrix}16&4\\13&11\end{pmatrix}\bmod 27$, and its ordered characteristic factor $T^2+16$ divides $T^4-2T^2+9$ modulo $27$. The reverse product gives a different matrix. The example includes a noncommuting genus-three calculation, a genus-two certificate modulo $3^7$, and the Legendre unit $39\bmod 49$. [Hyperelliptic Hasse–Witt example](examples/hyperelliptic-hasse-witt.asl)

### Lang maps and exact CM recognition

`rvECFrobenius[E]` gives the $q$-power endomorphism $(x,y)\mapsto(x^q,y^q)$ of a curve over the named field of size $q$. `rvECLang[E]` gives $1-\pi$. Its kernel is exactly the rational points over that field. `rvECSum` and `rvECNegate` use the group law on rational functions, including doubling, opposite maps, and the constant map to infinity. The two summands must have the same named source and target. `rvECChart[E,F](u)` gives the scaling $(x,y)\mapsto(u^2x,u^3y)$ when it takes the source equation to the target equation.

`rvECExtend[E',F'](map)` transports the coefficients of a rational map along the deterministic residue-field embedding and checks both extended curve equations. This uses the same rule for embedding defining roots as the closed-point Hasse–Witt calculation.

For an ordinary curve, `rvECCMRecognize[K](map)` returns the actual `RVQpQuad(K)` element $a+b\pi$ corresponding to a rational endomorphism. The field `K` must have the curve's Frobenius data. The calculation measures the degree $n$ and trace $1+n-\deg(1-\alpha)$ of the map. Degree polarization with $\pi$ gives $\operatorname{tr}(\alpha\pi)$, and the two traces determine $a,b$ by the trace pairing in the basis $(1,\pi)$. The operation checks the actual normalized rational functions in $[D]\alpha=[A]+[B]\pi$ for a common positive denominator $D$. `rvECCMCertify[K](map);` reports the two traces, degree, and checked integers $A,B,D$.

On $y^2=x^3+x+7$ over $\mathbb F_{11}$, the composite of the explicit degree-three and degree-five Vélu maps is $1-\pi$ of degree $15$. The same map over $\mathbb F_{121}$ has coordinates $(14+\pi_{121})/3$. On $y^2=x^3+2$ over $\mathbb F_{13}$ the scaling by $u=3$ gives the automorphism $(1+\pi)/3$ with trace $-1$ and degree $1$. [Lang and CM example](examples/elliptic-cm-recognition.asl)

### Hasse–Witt matrices over closed points

`rvhyperelliptic C[W]{f0,...,fd};` defines $y^2=f(x)$ over the residue field of a named `RVUnram(W)`. The coefficient expressions supply integral lifts in that same ring, in ascending degree order. The residue characteristic is odd, the leading coefficient is a unit, and the reduction of $f$ is squarefree. Both odd and even degrees are supported, with genus $g=\lfloor(d-1)/2\rfloor$.

`rvHWW[C,r]` gives the higher Hasse–Witt matrix with entry $(i,j)$ equal to $[x^{p^r j-i}]f^{(p^r-1)/2}$ at the precision of `W`, with $A_0=I$. If the curve is ordinary and `W` has precision $s$, `rvHWWRatio[C]` gives $M_s=A_s\sigma(A_{s-1})^{-1}$ using the ring's lifted Frobenius $\sigma$. With residue degree $h$, `rvHWWUnit[C]` computes the ordered product $P=M_s\sigma(M_s)\cdots\sigma^{h-1}(M_s)$. It satisfies $\sigma(P)=M_s^{-1}PM_s$, so its characteristic coefficients descend to integer scalars. Raw matrices and point counts are also available on nonordinary curves.

These values have type `RVHWWMatrix(W,g)` and support `zero`, `add`, `neg`, equality, and `rvHWWCompose(A,B)` for the product $AB$. Use `rvHWWInv` on matrices invertible modulo $p$, `rvHWWFrob[k]` for entrywise $\sigma^k$, `rvHWWEntry[i,j]` with indices starting at one, and `rvHWWChar[k]` for descending characteristic coefficient $k$ starting at zero. The entries and coefficients are `RVUnram(W)` elements. `rvHWWReduce[V]` reduces to a ring with the same defining polynomial at a lower or equal precision. `rvHWWDescend` sends a matrix with scalar entries to `RVHWMatrix(p,g,s)`, in particular for comparison with the prime-field operators at residue degree one.

`rvHWWPoints[C,r];` counts the smooth projective curve over $\mathbb F_{p^{hr}}$. It embeds the coefficient field by finding a root of its defining residue polynomial in the counted extension. The count includes the quadratic fiber at infinity for even degree. `rvHWWWeil[C];` reconstructs the descending Weil polynomial from the first $g$ extension traces. `rvHWWVerify[C];` checks that the descended characteristic polynomial of the ordered product divides this independent Weil polynomial at the precision of `W`. `rvHWWVerify[C](A);` checks a supplied matrix of the same shape, including unit determinant and descent of its characteristic coefficients.

The genus-two example over $\mathbb F_{27}$ gives Weil coefficients $[1,-2,34,-54,729]$ and ordered unit coefficients $[1,25,7]\bmod27$. Reversing the factors gives $[1,10,7]$, which fails the division. Changing one integral lift by $3$ changes the semilinear matrix and ordered product and preserves the certified characteristic coefficients. [Closed-point Hasse–Witt example](examples/closed-hasse-witt.asl)

### Rational cyclic kernels and normalized Vélu maps

`rvec E[W](a,b);` declares the short Weierstrass curve $E:y^2=x^3+ax+b$ with coefficients of type `RVUnram(W)`. The named ring has precision one and odd characteristic $p\ge5$. The interpreter checks that $4a^3+27b^2\ne0$. Use `rvECCurve[E]` for the value of type `RVECCurve(E)`, and `rvECA[E]` and `rvECB[E]` for its coefficients. `rvECPointCount[E];` counts the rational points, including infinity, and reports the Frobenius trace.

For an odd prime $\ell\ne p$, `rveckernel[E,ell]{[c0 coordinates],...,[1]}` gives a value of type `RVECKernel(E,ell)`. The lists give the field coordinates of the coefficients of a monic polynomial $K$ in ascending powers of $x$. Empty coordinate lists give zero. The degree is $(\ell-1)/2$. The interpreter checks squarefreeness and divisibility by the $\ell$-division polynomial, then evaluates the multiplication-by-$g$ rational function modulo $K$ for a generator $g$ of $\mathbb F_\ell^*/\{\pm1\}$. Closure under this action identifies one cyclic subgroup, since every geometric root is the $x$-coordinate of a nonzero point of order $\ell$. Coefficients in the declared field make the subgroup Frobenius-stable.

`rvECIsogenies[E,ell];` enumerates monic candidates within the work limit and lists the actual cyclic kernels. `rvECKCoeff[k](K)` returns the coefficient of $x^k$. An explicit kernel literal is checked with division-polynomial recurrences in the quotient by $K$. It can therefore succeed when enumerating the field or all monic candidates exceeds the work limit. A bare kernel type determines whether its domain is empty from the Frobenius characteristic polynomial. A certified kernel also supplies a witness that the domain is inhabited.

`rvecisogeny I[E,Q](K);` constructs the normalized Vélu map and registers its quotient curve under the fresh name `Q`. Use `rvECIso[I]` for the value of type `RVECIso(I)` and `rvECGetKernel` to recover its kernel. The map has $x'=R(x)$ and $y'=yR'(x)$, with
$$
R(x)=x+\sum_{r:K(r)=0}\left(\frac{6r^2+2a}{x-r}+\frac{4(r^3+ar+b)}{(x-r)^2}\right).
$$
The calculation uses Newton power sums and polynomial identities over the declared field. It expresses $R=N/K^2$ and $R'=(N'K-2NK')/K^3$. The coefficient queries `rvECXNum[k]`, `rvECXDen[k]`, `rvECDNum[k]`, and `rvECDDen[k]` take the isogeny as input and return the indicated coefficient. The constructor checks the cleared identity $f(N'K-2NK')^2=N^3+a'NK^4+b'K^6$ for the quotient coefficients $a',b'$. Use `rvECCertify[I];` to repeat the check.

`rvECPoint[E](x,y)` validates an affine point and `rvECInfinity[E]` is its group identity, both of type `RVECPoint(E)`. Use `rvECX` and `rvECY` for affine coordinates, `add` and `neg` for the group law, and `rvECScale[n]` for signed multiplication. `rvECApply(i,P)` gives the actual image point on the quotient, including infinity for the kernel points and for infinity itself.

The two curves over $\mathbb F_{13}$ in the example both have $19$ rational points. They have one and four rational degree-three kernels respectively. The example also computes the quotient for the irreducible kernel polynomial of the degree-five isogeny $x^2+3x+1$ over $\mathbb F_7$ and a kernel with nonconstant coefficients over $\mathbb F_{25}$. [Cyclic Vélu example](examples/cyclic-velu.asl)

### Local endomorphism orders and Vélu edge directions

`rvEllFrobCert[E,ell,k,a];` tests the arithmetic Frobenius of a named short Weierstrass curve on the full geometric $m$-torsion scheme, where $m=\ell^k$ for an odd prime $\ell\ne p$. It uses the monic $x$-division polynomial $\psi_m$ and reduces $x^q$ and $f(x)^{(q-1)/2}$ modulo it for $E:y^2=f(x)$ over $\mathbb F_q$. The candidate $a$ is an integer prime to $m$. Omitting it uses half the Frobenius trace modulo $m$. The interpreter chooses a congruent signed representative prime to $p$.

If the unreduced multiplication $x$-map is $N/D$, the invariant differential gives the $y$-multiplier $(N'D-ND')/(aD^2)$. The certificate prints the actual reduced Frobenius coordinates and the two remainders obtained by clearing these denominators. The denominators are units on the nonzero $m$-torsion since $a$ and $m$ are coprime. Polynomials appear as lists of exponent and coefficient pairs, with each coefficient in the named field basis. `rvEllFrobCoordinates[E,ell,k,a]` returns the pair of zero tests. This checks the sign even when the two candidate $x$-maps agree.

`rvEllLocalOrder[E,ell];` counts the trace of an ordinary curve and finds the largest $r$ for which Frobenius is scalar on $E[\ell^r]$. It searches through $\lfloor v_\ell(t^2-4q)/2\rfloor$, using the containment of successively larger torsion schemes. The resulting index exponent is $v_\ell([\operatorname{End}(E):\mathbb Z[\pi]])$. If $a$ is the scalar at that level, $\pi-a$ factors through multiplication by $\ell^r$, giving an endomorphism $\omega=(\pi-a)/\ell^r$. The report gives its integral trace, norm and discriminant and determines the order at $\ell$.

`rvellipticcm K[E];` declares the split Frobenius field of an ordinary Weierstrass curve for the existing quadratic and Frobenius-order operations. `rvEllLocalGen[K,E,ell]` gives $\omega$ in the compatible `RVQpQuad(K)` field. `rvEllEdge[I];` computes the two local indices of the named prime-degree Vélu isogeny and reports ascending, horizontal or descending according to the target index minus the source index.

Over $\mathbb F_{61}$ the curve $y^2=x^3+2$ has scalar $5$ on its complete degree-$40$ $9$-division polynomial. The negative candidate agrees on $x$ and fails on $y$. The example verifies $\pi=5+9\omega$ with trace $-1$ and norm $1$ for $\omega$. It also classifies five actual degree-three Vélu edges over $\mathbb F_{13}$. [Ordinary elliptic order example](examples/ordinary-tate-orders.asl)

### Dual isogenies and coprime Vélu compositions

`rvECMap(i)` turns a named Vélu isogeny into a value of type `RVECMap(E,F)` with nominal source and target curves in the same field. Such a map has coordinates $(R(x),yS(x))$ and carries infinity to infinity. Both rational functions are reduced with monic denominators. Equality compares both functions exactly. `zero[RVECMap(E,F)]` is the constant map to infinity. `rvECMap{XNum,XDen,YNum,YDen}[k](f)` denotes the four queries `rvECMapXNum[k](f)`, `rvECMapXDen[k](f)`, `rvECMapYNum[k](f)`, and `rvECMapYDen[k](f)`. They give the coefficient of $x^k$ in the indicated polynomial. `rvECMapDegree(f);` gives the degree, including inseparable degree. `rvECMapCertify(f);` verifies the exact function-field equation $f_E(x)S(x)^2=R(x)^3+a_F R(x)+b_F$.

`rvECMapCompose(f,g)` applies $g$ first and then $f$, and `rvECMapApply(f,P)` acts on rational points, with poles mapped to infinity. `rvECId[E]` is the identity and `rvECMultMap[E,n]` constructs the signed multiplication map $[n]$ by division polynomials, with $[0]$ the constant map. This includes multiplication by the field characteristic. `rvECChart[E,F](u)` gives the isomorphism $(x,y)\mapsto(u^2x,u^3y)$ when the named curves have the corresponding coefficients and $u\ne0$. Use the chart with $u=1$ to compare maps into two names for the same equation.

`rvECDual(f)` takes a map with a certified prime-degree cyclic kernel, including an already constructed dual. It computes the dual kernel as the minimal polynomial of $R$ in the finite algebra $\mathbb F_q[x]/(\psi_\ell/K)$, where $K$ is the kernel polynomial. It constructs the normalized Vélu quotient and then scales back to the original equation. For a normalized input the scale is $1/\ell$ and the resulting dual pulls the invariant differential back by $\ell$. `rvECMapKCoeff[k](f)` reads coefficients of the certified kernel. The example checks both composites against $[\ell]$ and takes the dual twice, comparing the complete rational functions each time.

`rvECPush(i,K)` takes a named Vélu isogeny and a kernel of a different prime degree on its source. It returns the kernel on the quotient by computing the minimal polynomial of the image $R$ in $\mathbb F_q[x]/K$. This uses the geometric kernel even if it has no rational generator. Over $\mathbb F_{11}$, $E:y^2=x^3+x+7$ has kernels $x+3$ and $x^2+3x+7$. Quotienting by them in either order reaches $y^2=x^3+3x+6$, with identical degree-$15$ maps after identifying the two target names. Using $x+7$ in place of $x+3$ returns to the original curve and gives an equal pair of degree-$15$ endomorphisms. Polynomial lengths and computation use the interpreter's resource limits. [Dual and coprime Vélu example](examples/dual-velu.asl)

### Jacobians of hyperelliptic curves over finite fields

`RVJac(C)` is the group of degree-zero divisor classes of a smooth monic odd-degree curve $C:y^2=f(x)$ with one point at infinity. Declare the curve with `rvhyperelliptic C[W]{coefficients};` for a named precision-one unramified field `W` of odd characteristic. A class is stored in its unique reduced Mumford form $(u,v)$, with $u$ monic, $\deg v<\deg u\le g$, and $u\mid f-v^2$. The identity is $(1,0)$.

`rvjac[C]{[u0 coordinates],...,[1]}{[v0 coordinates],...}` lists both polynomials in ascending powers of $x$. It checks the norm divisibility and applies Cantor reduction to semireduced input of arbitrary feasible degree. An empty coordinate list gives the zero field element and an empty polynomial list gives the zero polynomial. `rvJacPoint[C](x,y)` represents $[(x,y)-\infty]$ as $(X-x,y)$. Use `rvJacU[k]` and `rvJacV[k]` to read the two polynomials. `add` composes using the extended gcd of $u_1,u_2,v_1+v_2$, including common or repeated places. It then replaces an overdegree pair by the monic normalization of $(f-v^2)/u$ and reduces $-v$ modulo the new $u$. Use `neg`, `zero[RVJac(C)]`, or `rvJacScale[n](d)` for signed scalar multiplication. `rvJacCompose(a,b);` shows the actual semireduced pair and its exact norm quotient $(f-v^2)/u$ before reduction.

`rvjacchart N[C](a);` sends the finite branch point $(a,0)$ of a smooth odd- or even-degree hyperelliptic model `C` to infinity. The value `a` belongs to the named precision-one field of `C` and satisfies $f(a)=0$ with $d=f'(a)\ne0$. ASL registers `N` as a monic odd-degree curve over the same field in the coordinates $z=d/(x-a)$, $w=yz^{g+1}/d$. Its equation is $w^2=z^{2g+2}f(a+d/z)/d^2$. Use `RVJac(N)` for its Picard group of degree zero. If `C` is already monic of odd degree, `rvJacRechart[N](D)` transports a class from `RVJac(C)` to `RVJac(N)` and `rvJacUnchart[N](D)` transports it back. Both operations carry the finite support and account for the change of base point at infinity. [Runnable example](examples/hyperelliptic-chart.asl)

`rvJacHalfPoint[C](a,b){[r0 coordinates],...}` returns a class $H$ with $2H=[(a,b)-\infty]$. Supply the coefficients of $r$ in ascending powers, as for `rvjac`. The operation checks that $(a,b)$ lies on $C$, that $r^2=a-X$ in $K[X]/(f)$, and that the norm of $r$ to $K$ is $-b$. It computes the monic gcd $R(z)=\gcd(f(a-z^2),z-r(a-z^2))$ over $K$. Writing $R(z)=E(z^2)+zO(z^2)$ gives the Mumford pair $(\operatorname{monic}O(a-X),-E(a-X)\bmod u)$. The interpreter verifies the result by Cantor doubling. The example over $\mathbb F_3$ uses an irreducible quartic branch factor, two sign choices that give distinct halves, and a ramified input with $b=0$. [Jacobian point halving example](examples/hyperelliptic-jacobian-halving.asl)

`rvJacHalfFiber[C](a,b);` reports the empty fiber or the dimension $r$ and size $2^r$ of the rational halves of $[(a,b)-\infty]$. `rvJacHalfAt[C,k](a,b)` constructs the half at zero-based decimal index $k$ without a square-root witness. The algorithm factors $-f(a-z^2)$ over the named field and pairs each monic factor with its monic transform under $z\mapsto-z$. At a ramified input it removes $z^2$ before factoring and inserts $z$ into the selected product. A self-paired factor gives an empty fiber. For $b\ne0$, the first odd-degree pair is reserved to enforce $R(0)=b$. Factors are ordered lexicographically by their ascending coefficients in the named field basis, with the smaller member first in each pair. The bits of $k$, starting at the least significant bit, choose among the remaining pairs in that order. Each selected product gives a reduced Mumford pair whose doubling is checked. [Rational half-fiber example](examples/hyperelliptic-half-fiber.asl)

`rvJacCount[C];` counts $J(C)(\mathbb F_q)$ by evaluating the independently point-counted Weil polynomial $P(T)$ at $T=1$. `rvJacEnumerate[C];` counts all reduced pairs directly by the same divisibility condition, without Weil coefficients. The ordinary work limit restricts direct enumeration to small fields and genera. `rvJacOrder(d);` determines the exact order by factoring the group order within the work limit and testing multiples of the actual class.

`rvjacbase D[C,W];` registers the base extension of `C` to `W` with the deterministic residue-field embedding used by closed Hasse–Witt point counts. `rvJacExtend[D](d)` applies this embedding to the Mumford coefficients. `rvJacDescend[C](d)` solves for their preimages, failing if a coefficient is outside the subfield. The curves must match after coefficient extension. On `D`, `rvJacFrob[C,k](d)` applies the $k$th power of the $q$-power Frobenius for the base field of `C`. `rvJacWeilApply[C](d)` evaluates the full characteristic polynomial $P$ by composing Frobenius and adding signed multiples of the class.

The example gives $P(T)=T^4-2T^2+9$ and independently counts eight classes over $\mathbb F_3$. Over $\mathbb F_{27}$ the Jacobian has 776 rational classes. The class of $(z,z+z^2)$ has order 194 and a Frobenius orbit of length three. The example checks its doubling at a repeated place and shows that $F^4d-2F^2d+9d=0$ by actual arithmetic on that class. [Hyperelliptic Jacobian example](examples/hyperelliptic-jacobian.asl)

`rvJacQuadTrace[C](k){A}{B}{P}{Q}` traces a Mumford pair $U=A+\omega B$, $V=P+\omega Q$ to `RVJac(C)`. The nonzero radicand `k` belongs to the named field $K$ of `C`, with $\omega^2=k$. List the four polynomials in ascending powers using bracketed field coordinates as in `rvjac`. ASL checks that $U$ is monic and divides $f-V^2$, and reduces $V$ modulo $U$. For nonsquare $k$ the result adds the class and its quadratic Frobenius conjugate and descends to $K$. For square $k$ it adds the classes from the two substitutions for $\omega$. The calculation uses polynomial arithmetic over $K$, including at common and repeated places. [Quadratic Jacobian trace example](examples/hyperelliptic-quadratic-trace.asl)

`rvJacWeil2(T,U)` computes the Weil pairing on two classes killed by 2. `rvJacTate2(T,D)` computes the reduced quadratic Tate pairing with any divisor class `D` and a class `T` killed by 2. The two arguments must be on the same named curve. Both expressions return $1$ or $-1$ in its named `RVUnram` field. The Weil sign comes from the parity of the intersection of the even subsets of branch points, including infinity. The Tate pairing evaluates the function $u_T$ using polynomial resultants and its local leading unit at shared branch points, then takes the quadratic character of the field. The example over $\mathbb F_7$ checks both pairings and includes overlapping support in genus two. [Quadratic Jacobian pairings example](examples/hyperelliptic-pairings.asl)

`rvJacTate[n](T,D)` computes the reduced Tate pairing for $n\ge2$ dividing $q-1$, where $q$ is the order of the named field. The arguments must lie on the same named curve and $nT$ must be zero. The result is an $n$th root of unity in that field. Miller's algorithm uses the factors from Cantor composition and reduction, evaluating their local leading units at each place of `D` and at infinity. This accommodates shared support, intermediate zeros and poles, and repeated points. The example checks third and fifth roots of unity using explicit functions on curves of genus one and two. [Miller Tate pairing example](examples/hyperelliptic-miller-tate.asl)

`rvJacWeil[n](T,U)` computes the Weil pairing of two classes killed by $n$ on the same named curve, for $n\ge2$ coprime to the field characteristic. It returns an $n$th root of unity in the named field, including when only some of the $n$th roots are rational. It divides the two Miller evaluations and applies the tame reciprocity sign at shared places, counting residue degrees and multiplicities. Explicit functions give the third root $4$ over $\mathbb F_7$ in genus one and the fifth root $5$ over $\mathbb F_{11}$ in genus two. [Miller Weil pairing example](examples/hyperelliptic-miller-weil.asl)

`rvJacDivisible2(d)` is `S.hi` exactly when $d\in 2J(C)(\mathbb F_q)$. The étale Kummer map sends $d=(u,v)$ to the square class of $(-1)^{\deg u}u$ modulo $f/c$, where $c=\gcd(u,f)$ is monic. On the common branch factors $c$, it uses $(-1)^{\deg u-1}(u/c)(f/c)$, the square class obtained by evaluating at the ramification points. Distinct-degree factorization and modular exponentiation test the coordinates in the finite fields attached to the irreducible factors of $f$.

`rvJacKummer(d);` displays the monic product $B_d$ of the branch factors with nonsquare coordinates and counts these factors by degree. `rvJacKummerSupport[k](d)` gives the coefficient of $X^k$ in $B_d$ in the curve's named field, with zero beyond its degree. A divisible class has $B_d=1$. Adding classes takes the symmetric difference of their nonsquare branch factors. The example computes this cancellation for classes supported at ramification points and finds a class with nonsquare coordinates in both a quadratic and a cubic branch orbit. [Hyperelliptic Kummer example](examples/hyperelliptic-kummer.asl)

#### Certified torsion bases and Frobenius

`rvjactorsionbasis B[C,ell](zeta){e1,...,eg,f1,...,fg};` registers a symplectic basis of the full geometric $J[\ell]$ on `C`. The level $\ell$ is prime and differs from the residue characteristic. The root `zeta` has order $\ell$ in the curve's named field. The declaration checks that $\ell$ kills every listed class and that the Weil pairing is trivial within each half and equals $\zeta^{\delta_{ij}}$ on $(e_i,f_j)$. The resulting nondegenerate pairing matrix certifies that the $2g$ classes span $J[\ell]$.

`rvJacTorsionCoord[B,i](T)` returns coordinate `i` in the listed order, starting at zero, as `PrimeScalar(ell)`. It checks $\ell T=0$ and takes a discrete logarithm to base `zeta` of the appropriate Weil pairing by baby-step giant-step. `rvJacTorsionFrob[B,Base,k];` reports the matrix of the $k$th power of `Base` Frobenius, with images in columns. It also gives the multiplier $|K_{\mathrm{Base}}|^k\bmod\ell$ and the characteristic polynomial in descending order. The interpreter checks that each column reconstructs the actual divisor-class image and that $M^tJM=|K_{\mathrm{Base}}|^k J$ over $\mathbb F_\ell$. The example gives a nontrivial shear on genus-two 2-torsion and the matrix $\operatorname{diag}(1,2)$ on elliptic 3-torsion. [Certified Jacobian torsion example](examples/hyperelliptic-jacobian-level.asl)

#### Projective Kummer surface

`RVKummer(C)` is the projective quartic $J(C)/\{\pm1\}$ of a smooth monic quintic over a named odd finite field. `rvKummer[C](k0,k1,k2,k3)` checks the quartic equation and normalizes the first nonzero coordinate to one. `rvKummerCoord[i](K)` returns coordinate $i$ for $0\le i\le3$. The identity image is $(0:0:0:1)$ and `zero[RVKummer(C)]` gives this point.

`rvJacKummerPoint(D)` sends a Mumford class to the quartic and gives the same point for $D$ and $-D$. In degree one the coordinates are $(0:1:x:x^2)$. For $u=X^2-sX+p$ in degree two, write $q=(f-v^2)/u$. The coordinates are $(1:s:p:-q_0-pq_2)$, including at repeated support. `rvKummerLiftable(K)` tests for a rational lift to $J(C)$ and `rvKummerLift(K)` chooses a lift using the lexicographically smaller square root in the named field basis. Point enumeration includes rational points whose lifts lie on the quadratic twist. [Projective Kummer surface example](examples/hyperelliptic-kummer-surface.asl)

`rvKummerScale[n](K)` computes the image of multiplication by the signed integer `n` on any rational Kummer point. `rvKummerDiffAdd(A,B,D)` returns $\kappa(P+Q)$ from $A=\kappa(P)$, $B=\kappa(Q)$ and $D=\kappa(P-Q)$, where $\kappa$ is the quotient map. All three inputs must have the same named Kummer type. The interpreter uses Cantor arithmetic on a common monic model of the curve or its quadratic twist. It checks `D` against the two addition candidates. The identity and nodes work with either twist class. Two inputs away from the nodes with different twist classes have no rational difference on the quartic. [Kummer arithmetic example](examples/hyperelliptic-kummer-arithmetic.asl)

`rvKummerExtend[D](K)` embeds the normalized coordinates of `K` into the field of the coefficient extension `D`. `rvKummerDescend[C](K)` recovers those coordinates over `C` when they all lie in its field. This also works for twist-only points whose extension divisor lifts change sign under Frobenius. `rvKummerFrob[C,k](K)` applies the $k$th power of the relative Frobenius on the extension of `C` that contains `K`. The equation checks and residue embeddings are the same as for `RVJac`.

#### The marked Kummer theta group

For a certified genus-two level-2 basis `B`, a mask `a` from 0 to 15 labels the two-torsion class $T_a$ by its bits in the listed basis order, starting with the least significant bit. `rvKummerTranslate[B,a](K)` sends $\kappa(P)$ to $\kappa(P+T_a)$ on the named Kummer surface, including points whose lifts lie on the twist. The interpreter recovers the translating projectivities from all sixteen marked nodes and reuses them for the basis.

`rvKummerThetaEntry[B,a,i,j]` gives entry $(i,j)$ of the translating matrix $M_a$, with rows and columns numbered from zero and the first nonzero entry in row order equal to one. `rvKummerThetaFactor[B,a,b]` gives the scalar $c(a,b)$ in $M_aM_b=c(a,b)M_{a\mathbin{\mathrm{xor}}b}$. The quotient $c(a,b)/c(b,a)$ is the Weil pairing of $T_a$ and $T_b$. `rvKummerTropeCoeff[B,a,i]` reads a coefficient of the translated theta plane, normalizing its first nonzero coefficient to one. The plane at mask zero is $w=0$. [Theta group example over an extension field](examples/hyperelliptic-kummer-theta.asl)

#### Arithmetic of the resolved Kummer K3

`rvKummerK3Zeta[C];`, `rvKummerK3Points[C,r];`, and `rvKummerK3Picard[C,r];` use the same monic quintic, with $r\ge1$. The zeta function of the resolution has denominator $(1-T)(1-q^2T)\det(1-TF\mid\wedge^2H^1(J))\prod_d(1-(qT)^d)$ over $\mathbb F_q$. The displayed orbit lengths $d$ come from the zero class and the unordered pairs of branch points, including infinity. The two factors from eigenvalue $q$ in the abelian part are displayed separately from its cross quartic. Over $\mathbb F_{q^r}$ the point report gives the singular quartic count, the rational node count and the resolved count. Each rational node adds $q^r$ points on resolution. The calculation uses the base Weil polynomial and branch factor degrees.

`RVKummerK3(C)` is the finite type of rational points on the smooth resolution. `rvKummerSmooth(K)` takes any nonsingular Kummer point, including a point whose lifts lie on the quadratic twist. At a node, `rvKummerExceptional(K,h0,h1,h2,h3)` takes a projective tangent direction on the quadratic tangent cone. The interpreter removes the component along `K`, normalizes the first nonzero coordinate, and tests the coefficient of $\epsilon^2$ in the homogeneous quartic. Use `rvKummerBlowdown(X)` to recover the base Kummer point and `rvKummerIsExceptional(X)` to test whether it lies over a node. `rvKummerDirection[i](X)` reads coordinate $0\le i\le3$ of an exceptional direction. Over a field of size $q$, finite enumeration replaces each rational node with its $q+1$ tangent directions. [Kummer resolution example](examples/hyperelliptic-kummer-resolution.asl)

`rvKummerK3Extend[D](X)`, `rvKummerK3Descend[C](X)`, and `rvKummerK3Frob[C,k](X)` transport resolved points under the same field operations as the Kummer quartic. They act on both the normalized base coordinates and the exceptional direction. Descending an exceptional point requires both to lie over the named subfield.

`rvKummerK3Translate[B,a](X)` lifts translation by mask $a$ in the certified genus-two level-2 basis `B` to the resolution. It applies the same theta matrix to the base point and its exceptional direction. `rvKummerTropeParam[B,a](u,v)` parametrizes the entire strict trope by a nonzero projective pair $(u:v)$, including its six intersections with the exceptional curves. `rvKummerOnTrope[B,a](X)` tests the trope plane equation on the base and, at an exceptional point, also on the tangent direction. [Resolved trope example](examples/hyperelliptic-kummer-resolved-tropes.asl)

`rvKummerSwitch(X)` is the Gauss involution on `RVKummerK3(C)` for the odd-theta polarity $(w,s,p,t)\mapsto(t,-p,s,-w)$. It uses the quartic gradient at a smooth point and the Hessian applied to the direction at an exceptional point. When the image lies over a node, it solves the restricted Hessian equation to recover the image direction. The switch is defined on the full resolution over the named field and exchanges each exceptional curve with its paired strict trope. The example switches two distinct smooth points over $\mathbb F_3$ and exchanges a twist-only point and an exceptional direction over $\mathbb F_9$. [Gauss switch example](examples/hyperelliptic-kummer-switch.asl)

`rvKummerEnriquesPoints[B,C,a,r];` counts the Enriques quotient of the smooth Kummer K3 by the even theta switch of mask $a$. The certified genus-two level-2 basis `B` is on `C` or an extension of it. The mask must be fixed by Frobenius of `C` and lie outside the six masks $S_0$ on the trope defined by Weierstrass infinity. The report over the degree-$r$ field counts split fibers with two rational lifts and nonsplit fibers with Frobenius-conjugate lifts. The twisted cover count is twice the nonsplit count. [Enriques quotient example](examples/hyperelliptic-kummer-enriques.asl) [Ohashi, §5, Proposition 5.1](https://www.kurims.kyoto-u.ac.jp/preprint/file/RIMS1651.pdf) [Catanese, §§5 and 7](https://arxiv.org/abs/2101.10501)

If $M$ is normalized Frobenius on the node masks, let $f$ count the fixed nodes of $M^r$ and $m$ count masks $v$ with $v\mathbin{\oplus}M^r(v)\mathbin{\oplus}a\in S_0$. The Frobenius trace on the invariant divisor space is $t_r=(8+2f-m)/4$. On the span of the hyperplane class and the sixteen exceptional classes, the projector $(1+\sigma_a)/2$ averages the traces $1+f$ and $3-m/2$. In the switched trace the hyperplane contributes $3$, and each node incidence counted by $m$ contributes $-1/2$ to the diagonal. The invariant space has rank ten and gives the entire middle cohomology of the quotient. For $Q=q^r$ its point count is $1+Q^2+Qt_r$. `rvKummerEnriquesZeta[B,C,a];` uses the first ten normalized traces to compute $R$ by exact Newton division and reports the denominator $(1-T)R(qT)(1-q^2T)$.

`rvenriques Y[C,Q,B,a];` declares the finite rational-point type `RVEnriques(Y)`. The named curve `Q` is the coefficient extension of `C` to its quadratic field. The certified marking `B` may live over a different compatible extension. The interpreter descends the normalized translation matrix for the even mask $a$ to the base and then embeds its entries into `Q`. The resulting involution $\sigma_a$ uses the resolved switch and the projective translation on points and tangent directions. [Enriques rational-point example](examples/hyperelliptic-enriques-points.asl)

`rvEnriquesOrbit[Y](X)` accepts a resolved point on the named base `C` or cover `Q`. Over the cover it requires base Frobenius to fix `X` or exchange it with $\sigma_a(X)$. The value stores the smaller of the two resolved points in the canonical coordinate order on the cover. `rvEnriquesSplit(y)` detects the first case in `S`. `rvEnriquesRepresentative(y)` returns the stored cover point, and `rvEnriquesLift(y)` recovers a base point for a split orbit. Enumeration uses the actual resolved cover points within the finite-space bounds and checks its count against the arithmetic formula. The example recovers all 22 points, including six nonsplit orbits, and lifts an exceptional point with its tangent direction.

`rvKummerNodeInvolution(N,X)` is the involution of the double plane obtained by projection from a rational node `N`. Give `N` as an `RVKummer(C)` point and `X` as an `RVKummerK3(C)` point on the same named curve. It exchanges the two residual intersections of a line through `N` with the quartic. At other nodes it applies the differential to the exceptional directions, and over `N` it uses the residual point of the tangent ray or fixes the direction when that point is again `N`. On divisor classes it sends the hyperplane class $H$ to $3H-4E_N$ and the central exceptional class $E_N$ to $2H-3E_N$. Each other exceptional class is fixed. [Node projection example](examples/hyperelliptic-kummer-node-projection.asl)

By Tate's theorem for divisors on abelian varieties and the sixteen exceptional curves of the Kummer resolution, the Picard report counts the $q^r$-eigenvalues on $H^2$ over $\mathbb F_{q^r}$ and gives the geometric rank. A node orbit of length $d$ contributes $\gcd(r,d)$, and an abelian eigenvalue contributes when its ratio to $q$ is a root of unity of order dividing $r$. The example has ranks 10 and 20 over degrees one and three on the first curve. The supersingular example over $\mathbb F_{25}$ gives rank 22 and 1,176 resolved points. [Resolved Kummer example](examples/hyperelliptic-kummer-k3.asl)

`rvKummerK3ArtinTate[C,r];` gives the signed product $\operatorname{disc}\operatorname{NS}(X_{\mathbb F_Q})\,|\operatorname{Br}(X_{\mathbb F_Q})|$ for the resolved Kummer surface $X$ over $Q=q^r$. The displayed square class in $\mathbb Q^*/\mathbb Q^{*2}$ is that of the finite-field Néron–Severi discriminant, represented by the signed product. Over $\mathbb F_Q$ each node orbit of length $d$ contributes $(d/g)^g$ with $g=\gcd(d,r)$. For the Jacobian Weil coefficients $a,b$ in $t^4-at^3+bt^2-aQt+Q^2$, the abelian factor is $a^2-4b+8Q$ when nonzero. When it vanishes, use $4Q-(a/2)^2$ if nonzero and $Q$ otherwise. The sign is $(-1)^{\rho-1}$ for the finite-field Picard rank $\rho$. The calculation uses checked 128-bit arithmetic. [Artin–Tate formula, §4](https://wwwuser.gwdguser.de/~jjahnel/Arbeiten/pic_dst_cam2a.pdf)

When the geometric rank is 22, the report also gives the Artin invariant $\sigma$ and the geometric Néron–Severi discriminant $-p^{2\sigma}$. The invariant is 1 for a superspecial Jacobian, detected by the zero Cartier–Manin matrix in this rank, and 2 otherwise. When the finite-field rank is also 22, the two discriminants agree and the calculation verifies that the resulting Brauer order is a positive integer square. The example gives geometric discriminants $-25$ and $-81$, and exhibits a rank jump and a nontrivial Brauer order. [Supersingular Kummer surfaces, §6](https://arxiv.org/abs/1403.7894) · [Artin–Tate example](examples/hyperelliptic-kummer-artin-tate.asl)

#### Hilbert schemes of points on the resolved Kummer K3

For the resolution $S$ associated to `C`, `rvKummerHilbPoints[C,n,r];` counts points of the Hilbert scheme $S^{[n]}$ over $\mathbb F_{q^r}$. `rvKummerHilbTrace[C,n,r,i];` gives the signed trace of $F^r$ on $H^i(S^{[n]})$. `rvKummerHilbBetti[n,i];` gives the universal Betti number for a K3 Hilbert scheme. The length satisfies $n\ge0$ and the extension degree satisfies $r\ge1$. Length zero gives a point, $H^0$ has trace and Betti number one, and odd degrees and degrees above $4n$ give zero.

The arithmetic Göttsche–Nakajima decomposition uses surface cohomology in degree $b$ for each length $m\ge1$. Its Hilbert-scheme degree is $b+2(m-1)$ and its Frobenius eigenvalues acquire the factor $q^{m-1}$. ASL computes exact Newton recurrences from these data and the resolved Kummer Frobenius traces. Point counts and Betti numbers use the unsigned 128-bit range, traces use the signed 128-bit range, and exact intermediates are checked at 512 bits with the usual work and storage limits. The example gives point counts 327 and 4,480 at lengths two and three over $\mathbb F_3$. [Kummer Hilbert-scheme example](examples/hyperelliptic-kummer-hilbert.asl)

#### K3 elliptic genera and the Igusa cusp form

`rvK3EllipticGenus[n,d,r];` gives the coefficient of $P^n q^d y^r$ in the DMVV product for Hilbert schemes of a complex K3 surface. These coefficients are universal over complex K3 surfaces. The Frobenius reports above use a specified surface and finite field. Here $n,d\ge0$ and $r$ is signed. ASL computes the primitive $c(d,r)$ from $8\sum_{j=2}^4(\theta_j(z)/\theta_j(0))^2$ and expands $\prod_{m\ge1,d\ge0,r}(1-P^m q^d y^r)^{-c(md,r)}$. The length-zero term is one. [Dijkgraaf–Moore–Verlinde–Verlinde](https://arxiv.org/abs/hep-th/9608096)

`rvGenKummerElliptic[N,d,r];` gives $[q^d y^r]$ in the complex elliptic genus of the generalized Kummer variety of dimension $2(N-1)$, for $N\ge1$, $d\ge0$, and signed $r$. Let $A=\phi_{-2,1}=(y+y^{-1}-2)\prod_{j\ge1}(1-q^jy)^2(1-q^jy^{-1})^2/(1-q^j)^4$. The genus is the weight-zero, index-$(N-1)$ Jacobi form $N^4(A|_{-2}V_N)/A$. The interpreter forms the integral Hecke numerator by adding $N(N/a)^3 A_{Nd/a^2}(y^a)$ to row $d$ for each $a\mid\gcd(N,d)$. It subtracts the previously recovered rows convolved with the positive $q$-degrees of $A$, then divides the residual Laurent polynomial exactly by $A_0=(y-1)^2/y$, checking both remainders. At $q=0,y=1$ this gives the Euler number $N^3\sigma_1(N)$. [Nieper-Wißkirchen, Theorem 3](https://arxiv.org/abs/math/0208077) · [Generalized Kummer elliptic genus example](examples/generalized-kummer-elliptic.asl)

`rvIgusaPhi10[m,n,r];` gives the coefficient of $P^m q^n y^r$ in the genus-two weight-10 cusp form. The normalization uses the Weyl factor $Pq(y-2+y^{-1})$ and multiplies it by $(1-P^a q^b y^s)^{c(ab,s)}$ over $a,b\ge0$ with $a+b>0$, including both coordinate axes. Multiplying by the DMVV series cancels the factors with $a>0$, leaving only the $P^1$ Fourier–Jacobi term. Both reports use exact integers with a 512-bit intermediate limit and retain both signs of the Laurent exponents at each truncated $P,q$ degree. [Gritsenko–Nikulin, equation (0.2) squared](https://arxiv.org/abs/alg-geom/9603010) · [Runnable example](examples/k3-elliptic-genus.asl)

`rvK3SymplecticElliptic[h,n,d,r];` gives the $P^nq^dy^r$ coefficient of the equivariant DMVV lift of a complex K3 symplectic automorphism of exact order $1\le h\le8$. Order one gives the ordinary genus. For each higher order ASL computes the primitive elliptic trace by summing the local Lefschetz contributions over the real cyclotomic conjugates of $\zeta_h+\zeta_h^{-1}$. For each oscillator $X=P^m q^j y^r$, it evaluates the primitive traces $T_s$ of powers of the action at $(mj,r)$ and solves $T_s=\sum_{e\mid s}e k_e$ for $s\mid h$. The product contributes $(1-X^e)^{-k_e}$ for each divisor $e\mid h$. The calculation uses exact integers and depends only on the complex action type. [Creutzig–Höhn, Table 1 and §4.2](https://arxiv.org/abs/1309.2671) · [Höhn–Mason, §9.2](https://www.math.sinica.edu.tw/bulletins/20192/2019204.pdf) · [Runnable example](examples/k3-symplectic-elliptic-genus.asl)

At length one the lift gives the geometric elliptic trace. At $y=1$ it gives the Euler traces on the Hilbert schemes at all lengths. For the full two-variable equivariant elliptic trace, Höhn and Mason proved agreement at length two for orders 3, 4, 5, 7, and 8 and checked through $q^4$ for orders 2 and 6. The reports at general length compute the lift. [Höhn–Mason, Theorem 9.3](https://www.math.sinica.edu.tw/bulletins/20192/2019204.pdf) · [Bryan–Gyenge, fixed Hilbert-scheme Euler characteristics](https://arxiv.org/abs/1907.01535)

`rvKummerTranslateElliptic[B,a,n,d,r];` gives the same equivariant DMVV coefficient from a translation by mask `a` in the certified genus-two level-2 basis `B`. Mask zero gives the order-one coefficient. Each nonzero mask gives the order-two coefficient, whose eight fixed points contribute the primitive trace $8(\theta_2(z)/\theta_2(0))^2$ by holomorphic Lefschetz. [Nikulin class elliptic genus, equation (2.6)](https://arxiv.org/abs/1008.4924)

`rvKummerEnriquesElliptic[B,C,a,n,d,r];` gives the coefficient for the free even switch on the covering K3 Hilbert scheme. It uses the same marking, base extension, and Frobenius-stability checks as the Enriques quotient reports. The primitive trace of the free switch vanishes by holomorphic Lefschetz. In the equivariant DMVV product, the virtual $+1$ and $-1$ multiplicities $c_\pm=(c\pm c_g)/2$ contribute $(1-X)^{-c_+}(1+X)^{-c_-}$ at $X=P^m q^d y^r$, with $c$ and $c_g$ evaluated at $(md,r)$. For the free switch this reduces to $(1-X^2)^{-c/2}$, so coefficients vanish unless all three indices are even. [DMVV product](https://arxiv.org/abs/hep-th/9608096) · [Twined elliptic genus example](examples/k3-twined-elliptic-genus.asl)

For the free switch, the geometric identification holds at every length. An invariant subscheme descends along the étale double cover to a subscheme of the Enriques quotient $Y$, so the fixed scheme at length $2m$ is $Y^{[m]}$. Along it the invariant tangent is $T(Y^{[m]})$. The anti-invariant tangent comes from the trace-free ideal-sheaf deformation complex with coefficients in $K_Y$, the anti-invariant summand of the pushforward of the structure sheaf of the cover. Since $K_Y$ is torsion, Grothendieck–Riemann–Roch gives the two bundles the same rational Chern character. The two Lefschetz factors for each Chern root $x$ combine into half the ordinary elliptic factor at $(2x,2\tau,2z)$. Doubling the roots cancels the powers of one half in the top degree. Thus the geometric trace is $\operatorname{Ell}(Y^{[m]};2\tau,2z)$, whose ordinary DMVV series uses $\operatorname{Ell}(Y)=\operatorname{Ell}(K3)/2$.

#### Reciprocal Igusa coefficients in the two chambers at $y=1$

`rvIgusaInverse[m,n,r,inside];` and `rvIgusaInverse[m,n,r,outside];` give the coefficient of $P^m q^n y^r$ in $1/\Phi_{10}$. The signed exponents $m,n$ start at $-1$. ASL reverses the exponents of the Borcherds factors with $a,b\ge0$ and $a+b>0$. It keeps the full Laurent row of this product at $(a,b)=(m+1,n+1)$ and convolves it with the pole factor $y/(1-y)^2$. The `inside` chamber uses $\sum_{j\ge1}j y^j$ and the `outside` chamber uses $\sum_{j\ge1}j y^{-j}$. These are the adjacent expansions at $y=1$ near the $P,q$ cusp. They use exact integers with checked 512-bit intermediates. At $m=n=-1$, the respective coefficients are $\max(r,0)$ and $\max(-r,0)$. [Gritsenko–Nikulin, equation (0.2) squared](https://arxiv.org/abs/alg-geom/9603010) · [Dabholkar–Murthy–Zagier, wall crossing](https://arxiv.org/abs/1208.4074)

`rvIgusaResidue[m,n];` gives the double-pole coefficient $d(m)d(n)$, with $d(u)=[q^u]\eta^{-24}$. Specializing the K3 elliptic genus at $y=1$ gives $24$, so only the two axes of the product contribute to this value. The Laurent row is symmetric under $y\mapsto y^{-1}$. Its first moment therefore vanishes, and `rvIgusaWallJump[m,n,r];` gives the inside coefficient minus the outside coefficient as $r\,d(m)d(n)$. The same product appears in the primitive $K3\times E$ reduced Gromov–Witten formula, which uses the negative reciprocal, the `outside` chamber, and the substitution $y=e^{iu}$. [Oberdieck–Pixton, equation (3) and Theorem 1](https://arxiv.org/abs/1706.10100) · [Runnable example](examples/igusa-meromorphic-chambers.asl)

#### Polar and finite parts of the Igusa Fourier–Jacobi coefficients

For $m>0$ write $\psi_m=[P^m]\Phi_{10}^{-1}$. `rvIgusaPolar[m,n,r,inside];` and `rvIgusaPolar[m,n,r,outside];` report the coefficient of $q^ny^r$ in its polar part

$$
\psi_m^{\mathrm P}=\frac{d(m)}{\Delta(q)}\sum_{s\in\mathbb Z}q^{ms^2}y^{2ms}\frac{q^sy}{(1-q^sy)^2}.
$$

The $s=0$ term uses the named chamber. The other terms use positive powers of $q^sy$ for $s>0$ and negative powers for $s<0$, as appropriate in the two strips adjacent to $y=1$. This includes every elliptic translate of the pole that can contribute to the requested coefficient. [Dabholkar–Murthy–Zagier, Theorem 8.3 and equation (403)](https://arxiv.org/abs/1208.4074)

`rvIgusaFinite[m,n,r];` gives the chamber-independent coefficient of the holomorphic-in-$z$ part $\psi_m^{\mathrm F}=\psi_m-\psi_m^{\mathrm P}$. It uses the elliptic law $(n,r)\mapsto(n+\lambda r+m\lambda^2,r+2m\lambda)$ and evenness in $r$ to reduce to $0\le r\le m$. The polar coefficient vanishes there in the `outside` expansion. `rvIgusaMockWeight2[m,n,r];` gives the corresponding coefficient of $\Delta\psi_m^{\mathrm F}$, of mock Jacobi weight two. The commands compute these Fourier coefficients with checked 512-bit integer arithmetic. In the example, the full coefficients 315255 and 327 reduce to the same finite coefficient 3 after the polar subtraction. [Dabholkar–Murthy–Zagier, equations (128), (249), and (404)](https://arxiv.org/abs/1208.4074) · [Runnable example](examples/igusa-mock-jacobi.asl)

#### Hurwitz class numbers and level-one Hecke traces

`rvHurwitzClass[N];` gives the exact rational Hurwitz number $H(N)$. For $N>0$ it sums $2h(-N/d^2)/w(-N/d^2)$ over square divisors $d^2\mid N$ for which $-N/d^2$ is a quadratic-order discriminant. Here $h$ counts proper primitive classes and $w$ is the order of the unit group, equal to 6 at $-3$, 4 at $-4$, and 2 otherwise. The square divisors include the imprimitive forms. The convention is $H(0)=-1/12$, with zero for negative $N$ or $N\equiv1,2\pmod4$. `rvHurwitzJacobi[n,r];` evaluates $H(4n-r^2)$, the $q^ny^r$ coefficient of Zagier's index-one mock Jacobi series $\mathcal H$. The discriminant is computed exactly before enumerating the classes. This series enters $\Delta\psi_1^{\mathrm F}=3E_4\phi_{-2,1}-648\mathcal H$. [Zagier, Utrecht lectures, §6.1](https://people.mpim-bonn.mpg.de/zagier/files/tex/UtrechtLectures/UtBook.pdf) · [Dabholkar–Murthy–Zagier, equations (128), (168), and (249)](https://arxiv.org/abs/1208.4074)

`rvHeckeTrace[k,n];` gives the trace of the classical $T_n$ on $S_k(\operatorname{SL}_2(\mathbb Z))$ for even $k\ge2$ and $n\ge1$. The normalization sends the $j$th Fourier coefficient $a(j)$ to $\sum_{d\mid\gcd(n,j)}d^{k-1}a(nj/d^2)$. ASL evaluates the Eichler–Selberg formula

$$
-\frac12\sum_{t^2\le4n}P_{k-2}(t,n)H(4n-t^2)-\frac12\sum_{dd'=n}\min(d,d')^{k-1}+\mathbf1_{k=2}\sigma_1(n),
$$

where $P_0=1$, $P_1=t$, and $P_j=tP_{j-1}-nP_{j-2}$. It uses binary doubling for the recurrence and checked 512-bit integers for the trace sum. The example gives $\operatorname{Tr}(T_5\mid S_{12})=4830$ and $\operatorname{Tr}(T_2\mid S_{24})=1080$. [Zagier, Utrecht lectures, §6.3](https://people.mpim-bonn.mpg.de/zagier/files/tex/UtrechtLectures/UtBook.pdf) · [Runnable example](examples/hurwitz-eichler-selberg.asl)

#### Formal level-one Hecke operators

`RVHecke(k)` is the integral Hecke algebra in even weight $k\ge2$. Construct $T_n$ with `rvHeckeT[k,n]` for $n\ge1$, or its unit $T_1$ with `rvHeckeUnit[k]`. Use `zero[RVHecke(k)]`, `add`, `neg`, and `rvHeckeScale[c](A)` for finite integer combinations. The product `rvHeckeMul(A,B)` extends $T_m T_n=\sum_{d\mid\gcd(m,n)}d^{k-1}T_{mn/d^2}$ bilinearly. These operations preserve the weight, and `==` compares the formal coefficients. [Zagier, Utrecht lectures, §2.1](https://people.mpim-bonn.mpg.de/zagier/files/tex/UtrechtLectures/UtBook.pdf)

`rvHeckeActionTrace(A);` reports the trace of the action on $S_k(\operatorname{SL}_2(\mathbb Z))$. `rvHeckeChar(A);` gives the characteristic polynomial as descending coefficients $[1,c_1,\ldots,c_d]$, where $d=\dim S_k$. `rvHeckeCharCoeff[j](A);` gives $c_j$, with $c_0=1$ and zero for $j>d$. The calculation combines Eichler–Selberg traces with exact Newton identities, and uses $\det(X-cT_1)=(X-c)^d$ for scalar operators. The commands use checked 512-bit integer intermediates and bounded computation.

`rvHeckeSameAction(A,B)` is true when $\operatorname{Tr}((A-B)^2)=0$. Each $T_n$ is self-adjoint for the Petersson inner product, so this trace is the sum of squares of the real eigenvalues of $A-B$. In weight 24, $T_2^2-1080T_2-20468736T_1$ has zero action on the two-dimensional cusp space and a nonzero formal $T_4$ coefficient. [Runnable example](examples/graded-hecke-operators.asl)

#### Simultaneous Hecke eigenpackets

`rvheckepacket E(A);` names the rational algebra $\mathbb Q[X]/(f)$ where $f$ is the characteristic polynomial of $A$ acting on the nonzero cusp space $S_k$. It requires $f$ to be squarefree. The certificate uses the trace matrix $H_{ij}=\operatorname{Tr}(A^{i+j})$ for $0\le i,j<\dim S_k$, which is the Gram matrix of the Vandermonde matrix of the real eigenvalues. Invertibility certifies that $A$ separates the simultaneous eigenforms. The resulting algebra is étale over $\mathbb Q$ and may have more than one field factor. `rvHeckeEigenData[E];` reports the weight and the descending coefficients of $f$. [Zagier, Utrecht lectures, §2.1](https://people.mpim-bonn.mpg.de/zagier/files/tex/UtrechtLectures/UtBook.pdf)

`rvHeckeEigenImage[E](B)` identifies the action of any same-weight operator $B$ with the unique polynomial $G_B(X)$ of degree less than $\deg f$. Since $B$ commutes with $A$, it preserves each eigenspace, and the coordinates solve $Hc=(\operatorname{Tr}(A^iB))_i$ by exact fraction-free elimination. `rvHeckeEigenCoeff[E,n]` is the image of $T_n$. At each embedding it gives the $n$th Fourier coefficient of the corresponding normalized eigenform. The calculation uses checked 512-bit stored integers and bounded exact arithmetic.

Elements have type `RVHeckeEigen(E)` and support `zero`, `add`, `neg`, `rvHeckeEigenMul(P,Q)`, and rational constants `rvHeckeEigenScalar[E,num,den]` with positive denominator. `rvHeckeEigenGenerator[E]` is the class of $X$. `show` prints the ascending numerator coordinates with a common denominator. `rvHeckeEigenCoord[j](P);` gives the reduced rational coefficient of $X^j$. `rvHeckeEigenTrace(P);` and `rvHeckeEigenNorm(P);` give the trace and determinant of multiplication by $P$ in the algebra. In weight 24, the $T_2$ generator gives $a_3=195660-48X$. With generator $5T_2+7T_1$ the same coefficient becomes $(978636-48X)/5$. [Runnable example](examples/hecke-eigenpackets.asl)

#### Integral period chambers and Farey wall crossing

`rvsiegelchamber H[A,B,C];` names the open chamber containing the positive definite imaginary-period ray with matrix $Y=\left(\begin{smallmatrix}A&B\\B&C\end{smallmatrix}\right)$. The coordinates are signed integers. An integral change of basis that makes $Y$ diagonal places the ray on a separating wall, so the declaration requires an open chamber. The reference chamber has $0<B<\min(A,C)$. Projectivizing the cone identifies it with the upper half-plane by $(B+i\sqrt{AC-B^2})/A$, where the reference chamber is the ideal Farey triangle with vertices $0,1,\infty$. [Cheng–Verlinde, §2 and Appendix A](https://arxiv.org/abs/0806.2337)

`rvIgusaChamberCoeff[H,m,n,r];` computes the coefficient of $1/\Phi_{10}$ in the named chamber at the Fourier index $T=\left(\begin{smallmatrix}m&r/2\\r/2&n\end{smallmatrix}\right)$. Integer Gauss reduction chooses a basis $U$ for which $U^tYU$ has $0<2B\le A\le C$ in its reduced coordinates. `rvSiegelChamberData[H];` reports those coordinates and the columns of $U$ as a matrix. The coefficient uses the `inside` expansion at $U^{-1}TU^{-t}$. Positive rescaling of $Y$ leaves the coefficient unchanged, and $\det U=\pm1$ gives the same multiplier because the form has even weight. These are meromorphic Fourier coefficients for open cusp contours specified by real imaginary-period rays. [Cheng–Verlinde, §§2 and 5](https://arxiv.org/abs/0806.2337)

`rvIgusaChamberPath[H0,H1,m,n,r];` reports how many separating walls lie between the two open chambers, and the coefficient in `H0` minus the coefficient in `H1`. Each crossing changes to an integral frame in which its wall is $B=0$ and the current chamber is on the positive side. It adds the local jump $r' d(m')d(n')$ using the transformed Fourier indices and the eta coefficients above. The sum is calculated from the residues. A path between points of the same chamber has zero walls. Explicit paths have a 256-wall limit. Direct coefficient lookup reduces the period using whole integer quotients. [Cheng–Verlinde, §5.3](https://arxiv.org/abs/0806.2337) · [Runnable example](examples/igusa-farey-chambers.asl)

#### Rank-three Borcherds denominator and Weyl weights

`rvSiegelWeight[A,B,R]` constructs an `RVSiegelWeight` with monomial exponents $(a,b,r)=(A,B,R)/2$. The doubled coordinates are signed integers and support `zero[RVSiegelWeight]`, `add`, and `neg`. The report `rvSiegelPair(w,v);` gives $B(w,v)=2rr'-4(ab'+ba')$ as an exact integer or half-integer. `rvSiegelReflect[i](w)` acts by $s_i(w)=w-B(w,e_i)e_i$, with indices $0,1,2$ in the order $(0,0,-1),(1,0,1),(0,1,1)$. These roots have square $2$ and pairwise pairing $-2$. The Weyl vector $\rho=(1/2,1/2,1/2)$ pairs to $-1$ with each. [Gritsenko–Nikulin, the rank-three matrix $A_{1,\mathrm{II}}$](https://arxiv.org/abs/alg-geom/9603010)

`rvSiegelRoot(w);` requires integral $(a,b,r)$ and gives the root class, its orientation, its coordinates $(a+b-r,a,b)$ in the simple basis, its norm, and its signed supermultiplicity. Positive and negative roots use the same multiplicity. An out-of-support weight or the zero vector gives zero multiplicity. For a root the multiplicity is the coefficient of $q^{ab}y^r$ in the weak Jacobi form $\phi_{0,1}=\operatorname{Ell}(K3)/2$. ASL uses the discriminant $4ab-r^2$ and the parity of $r$ to reduce the Jacobi coefficient before computing it. These are supermultiplicities of the root spaces in the denominator product. [Gritsenko–Nikulin, equation (0.2)](https://arxiv.org/abs/alg-geom/9603010)

`rvSiegelDelta5[A,B,R];` and `rvSiegelDelta5At(w);` give the coefficient of $P^{A/2}q^{B/2}y^{R/2}$ in the weight-five denominator. Its leading factor is $(Pq)^{1/2}(y^{1/2}-y^{-1/2})$. This factor absorbs the positive root $(0,0,-1)$ of multiplicity one. The remaining factors are $(1-P^a q^b y^r)^{[q^{ab}y^r]\phi_{0,1}}$ for $a,b\ge0$ with $a+b>0$, including both axes. The normalized denominator satisfies $\Delta_5^2=\Phi_{10}$ and changes sign under each simple reflection. [Gritsenko–Nikulin, equation (0.2)](https://arxiv.org/abs/alg-geom/9603010) · [Runnable example](examples/siegel-borcherds-denominator.asl)

#### Jacobi cusp forms and Maaß additive lifts

Declare an index-one Jacobi cusp form with `rvsiegelmaass F[k]{minus(a,b,d,c),zero(a,b,d,c),...};`. A term is $c E_4^a E_6^b\Delta^d\phi_{j,1}$, with $j=-2$ for `minus` and $j=0$ for `zero`. Each term must have the declared weight $k=4a+6b+12d+j$. The exponents $a,b$ are nonnegative, $d\ge1$, and the scalar $c$ is signed. The normalized Eisenstein series have constant one, $\Delta=q\prod_{i\ge1}(1-q^i)^{24}$, and the leading Laurent polynomials of $\phi_{-2,1}$ and $\phi_{0,1}$ are $y^{-1}-2+y$ and $y^{-1}+10+y$. `rvSiegelMaassJacobi[F,n,r];` gives the seed coefficient $c_F(n,r)$ of $q^n y^r$. [Eichler–Zagier, *The Theory of Jacobi Forms*](https://people.mpim-bonn.mpg.de/zagier/files-restricted/scanned/TheoryOfJacobiForms/fulltext.pdf)

`rvSiegelMaassCoeff[F,m,n,r];` gives the coefficient of $P^m q^n y^r$ in the degree-two Saito–Kurokawa lift of the declared seed. For $D=4mn-r^2>0$, it evaluates $\sum_{e\mid\gcd(m,n,r)}e^{k-1}c_F(mn/e^2,r/e)$. A nonpositive discriminant gives zero. The seed coefficients use the index-one elliptic law to reduce to Laurent exponent $0$ or $1$ before expansion. The discriminant is computed exactly from the original coordinates, so a binary quadratic form with large coefficients can still require only a short expansion. [Maaß, *Über eine Spezialschar von Modulformen zweiten Grades*](https://doi.org/10.1007/BF01389857)

`rvIgusaChi12[m,n,r];` uses the seed $\Delta\phi_{0,1}$. Its first term is $Pq(y^{-1}+10+y)$, and its restriction at $y=1$ is $12\Delta(P)\Delta(q)$. The seed $\Delta\phi_{-2,1}$ gives the same normalized $\Phi_{10}$ as the Borcherds product above. [Igusa, *On Siegel Modular Forms of Genus Two*](https://doi.org/10.2307/2372812) · [Runnable example](examples/siegel-maass-lifts.asl)

#### Hilbert schemes of Jacobians and generalized Kummer fibers

For a genus-two Jacobian $A=J(C)$, `rvJacHilbPoints[C,n,r];` counts $\operatorname{Hilb}^n(A)$ over $\mathbb F_{q^r}$, including the single length-zero point. The generating function over $\mathbb F_Q$ is $\prod_{j\ge1}Z(A,Q^{j-1}t^j)$. The interpreter evaluates $\#A(\mathbb F_{Q^s})=\det(1-F^s)$ from the base Weil polynomial and uses the logarithmic derivative recurrence to recover the requested coefficient, with checked 128-bit arithmetic and exact division.

`rvGenKummerPoints[C,n];` counts the zero fiber of the summation map over the named field. `rvGenKummerFiber[C,n](D);` counts the fiber over the supplied `RVJac(C)` class, with $n\ge1$. For $p\nmid n$ these are smooth generalized Kummer fibers. For $p\mid n$ the report gives the summation fiber count. Put $G=A(\mathbb F_q)$ and $H_m=\#\operatorname{Hilb}^m(A)(\mathbb F_q)$. The fiber count is $|G|^{-1}\sum_{d\mid n}w_d(D)q^{n-n/d}H_{n/d}$, where $w_d(D)$ sums the characters of order $d$ evaluated at $-D$. The trace $A(\mathbb F_{q^s})\to G$ is onto: its norm operator is an isogeny and $(F-1)(1+F+\cdots+F^{s-1})=F^s-1$. In the character-weighted Euler product, trace surjectivity leaves precisely the factors with punctual index divisible by $d$. Möbius inversion recovers $w_d(D)$ from the sums for orders dividing $e$, which equal $|G[e]|$ when $D\in eG$ and zero otherwise. This uses rational branch orbits and Kummer descent for $e=2$. The other required kernels and images use a single bounded pass over reduced Mumford pairs.

`rvGenKummerLevel[B,C,r];` computes the prime-length zero fiber over $\mathbb F_{q^r}$ using the certified full torsion basis `B` on a base extension of `C`. The nullity of $F^r-I$ on the basis determines the rational torsion order $\ell^k$ and the contribution $(\ell^k-1)Q^{\ell-1}$ to $H_\ell/|A(\mathbb F_Q)|$. The example compares length two with the resolved Kummer K3 over two fields and checks nonuniform fibers at lengths three and four. [Generalized Kummer example](examples/hyperelliptic-generalized-kummer.asl)

`rvGenKummerBetti[n,i];` gives the geometric Betti number in degree $i$ of the smooth generalized Kummer variety of length $n\ge1$, in characteristic zero or characteristic prime to $n$. Its dimension is $2(n-1)$. `rvGenKummerTrace[B,C,r,i];` gives the actual degree-$i$ trace for prime length $\ell$ from a certified full torsion basis on a base extension of `C`. `rvJacHilbBetti[n,i];` and `rvJacHilbTrace[C,n,r,i];` give the corresponding invariants of the full abelian-surface Hilbert scheme. Degrees outside the geometric range give zero. The signed calculations retain odd degrees and remove the alternating sign when reporting a trace or Betti number.

The graded recurrence uses the Göttsche–Nakajima cohomology decomposition. A copy of $H^b(A)$ at length $m$ lies in degree $b+2(m-1)$ with the Tate factor $Q^{m-1}$. Odd-degree copies contribute exterior powers, so the signed logarithmic derivative has sign $(-1)^b$ at every oscillator power. With $G$ denoting the signed, graded Frobenius polynomial, the generalized Kummer polynomial satisfies $G_A(t)G_K(t)=\sum_{d\mid n} E_d(F^r)(Qt^2)^{n-n/d}G_{\operatorname{Hilb}^{n/d}(A)}(t)$. Here $E_d$ is the trace on the finite set of exact-order-$d$ torsion points. The identity follows from the point counts over all extensions and separation by cohomological weights. For Betti numbers, $Q=1$ and the torsion contribution is the Jordan totient $J_4(d)=\sum_{e\mid d}\mu(d/e)e^4$. The calculation uses checked 512-bit intermediates and Poincaré duality, with unsigned 128-bit Betti numbers and signed 128-bit traces. The example gives $b_3=8$ and a degree-three trace of $-1000$ for a generalized Kummer fourfold over $\mathbb F_{25}$. [Graded generalized Kummer example](examples/hyperelliptic-generalized-kummer-cohomology.asl)

#### Torsion translations on generalized Kummer varieties

`rvGenKummerTranslateChar[n,h,i];` gives the trace on complex $H^i(K_{n-1})$ of translation by a point of exact order $h$, for $h\mid n$ and $h\ge1$. In the Göttsche–Soergel decomposition the sum over characters of order dividing $e$ is $e^4$ if $eh\mid n$, and zero otherwise. Möbius inversion gives the exact-order contribution in the graded formula. Order one gives the Betti numbers. `rvGenKummerTranslateFixed[n,h];` reports the fixed locus, with $h^3$ components of deformation type $K_{n/h-1}$ and total Euler characteristic $n^3\sigma_1(n/h)$. A translation-invariant subscheme descends to a subscheme of length $n/h$ on the abelian surface $A/\langle t\rangle$. The induced isogeny $A/\langle t\rangle\to A$ factors $[h]$ and has degree $h^3$. The sum of the lifted subscheme also has the constant shift $(n/h)\sum_{j=0}^{h-1}jt$, so the summation condition selects fibers over a translate of that kernel. [Göttsche–Soergel](https://doi.org/10.1007/BF01445104)

`rvGenKummerTranslateFrob[B,C,r,i](t);` gives the trace on degree $i$ of translation by $t$ composed after the $r$-th Frobenius of `C`. Here $r\ge1$ and `B` is a certified genus-two prime-level $\ell$ marking on a compatible extension of `C`. The input $t$ is an actual $\ell$-torsion Jacobian class on the curve of `B`. `rvGenKummerTranslatePoints[B,C,r](t);` counts the geometric fixed points of the same composition. Frobenius acts on torsion with matrix $M_r$, and the affine fixed-point count is the number of solutions of $(I-M_r)x=v$, where $v$ is the coordinate vector of $t$. Linear elimination gives either zero solutions or the size of the kernel. This count replaces the rational torsion order in the middle-degree correction of weight $Q^{\ell-1}$ for $Q=q^r$. Zero translation recovers the ordinary Frobenius reports. Traces use signed 128-bit results with checked 512-bit intermediates. Counts and fixed-locus invariants use checked unsigned 128-bit arithmetic and the usual work limits. [Torsion translation example](examples/hyperelliptic-generalized-kummer-translation.asl)

#### Hodge numbers and torsion spectra

`rvK3HilbHodge[n,p,q];` and `rvJacHilbHodge[n,p,q];` give $h^{p,q}$ over $\mathbb C$ for length-$n$ Hilbert schemes of a K3 surface and an abelian surface, respectively. Length zero gives a point. `rvGenKummerHodge[N,p,q];` gives the Hodge number of the length-$N$ generalized Kummer fiber $K_{N-1}$, of dimension $2(N-1)$ for $N\ge1$. Bidegrees outside the dimension give zero. A length-$j$ generator from $H^{a,b}$ of the surface contributes in bidegree $(a+j-1,b+j-1)$. Generators of even total degree contribute symmetric powers and those of odd total degree contribute exterior powers. [Göttsche–Soergel, §5](https://doi.org/10.1007/BF01445104)

`rvGenKummerTranslateHodge[N,h,p,q];` gives the trace of an exact-order $h$ translation on $H^{p,q}(K_{N-1})$, with $h\ge1$ dividing $N$. In the cohomology decomposition a partition with parts of greatest common divisor $g$ has copies indexed by $A[g]$. The translation permutes these copies by adding $(N/g)t$. Grouping partitions with all parts divisible by $d$ gives the character sum $C_d(t)$ from the graded trace above and shifts both Hodge indices by $N-N/d$. Division by $E(A)=(1-x)^2(1-y)^2$ gives the signed polynomial of the fiber. The report removes the sign $(-1)^{p+q}$. [Göttsche–Soergel, Theorem 7](https://doi.org/10.1007/BF01445104) · [Morrison–Shen, relative power structures](https://www.mathjournals.org/mmj/2021-021-004/2021-021-004-006.html)

`rvGenKummerHodgeSpectrum[N,h,p,q];` returns the pairs $(d,m_d)$ with $m_d>0$. The characteristic polynomial of translation on that Hodge piece is $\prod_{d\mid h}\Phi_d^{m_d}$. Inverting the traces of all orders of powers of the translation by Ramanujan sums gives the multiplicity $m_d$ of each primitive $d$-th root. The interpreter checks integrality, nonnegativity, and $\sum_{d\mid h}\varphi(d)m_d=h^{p,q}$. The calculation uses checked 512-bit intermediates and the usual work and storage limits. Hodge numbers and multiplicities use unsigned 128-bit results, and traces use signed 128-bit results. [Hodge and translation spectrum example](examples/hyperelliptic-generalized-kummer-hodge.asl)

#### Affine involutions on generalized Kummer varieties

`rvGenKummerInversionHodge[N,h,p,q];` gives the complex Hodge trace of the natural action $x\mapsto t-x$ on $K_{N-1}$, for a translation parameter $t$ of exact order $h\mid N$. The parameter of order one is zero. Inversion interchanges each torsion character with its inverse, so only the orders $1$ and $2$ contribute to the Hodge trace. For even $N$, the order-two term has character sum $C_2(t)=16[2h\mid N]-1$ and shifts both indices by $N/2$. Ordinary abelian Hilbert cohomology contributes the inversion sign $(-1)^{p+q}$. [Boissière–Nieper-Wißkirchen–Sarti](https://arxiv.org/abs/1001.4728) · [Göttsche–Soergel, Theorem 7](https://doi.org/10.1007/BF01445104)

`rvGenKummerInversionEuler[N,h];` gives the Euler characteristic of the fixed locus. Put $H_m=[z^m]\prod_{j\ge1}((1+z^j)/(1-z^j))^8$. The Lefschetz sum is $H_N/16$ plus $C_2(t)H_{N/2}/16$ for even $N$. At length two the Euler characteristic is 24 for zero translation and 8 for nonzero translation. At length three it is 60 for every $t$, from a Kummer K3 component and 36 isolated points. [Boissière–Nieper-Wißkirchen–Tari, generalized Kummer fourfolds](https://arxiv.org/abs/1902.01685)

`rvGenKummerInversionFrob[B,C,r,i](t);` and `rvGenKummerInversionPoints[B,C,r](t);` give the trace in degree $i$ and the geometric fixed-point count for $((t-\mathrm{id})\circ F_q)^r$, with $r\ge1$. The level $\ell$ of the certified genus-two marking `B` is the prime length. The class $t$ belongs to the marked curve's $\ell$-torsion. The homogeneous torsion matrix has linear block $-M$ for Frobenius matrix $M$ and translation column $v$ for $t$. Raising that matrix to the $r$-th power evolves the translation along with the linear action. The resulting fixed torsion count determines the middle-degree correction. The other degrees have the ordinary Frobenius trace multiplied by $(-1)^{ri}$. The fixed count uses the abelian Hilbert Lefschetz sum with Frobenius and inversion, with exact division and the same middle correction. Counts use unsigned 128-bit results, traces use signed 128-bit results, and the Hilbert recurrence uses checked 512-bit intermediates. [Affine inversion example](examples/hyperelliptic-generalized-kummer-inversion.asl)

#### Hilbert–Chow partition summands

`rvhilbstratum S{(m,c),...};` declares a multiset with $c$ parts of positive size $m$. Repeated sizes are combined and zero multiplicities are omitted. Write $n=\sum mc$, $L=\sum c$, and $e=\gcd(m)$ over the occupied sizes. The empty partition has $n=L=e=0$. `rvHilbStratumData[S];` reports the fiber dimension $n-L$, the stratum codimension $2(n-L)$, and the signed local intersection multiplier $(-1)^{n-L}\prod m^c$. The report retains a factorization if the magnitude exceeds 128 bits. [de Cataldo–Migliorini, Lemma 5.1.2](https://arxiv.org/abs/math/0005249)

`rvK3StratumHodge[S,p,q];` and `rvJacStratumHodge[S,p,q];` give the Hodge number of this partition summand for the two surfaces. Each size contributes the signed Hodge polynomial of the symmetric product $\operatorname{Sym}^c$ of the surface. The product shifts both indices by $n-L$. The symmetric products use symmetric powers of even cohomology and exterior powers of odd cohomology. `rvGenKummerStratumHodge[S,h,p,q];` gives the corresponding trace on the generalized Kummer fiber for a translation of exact order $h\mid n$, with $n\ge1$. Its summands have labels in $A[e]$, shifted by $(n/e)t$ under translation by $t$. The trace is zero if $h\nmid n/e$. Otherwise it is $e^4$ times the shifted abelian product divided by $E(A)=(1-x)^2(1-y)^2$. [Xu, generalized Kummer motives](https://arxiv.org/abs/1506.04297) · [Göttsche–Soergel, Theorem 7](https://doi.org/10.1007/BF01445104)

`rvNakajimaStratum[S](x)` keeps the monomials of the integral Nakajima state whose mode-length multiplicities give the partition `S`, summing multiplicities across colors within a length. It preserves the input type and returns the typed zero on other weights or formal zero grades. The projectors are idempotent and orthogonal for the Wick pairing. The Hodge calculations use checked 512-bit intermediates, unsigned 128-bit results, and the usual work and storage limits.

`rvGenKummerStratumFrob[B,C,S,r,i](t);` gives the degree-$i$ trace on this partition summand for translation by $t$ after the $r$-th Frobenius of `C`. The partition weight equals the prime level $\ell$ of the compatible marking `B`, and $t$ is an $\ell$-torsion class on the marked curve. The full diagonal is indexed by $A[\ell]$ and lies in degree $2(\ell-1)$. Its trace is the affine torsion fixed count times $q^{r(\ell-1)}$. Every other prime-weight partition has connected weighted summation fiber and uses the symmetric-power Frobenius traces and division by the signed graded trace polynomial of $A$. The report uses signed 128-bit traces. [Hilbert–Chow partition example](examples/hyperelliptic-hilbert-chow-strata.asl)

`rvK3StratumFrobTrace[W,C,S,r,i];` traces $\Psi^r$ on degree $i$ of a resolved Kummer K3 partition summand, where $\Psi=W\circ F_q$ and $r\ge1$. For each size $m$ with multiplicity $c$, it uses the graded symmetric power $\operatorname{Sym}^c H^*(X)$ of the resolved surface $X$. The product shifts degree by $2(n-L)$ and acquires the Tate factor $q^{r(n-L)}$. The surface power sums come from powers of the single composite $\Psi$, including when $W$ moves under Frobenius. The sum over partitions is `rvK3FrobHilbTrace[W,C,n,r,i];`. [de Cataldo–Migliorini](https://arxiv.org/abs/math/0005249) · [Nakajima](https://arxiv.org/abs/alg-geom/9507012)

`rvK3StratumEulerFactor[W,C,S,r,i,j];` returns the coefficient of $T^j$ in $\det(1-T\Psi^r)$ on that summand. It uses Newton identities and the traces at powers $r,2r,\ldots,jr$, with exact division. The constant coefficient is one, and indices above the cohomological rank give zero. These determinants give formal Lefschetz factors for general words, including words of infinite order. With the identity word they are finite-field Frobenius factors. Both commands use signed 128-bit results with checked 512-bit intermediates and the usual work and storage limits. [Equivariant Hilbert–Chow example](examples/hyperelliptic-hilbert-chow-factors.asl)

#### Specht local systems on Hilbert–Chow strata

`rvspechtstratum D[S]{(k,[mu1,mu2,...]),...};` decorates each occupied size $k$ of a declared Hilbert–Chow partition `S` by a Young diagram $\mu_k$ whose positive, weakly decreasing rows sum to the multiplicity $a_k$ of $k$. Omitted sizes use the one-row diagram and hence the trivial representation. Over the open weighted stratum, the cover that orders equal-weight points has group $G=\prod_k\mathfrak S_{a_k}$. The finite small map from the ordered surface product extends that cover across the closure. Taking the multiplicity space of $\bigotimes_k\operatorname{Specht}_{\mu_k}$ gives the intersection cohomology with its local system. Its graded realization is $\bigotimes_k\mathbb S_{\mu_k}H^*(X)$ for the resolved Kummer K3 surface $X$. The ordinary Hilbert–Chow decomposition uses only the trivial local systems. [de Cataldo–Migliorini, small maps and intersection cohomology](https://arxiv.org/abs/0712.0349) · [Hilbert–Chow motives](https://arxiv.org/abs/math/0005249)

`rvK3SpechtFrobTrace[W,C,D,r,i];` traces $(W\circ F_q)^r$ in intrinsic intersection-cohomology degree $i$, with $r\ge1$. `rvK3SpechtEulerFactor[W,C,D,r,i,j];` gives its degree-$j$ coefficient in $\det(1-T(W\circ F_q)^r)$. Intrinsic grading counts ordinary degrees on the ordered product. For a trivial decoration, passing to the Hilbert-scheme grading in `rvK3StratumFrobTrace` adds $2(n-L)$ to the degree and multiplies the trace by $q^{r(n-L)}$. The Schur traces use Jacobi–Trudi on graded symmetric powers with actual Adams powers of the composite. A diagram with more than 24 rows gives zero cohomology. Factor coefficients above the cohomological rank vanish and the constant coefficient is one. Both queries use signed 128-bit results and checked 512-bit intermediates. [Macdonald, Chapter I, §3](https://math.berkeley.edu/~corteel/MATH249/macdonald.pdf) · [Specht local system example](examples/hyperelliptic-specht-strata.asl)

#### Rational Witt spectra

`RVWittSpectrum` represents rational determinant series $P(T)/Q(T)$ with integral coefficients and $P(0)=Q(0)=1$. Use `rvWittSpectrum{[1,p1,...]}` for an effective polynomial and `rvWittSpectrum{[1,p1,...],[1,q1,...]}` for a virtual difference. Trailing zeros and common factors are removed exactly by primitive polynomial gcd and monic division on the reciprocal polynomials. The factors encode virtual multisets of nonzero algebraic-integer eigenvalues, or spectra of endomorphisms invertible over $\mathbb Q$. The spectrum records eigenvalue multiplicities and forgets Jordan data and the underlying motive. [Almkvist’s theorem on rational Witt vectors](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/CLMPZ.pdf)

`rvWittSpectrumLine[a]` is $1-aT$ for a nonzero integer $a$. The line of value one is the tensor unit, and `zero[RVWittSpectrum]` is the additive zero with series $1$. `add` multiplies series, `neg` swaps the two factors, and `rvWittSpectrumTensor(x,y)` multiplies the ghosts. The convention is $w_k=-[T^k]T\frac{\mathrm d}{\mathrm dT}\log(P/Q)$ for $k\ge1$. `rvWittSpectrumAdams[r](x)` sends $w_k$ to $w_{rk}$ for $r\ge1$. `rvWittSpectrumExterior[m](x)` takes the $m$-th exterior power of an effective spectrum, with the tensor unit at $m=0$. The operations compute with power sums, exact Newton identities, and reduction modulo monic companion polynomials. `rvWittSpectrumData(x);` reports both factor degrees and their difference, the virtual rank. `rvWittSpectrumGhost[k](x);` reports one power sum and `rvWittSpectrumCoeff[j](x);` gives the degree-$j$ coefficients of the two reduced factors. [Rational Witt vectors and characteristic polynomials](https://doi.org/10.1112/S0010437X22007254)

`rvK3SpechtSpectrum[W,C,D,r,i]` produces the full effective determinant in intrinsic degree $i$ of a Specht decoration. `rvK3StratumSpectrum[W,C,S,r,i]` gives the corresponding Hilbert-graded partition spectrum. They reuse the actual powers of the semilinear composite and compute the determinant up to the cohomological rank, subject to the same arithmetic and work bounds as the ring operations. Trimmed literal factors, reduced factors, and reconstructed effective factors each have degree at most 256. Stored coefficients and ghost reports use signed 128-bit integers, and intermediates use checked 512-bit integers. [Witt spectrum example](examples/hyperelliptic-witt-spectrum.asl)

`rvWittSpectrumVerschiebung[n](x)` substitutes $T^n$ for $T$ in both factors for $n≥1$. Its $k$-th ghost is $n w_{k/n}$ when $n$ divides $k$, and zero otherwise. This is the arithmetic Verschiebung on rational Witt series. `rvWittSpectrumNecklace[n](x);` computes the signed integer $b_n=\frac1n\sum_{d\mid n}\mu(d)w_{n/d}$. These are the exponents in the formal factorization $P(T)/Q(T)=\prod_{n\ge1}(1-T^n)^{b_n}$. The calculation uses exact division and bounded trial factorization of the period. [Dold sequences and periodic points](https://arxiv.org/abs/2007.04031)

`rvK3FrobHilbOrbits[W,C,m,n];` counts geometric cycles of least period $n$ for $\Psi=W^{[m]}\circ F_q$ on the length-$m$ Hilbert scheme. It computes the fixed-point counts at the divisor powers of the actual composite and applies Möbius inversion before narrowing to an unsigned 128-bit result. The powers use the semilinear marked action, including when Frobenius moves the generators of $W$. With the identity word, these cycles are the degree-$n$ closed points. Length zero is a single fixed point. `rvK3FrobHilbSpectrum[W,C,m,r]` produces the effective Lefschetz determinant of $\Psi^r$ on the full even cohomology. The Hilbert scheme has only even cohomology, so its series is the inverse of the dynamical zeta function for $\Psi^r$. The bridge uses the same degree and integer bounds as other typed spectra. Under the degree limit of 256, the full determinant fits in length zero or one, and orbit queries can also use larger lengths. [Witt and Hilbert-scheme orbit example](examples/hyperelliptic-witt-orbits.asl)

#### Crystalline Newton polygons

`rvKummerHilbNewton[C,n,i];`, `rvJacHilbNewton[C,n,i];`, and `rvGenKummerNewton[C,n,i];` list the normalized crystalline slopes with multiplicity in degree $i$ of the resolved Kummer K3 Hilbert scheme, the Jacobian Hilbert scheme, and the generalized Kummer variety. The first two accept length zero. The last requires $n\ge1$ prime to the residue characteristic so the summation fiber is smooth. Degrees outside the cohomological range give an empty polygon. `rvKummerK3Height[C];` reports the formal Brauer height, which is 1, 2, or infinite for these Kummer surfaces.

The lower hull of the actual Jacobian Weil coefficients $(j,v_p(c_j)/v_p(q))$ gives the four $H^1$ slopes. Exterior powers give the Jacobian slopes in each degree. On the resolved Kummer K3 the sixteen exceptional classes add slope 1 to degree two, regardless of their Frobenius permutation. That degree's slopes below 1 determine the finite formal Brauer height. For Hilbert schemes, a length-$m$ copy of degree $b$ and slope $s$ contributes degree $b+2(m-1)$ and slope $s+m-1$. The calculation uses symmetric powers for even $b$ and exterior powers for odd $b$.

For a generalized Kummer variety the cohomological identity described above also holds with a variable $v$ for twice the slope. The shift for divisor $d$ of $n$ multiplies by $(t^2v^2)^{n-n/d}$, and the torsion permutation has slope zero and geometric multiplicity $J_4(d)$. Ascending polynomial division by the Jacobian's cohomology recovers each slope bin. Counts are checked against the graded Betti numbers within the work, storage, and exact integer limits. The example compares the ordinary, height-two, and supersingular cases and displays the odd degree-three slopes of two generalized Kummer fourfolds. [Crystalline Newton example](examples/hyperelliptic-crystalline-newton.asl)

#### Exceptional Kummer lattice and its discriminant

A certified genus-two level-2 basis `B` labels the 16 exceptional curves $E_v$ by $v\in\mathbb F_2^4$ in basis order. Their integral intersection pairing is $E_v\cdot E_w=-2\delta_{vw}$. The primitive exceptional lattice $K$ consists of $(\sum_v a_vE_v)/2$ with $a\bmod2$ in the affine code $\mathrm{RM}(1,4)$. It has rank 16 and determinant $+64$. The full Néron–Severi lattice of the resolution also has classes from the abelian surface.

`RVKummerDisc(B)` is the 64-element additive discriminant group $K^*/K=\mathrm{RM}(2,4)/\mathrm{RM}(1,4)$. `rvKummerDisc[B,mask]` uses the six squarefree monomials $x_ix_j$ for $i<j$, with bit zero on the first lexicographic pair. `rvKummerDiscQuadratic(x)` gives $q(x)=-\operatorname{wt}(x)/2\bmod2$ from its truth table, and `rvKummerDiscPair(x,y)` gives $2b(x,y)\bmod2$ by polarization. Both return `PrimeScalar(2)`. `rvKummerDiscFrob[Base,k](x)` pulls back by the inverse of the certified node permutation for $\mathrm{Frob}_{\mathrm{Base}}^k$ and discards the affine terms.

`rvKummerNodeLattice[B,Base,k];` lists the node orbits $O_i$ of that Frobenius power. The fixed exceptional lattice has rank $m$ equal to the number of orbits. If the invariant affine code has dimension $d$, the orbit sums have index $2^d$ in the fixed lattice. The signed determinant is $(-2)^m\prod_i|O_i|/2^{2d}$. The example compares the split lattice with a nontrivial Frobenius shear. [Exceptional Kummer lattice example](examples/hyperelliptic-kummer-nikulin.asl)

#### Finite Weil representation and theta shells

`RVKummerWeil(B)` is the rational group algebra of `RVKummerDisc(B)`, with basis vector `rvKummerWeilBasis(x)` for each discriminant class `x`. Use `zero`, `add`, `neg`, and `rvKummerWeilScale[a,b]` for its rational vector operations, with positive denominator `b`. `rvKummerWeilCoeff(v,x);` reports the coefficient at `x`. `rvKummerWeilT` multiplies $e_x$ by $(-1)^{q(x)}$, and `rvKummerWeilS` sends it to $\sum_y(-1)^{2b(x,y)}e_y/8$. The signature is a multiple of eight, so the Fourier transform has phase one. These actions satisfy $S^2=T^2=(ST)^3=1$. Their $S_3$ character has 15 trivial, 7 sign, and 21 standard constituents. `rvKummerWeilFrob[Base,r]` permutes the basis by certified relative Frobenius on the discriminant group and commutes with both operators.

`rvKummerThetaShell[B,m]` gives the coefficient vector at $q^{m/2}$ of the discriminant-coset theta series of the positive lattice $K(-1)$. Its component at $x$ counts the vectors of squared length $m$ in the corresponding coset. The vector-valued series is modular of weight eight for the stated Weil action. For each affine word, the interpreter adds it to the quadratic representative of the coset and groups the results by Hamming weight. It then computes the degree-$2m$ coefficient of the corresponding product of even- and odd-coordinate square series. The example gives zero-coset counts $1$, $32$, and $8160$ in norms zero, two, and four. [Finite Weil and theta shell example](examples/hyperelliptic-kummer-weil-theta.asl)

#### Marked polarized Kummer divisors and Hilbert schemes

`RVKummerDiv(B)` is the rank-17 lattice generated by the quartic hyperplane $H$, the sixteen exceptional curves $E_v$, and the sixteen tropes $T_a$. The marking `B` must be a certified genus-two level-two basis, and the node and trope indices are the integers 0 through 15 encoding its four basis coordinates. `rvKummerHyperplane[B]`, `rvKummerNodeDiv[B,v]`, and `rvKummerTropeDiv[B,a]` construct these generators. The six nodes on the zero trope are the classes of Mumford degree at most one. The support $S_a$ is their translate by $a$, and $T_a=(H-\sum_{v\in S_a}E_v)/2$. The full geometric Néron–Severi rank is reported by `rvKummerK3Picard`.

`rvKummerDiv[B]{h,b0,...,b15}` writes $(hH+\sum_v b_vE_v)/2$. For integrality, the parity vector of $b$ plus $(h\bmod 2)1_{S_0}$ must give an affine function on the four binary coordinates. Use `zero`, `add`, `neg`, and `rvKummerDivScale[n]` for the lattice operations. `rvKummerDivPair(x,y);` reports the integer intersection pairing, with $H^2=4$ and $E_v^2=-2$. `rvKummerDivFrob[Base,k](x)` permutes exceptional coefficients by the certified relative Frobenius action and fixes $H$.

`rvKummerFujiki[n,j,k,l](x,y);` computes $(x^{[n]}+k\delta)^{2n-j}(y^{[n]}+l\delta)^j$ on the Hilbert scheme of $n$ points. Here $\delta$ is half the Hilbert–Chow exceptional divisor, orthogonal to the surface classes and of Beauville–Bogomolov square $-2(n-1)$. The computation sums over perfect pairings, with coefficient $(2n-1)!!$ on a single class. It requires $n\ge1$ and $0\le j\le2n$, with zero shifts for $n=1$. Coefficients and results use checked signed 128-bit integers. [Marked divisor and Fujiki example](examples/hyperelliptic-kummer-polarized-divisors.asl)

`rvKummerHilbChar[C,n,r,i,ell,k];` gives the coefficients of $\det(1-TF^r\mid H^i(S^{[n]}))$ from degree zero through `k`, reduced modulo the prime `ell`. `rvKummerHilbCharCoeff[C,n,r,i,ell,k]` returns the degree-`k` coefficient as `PrimeScalar(ell)`. These are reductions of the integer Frobenius polynomial, including when `ell` is the residue characteristic. Coefficients above the degree of the polynomial are zero.

The calculation uses division-free symmetric power recurrences to compute the graded power traces modulo $\ell^{1+\lfloor\log_\ell k\rfloor}$ for $k>0$. Big Witt ghost inversion then recovers the coefficients. It divides out the exact power of $\ell$ in each index before inverting the unit part at the remaining precision. This works for small primes at every index within the usual work, storage, and report limits. Frobenius powering works with extension degrees up to unsigned 128 bits, even when the corresponding exact traces exceed their output range. [Modular Hilbert Frobenius example](examples/hyperelliptic-kummer-hilbert-characteristic.asl)

#### Marked unimodular K3 lattice

`RVK3Lattice(B)` realizes the abstract even unimodular lattice of signature $(3,19)$ by gluing the exceptional lattice $K$ to $H^2(A,\mathbb Z)(2)$. The marking `B` labels the exceptional nodes and supplies the reduction modulo 2 of a chosen standard integral symplectic frame. This construction includes the transcendental lattice. An integral Frobenius action requires additional data.

`rvK3Lattice[B]{[w],[b]}` enters the doubled coordinates of $(w,b)/2$. Give the six wedge coordinates in order $(01,02,03,12,13,23)$ and the sixteen exceptional coordinates in node-mask order. Membership requires $b\bmod2$ plus the truth table of $\sum_{i<j}(w_{ij}\bmod2)x_ix_j$ to lie in $\mathrm{RM}(1,4)$. The pairing is $(W(w,w')-\sum_v b_vb'_v)/2$, where the complementary wedge pairs $(01,23),(02,13),(03,12)$ have signs $-,+,-$. The class $e_{02}+e_{13}$ has square $4$ in this normalization.

`rvK3LatticeBasis[B,i]` reads basis vector $0\le i<22$. The basis uses the binary reduced-echelon gluing rows and twice the standard vectors at the nonpivot columns. `rvK3LatticeGram[B];` gives its Gram matrix, determinant $-1$, and index $2^{11}$ over the lattice of the unglued wedge and exceptional generators. `rvK3LatticePair(x,y);` evaluates an integral pairing. Use `add`, `neg`, and `zero[RVK3Lattice(B)]` for the additive group. `rvK3Reflect(r,x)` checks $r^2=-2$ and returns $x+(x,r)r$. `rvK3Transvection[a,k](x)` changes the integral marking by $M=I+k u(Ju)^t$ for the four-bit mask $u=a$, with standard symplectic pairs $(0,2)$ and $(1,3)$. It acts on wedge coordinates by $\wedge^2M^{-t}$ and carries the exceptional coefficient at $v$ to $Mv\bmod2$. [Marked K3 lattice example](examples/hyperelliptic-kummer-unimodular.asl)

#### Marked rational K3 Clifford algebra

`RVK3Clifford(B,p)` is parity `p` of the rational Clifford algebra of the full marked K3 lattice, with `p=0` for even and `p=1` for odd. `rvCliffordVector(x)` embeds a lattice vector in odd parity and satisfies $v^2=(v,v)/2$, so $uv+vu=(u,v)$. Use `rvCliffordOne[B]`, `zero`, `add`, `neg`, and `rvCliffordScale[a,b](x)` for scalar and additive operations, with $b>0$. `rvCliffordMul(x,y)` adds the parities modulo two. `rvCliffordReverse(x)` reverses each monomial's order.

The ordered generators are $e_I+e_J,e_I-e_J$ for each complementary wedge pair in the order above, followed by the sixteen $E_v$ in node-mask order. The six wedge generators have Clifford squares $-2,2,2,-2,-2,2$, and each exceptional generator squares to $-1$. `rvCliffordCoeff[mask](x);` reports the rational coefficient of the ordered product selected by its 22-bit mask. The computation uses sparse monomials. Each numerator and the shared reduced denominator are bounded to 256 bits, with wider intermediate arithmetic. `rvCliffordVolume[B]` is the ordered product of all 22 generators divided by eight. It squares to one, anticommutes with vectors, and gives complementary central idempotents $(1\pm\omega)/2$ of the even algebra.

`rvCliffordAct(g,x)` uses twisted conjugation $(-1)^p g x\widetilde g/(g\widetilde g)$ for a homogeneous actor of parity $p$. It requires $g\widetilde g$ to be a nonzero scalar and the resulting vector to satisfy the coordinate integrality and affine gluing conditions of the marked lattice. A norm-$-2$ vector acts as `rvK3Reflect`, and the normalized volume acts by negation. [Marked Clifford algebra example](examples/hyperelliptic-kummer-clifford.asl)

#### Rational K3 half-spin modules

`RVK3Spinor(B,0)` and `RVK3Spinor(B,1)` are the $+1$ and $-1$ eigenspaces of `rvCliffordVolume[B]`, respectively. Each has dimension $1024$ over $\mathbb Q$. `rvK3SpinorBasis[B,c,j]` gives its basis vector with $0\le j<1024$. Use `zero`, `add`, `neg`, and `rvK3SpinorScale[a,b](s)` with $b>0$. `rvK3SpinCoeff[j](s);` gives an exact rational coordinate. `rvK3SpinAct(g,s)` acts by a homogeneous Clifford element on a spinor with the same marking, changing chirality by the parity of `g` modulo two.

The basis comes from the graded tensor product of three two-dimensional wedge blocks and two sixteen-dimensional octonion blocks. The three low bits of `j` encode the wedge-block grades, the next three the first octonion basis index, the next its block grade, and the last three the second octonion basis index. The volume eigenvalue fixes the second octonion block grade. Each Clifford generator permutes the basis with multipliers in $\{-2,-1,1,2\}$, and sparse products apply the rightmost monomial factor first. For the connection between octonionic multiplication operators and Clifford modules, see [Baez, *The Octonions*](https://arxiv.org/abs/math/0105155).

The nondegenerate full Clifford algebra acts faithfully on the $2048$-dimensional combined space because it is central simple. Its dimension $2^{22}$ equals the dimension of the full endomorphism algebra. Its even part acts on the two volume eigenspaces, and dimensional equality gives the full endomorphism algebra on each, proving irreducibility. [Rational half-spin example](examples/hyperelliptic-k3-rational-half-spin.asl)

#### Rational K3 periods and the Kuga–Satake action

`rvk3period P(u,v);` registers an oriented rational period plane in the marked K3 lattice. Its two ordered generators must have the same marking and a positive definite Gram matrix. Write $D=(u,u)(v,v)-(u,v)^2$. `rvK3PeriodBlade[P]` gives the even Clifford element $\Omega=2uv-(u,v)1$ with $\Omega^2=-D$. `rvK3PeriodCM[P](s)` acts by $\Omega$ on either half-spin module with that marking. The action of $\mathbb Q(\sqrt{-D})$ is faithful and each half-spin module has rank $512$ over this field. The real complex structure is $J=\Omega/\sqrt D$. A change of generators with determinant $k$ scales $\Omega$ by $k$ and $D$ by $k^2$, so $J$ depends on the orientation of the plane.

`rvK3PeriodReport[P];` gives $D$ and the saturation index $I$, the gcd of the two-by-two minors in the integral unimodular basis. The primitive transcendental lattice has discriminant $D/I^2$. Its orthogonal Néron–Severi lattice has rank $20$, signature $(1,19)$ and determinant $-D/I^2$. `rvK3PeriodIs11[P](x);` tests orthogonality to both generators. The report also gives the half-spin dimension and the dimension $2^{21}$ of the full regular even Clifford realization of weight one.

The supplied plane defines an abstract marked Hodge structure. The finite-field curve supplies the level-two marking. The commands report lattice invariants and construct the rational Clifford action, with the existing exact coefficient bounds for blades and spinors. See Kuga and Satake, *Abelian varieties attached to polarized K3-surfaces* (1967), Deligne, *La conjecture de Weil pour les surfaces K3* (1972), and Huybrechts, *Lectures on K3 Surfaces* (2016). [Rational period example](examples/hyperelliptic-k3-kuga-satake-period.asl)

#### Shioda–Inose CM data

`rvK3CMForm[P];` saturates the period plane in the integral K3 lattice and reduces its oriented binary form $[a,b,c]$ by $\mathrm{SL}_2(\mathbb Z)$. The basis has Gram matrix $\begin{pmatrix}2a&b\\b&2c\end{pmatrix}$. The reduced form satisfies $|b|\le a\le c$ with $b\ge0$ when $|b|=a$ or $a=c$. The report gives the content $m=\gcd(a,b,c)$, the discriminant $d=b^2-4ac$, and the CM order discriminant $d/m^2$. `rvK3CMBasis[P,j]` gives the corresponding lattice vector for $j=0$ or $1$, within the lattice coordinate bounds. `rvK3CMEquivalent[P,Q];` compares the oriented transcendental lattices by their proper reduced forms. These forms classify the associated singular complex K3 surfaces up to isomorphism, independently of their lattice markings.

`rvK3CMPeriods[P];` gives the Shioda–Mitani periods $\tau_1=(-b+\sqrt d)/(2a)$ and $\tau_2=(b+\sqrt d)/2$ with positive imaginary part. For their product abelian surface $A$, the Shioda–Inose K3 has transcendental lattice $T_A$, and the Kummer resolution of $A$ has $T_A(2)$. `rvK3CMClassNumber[P];` counts the proper classes of the primitive order by enumerating primitive reduced forms up to $a\le\sqrt{|d/m^2|/3}$ within the computation budget. Form and period queries use direct reduction.

See Shioda and Mitani, *Singular abelian surfaces and binary quadratic forms* (1974), and Shioda and Inose, *On singular K3 surfaces* (1977). [Shioda–Inose CM example](examples/hyperelliptic-k3-shioda-inose.asl)

#### Imaginary quadratic ring class groups

`rvcmorder O[D];` names the imaginary quadratic order of negative discriminant $D\equiv 0,1\pmod 4$. `RVCMClass(O)` is its proper invertible ideal class group. `rvCMClass[O,a,b,c]` reduces a primitive positive definite form of discriminant $D$ to its canonical proper representative. `rvCMPrincipal[O]` gives the identity. Use `rvCMCompose(x,y)`, `rvCMInverse(x)`, and `rvCMPow[n](x)` for composition, inversion, and signed powers within the same named order.

The form $[a,b,c]$ represents the ideal generated by $a$ and $(-b+\sqrt D)/2$. Composition multiplies their generators in the order, removes a common integral scale, and recovers the primitive form by Bézout reduction. This works with shared leading factors and nonmaximal orders. `rvCMClassOrder(x);` gives the element order. `rvCMClassGroup[O];` gives the class number, the finite abelian invariant factors, and the number of genera. `rvCMSameGenus(x,y);` tests whether $xy^{-1}$ belongs to the subgroup of squares. Reports that enumerate forms use the computation and finite-value limits.

`rvK3CMClass[O,P]` gives the proper class of the primitive part of the oriented transcendental form of a rational K3 period. Its primitive discriminant must agree with $D$. For the correspondence between forms and ideals, see Gauss, *Disquisitiones Arithmeticae*, Section V, and Dedekind, *Über die Theorie der ganzen algebraischen Zahlen*. [Ring class group example](examples/hyperelliptic-k3-cm-class-group.asl)

#### Cyclic CM isogeny neighbors

`rvCMNeighbor[Otarget,ell,j](x)` returns the target proper class of the cyclic isogeny of prime degree $\ell$ at address $j$. Its discriminant must agree with the named target order. `rvCMNeighborData[ell,j](x);` gives one address, including the reduced form, discriminant, direction, and content. `rvCMNeighbors[ell](x);` groups reduced targets with their multiplicities. `rvCMNeighborCounts[ell](x);` counts ascending, horizontal, and descending cyclic kernels. These counts retain distinct kernels with the same target, including those related by extra Gaussian or Eisenstein units.

For the source form $[a,b,c]$ of discriminant $D$, use the positively oriented periods $z_1=a$ and $z_2=(b+\sqrt D)/2$ with positive imaginary square root. Their norm form divided by $a$ is $[a,b,c]$. This convention identifies the torus with the inverse of the ideal class represented by $a,(-b+\sqrt D)/2$ in Gauss composition. For $0\le j<\ell$, the index-$\ell$ sublattice is spanned by $\ell z_1,jz_1+z_2$, and the kernel of the isogeny from the source is generated by $(jz_1+z_2)/\ell$ modulo its lattice. The address $j=\ell$ is infinity, with sublattice spanned by $z_1,\ell z_2$ and kernel generated by $z_1/\ell$.

Restrict the norm form to the sublattice, divide by its content $g$, and reduce properly. The target has discriminant $D'=\ell^2D/g^2$. The cases $g=1,\ell,\ell^2$ give descending, horizontal, and ascending edges, respectively. At the top of the $\ell$-conductor tower there are $1+\left(\frac{D}{\ell}\right)$ horizontal kernels, using the Kronecker symbol. Below the top, one kernel ascends and $\ell$ descend. For background, see Kohel, *Endomorphism rings of elliptic curves over finite fields* (1996), and Sutherland, *Isogeny volcanoes* (2013). [Cyclic CM isogeny example](examples/hyperelliptic-cm-cyclic-isogeny-volcanoes.asl) · [Sutherland, *Isogeny volcanoes*](https://arxiv.org/abs/1208.5370)

#### Classical elliptic modular polynomials

`rvModularPolyDegree[ell];` gives the bidegree $(\ell+1,\ell+1)$ of the classical modular polynomial $\Phi_\ell$ at prime level. `rvModularPolyCoeff[ell,i,j];` gives the exact integer coefficient of $X^iY^j$. `rvModularPolyEval[ell](x,y)` evaluates the polynomial in a named precision-one `RVUnram` field. `rvModularPolySpecCoeff[ell,j](x)` gives the coefficient of $Y^j$ in $\Phi_\ell(x,Y)$ in the same field. Coefficients outside the bidegree are zero.

The construction starts with $j=E_4^3/\Delta$. The power sum of the $\ell+1$ cyclic quotients is $j(q^\ell)^k+\ell\sum_n [q^{\ell n}]j(q)^k q^n$. Matching the nonpositive Laurent coefficients in the triangular basis $1,j,j^2,\ldots$ recovers that power sum as a polynomial in $j$. Newton identities recover the integral coefficients by exact division before any field reduction. Exact construction is available for levels $2,3,5,7$ within the work budget. In characteristic $\ell$, evaluation uses the Kronecker congruence $(X^\ell-Y)(X-Y^\ell)$, also at larger prime levels.

At level two the fibers at $j=0$ and $1728$ are $(Y-54000)^3$ and $(Y-1728)(Y-287496)^2$, respectively. Their multiplicities come from the three cyclic kernels at the Eisenstein and Gaussian CM vertices. For the product construction and the Kronecker congruence, see Cox, *Primes of the Form $x^2+ny^2$*, and Bröker, Lauter and Sutherland, *Modular Polynomials via Isogeny Volcanoes*. [Classical modular polynomial example](examples/elliptic-modular-polynomials.asl)

#### Ring class polynomials by modular norms

`rvcmclasspoly H[D];` registers the integral ring class polynomial at a class-number-one fundamental discriminant. The available seeds are $-3,-4,-7,-8,-11,-19,-43,-67,-163$. `rvcmclasspolylift T[H,ell];` constructs the polynomial for discriminant $\ell^2D$ from $H$ at prime level $\ell$. If the source has an ascending $\ell$-isogeny, register the polynomial for $D/\ell^2$ before lifting. Both conductor paths can register the same discriminant when their exact coefficients agree.

The norm $\operatorname{Res}_X(H_D(X),\Phi_\ell(X,Y))$ is the determinant of multiplication by $\Phi_\ell$ in $\mathbb Z[Y][X]/(H_D)$. The implementation uses fraction-free elimination and monic polynomial division. Each direction with $c$ cyclic kernels per source class contributes $H_{D'}^{h(D)c/h(D')}$, where $h$ is the proper class number. The lift divides out the horizontal and ascending factors, extracts the integral monic root for the descending multiplicity, and checks the reconstructed norm and degree. This keeps the extra unit multiplicities at discriminants $-3$ and $-4$. All constructed descendants come from the norm calculation.

`rvCMClassPolyData[H];` reports the discriminant and degree. `rvCMClassPolyCoeff[H,i];` gives the exact coefficient of $Y^i$. `rvCMClassPolyModCoeff[H,F,i]` reduces that coefficient in a named precision-one `RVUnram` field, and `rvCMClassPolyEval[H](x)` evaluates there. Reduction retains root multiplicities. For the singular moduli and class polynomials, see Cox, *Primes of the Form $x^2+ny^2$*. For the kernel counts, see Bröker, Lauter and Sutherland, *Modular Polynomials via Isogeny Volcanoes*. [Ring class polynomial example](examples/ring-class-polynomial-norms.asl)

#### Ordinary elliptic endomorphism orders

`rvecendorder R[E,H];` certifies the endomorphism order of a named short Weierstrass curve $E/\mathbb F_q$ using a certified ring class polynomial $H_D$. It computes $j(E)$ and the Frobenius trace $a=q+1-\#E(\mathbb F_q)$. It requires $p\nmid a$, $H_D(j(E))=0$, a split characteristic $p$ prime to $D$, and $a^2-4q=Dm^2$ with $m>0$. The ordinary reduction and canonical lift correspondence identifies the endomorphism order with $\mathcal O_D$. In the ordinary case the geometric endomorphisms commute with Frobenius and are defined over the field of the curve. `rvECEndCheck[E,H];` reports the evidence and the certification result, including the polynomial residual. `rvECEndJ[E]` returns $j(E)$ in the curve's named field.

`rvECEndData[R];` reports the discriminant, Frobenius trace and index $[\mathcal O_D:\mathbb Z[\pi]]=m$. The type `RVECEnd(R)` gives the abstract integral order in the basis $1,\omega=(b+\sqrt D)/2$, where $b\in\{0,1\}$ is the parity of $D$. Enter coordinates with `rvECEndElem[R,x,y]` and use `zero`, `add`, `neg`, `rvECEndMul(x,y)` and `rvECEndConj(x)`. `rvECEndTrace(x);` and `rvECEndNorm(x);` report the quadratic invariants. `rvECEndFrobenius[R]` gives $\pi=(a+m\sqrt D)/2$ up to conjugation. The certificate verifies $\pi^2-a\pi+q=0$ and gives the order and Frobenius coordinates in this abstract basis.

See Deuring, *Die Typen der Multiplikatorenringe elliptischer Funktionenkörper*, Abh. Math. Sem. Hansischen Univ. **14** (1941), 197–272, and Cox, *Primes of the Form $x^2+ny^2$*, for ordinary reduction and its relation to ring class polynomials. [Ordinary endomorphism order example](examples/elliptic-deuring-endomorphism-orders.asl)

#### Serre–Tate canonical CM lifts

`rvECEndCanonicalJ[M,W]` returns the canonical $j$-invariant of the certified ordinary curve in the named truncated unramified ring `W`. Its residue presentation must agree with the curve's field, including the defining polynomial modulo $p$. The operation evaluates every integral coefficient of $H_D$ modulo the full $p^s$, checks that the derivative at the residue root is a unit, and applies Hensel–Newton lifting. Uniqueness gives compatibility under truncation. The computation checks that the Frobenius conjugate still satisfies $H_D$ and reduces to the conjugate of the original $j$.

For certified Frobenius $\pi=c+b\omega$ with trace $a$ and norm $q$, `rvECEndUnitRoot[M,W]` gives the root $\alpha$ of $X^2-aX+q$ reducing to $a$ modulo $p$. `rvECEndNonunitRoot[M,W]` gives $\beta=q/\alpha$. The two maps `rvECEndEmbed[M,W](x)` and `rvECEndEmbedConj[M,W](x)` evaluate an element of `RVECEnd(M)` using $\omega=(\alpha-c)/b$ or $(\beta-c)/b$. The coefficient $b$ is a unit for an ordinary certificate. The implementation checks the quadratic relation and the root sum and product.

`rvECEndCanonicalPhi[M,W]` evaluates $\Phi_p(J,\sigma J)$ in the full ring using integral coefficients at available levels. The short Weierstrass certificates start at $p=5$, so this check applies there at $p=5,7$. `rvModularPolyWittEval[ell](x,y)` also evaluates the integral polynomial in a shared named unramified ring at levels $2,3,5,7$. See Serre and Tate, *Good reduction of abelian varieties*, Katz, *Serre–Tate local moduli*, and Satoh, *The canonical lift of an ordinary elliptic curve over a finite field and its point counting*. [Canonical CM lift example](examples/elliptic-serre-tate-canonical-lifts.asl)

#### Marked Mukai lattice

`rvK3FromKummerDiv(d)` embeds `RVKummerDiv(B)` into `RVK3Lattice(B)` by sending the hyperplane coefficient $h$ to $h(e_{02}+e_{13})$ in the doubled wedge coordinates and preserving the exceptional coefficients. `rvK3ToKummerDiv(c)` inverts this on that sublattice. The zero trope incidence vector has the same quadratic polarization as the Weil form in the marking, so the two gluing conditions agree up to an affine term.

`RVMukai(B)` is the marked rank-24 lattice of signature $(4,20)$, with all of `RVK3Lattice(B)` in its middle coordinate. `rvMukai[r,s](c)` forms $(r,c,s)$ and `rvMukaiC(x)` returns $c$. Use `zero`, `add`, and `neg` for its additive group. `rvMukaiPair(x,y);` gives $c\cdot d-rt-r's$ for $y=(r',d,t)$. For sheaf vectors its negative is the Euler pairing. The full marked lattice also includes nonalgebraic directions.

`rvMukaiLine(d)` gives $(1,d,1+d^2/2)$ using the embedded divisor, and `rvMukaiIdeal[B,n]` gives $(1,0,1-n)$ for $n\ge0$. `rvMukaiTensor(x,d)` acts by tensoring with the line bundle of `d`. `rvMukaiTwist(v,x)` checks $v^2=-2$ and returns $x+(x,v)v$, the cohomological action of a spherical twist when `v` is represented by a spherical object. `rvMukaiHilbFujiki[n,j](x,y);` gives the mixed top intersection in degrees $2n-j$ and $j$ on the orthogonal lattice of $(1,0,1-n)$ for $n\ge2$. It uses the Mukai pairing as the Beauville–Bogomolov form and the same perfect-pairing sum as `rvKummerFujiki`. [Marked Mukai example](examples/hyperelliptic-kummer-mukai.asl)

##### Perpendicular lattices of Mukai vectors

`rvmukaimoduli M(v);` accepts a primitive nonzero marked Mukai vector of nonnegative square $s=(v,v)$. For positive $s$ it registers the integral lattice $L=v^\perp$. Primitivity is tested by the gcd of pairings with the actual unimodular ambient basis. Unimodular column operations give an integral basis of $L$ and a vector $w$ with $(w,v)=1$. In that ambient basis the component of $w$ orthogonal to $L$ is $v/s$, of square $1/s$, so unimodularity gives $\det L=s$. `rvMukaiModuliBasis[M,j]` returns basis element $j$ for $0\le j<23$, and `rvMukaiModuliGlue[M]` returns $g=sw-v$ in the same `RVMukai` type. The class of $g/s$ generates $L^*/L$ because $v$ is primitive, and $(g/s)^2=w^2-1/s\equiv-1/s\pmod{2\mathbb Z}$.

`rvMukaiModuliLattice[M];` gives $n=1+s/2$, dimension $2n$, rank 23, signature $(3,20)$, determinant $s$, and the cyclic discriminant form. `rvMukaiModuliDiv[M](x);` gives the nonnegative gcd of the pairings of $x$ with the integral basis. `rvMukaiModuliFujiki[M,j](x,y);` gives the mixed top intersection $x^{2n-j}y^j$ by summing products of pairings over perfect matchings. `rvMukaiModuliRR[M](x);` evaluates $\binom{q(x)/2+n+1}{n}$ as a polynomial in the upper entry, including negative upper integers. These commands require perpendicular inputs on the registered marking. The defining square, returned Mukai coordinates, and arithmetic results use checked signed 128-bit integers. The example compares the half exceptional class for an ideal-sheaf vector with a primitive vector whose rank and last coordinate are both even. [Mukai perpendicular example](examples/hyperelliptic-mukai-moduli.asl)

For a primitive effective algebraic Mukai vector on a projective complex K3 with generic stability and a nonempty smooth projective moduli space of stable sheaves, the Mukai–Yoshioka isometry identifies $L$ with its integral second cohomology and the space has K3 Hilbert-scheme deformation type. The commands also apply to abstract marked vectors in the ambient lattice. [Yoshioka, *Moduli spaces of stable sheaves on abelian surfaces*, §8](https://arxiv.org/abs/math/0009001) [Ellingsrud, Göttsche and Lehn, *On the cobordism class of the Hilbert scheme of a surface*](https://arxiv.org/abs/math/9904095)

##### Isotropic Mukai quotients

When $v^2=0$, the registered lattice is $v^\perp/\mathbb Z v$. The integral Bézout vector $w_0$ has $(v,w_0)=1$, so $w=w_0-(w_0^2/2)v$ is integral, isotropic, and still pairs to 1 with $v$. `rvMukaiModuliDual[M]` returns $w$. The pair $v,w$ spans a unimodular hyperbolic plane. Its perpendicular is the saturated rank-22 lattice of signature $(3,19)$ and determinant $-1$ returned by `rvMukaiModuliBasis[M,j]`. The report has $n=1$, dimension 2, and trivial discriminant, with `rvMukaiModuliGlue[M]` equal to zero.

`rvMukaiQuotient[M](x)` sends $x\in v^\perp$ to `RVMukaiQuotient(M)`. The canonical representative is $x-(x,w)v$, which is invariant under integral shifts by $v$. Use `rvMukaiQuotientLift(q)` to obtain that representative as an ordinary `RVMukai` vector and `rvMukaiQuotientPair(q_1,q_2);` to pair two classes of the same quotient. The quotient supports `zero`, `add`, `neg`, and exact equality. Divisibility, Fujiki intersections, and Riemann–Roch accept either quotient classes or perpendicular Mukai representatives for the same declaration. The Fujiki pairing in dimension two is the intersection form, and Riemann–Roch evaluates $2+q(x)/2$. [Isotropic Mukai example](examples/hyperelliptic-mukai-isotropic.asl)

For a primitive effective algebraic isotropic vector with a generic polarization and a nonempty smooth projective moduli space of stable sheaves on a complex K3 surface, Mukai identifies the quotient with the second integral cohomology of the resulting K3 surface. An untwisted universal sheaf and an untwisted Fourier–Mukai equivalence require additional hypotheses. [Mukai, *On the moduli space of bundles on K3 surfaces, I*](https://www.kurims.kyoto-u.ac.jp/~mukai/paper/Tata.pdf)

##### The O’Grady tenfold lattice

`RVOG10(M)` uses a registered primitive Mukai vector $w$ of square 2. With $W=w^\perp$ and $g=2b-w$ from the integral Bézout vector $b$, introduce an exceptional class $\sigma$ orthogonal to $W$ with $\sigma^2=-6$. The lattice adjoins $G=(g+\sigma)/2$ to $W\oplus\mathbb Z\sigma$. Since $g^2\equiv-2\pmod 8$ and $(g,W)$ is even, the extension is even integral. It has index 2, signature $(3,21)$, and determinant $-3$. The class of $\sigma/3$ generates the discriminant group $\mathbb Z/3$ with quadratic value $-2/3$ modulo $2\mathbb Z$.

Use `rvOG10[M,k](u)` for $u+kG$ with $u\in W$, `rvOG10Exc[M]` for $\sigma=2G-g$, and `rvOG10Basis[M,j]` for the 23 marked perpendicular basis vectors followed by $G$. The type supports `zero`, `add`, `neg`, and exact equality on the same declaration. `rvOG10Pair(x,y);` gives the Beauville pairing and `rvOG10Div(x);` gives the nonnegative gcd of the pairings with the whole integral basis. `rvOG10Lattice[M];` reports the lattice invariants. `rvOG10Fujiki[j](x,y);` gives the intersection $x^{10-j}y^j$ with Fujiki constant 945. `rvOG10RR(x);` evaluates $\binom{q(x)/2+6}{5}$, including negative upper entries. The example gives $(\sigma,G)=-3$, $\operatorname{div}(\sigma)=3$, $\int\sigma^{10}=-7348320$, and $\chi(\mathcal O)=6$. [O’Grady tenfold example](examples/hyperelliptic-ogradys-tenfold.asl) [Rapagnetta, *On the Beauville form of the known irreducible symplectic varieties*](https://arxiv.org/abs/math/0606409) [Ríos Ortiz, *Riemann–Roch polynomials of the known hyperkähler manifolds*](https://arxiv.org/abs/2006.09307)

Geometrically, for an OLS triple on a projective complex K3 surface with the effective algebraic vector $2w$ and a $2w$-generic polarization, the lattice is the second integral cohomology of the symplectic resolution of the moduli space of semistable sheaves of vector $2w$. The same lattice construction applies to any primitive square-two vector in the abstract marking. [Perego and Rapagnetta, *Deformation of the O’Grady moduli spaces*](https://arxiv.org/abs/1008.0190) [*O’Grady tenfolds as moduli spaces of sheaves*, §3](https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/ogrady-tenfolds-as-moduli-spaces-of-sheaves/390DF9D375D2DDED972F3F0D0FEE38F0)

##### Numerical walls in an O’Grady plane

`rvOG10Wall(d);` divides a nonzero integral class by its content and reports the primitive square, divisibility, and wall type. The pairs $(-2,1)$ and $(-6,3)$ are divisorial, $(-4,1)$ covers a codimension-five locus, and $(-24,3)$ covers a codimension-three locus. The remaining pairs are regular. These numerical types give geometric wall divisors when the class is of type $(1,1)$ on an O’Grady tenfold. `rvOG10Reflect(d,x)` applies the integral reflection $x-2(x,d)d/d^2$ for a divisorial ray and preserves the positive cone component. [Mongardi and Onorati, Proposition 5.4](https://arxiv.org/abs/2010.12511)

`rvOG10Walls(h0,h1);` reports crossings in order along the closed segment from `h0` to `h1`. Both endpoints must be positive, in the same cone component, and span a plane of signature $(1,1)$. The search includes primitive walls in their rational plane and tests integrality in the full O’Grady lattice. This recovers walls from non-saturated endpoint presentations. `rvOG10WallAt[j](h0,h1)` returns the zero-based wall vector oriented to pair nonnegatively with `h0` and nonpositively with `h1`.

Writing $A=h_0^2$, $B=(h_0,h_1)$, $C=h_1^2$, and $\Delta=B^2-AC$, a wall of square $-s$ gives pairings $u=(d,h_0)\ge0$ and $v=-(d,h_1)\ge0$ satisfying $Cu^2+2Buv+Av^2=\Delta s$. The search bounds one pairing and recovers the other from $(Av+Bu)^2=\Delta(u^2+As)$ by exact integer square root. It reports $u/(u+v)$ in lowest terms, including endpoint values 0 and 1. Work and output limits apply. The example recovers the two small wall types associated to $D=(3,2H,3)$ when $H^2=2$, with $D$ itself of square $-10$ and regular. [Wall example](examples/hyperelliptic-og10-walls.asl) [Mongardi, 2023 erratum, Propositions 1–2](https://content.algebraicgeometry.nl/2023-2/2023-2-008.pdf)

##### The abelian Mukai lattice and O’Grady sixfolds

`RVAbMukai(B)` is the marked abelian Mukai lattice $H^{\mathrm{even}}(A,\mathbb Z)\cong U^4$. The literal `rvAbMukai[B,r,s]{u01,u02,u03,u12,u13,u23}` gives the rank, last Mukai coordinate, and six middle wedge coordinates. `rvAbMukaiBasis[B,j]` orders its eight axes by rank, wedge coordinates, then last coordinate. The group supports `zero`, `add`, and `neg`. `rvAbMukaiPair(x,y);` pairs complementary wedge coordinates with signs $-,+,-$ on $(01,23),(02,13),(03,12)$ and subtracts $rs'+sr'$. `rvAbMukaiKummerH2(x)` sends the middle class into `RVK3Lattice(B)` using doubled wedge coordinates $2u$ and zero exceptional coordinates. This doubles the middle pairing and takes the marked principal polarization to the Kummer hyperplane.

`rvog6 M(w);` takes a primitive abelian Mukai vector $w$ of square 2 and forms $w^\perp\oplus\mathbb Z e$. Here $e$ is half the exceptional divisor and has square $-2$ and divisibility $2$. Use `rvOG6[M,k](u)` for $u+ke$, `rvOG6ExcHalf[M]` for $e$, and `rvOG6Exc[M]` for $2e$. `rvOG6Basis[M,j]` returns the seven integral perpendicular basis vectors followed by $e$. The resulting `RVOG6(M)` supports `zero`, `add`, `neg`, and exact equality. `rvOG6Pair(x,y);` gives the Beauville pairing and `rvOG6Div(x);` computes divisibility in the full basis. `rvOG6Lattice[M];` reports rank 8, signature $(3,5)$, determinant $-4$, and discriminant form $(\mathbb Z/2)^2$ with diagonal $(-1/2,-1/2)$ modulo $2\mathbb Z$.

`rvOG6Fujiki[j](x,y);` gives the mixed degree-six intersection $x^{6-j}y^j$ by weighting each perfect matching by 4, so the Fujiki constant is $60$. `rvOG6RR(x);` evaluates $4\binom{q(x)/2+3}{3}$, including negative upper entries. The example computes $\int e^6=-480$, $\int(2e)^6=-30720$, and $\chi(\mathcal O)=4$. [O’Grady sixfold example](examples/hyperelliptic-ogrady-sixfold.asl) [Rapagnetta, *Topological invariants of O’Grady’s six dimensional irreducible symplectic variety*, Theorem 3.3.1 and §3.5](https://arxiv.org/abs/math/0406026) [Ríos Ortiz, *Riemann–Roch polynomials of the known hyperkähler manifolds*](https://arxiv.org/abs/2006.09307)

For effective algebraic $2w$ on a complex projective abelian surface with a $2w$-generic polarization under the OLS hypotheses, this lattice is the integral second cohomology of the symplectic resolution of the Albanese fiber of the sheaf moduli space. The marking supplies the underlying abstract lattice. Hodge classes and universal sheaves require geometric data. [Perego and Rapagnetta, *Deformation of the O’Grady moduli spaces*](https://arxiv.org/abs/1008.0190) [Perego and Rapagnetta, *The second integral cohomology of moduli spaces of sheaves on K3 and Abelian surfaces*](https://arxiv.org/abs/2012.10649)

##### O’Grady sixfold monodromy and walls

`rvOG6Orbit(x);` reports the content and the primitive square, divisibility, and marked discriminant residue of an `RVOG6` class. `rvOG6SameOrbit(x,y);` compares two classes of the same declaration and reports both stable and full oriented monodromy equivalence. Content is the positive gcd of the abelian Mukai coordinates and the half-exceptional coefficient, since $w^\perp$ is saturated. The zero class has its own orbit.

For a primitive class $x=u+ke$, divisibility one gives residue $(0,0)$. At divisibility two, $u\bmod 2=bw\bmod 2$ and the residue of $x/2$ is $(b,k\bmod 2)$ in the ordered perpendicular-glue and half-exceptional basis. Its quadratic value is $-(b+(k\bmod 2))/2$ modulo $2\mathbb Z$. This follows by reducing the unimodular pairing modulo two and taking the orthogonal complement of $w^\perp\bmod 2$. Every nonzero residue has nonzero quadratic value, so a primitive isotropic vector has divisibility one. Eichler’s criterion gives the stable orbits from the primitive square and marked residue because the lattice contains two hyperbolic planes. Full monodromy is $O^+$ and can exchange the two residues of square $-1/2$. The reflection in the sum of the two orthogonal square-$-2$ generators realizes that exchange. [Mongardi and Rapagnetta, Lemma 2.6 and Theorem 5.5](https://arxiv.org/abs/1909.07173)

`rvOG6Wall(d);` reports the primitive wall type. Square and divisibility pairs $(-2,2)$ and $(-4,2)$ are divisorial, $(-2,1)$ is a small wall, and the remaining pairs are regular. These are geometric wall divisors when the class is of type $(1,1)$. `rvOG6Reflect(d,x)` applies $x-2(x,d)d/d^2$ for a divisorial ray and preserves the positive orientation. The example compares the distinct marked residues in one full orbit and applies the classifier to the geometric exceptional class $2e$. [Monodromy example](examples/hyperelliptic-og6-monodromy.asl) [Mongardi and Rapagnetta, Proposition 6.8](https://arxiv.org/abs/1909.07173)

##### Polarized O’Grady sixfold cusps

`rvog6cusp P(h,e,f);` registers a marked isotropic plane in the type-IV lattice $h^\perp$. The three classes belong to one `RVOG6` declaration. The polarization $h$ is primitive and positive, and $e,f$ are orthogonal null vectors perpendicular to $h$. Their order orients $I=\mathbb Z e+\mathbb Z f$. The gcd of their Plücker minors in the saturated ambient lattice must be 1. The declaration computes the integral quotient $K=(h^\perp\cap I^\perp)/I$ by taking saturated kernels and splitting off $I$ with integral dual coordinate functionals. `rvOG6CuspGram[P];` reports the rank-three negative-definite Gram matrix and its determinant. `rvOG6CuspBasis[P,j]` for $0\le j\le 2$ returns lifts of the three quotient basis classes. This describes the specified marked rational boundary plane. [Gritsenko, Hulek and Sankaran, *Moduli of K3 surfaces and irreducible symplectic manifolds*, Theorem 5.5](https://arxiv.org/abs/1012.4155)

`rvOG6CuspTransvection[P,j](u,x)` uses the axis $a=e$ for $j=0$ and $a=f$ for $j=1$. It requires $u\in h^\perp\cap I^\perp$ and sends $x$ to $x+(a,x)u-(u,x)a-(u,u)(a,x)a/2$. The full lift $u$ matters. For instance, $T_{e,f}$ is the central transformation $C_1$. `rvOG6CuspCentral[P,k](x)` gives $C_k(x)=x+k((e,x)f-(f,x)e)$. Both are integral isometries of the marked lattice fixing $h$ and $I$ pointwise. Scaling the formulas continuously from zero preserves the positive orientation. With rightmost operation applied first, $T_{e,u}T_{f,v}T_{e,-u}T_{f,-v}=C_{(u,v)}$. The example has $K\cong\langle-2\rangle^3$ and demonstrates the commutator on a vector pairing to 1 with the first axis. [Type-IV cusp example](examples/hyperelliptic-og6-type-iv-cusp.asl) [Gritsenko, Hulek and Sankaran, *Abelianisation of orthogonal groups and the fundamental group of modular varieties*, §3](https://arxiv.org/abs/0810.1614)

##### Integral nilpotent limits of polarized OG6 lattices

`rvog6nilpotentii N[P,k];` takes the logarithm $N(x)=k((e,x)f-(f,x)e)$ with $k\ne0$ at the oriented cusp `P`. `rvog6nilpotentiii N[P,j](u);` uses $a=e$ for $j=0$ or $a=f$ for $j=1$. It requires $h\cdot u=a\cdot u=0$ and $u^2>0$, and sets $N(x)=(a,x)u-(u,x)a$. This has $N^2(x)=-u^2(a,x)a$ and $N^3=0$. The direction may pair nontrivially with the other cusp axis. `rvOG6NilpotentPower[N,j](x)` evaluates $N^j(x)$ and `rvOG6NilpotentExp[N,m](x)` gives the integral isometry $x+mN(x)+m^2N^2(x)/2$ for signed $m$. Both act on the marked lattice and preserve the polarization line and its perpendicular.

`rvOG6NilpotentData[N];` gives the weight-graded ranks in degrees 0 through 4, and the Smith invariant factors and saturation indices of the images of $N$ and $N^2$ on the integral lattice $H=h^\perp$. The images are computed on a saturated integral basis of $H$. In Type II the saturated filtration is $W_0=0$, $W_1=I$, $W_2=H\cap I^\perp$, $W_3=W_4=H$. In Type III it is $W_0=W_1=\mathbb Z a$, $W_2=W_3=H\cap a^\perp$, $W_4=H$. Over $\mathbb Q$, $N$ in Type II and $N^2$ in Type III identify the highest and lowest nonzero graded pieces, as required by the centered weight filtration. `rvOG6NilpotentInWeight[N,j](x);` tests membership for a class in $H$, also accepting indices below 0 or above 4. The example shows an image of index 2 in $W_0$ in Type III, and the tests give index 4 when the polarization forces $(a,H)=2\mathbb Z$. A limiting Hodge filtration or a geometric degeneration requires additional data. [Nilpotent limit example](examples/hyperelliptic-og6-nilpotent-limits.asl) [Schmid, *Variation of Hodge structure*, §6](https://doi.org/10.1007/BF01389674) [Cattani and Kaplan, *Degenerating variations of Hodge structure*](https://www.numdam.org/item/AST_1989__179-180__67_0.pdf)

##### Kuga–Satake spin monodromy of polarized OG6 lattices

`RVOG6Clifford(P,p)` is the parity-$p$ part, for $p=0,1$, of the integral Clifford algebra of the polarized lattice $H=h^\perp$ of cusp `P`. The ordered seven-generator blade basis uses $v^2=(v,v)/2$ and $vw+wv=(v,w)$. `rvOG6PolarizedBasis[P,j]` gives the $j$th generator as a lattice vector for $0\le j<7$. `rvOG6CliffordVector[P](x)` embeds any vector of $H$ in the odd part, and `rvOG6CliffordBasis[P,mask]` gives the ordered blade with the selected generator bits for $0\le\mathrm{mask}<128$. The empty blade is the unit. Use `rvOG6CliffordMul(x,y)`, `rvOG6CliffordScale[k](x)`, `rvOG6CliffordCoeff[mask](x);`, and `zero`, `add`, `neg`, and exact equality.

`rvOG6SpinLog[N]` gives the even Clifford element $s=kfe$ in Type II, with the oriented cusp order $(e,f)$, and $s=ua$ in Type III for direction $u$ and selected axis $a$. In both cases $s^2=0$ and $[s,x]=N(x)$ for $x\in H$. `rvOG6SpinExp[N,m]` gives $1+ms$ for signed $m$. Left multiplication acts on the even regular module, and conjugation by this unit sends $x$ to $\exp(mN)x$ in the odd part. The example applies both operations and checks the conjugation identity on an integral generator. [Spin monodromy example](examples/hyperelliptic-og6-spin-monodromy.asl)

`rvOG6SpinData[N];` computes the Smith factors and rank of left multiplication $S=L_s$ on the 64-dimensional even regular module. For its weight-one filtration, $W_0$ is the saturation of $\operatorname{im}S$, $W_1=\ker S$, and $W_2=C^+(H)$. The graded ranks in degrees 0, 1, 2 are $(16,32,16)$ in Type II and $(32,0,32)$ in Type III. The integral Clifford marking has an associated Kuga–Satake abelian dimension of 32. A geometric degeneration and its limiting Hodge filtration require further data. [Schreieder and Soldatenkov, *The Kuga–Satake construction under degeneration*](https://arxiv.org/abs/1710.02102)

##### Kuga–Satake exterior cohomology monodromy

`RVKSExterior(P,k)` is cohomological degree $k$ of the exterior algebra of the rank-64 even regular spin lattice of `P`. This is the marked cohomology representation of a Kuga–Satake abelian variety of dimension 32. Degrees outside $0\le k\le64$ contain only a formal zero. In `rvKSExteriorBasis[P,mask]` the unsigned 64-bit mask selects spin basis vectors in increasing even PBW mask order. The empty mask is the cup unit. `rvKSExteriorEmbed(x)` identifies the even Clifford lattice with degree one, and `rvKSExteriorToSpin(x)` is its inverse there. `rvKSExteriorCup(x,y)` is the alternating product. `rvKSExteriorScale[m](x)`, `rvKSExteriorCoeff[mask](x);`, `zero`, `add`, and `neg` operate on sparse integral coefficients.

For $S=L_s$, `rvKSExteriorPower[N,j](x)` applies the $j$th power of the derivation $D$ that applies $S$ to one factor at a time. `rvKSExteriorExp[N,m](x)` applies $1+mS$ to every factor, giving $\exp(mD)$ integrally. `rvKSExteriorData[N,k];` gives the dimension $\binom{64}{k}$, nilpotency index $1+\min(k,64-k,r)$, weight-graded dimensions centered at $k$, and Jordan block sizes and multiplicities. For the Smith rank $r$ of $S$, the relative character is the coefficient of $z^k$ in $(1+zy)^r(1+z/y)^r(1+z)^{64-2r}$. A Jordan block of size $j+1$ has relative weights $j,j-2,\ldots,-j$, so its multiplicity is $c_j-c_{j+2}$ for character coefficients $c_j$. The character calculation works in every degree without listing exterior basis vectors. Individual operations are subject to the usual sparse arithmetic bounds. A geometric family and its Hodge filtration and polarization require additional data. [Exterior monodromy example](examples/hyperelliptic-ks-exterior-monodromy.asl) [Schreieder and Soldatenkov, *The Kuga–Satake construction under degeneration*](https://arxiv.org/abs/1710.02102)

#### Integral Nakajima states

`RVNakajima(B,n,i)` is the integral polynomial form generated by the unnormalized Nakajima creation operators in length $n$ and cohomological degree $i$. It uses a certified genus-two level-2 marking `B`. Tensoring with $\mathbb Q$ gives the even cohomology of $S^{[n]}$ under the Nakajima identification. The full integral cohomology lattice can require divided powers. The signed grades outside $n\ge0$ and even $0\le i\le4n$ contain only zero.

`rvNakajimaVac[B]` is the vacuum. Color 0 denotes the unit, colors 1 through 22 denote `rvK3LatticeBasis[B,c-1]`, and color 23 denotes the point class. The pairing $G$ extends the marked Gram matrix by pairing the unit and point with value 1. For length $m>0$, `rvNakajimaCreate[m,c](x)` multiplies by the corresponding variable. This adds $m$ to the length and $d_c+2(m-1)$ to the degree, with $d_c$ the surface degree of the color. `rvNakajimaAnn[m,c](x)` contracts the variables of length $m$ by $mG_{cj}$ times their multiplicities. Its grade shifts are $-m$ and $d_c-2m-2$. Use `add`, `neg`, `zero[RVNakajima(B,n,i)]`, and `rvNakajimaScale[k](x)` for the additive group and scalar action. `rvNakajimaMonomial[B]{(m,c,e),...}` supplies a monomial with multiplicity $e$ on each listed variable.

`rvNakajimaPair(x,y);` gives the Wick–Poincaré pairing on the common marking. It vanishes unless both states have length $n$ and their degrees sum to $4n$. A contraction at length $m$ has weight $(-1)^{m-1}mG_{cj}$. The computation includes every pairing of repeated factors and keeps powers in sparse multiplicity form. Coefficients and pairings use signed 128-bit integers, and Wick intermediates are checked at 512 bits with the usual work and storage limits. The example checks the Heisenberg bracket and gives pairings 1 and $-2$ for the unit and point at lengths one and two. [Nakajima example](examples/hyperelliptic-kummer-nakajima.asl)

`rvNakajimaBoundary(x)` cups with the effective half-boundary $\delta=-c_1(\mathcal O^{[n]})$, taking `RVNakajima(B,n,i)` to `RVNakajima(B,n,i+2)`. In the Fock convention above, Lehn's operator is half of the ordered join sum minus half of the ordered split sum. The canonical class term vanishes on a K3 surface. A join of lengths $a,b$ carries the weight $ab$ and uses the surface cup product on the colors. A split of length $a+b$ carries the weight $a+b$ and uses the surface diagonal from $G^{-1}$. The computation combines swapped terms and accounts for repeated copies to give integral coefficients. Sparse sums use wide intermediates with checked signed 128-bit output. Large splits are subject to the usual work and storage limits. The example checks Hilbert-square cup product and the Hilbert-cube top intersection $\int\delta^6=-960$. [Half-boundary example](examples/hyperelliptic-kummer-lehn-boundary.asl) · [Lehn, *Chern classes of tautological sheaves on Hilbert schemes of points on surfaces*](https://arxiv.org/abs/math/9803091)

`rvNakajimaCup(x,y)` multiplies classes from the same marking and Hilbert length $n$, adding their cohomological degrees. It computes Lehn–Sorger cup product on the integral Fock coordinate lattice, which spans the ring over $\mathbb Q$. A colored Fock monomial represents its full conjugator sum, so $q_1(1)^n$ acts by $n!$. Each pair of monomials uses a single sum over relative conjugators in $S_n$. On a permutation orbit the product uses the surface multiplication, the Euler class $24[\mathrm{pt}]$ to the excess power, and the surface diagonal to distribute colors to the output cycles. The resolution sign is $(-1)^{(|\sigma|+|\tau|-|\sigma\tau|)/2}$, where $|\sigma|=n-\#\operatorname{cycles}(\sigma)$. Coefficients are checked in the signed 128-bit range after wide accumulation. The factorial conjugator count, diagonal expansion, and sparse output storage are checked against the computation limits. The example compares the boundary multiplication and computes a nonzero excess contribution at length three. [Lehn–Sorger example](examples/hyperelliptic-kummer-lehn-sorger.asl) · [Lehn and Sorger, *The cup product of the Hilbert scheme for K3 surfaces*](https://arxiv.org/abs/math/0012166)

`rvHilbTautChern[n,j,k](a)` evaluates $k c_j(L^{[n]})$ in `RVNakajima(B,n,2j)`, with $a=c_1(L)$ in `RVK3Lattice(B)`. For an integral class of type $(1,1)$ this is the geometric tautological bundle. The same formula defines a universal polynomial for every marked lattice vector. Lehn's total Chern class is $\exp(\sum_{m\ge1}(-1)^{m-1}q_m(1+a)/m)|0\rangle$. In codimension $j$ there are exactly $n-j$ unit-colored cycles. A factor of length $m$ repeated $f$ times at one color contributes denominator $m^f f!$. `rvHilbTautChernDen[n,j](a);` gives the least positive integer clearing every coefficient, and the expression checks that `k` clears them. Classes above rank $n$ and classes scaled by zero are zero. Enumeration and coefficient storage are checked, with unsigned 128-bit denominators and signed 128-bit output coordinates.

`rvHilbTautSegre[n](a);` computes the integral of $s_{2n}(L^{[n]})$ as $2^n\binom{a^2/2-2n+2}{n}$ with the generalized integer binomial. It evaluates the square before checking the signed 128-bit result. The example computes the Chern classes of the trivial bundle at length three and verifies its top Segre number using the cup product. [Tautological Chern and Segre example](examples/hyperelliptic-kummer-tautological.asl) · [Lehn, *Chern classes of tautological sheaves on Hilbert schemes of points on surfaces*](https://arxiv.org/abs/math/9803091) · [Marian, Oprea and Pandharipande, *Segre classes and Hilbert schemes of points*](https://arxiv.org/abs/1507.00688)

`RVHilbCoh(B,n,i)` is the homogeneous rational cohomology of the marked Hilbert scheme in degree $i$. Use `rvHilbCohQ[k,d](x)` to give $kx/d$ from an integral Nakajima state, and `rvHilbCohScale[k,d]`, `add`, `neg`, and `zero` for the rational vector space operations. The fraction denominator `d` is positive. `rvHilbCohDen(x);` reports the least positive multiplier clearing the Fock coordinates, and `rvHilbCohFock[k](x)` returns the integral state for $kx$ when this multiplier divides `k`.

`rvHilbCohUnit[B,n]` is $q_1(1)^n/n!$ and `rvHilbCohPoint[B,n]` is $q_1(\mathrm{pt})^n$. They coincide at $n=0$. `rvHilbCohCup(x,y)` uses the Lehn–Sorger product at a common marking and length, adding cohomological degrees. `rvHilbCohIntegral(x);` gives the rational coefficient of the point class in top degree. Signed formal grades outside the cohomological range contain only zero. The common denominators and numerators are reduced before checking unsigned and signed 128-bit output ranges, respectively. The unit uses a single monomial and its checked factorial denominator.

`rvHilbCohChern[n,j](a)` gives the unscaled class $c_j(L^{[n]})$ on the marking of $a\in\mathrm{RVK3Lattice}$. `rvHilbCohCh[n,j](a)` gives $\operatorname{ch}_j(L^{[n]})$ by Newton's identities with cup powers. The degree-zero component is the rank $n$ times the unit. Chern characters can be nonzero above the rank, and both constructions vanish above the cohomological dimension. They use the same formal interpretation for marked classes as the scaled Chern operations. [Rational Hilbert cohomology example](examples/hyperelliptic-kummer-rational-hilbert.asl)

`rvHilbCohSchur[n,j]{lambda1,...}(a)` gives the codimension-$j$ character of $\mathbb S_\lambda L^{[n]}$. The rows of $\lambda$ are weakly decreasing, with trailing zeros allowed. One row gives a symmetric power, one column gives an exterior power, and the empty partition gives the trivial line. A column of height $n$ gives the determinant, and a diagram with more rows gives zero. The character uses Adams operations on the cohomological components, multiplication truncated at codimension `j`, and Jacobi–Trudi in the smaller of the row and column descriptions. The rank uses the Weyl dimension formula. Shape, degree, and exact coefficient bounds apply to the truncated computation.

`rvHilbExteriorEuler[n,k](a,b);` gives the Euler characteristic of $\bigwedge^k L^{[n]}\otimes D_A$, where $a=c_1(L)$ and $b=c_1(A)$ are marked lattice classes and $D_A$ descends $A^{\boxtimes n}$ along Hilbert–Chow. Scala's exterior cohomology formula gives $\binom{\chi(L\otimes A)}{k}\binom{\chi(A)+n-k-1}{n-k}$ with $\chi(A)=2+b^2/2$ on the K3 and generalized integer binomials. The holomorphic interpretations require integral classes of type $(1,1)$. The same formulas evaluate formal marked classes. [Schur and exterior Euler example](examples/hyperelliptic-kummer-schur-hilbert.asl) · [Scala, *Cohomology of the Hilbert scheme of points on a surface with values in representations of tautological bundles*, Theorem 5.2.1](https://arxiv.org/abs/0710.3072)

`rvHilbCohLefschetzInverse(h,x)` inverts the cup map $L_h^{D-d}:H^d\to H^{2D-d}$ for $D=2n$ and $d\le D$. Both classes have the same marking and Hilbert length, and $h$ has degree two. The middle-degree map is the identity. The other maps use the exact complementary Nakajima bases and rational cup product, with a scalar calculation for the top degree. A singular slice has no inverse. Hard Lefschetz guarantees the isomorphisms for a Kähler class. The interpreter also tests the actual cup maps for rational divisors outside the Kähler cone.

`rvHilbCohLower(h,x)` gives the normalized Lefschetz operator $\Lambda_h$ in degree two less than $x$. For a primitive class $p$ of degree $d$, it sends $L_h^j p$ to $j(D-d-j+1)L_h^{j-1}p$. Consequently $[L_h,\Lambda_h]$ acts by $i-D$ on $H^i$ and scaling $h$ by $a\ne0$ scales $\Lambda_h$ by $a^{-1}$. In degree two the operator uses the polarized Fujiki pairing. Higher degrees extract primitive strings through the exact required Lefschetz inverses. The graded matrices use sparse elimination with wide intermediates and checked output fractions. Enumeration, cup products, and matrix work are subject to their usual limits. [Lefschetz example](examples/hyperelliptic-kummer-hilbert-lefschetz.asl) · [Looijenga and Lunts, *A Lie algebra attached to a projective variety*](https://arxiv.org/abs/alg-geom/9604014)

`rvHilbCohBeauville(a,b);` evaluates the Beauville–Bogomolov pairing on two rational degree-two classes at a common marking and length. The surface summand uses the marked K3 Gram form. For $n\ge2$, the raw monomials $q_1(1)^{n-1}q_1(v)$ and $q_1(1)^{n-2}q_2(1)$ represent $(n-1)!\,\mu(v)$ and $2(n-2)!\,\delta$, respectively, where $q(\delta)=-2(n-1)$ and $\delta$ is orthogonal to the surface summand. At length one the pairing is the surface intersection form.

`rvHilbCohLLVLower(b,x)` gives the linear LLV lowering $F_b=q(b)\Lambda_b$ in degree two less than $x$. For isotropic $b$ it uses $F_{b+t}-F_t$, with both nonisotropic squares certified from the pairing. It accepts the same graded rational classes as normalized lowering. `rvHilbCohRotate(a,b,x)` acts in the degree of $x$ by $([L_a,F_b]-q(a,b)(i-2n)\mathrm{id})/2$ in degree $i$. On divisors this is $q(b,x)a-q(a,x)b$. The rotation is a derivation of the cup ring, and the lowering and rotation use the same exact cohomology and resource limits as the Lefschetz operations. At length three the example rotates $q_2(1)q_1(v)$ in degree four. In its two-cycle sector the products of divisors span the combinations $q_2(1)q_1(v)+2q_2(v)q_1(1)$. [LLV example](examples/hyperelliptic-kummer-hilbert-llv.asl) · [Looijenga and Lunts, *A Lie algebra attached to a projective variety*](https://arxiv.org/abs/alg-geom/9604014) · [Verbitsky, *Cohomology of compact hyperkaehler manifolds*](https://arxiv.org/abs/alg-geom/9501001)

`rvHilbLLVMultiplicity[n]{lambda1,...};` gives the multiplicity of a dominant integral highest weight in the full LLV module of K3$^{[n]}$. The coordinates are weakly decreasing in the orthogonal weight basis of type $B_{12}$, with up to twelve nonzero entries. The empty diagram gives the trivial representation. Length zero gives the vacuum. Length one has type $D_{12}$ and is outside these Hilbert-scheme $B_{12}$ queries.

`rvHilbLLVDimension{lambda1,...};` gives the Weyl dimension of an integral $B_{12}$ highest weight. `rvHilbLLVDecomposition[n];` lists the nonzero triples of highest weight, multiplicity, and dimension. It reports the total dimension and the even Betti numbers in degree order. The calculation extracts the dominant shifted coefficient of the Green–Kim–Laza–Robles character product with the Weyl denominator. Each one-variable oscillator coefficient enters a determinant, which is truncated by length and pruned by the minimum degree of each completion. The Betti vector comes from specializing each summand's Weyl character at $x_1=t$, $x_i=1$ for $i>1$, so exponent $a$ contributes in degree $2n+2a$. Polynomial work and storage are checked. Multiplicities use signed 128-bit output after 640-bit determinant arithmetic. Weyl dimensions and the specialization use checked 1024-bit arithmetic. [LLV decomposition example](examples/hyperelliptic-kummer-llv-decomposition.asl) · [Green, Kim, Laza, and Robles, *The LLV decomposition of hyper-Kähler cohomology*](https://arxiv.org/abs/1906.03432)

`rvVirasoro[k](x)` acts by $L_k$ on the same polynomial form, taking $(n,i)$ to $(n-k,i-2k)$, with $L_0x=nx$. For the oscillator convention $\alpha_{-m}=\operatorname{Create}_m$, $\alpha_m=\operatorname{Ann}_m$, and $\alpha_0=0$, it uses the normal-ordered quadratic expression $L_k=\tfrac12\sum_{\ell,a,b}(G^{-1})_{ab}:\alpha_{k-\ell,a}\alpha_{\ell,b}:$. The calculation verifies the integral inverse Gram matrix and uses its even diagonal for the equal-mode creation terms. The central charge is $\chi(S)=24$. The example verifies $L_2L_{-2}$ and $L_3L_{-3}$ on the vacuum and checks the central term on an occupied state. [Virasoro example](examples/hyperelliptic-kummer-virasoro.asl)

`rvNakajimaPrimaryDim[B,n,i];` gives the dimension of the Virasoro highest-weight space of conformal weight $n$ and Hilbert cohomological degree $i$. `rvNakajimaIsPrimary(x);` tests any Nakajima state, including zero, by evaluating $L_1$ and $L_2$. These suffice in characteristic zero because $[L_1,L_m]=(1-m)L_{m+1}$ for $m\ge2$.

`rvNakajimaPrimaryGen[B,n,i,j]` gives the zero-based $j$th generator as an `RVNakajima(B,n,i)` value. The columns are the monomials of the requested grade in lexicographic order of their sorted `(mode length, color, multiplicity)` tuples. Each free column of the joint kernel gives a basis vector with its denominators cleared, its coefficient gcd divided out, and its first nonzero coefficient positive. The resulting generators form a basis over $\mathbb Q$. Graded enumeration and exact sparse elimination are subject to the usual work and storage limits, and exact coefficients use the checked signed 128-bit range. The example compares a weight-two primary with a descendant of the vacuum. [Virasoro highest-weight example](examples/hyperelliptic-kummer-virasoro-primaries.asl)

#### Genus-one pencils through two nodes

`rvk3pencil P[B,i,j];` declares the pencil of planes through the two marked nodes $i,j$ on the Kummer quartic. Its divisor is $F=H-E_i-E_j$, with $F^2=0$, and `rvK3PencilDiv[P]` gives that class in `RVKummerDiv(B)`. The secant meets the quartic exactly twice at each node. The interpreter checks this by restricting the quartic to the secant. It uses the reduced row-echelon basis of the two linear forms vanishing on the secant, with the same coordinates for either ordering of the nodes.

`RVK3PencilBase(P)` is the pencil's rational projective line. `rvK3PencilParam[P](u,v)` normalizes a nonzero pair over the marking's field and `rvK3PencilMap[P](X)` evaluates the two forms on a resolved point. At either base node it evaluates on the exceptional tangent direction. A strict trope through exactly one base node is a section, and `rvK3PencilSection[P,a](p)` returns its unique point above $p$. The interpreter recovers the degree-one map on the trope from three projective parameters and inverts it without enumerating points. Each of the two tropes through both nodes gives the fiber class $2T_a+\sum_{k\in S_a\setminus\{i,j\}}E_k$.

The generic fiber is smooth in odd characteristic. Over the algebraic closure, translate the first node to zero and choose Weierstrass coordinates so the second is $[(0,0)-\infty]$. The planes then fix $x_1x_2=\lambda$ in the symmetric Kummer coordinates. The residual genus-one curve is birational to a double cover in $x_1+x_2$ branched at the four generically distinct values $\alpha+\lambda/\alpha$ for the other finite branch points $\alpha$. The eight rational trope sections make the pencil elliptic. `rvK3PencilFiber[P](p);` counts rational resolved points of a sampled fiber. `rvK3PencilCensus[P];` reports the $q+1$ base points, the total number of resolved points, and pairs of fiber size and multiplicity. Both reports use bounded enumeration on the resolved surface and include points with only quadratic-twist lifts. [Two-node pencil example](examples/hyperelliptic-k3-two-node-pencil.asl) [Kumar, *Elliptic fibrations on a generic Jacobian Kummer surface*](https://arxiv.org/abs/1105.1715)

##### Fiber quartics and Jacobians

`rvK3FiberQuartic[P](p);` restricts the Kummer equation to the selected plane. In canonical coordinates $(X:Y:Z)$ the larger node mask is $(1:0:0)$, the smaller is $(0:0:1)$, and the transverse vector is fixed by the two pivot columns of the pencil forms. The report gives the coefficients of $A_2Z^2+YB_2Z+Y^2C_2$ and of $F_4=B_2^2-4A_2C_2$, in descending powers of $X$. Over $Y\ne0$, completing the square using $2A_2(Z/Y)+B_2$ gives the double cover with binary quartic $F_4$.

For $F_4=aX^4+bX^3Y+cX^2Y^2+dXY^3+eY^4$, `rvK3FiberJacobian[P](p);` gives the coefficients of $y^2=x^3+cx^2+(bd-4ae)x+(ad^2+b^2e-4ace)$ and its discriminant. This formula and the integral Weierstrass invariants work in every odd characteristic. A zero discriminant is reported as singular, with the original plane coefficients retained in the quartic report. This includes the doubled conics where $F_4=0$. For smooth fibers, `rvK3FiberJ[P](p)` gives $c_4^3/\Delta$ and `rvK3FiberPoints[P,r](p);` counts points over the degree-$r$ extension. It sums the quadratic character over the marking's projective line and powers the genus-one Frobenius trace. [Fisher and Sills, *Local solubility and height bounds for coverings of elliptic curves*, Theorem 3.8](https://arxiv.org/abs/1103.4944) [Kumar, *Elliptic fibrations on a generic Jacobian Kummer surface*](https://arxiv.org/abs/1105.1715)

For a smooth fiber in characteristic at least five, `rvk3fiberjac E[P](p);` registers the Jacobian as a short Weierstrass curve for the ordinary `RVECPoint(E)` and `rvECPointCount[E];` operations. Characteristic three uses the long equation and the fiber point-count command. The example gives 12 points over $\mathbb F_{13}$, 192 over its quadratic extension, and the original doubled conic coefficients at a singular parameter. [Fiber Jacobian example](examples/hyperelliptic-k3-pencil-fiber-jacobian.asl)

##### The relative Jacobian and its discriminant

`rvK3PencilWeierstrass[P];` gives the relative model $y^2=x^3+a_2x^2+a_4x+a_6$ over the pencil base, with $a_2,a_4,a_6$ homogeneous of degrees $4,8,12$ in $[u:v]$. The report includes $c_4,c_6$ of degrees $8,12$ and the values at $[0:1]$. For a form of degree $d$, entry $i$ is the coefficient of $u^{d-i}v^i$. The interpreter expands the plane equation with a symbolic transverse vector, so the coefficient lists also apply over fields too small for evaluation and interpolation to recover these degrees. The construction checks $c_4^3-c_6^2=1728\Delta$ over the base ring. [Fisher and Sills, *Local solubility and height bounds for coverings of elliptic curves*](https://arxiv.org/abs/1103.4944) [Kumar, *Elliptic fibrations on a generic Jacobian Kummer surface*](https://arxiv.org/abs/1105.1715)

`rvK3PencilDiscriminant[P];` gives the homogeneous model discriminant of degree $24$, its order at infinity, and each finite closed place as its degree, order, and monic irreducible polynomial in $v/u$. These factor coefficients are in ascending order. The report counts each place with its degree times its order. Inseparable multiplicities are retained.

`rvK3FiberReduction[P](p);` gives the raw coefficient and invariant orders in the uniformizer $v/u-(v/u)(p)$ or $u/v$ at infinity. An identically zero coefficient has order `inf`. In characteristic at least five, the interpreter minimizes the short model by uniformizer scaling and gives its geometric Kodaira symbol from the minimal invariants, including the order of $j$. In characteristic three it reports the raw orders. The example exhibits good, multiplicative, and $I_0^*$ fibers. [Relative Jacobian example](examples/hyperelliptic-k3-pencil-relative-jacobian.asl) [Tate, *Algorithm for determining the type of a singular fiber in an elliptic pencil*](https://wstein.org/Tables/antwerp/tate/tate.pdf)

##### Euler factors and the function-field L-polynomial

`rvK3FiberEuler[P](p);` reports the local Euler polynomial in $U$ for a rational pencil parameter in characteristic at least five. A good fiber gives $1-a_pU+qU^2$. For multiplicative reduction, the quadratic character of the minimal $-c_6$ gives the split sign, hence $1-U$ or $1+U$. A split $I_n$ fixes all components. A nonsplit $I_n$ acts by $i\mapsto n-i$ on the component labels, with the zero section meeting component $0$. For additive reduction the Euler polynomial is $1$. At $I_0^*$ the first reported component cycle is the fixed central component. The other cycles come from the roots of the residual cubic with coefficients $A/\pi^2$ and $B/\pi^3$ in the minimal short model. Its number of rational roots plus one is the Tamagawa number. [Tate, *Algorithm for determining the type of a singular fiber in an elliptic pencil*](https://wstein.org/Tables/antwerp/tate/tate.pdf)

`rvK3PencilLFunction[P];` gives $L(T)$ in ascending coefficient order for the elliptic curve over the marking's $\mathbb F_q(\mathbb P^1)$. The report covers nonconstant $j$ in characteristic at least five with rational discriminant support and bad fibers of types $I_n$ and $I_0^*$. It takes the actual Frobenius polynomial on the Kummer surface's $H^2$, divides exactly by $(1-qT)^2$ for the fiber and section, and divides by $1-(qT)^d$ for each length-$d$ cycle on the nonidentity fiber components. The report checks that the resulting degree equals the conductor degree minus four. The order at $T=1/q$ is the arithmetic Mordell–Weil rank by the Tate theorem for Kummer surfaces and Shioda–Tate. [Shioda, *Some remarks on elliptic curves over function fields*](https://www.numdam.org/article/AST_1992__209__99_0.pdf) [Ulmer, *Elliptic curves over function fields*](https://math.arizona.edu/~ulmer/research/papers/2011.pdf) [Milne, *On a conjecture of Artin and Tate*](https://math.stanford.edu/~conrad/BSDseminar/refs/MilneAT.pdf)

The example over $\mathbb F_{13}$ reports conductor degree 10 and arithmetic Mordell–Weil rank two. The marks and heights below describe the explicit marked divisor lattice of the pencil. [Function-field L-polynomial example](examples/hyperelliptic-k3-pencil-l-function.asl)

##### Central special values and the BSD quotient

`rvK3PencilBSD[P];` computes the order $r$ and leading coefficient $L^*$ at $u=0$ in the exact expansion of $L((1-u)/q)$, where $u=1-qT$. It lists the Taylor coefficients in ascending powers of $u$. It verifies $L(T)=w(qT)^dL(1/(q^2T))$ for $d=\deg L$ and checks the root number $w=(-1)^r$. The report uses the same arithmetic scope as the L-polynomial command.

For this elliptic K3 over $\mathbb F_q(\mathbb P^1)$, the geometric height normalization gives the Birch–Swinnerton-Dyer quotient $|\Sha|\operatorname{Reg}/|E(\mathbb F_q(\mathbb P^1))_{\mathrm{tors}}|^2=qL^*/\prod_v c_v$. The factor $q$ comes from $g(\mathbb P^1)-1+\chi(\mathcal O_S)=1$. The interpreter computes the arithmetic trivial lattice determinant from the Gram matrices of the orbit sums of the nonidentity fiber components. It checks that its magnitude is $\prod_v c_v\prod_O |O|$. Independently, the Kummer Artin–Tate computation gives the signed product of the finite-field Néron–Severi discriminant and the Brauer order. Its magnitude is $qL^*\prod_O |O|$, whose ratio to the trivial determinant gives the same BSD quotient. [Shioda, *Some remarks on elliptic curves over function fields*, §2, Theorem 8](https://www.numdam.org/article/AST_1992__209__99_0.pdf) [Milne, *On a conjecture of Artin and Tate*](https://math.stanford.edu/~conrad/BSDseminar/refs/MilneAT.pdf)

The example over $\mathbb F_7$ has rank one, root number $-1$, leading coefficient $32/7$, and BSD quotient $1/16$. The independent marked height calculation gives height one for its chosen trope class. [Pencil BSD example](examples/hyperelliptic-k3-pencil-bsd.asl)

##### Marked Mordell–Weil coordinates and translations

`rvk3mordell M[P,z];` chooses the strict trope of mask $z$ as the zero section $O$ of a declared pencil `P`. The trope must meet exactly one base node. For the canonical node order $i<j$, put $S=E_i-E_j$. `rvK3MordellCoord[M](D);` gives the degree $d=D\cdot F$ and the integral free coordinate $\ell(D)=(D-dO)\cdot S/2$. The marked height `rvK3MordellHeight[M](D,E);` is $\ell(D)\ell(E)$. Indeed, $S$ is orthogonal to the fiber and its vertical components, and $S^2=-4$.

`rvK3MordellLattice[M];` verifies the trivial lattice $U\oplus D_4(-1)^2\oplus A_1(-1)^6$ and the marked divisor lattice by their Gram matrices. Their determinant magnitudes are 1024 and 64, respectively. Smith reduction of the inclusion gives the marked quotient $\mathbb Z\oplus(\mathbb Z/2)^2$. The four trope sections incident with the same base node as $O$ represent the four torsion classes. The other four have free coordinate $+1$ or $-1$ according to the chosen zero. For a very general complex Jacobian Kummer surface the marked lattice gives the full Mordell–Weil calculation for this pencil. Specializations can carry additional divisor classes and sections. [Schütt and Shioda, *Elliptic Surfaces*](https://arxiv.org/abs/0907.0298) [Kumar, *Elliptic fibrations on a generic Jacobian Kummer surface*](https://arxiv.org/abs/1105.1715)

`rvK3MordellShiftDiv[M,k](D)` applies the $k$th translation to a marked divisor by $D\mapsto D+k(D\cdot F)S-k(D\cdot S)F+2k^2(D\cdot F)F$. It uses wide checked intermediates and signed 128-bit input and output. A single translation shifts the free coordinate of each section by $-2$ and has height 4. On resolved points, `rvK3MordellShift[M,k](X)` uses the node projections first from $i$ and then from $j$, in reverse order for negative $k$. It preserves each pencil fiber, including at exceptional directions. For long iteration counts the point computation uses a bounded orbit walk. [Mordell–Weil pencil example](examples/hyperelliptic-k3-mordell-weil-pencil.asl)

##### Rational sections and full-split two-descent

`RVK3PencilSection(P)` is the group of points of the generic relative Jacobian over $K=\mathbb F_q(t)$ with $t=v/u$. Use `zero[RVK3PencilSection(P)]` for infinity, `neg` for the inverse, and `add` for the chord-and-tangent law on the long cubic. In particular, `add(x,x)` doubles a section. `rvK3PencilSectionPoint[P]{Xn}{Xd}{Yn}{Yd}` gives affine coordinates $(Xn/Xd,Yn/Yd)$. Each brace lists coefficients in ascending degree as bracketed residue coordinates, and an empty numerator gives zero. The interpreter cancels polynomial common factors, makes denominators monic, and checks the generic equation.

`rvK3PencilTorsion2[P];` reports the order of the rational two-torsion, including infinity, and the affine roots $e_i$ in lexicographic coefficient order. `rvK3PencilTorsion2At[P,i]` gives $(e_i,0)$ using zero-based indices. The monic cubic forces a rational root to be polynomial, and its leading degrees bound the root degree by four. The interpreter factors the constant term and tests bounded-degree monic divisors. For each divisor it finds the allowable scalar multiples using a polynomial gcd over the coefficient field. It peels off $x=0$ first when that is a root.

With three rational affine roots, `rvK3PencilDescent2(x)` gives the three square classes $[x-e_i]$ in `RVK3PencilKummer2(P)`. At $(e_i,0)$ the $i$th entry is $[(e_i-e_j)(e_i-e_k)]$, and infinity gives three trivial classes. Each coordinate records the leading constant's square character and the sorted monic irreducibles with odd order. The sum of their degrees modulo two gives the order parity at infinity. Use `add` to multiply square classes, `zero` for the trivial triple, and `neg` for inversion.

For each fixed $i$, evaluating the squared chord or tangent at $e_i$ gives the product of the three $x-e_i$ values at the intersections. The third intersection has the same $x$ coordinate as the sum. Taking the leading term at a torsion point gives the stated exceptional value. This proves the Kummer map is a homomorphism. When all classes are trivial, choose $r_i^2=x-e_i$ with $\prod_i r_i=-y$. With $s_1=\sum_i r_i$ and $s_2=\sum_{i<j}r_ir_j$, the section $(x+s_2,-y-s_1s_2)$ doubles to $(x,y)$. Thus `rvK3PencilDivisible2(x)` tests membership in $2E(K)$. Selmer computations additionally impose local solubility on classes outside the given rational-section image. [Silverman, *The Arithmetic of Elliptic Curves*, X.1](https://www.math.ens.psl.eu/~obenoist/refs/Silverman.pdf) [Zarhin, *The divisibility by 2 of rational points on elliptic curves*](https://arxiv.org/abs/1702.02255)

The example over $\mathbb F_{13}$ uses the section $(t+6t^3,7t^2+6t^3+3t^4)$, whose good specializations at $t=1$ and $t=5$ have orders three and four. These distinct orders prime to the characteristic certify infinite order. Its Kummer class is nontrivial and its double has trivial class. [Pencil two-descent example](examples/hyperelliptic-k3-pencil-two-descent.asl)

##### Global two-Selmer space

`rvK3PencilSelmer[P];` computes the global $\mathbb F_2$-Selmer space for the full-split generic Jacobian. The two independent square classes are supported on the constant nonsquare and the monic factors of the differences of the three torsion roots. At other finite places the Kummer image consists of the unramified classes. At each finite place of the support and at infinity, the calculation checks exact local solubility by expanding residue prefixes that contain at least two roots. Once a prefix contains exactly one root, the product-square condition forces its torsion class. The infinity chart uses $s=1/t$ and the integral coordinates $s^4x,s^6y$.

The local codes use two bits per square class, with the valuation parity first and the angular nonsquare bit second. The first two torsion-root coordinates supply the four bits and their sum gives the third square class. Every reported local image has dimension two. The computation currently requires rational places in the finite support. `rvK3PencilSelmerBasis[P,i]` returns the zero-based basis as `RVK3PencilKummer2(P)` values. `rvK3PencilSelmerContains(c)` tests a typed class against the local conditions, including classes obtained with `rvK3PencilDescent2`.

When the pencil's L-function calculation certifies the arithmetic Mordell–Weil rank $r$, the Selmer exact sequence gives $\dim\Sha(E/K)[2]=\dim\operatorname{Sel}_2(E/K)-r-2$. In characteristic three the report gives the Selmer dimension and local images. The example over $\mathbb F_{23}$ has Selmer dimension six and rank two, so $\Sha(E/K)[2]\cong (\mathbb Z/2)^2$. [Pencil Selmer example](examples/hyperelliptic-k3-pencil-selmer.asl) [Silverman, *The Arithmetic of Elliptic Curves*, X.1 and X.4](https://www.math.ens.psl.eu/~obenoist/refs/Silverman.pdf)

#### Geometric Kummer automorphism words

`rvk3autom W[B]{translate(a),switch(b),project(c),...};` registers a word on the marked Kummer resolution. Each index is a node mask between 0 and 15, and the word acts from left to right. Translation fixes $H$ and takes $E_v$ to $E_{v+a}$. The Gauss switch takes $E_v$ to $T_v$ and $H$ to $3H-\sum_v E_v$. `switch(b)` applies this switch and then translates by $b$. Projection from node $c$ takes $(H,E_c)$ to $(3H-4E_c,2H-3E_c)$ and fixes the other exceptional curves. Each generator squares to the identity. `rvK3AutoDiv[W](D)` computes the word's action on `RVKummerDiv(B)`.

`rvK3AutoLattice[W](x)` applies the same word to `RVK3Lattice(B)`, including its glued classes. On the rational orthogonal complement of the hyperplane and exceptional classes, translations act by $+1$ and each switch or projection acts by $-1$. The calculation uses signed 128-bit doubled coordinates with checked 416-bit intermediates and verifies the integral gluing condition on the result.

`rvK3AutoNikulin[W];` verifies that the integral rank-22 action squares to the identity and reports the invariant lattice of the involution. Its rank is $r$, its discriminant group is $(\mathbb Z/2)^a$, and $\delta=0$ means that the discriminant quadratic form is integer-valued. The reported fixed locus is geometric on the complex marked K3 model. The identity fixes the whole surface, a nontrivial symplectic involution fixes eight isolated points, and the antisymplectic cases use the Nikulin classification. The example includes a free Enriques switch, an odd Gauss switch, and a node projection with six rational fixed curves. [K3 involution example](examples/hyperelliptic-k3-involution-nikulin.asl) [Nikulin's classification, as recalled by Artebani, Sarti and Taki, §3](https://arxiv.org/abs/0903.3481)

`rvK3Entropy[W];` gives the exact monic characteristic polynomial in descending coefficient order and the multiplicities of its cyclotomic factors by order. The action on the rank-17 marked lattice suffices for the first dynamical degree because its orthogonal complement in Néron–Severi is negative definite. A reciprocal quadratic or Salem factor gives the degree as its unique real root greater than one. With a wholly cyclotomic spectrum the degree is one. The example gives $9+4\sqrt5$ for three distinct projections, a degree-four Salem factor for four, and a nonidentity word of dynamical degree one for two. [Kummer dynamics example](examples/hyperelliptic-kummer-lattice-dynamics.asl)

`rvK3AutoPoint[W](X)` applies the same word to a resolved point on the marking's named curve, including exceptional directions. `rvK3AutoIter[W,k](X)` accepts a signed 128-bit iteration count. Negative counts use the reversed word, and long iterations can reduce the count modulo the first return to `X`. `rvK3AutoPeriod[W](X);` gives that least positive return time without enumerating the surface. `rvK3AutoCycles[W];` enumerates the resolved points over the marking's finite field and reports the total, fixed points, and cycle multiplicities by length. These finite permutation statistics depend on the field. The example exhibits one word of geometric dynamical degree $9+4\sqrt5$ with finite cycles of lengths 1, 2, 3 and 5 over $\mathbb F_5$. [Resolved Kummer automorphism example](examples/hyperelliptic-kummer-automorphism-points.asl)

`rvK3FrobFixed[W,Base,r];` counts the geometric fixed points of $\Psi^r$ for $\Psi=W\circ F_q$, where $F_q$ raises coordinates to the $q$th power for the field of `Base`. Use a positive iteration count $r$. The equation on the marking must be the base extension of `Base`. The word can use generators whose node masks move under $F_q$. The calculation composes the word's marked action $G$ with the certified Frobenius permutation $P$ before taking powers of $A=GP$. The vanishing differential of $\Psi^r$ makes the fixed points isolated and reduced.

`rvK3FrobZeta[W,Base];` gives the dynamical zeta function $\exp(\sum_{r\ge1}\#\operatorname{Fix}(\Psi^r)T^r/r)$ in factored form. It lists $D(U)=\det(1-UA)$ and the abelian cross quartic $h(T)$ in ascending coefficient order. Translations act trivially on primitive cohomology. Each switch and projection contributes a minus sign, so the primitive factors are $(1-\epsilon qT)h(\epsilon T)$ for the resulting sign $\epsilon$. Counts and coefficients use checked integer arithmetic. The example has 12 fixed points for Frobenius followed by a free Enriques switch over $\mathbb F_3$, and 248 for its square. [Semilinear Kummer dynamics example](examples/hyperelliptic-k3-semilinear-dynamics.asl)

`rvK3FrobHilbFixed[W,C,n,r];` counts the geometric fixed points of $\Psi_n^r$ on the length-$n$ Hilbert scheme of the resolved Kummer K3, where $\Psi_n=W^{[n]}\circ F_q$ and $q$ is the field size of `C`. `rvK3FrobHilbTrace[W,C,n,r,i];` gives the trace of the same action in degree $i$. Length zero is a point, $r\ge1$, and traces are zero in odd degrees and outside $0\le i\le4n$. The Frobenius here acts on the Hilbert scheme over the base field, and the word acts by functoriality on subschemes. This defines the composition even when the word's marked generators move under Frobenius. Its differential vanishes, so the fixed points counted by the Lefschetz formula are isolated and reduced.

The Nakajima decomposition places a length-$m$ copy of surface $H^b$ in degree $b+2(m-1)$ with a Tate factor. At power $s$ of $\Psi_n^r$ the recurrence uses the surface traces $1$, $\operatorname{Tr}(\Psi^{rs}|H^2)$ and $q^{2rs}$ with mode factor $q^{rs(m-1)}$. The middle surface traces come from a single marked matrix $A=GP$ and the primitive sign formula above. The normalized middle surface action preserves the intersection pairing, so Poincaré duality gives the higher Hilbert traces from lower degrees with the corresponding Tate factor. If $C_r=\#\operatorname{Fix}(\Psi^r)$, the Hilbert square has $(C_r^2+C_{2r})/2+q^r C_r$ fixed points. The example gives 232 fixed points and even traces $1,5,100,45,81$ for the Enriques switch after Frobenius over $\mathbb F_3$. Counts use unsigned 128-bit results, traces use signed 128-bit results, and the Hilbert coefficients use checked 512-bit intermediates with work and storage limits. [Nakajima](https://arxiv.org/abs/alg-geom/9507012) · [Göttsche–Soergel, §5](https://doi.org/10.1007/BF01445104) · [Semilinear Hilbert-scheme example](examples/hyperelliptic-k3-semilinear-hilbert.asl)

#### Rational cohomology of the Hilbert square

`RVHilbSquare(B,d)` is the rational cohomology in degree $d$ of the Hilbert square $S^{[2]}$ of the marked Kummer K3 surface. Its dimensions in degrees $0,2,4,6,8$ are $1,23,276,23,1$. Degrees are signed. Every other degree has only the formal zero. `rvHilbSquareBasis[B,d,j]` selects coordinate $j$ starting at zero. Degrees two and six use the actual basis of the marked K3 lattice followed by $\delta$, half the exceptional divisor, with Beauville square $-2$. Degree six uses the Beauville form to identify this basis with its Poincaré dual space. Degree four uses the symmetric monomials $e_i e_j$ in lexicographic order with $i\le j$. `rvHilbSquareDiv[k](c)` places a class `c` of `RVK3Lattice(B)` in degree two and adds $k\delta$.

Use `rvHilbSquareCup(x,y)` to multiply classes on the same marking, `rvHilbSquareScale[n,d](x)` for a rational multiple, and `zero`, `add`, and `neg` in each degree. The cup product uses the polarized Beauville–Fujiki relation $\int xyzw=q(x,y)q(z,w)+q(x,z)q(y,w)+q(x,w)q(y,z)$ for degree-two classes. In degree four, `rvHilbSquareC2[B]` gives $(6/5)q^\vee$ using the inverse Beauville matrix. In degree eight, `rvHilbSquareC4[B]` gives $324[\mathrm{pt}]$ and `rvHilbSquareIntegral(x);` evaluates a class on the fundamental cycle. The example computes $\int\delta^4=12$, $\int c_2^2=828$, and $\int (3c_2^2-c_4)/720=3$. [Hilbert-square cohomology example](examples/hyperelliptic-hilbert-square-ring.asl) [Taelman, *Derived equivalences of hyperkähler varieties*, §§7.2–7.3](https://msp.org/gt/2023/27-7/gt-v27-n7-p02-s.pdf)

##### Lefschetz and LLV operators

For $a\in H^2$, multiplication by $a$ gives the raising operator $L_a$. Use `rvHilbSquareLower(a,x)` for the linear lowering operator $F_a$, including when $q(a,a)=0$. It maps $z\in H^2$ to $4q(a,z)\cdot 1$ and sends $H^0$ to the formal zero component in degree $-2$. For $q(a,a)\ne0$, `rvHilbSquareLambda(a,x)` gives $F_a(x)/q(a,a)$ and `rvHilbSquarePrimitive(a,x)` projects to the primitive summand in degrees $0,2,4$.

Use `rvHilbSquareRotate(a,b,x)` for the degree-preserving operator $R_{ab}$. It sends $z$ to $q(b,z)a-q(a,z)b$ on $H^2$ and $H^6$, acts on $H^4$ by the product rule, and kills the unit and point class. All inputs use the same marking. On degree $d$, $[L_a,F_b]=2R_{ab}+q(a,b)(d-4)$. These operators give the Looijenga–Lunts–Verbitsky action on the rational Hilbert-square ring. Its Lie algebra is $\mathfrak{so}(H^2\oplus U)$ for a hyperbolic plane $U$, with underlying quadratic space of dimension $25$ and real signature $(4,21)$. [LLV example](examples/hyperelliptic-hilbert-square-llv.asl) [Green, Kim, Laza and Robles, *The LLV decomposition of hyper-Kähler cohomology*](https://arxiv.org/abs/1906.03432)

##### The weight-two Nakajima correspondence

`rvNakajimaHilbSquare(x)` sends a length-two `RVNakajima(B,2,d)` state to `RVHilbSquare(B,d)`. `rvHilbSquareNakajima[k](x)` applies the inverse to the multiple $kx$, returning an integral Fock state when all the coefficients are integers. `rvHilbSquareNakajimaDen(x);` gives the least positive integer scale, or 1 for zero. Both maps preserve the marking and cohomological degree. The unit has denominator 2 because the unnormalized state $q_1(1)^2$ maps to $2\cdot 1$. [Weight-two Nakajima example](examples/hyperelliptic-hilbert-square-fock.asl)

Write $a,b$ for marked degree-two surface classes, $p$ for the point class on the surface, and $P=c_2/24+\delta^2/8$ for its incidence class on the Hilbert square. The map sends $q_1(1)q_1(a)$ to $\mu(a)$ and $q_1(1)q_1(p)$ to $P$. It takes $q_1(a)q_1(b)$ to $\mu(a)\mu(b)-q(a,b)P$, $q_1(a)q_1(p)$ to $a$ in the dual-vector model of $H^6$, and $q_1(p)^2$ to the point class in $H^8$. The double-point operators $q_2(1),q_2(a),q_2(p)$ give $2\delta$, $\delta\mu(a)$, and $\delta/2$ in $H^6$, respectively.

On the blowup of $S\times S$ along the diagonal, the double cover of $S^{[2]}$ pulls $\delta$ back to the exceptional divisor $E$. Pullbacks of the two-support states are symmetric tensors, and the double-point states pull back to twice the exceptional incidence classes. The diagonal relation $E^2=-[\Delta]$ uses $c_1(T_S)=0$. Its Künneth expansion yields $20 P=G^\vee+2\delta^2$, where $G^\vee$ is the inverse surface pairing in symmetric-degree-two coordinates. In particular, if $P_F=q_1(1)q_1(p)$ and $Q_{ij}=q_1(e_i)q_1(e_j)$, the inverse sends $e_i e_j$ to $Q_{ij}+G_{ij}P_F$ and $\delta^2$ to $-P_F-\frac12\sum_{i,j}G^{ij}Q_{ij}$. The example checks $\int P^2=1$, $\int P\delta^2=-1$, and the Wick pairing $-2$ for the unit and point double-mode states. [Lehn and Sorger, *The cup product of the Hilbert scheme for K3 surfaces*](https://arxiv.org/abs/math/0012166)

#### Richelot isogeny targets

`rvrichelot D[C]{G1}{G2}{G3};` registers a genus-two Richelot target over the same named odd finite field as `C`. The source curve has precision one and a smooth quintic or sextic equation $y^2=f(x)$. Each brace lists the coefficients of one factor as bracketed signed field coordinates in ascending degree, as in `rvjac`. The product must equal $f$ exactly. Each factor is quadratic, with one linear factor allowed for a quintic. The linear factor accounts for the branch point at infinity.

Let $\Delta$ be the determinant of the three ascending coefficient rows, with a zero quadratic coefficient for a linear factor. For cyclic $(i,j,k)$ put $B_i=G'_jG_k-G_jG'_k$. The target equation is $Y^2=B_1B_2B_3/\Delta$. The declaration rejects a zero determinant, for which the quotient is a product of elliptic curves over the algebraic closure. `rvRichelotData[D];` shows the source name, determinant, target polynomial, and stored splitting $(B_1/\Delta,B_2,B_3)$. `rvrichelotdual E[D];` applies the same construction to that splitting and gives $Y^2=4f(x)$, which identifies with the original curve by $(x,y)\mapsto(x,2y)$. The dual declaration also stores a splitting, so it can be iterated.

The example starts with a quintic over $\mathbb F_7$ and obtains a sextic target. The two curves have Frobenius trace $-2$ and Weil polynomial $T^4+2T^3-2T^2+14T+49$. The second target has the same Weil polynomial. [Richelot example](examples/hyperelliptic-richelot.asl)

`rvRichelotPoint[N](x,y)` gives the Richelot image of the class $[(x,y)-\infty]$ for an affine point on the monic quintic source. The target `N` can be the Richelot target itself if it is a monic quintic, or a branch chart registered directly on that target. The input coordinates belong to the source field and the output has type `RVJac(N)`.

The correspondence cyclically orders the splitting with the linear $G_3$ last and uses the unscaled brackets $B_i$. For $y\ne0$, the projective fiber has equation $L(z)=G_1(x)B_1(z)+G_2(x)B_2(z)=0$ and ordinate $w=G_1(x)B_1(z)(x-z)/y$. The fiber above the source infinity is linearly equivalent to twice the target base point. In a branch chart $t=d/(z-a)$, the interpreter transports the fiber using $t^2L(a+d/t)$ and $t^3w(a+d/t)/d$, then makes the resulting Mumford pair. This includes points at the old infinity and drops support at the new infinity. For $y=0$, a root of $G_1$ maps to the branch pair of $B_2$, a root of $G_2$ to the pair of $B_1$, and a root of $G_3$ to the identity. The example compares a split image and its dual with explicit Jacobian addition over $\mathbb F_9$. It also computes a repeated point at the old infinity on a charted sextic over $\mathbb F_7$. [Richelot point example](examples/hyperelliptic-richelot-point.asl)

`rvRichelotPush[N](D)` applies the same correspondence to any class `D` of the named source Jacobian and returns a class in `RVJac(N)`. For a squarefree quadratic Mumford support $u(X)=(X-h)^2-k$, it computes the fiber of $(h+\omega,v(h+\omega))$ in the quadratic algebra $\omega^2=k$ and takes its Jacobian trace. A repeated point gives twice its point image. Branch pairs and fibers meeting the chart infinity use their specialized images. The example verifies a nonsplit push by adjoining its two conjugate source points, then compares a repeated-point push with doubling. [Richelot push example](examples/hyperelliptic-richelot-push.asl)

`rvKummerRichelot[N](K)` applies the induced Richelot map to any rational point of the source Kummer quartic. The result belongs to `RVKummer(N)`, where `N` is the monic quintic target or a branch chart on the target. For points that lift only on a quadratic twist, the interpreter rescales the Richelot splitting and performs the divisor push over the named field. A chart on the target changes the twist parameter from $k$ to $k^3$ before the final Kummer projection. [Kummer Richelot example](examples/hyperelliptic-kummer-richelot.asl)

`rvJacRichelotKernels[B,C];` uses a certified genus-two level-2 basis `B` on `L` and one step of Frobenius for the base curve `C`. It reports the orbit lengths on all 15 maximal isotropic subspaces of $J(L)[2]$. Each stable plane appears as the sorted triple of its nonzero binary masks, with the bits in basis order $(e_1,e_2,f_1,f_2)$ and $e_1$ in the least significant bit. The report gives the base-field kernel size, including the identity, and indicates whether the Richelot determinant vanishes.

`rvJacRichelotKernel[B,C,i,j]` returns the actual class in `RVJac(L)` at mask index `j` in stable plane `i`. The planes are listed lexicographically and both indices start at zero. `rvrichelotfromlevel DL[B,C,i];` uses their three Mumford support polynomials to register a Richelot target over the field of `L`. A zero determinant leaves the kernel available and prevents this single-curve construction. In the example, one transposition of the branch points fixes three planes, and the smooth target kills all three nonzero classes of its kernel. [Richelot Hecke example](examples/hyperelliptic-richelot-hecke.asl)

`rvrichelotdescend D[C,DL];` descends the actual target equation of `DL` to the field of `C`. The source of `DL` must be the coefficient extension of the quintic or sextic `C`, and Frobenius over the base field must permute the three branch pairs up to scaling of their factor polynomials. `rvRichelotDescent[D];` gives the relative degree, the permutation in the stored factor order starting at zero, and the rational kernel order, which is one plus the number of fixed pairs. When both base equations are monic quintics, `rvRichelotDescPush[D](d)` extends a class on `C`, applies the stored correspondence, and descends the image to `RVJac(D)`. The example over $\mathbb F_3$ uses a transposition over $\mathbb F_9$, giving rational kernel order two and image order ten. [Richelot descent example](examples/hyperelliptic-richelot-descent.asl)

#### Richelot matrices on torsion

`rvJacRichelotMatrix[B,H];` computes the matrix $M$ of the Richelot push on the certified torsion bases `B` and `H` of the same prime level $\ell$. The curve of `H` determines the direct or charted Richelot target, and the curve of `B` must be its named source. The output writes matrices by rows with generator images in columns. Let $c=\log_{\zeta_H}\zeta_B$ and $J=\begin{pmatrix}0&I\\-I&0\end{pmatrix}$. The interpreter reconstructs each pushed class, verifies $M^tJM=2cJ$, and computes the polarized adjoint $A=c^{-1}J^{-1}M^tJ$. It checks $AM=MA=2I$ and that $M$ has rank two at level two and rank four at odd levels. The example uses two primitive-root choices at level three and exhibits the level-two kernel. [Richelot torsion matrix example](examples/hyperelliptic-richelot-matrix.asl)

`rvKummerRichelotDesc[D](K)` uses a Richelot descent registered on monic quintics. It extends the source point to the field of the stored correspondence, applies the Kummer map there, and descends the normalized image. The example starts with a twist-only point over $\mathbb F_3$ whose lift changes sign over $\mathbb F_9$. Its image is also twist-only over $\mathbb F_3$. [Kummer base change and descended Richelot example](examples/hyperelliptic-kummer-base-descent.asl)

#### Bielliptic quotients in a split chart

`rvbielliptic C[Plus,Minus](a,b,c,d);` registers two degree-two elliptic quotients of a smooth genus-two quintic or sextic. The four coefficients belong to the curve's named precision-one odd field and satisfy $ad-bc\ne0$. The change $x=(at+b)/(ct+d)$ and $Y=(ct+d)^3y$ must give $Y^2=H(t^2)$ for $H(u)=h_0+h_1u+h_2u^2+h_3u^3$ with $h_0h_3\ne0$. The three branch pairs exchanged by $t\mapsto-t$ give a split Richelot quotient with target $\mathrm{Plus}\times\mathrm{Minus}$. The declaration requires the explicit chart over the base field.

The quotient by $(t,Y)\mapsto(-t,Y)$ uses $(t^2,Y)$ and the quotient by $(t,Y)\mapsto(-t,-Y)$ uses $(t^{-2},Y/t^3)$. The interpreter scales both coordinates by the leading coefficient of the corresponding cubic, $h_3$ for $H(u)$ and $h_0$ for $u^3H(1/u)$. This gives monic cubic models that work with the existing `RVJac` and `rvHWW` operations.

`rvBiellipticPoint[E](x,y)` checks the affine source point $(x,y)$ and returns its image in `RVJac(E)` for either factor. A pole of the selected quotient maps to the elliptic point at infinity, represented by the zero class. The example computes both images and the Frobenius polynomials over $\mathbb F_7$. [Bielliptic quotient example](examples/hyperelliptic-bielliptic.asl)

`rvBiellipticPush[E](D)` pushes a class from `RVJac(S)` to `RVJac(E)`, and `rvBiellipticPull[S](Q)` pulls a class from `RVJac(E)` back to `RVJac(S)`. Here `E` is one of the two declared factors and `S` is the original source when it is a monic quintic, or a branch chart registered directly on that source with `rvjacchart`.

The pushforward traces squarefree quadratic support over its quadratic algebra and doubles the image of a repeated point. It sends support over a pole to elliptic infinity and subtracts the image of the source base point with the original multiplicity. In working coordinates $X=h(N/M)^2$ and $W=A y/M^3$, the pullback of $(X,W)$ has finite Mumford support $hN^2-XM^2$ and ordinate $WM^3/A$. Normalizing the support accounts for a fiber meeting the source infinity.

The example checks on all 64 source classes and all eight classes of each factor that push after pull doubles on the same factor and vanishes on the other. The sum of the two pull-after-push maps doubles on the source. [Split Picard isogeny example](examples/hyperelliptic-bielliptic-picard.asl)

#### Gluing elliptic curves along two-torsion

`rvellglue C[E,F](a,b,c,d);` glues the monic cubics $E:y^2=f(x)$ and $F:y^2=g(x)$ over the same named odd finite field. The Möbius map $(ax+b)/(cx+d)$ must have nonzero determinant and $c\ne0$. Set $A=-d/c$, $B=a/c$, $m=(bc-ad)/c^2$, and $k=g(B)$. The declaration checks $f(A)k\ne0$ and $u^3g(B+m/u)=k f(A+u)$ as a polynomial identity. This identifies the three finite branch points, including nonsplit orbits, and gives the genus-two sextic $C:Y^2=f(A+kt^2)$.

`rvEllGlueLeft[C](t,Y)` gives $(A+kt^2,Y)$ in `RVJac(E)`. `rvEllGlueRight[C](t,Y)` gives $(B+m/(kt^2),Y/(kt^3))$ in `RVJac(F)`, with $t=0$ sent to infinity. Both expressions check that the input lies on `C`. The example over $\mathbb F_5$ glues two irreducible cubics by $x\mapsto1/x$ and checks the maps and Frobenius polynomials of all three curves. [Elliptic gluing example](examples/hyperelliptic-elliptic-gluing.asl)

`rvEllGlueMap[S](P,Q)` adds the pullbacks of classes on `E` and `F` to give a class in `RVJac(S)`. Declare `S` directly on the glued sextic with `rvjacchart S[C](r);` for a rational branch point `r`. `rvEllGlueDual[S](D)` returns the two norms as a pair of type `RVJac(E) * RVJac(F)`. The factor names come from the original gluing declaration in the same order, including when both factors use one name. Both maps use the branch chart's coordinates and its point at infinity as the source base point for the norms.

Both compositions are multiplication by two on their respective groups. The kernel of the map from the product is the graph of the declared identification on two-torsion. The example checks the compositions on every class of the product and every class of the glued Jacobian, and verifies the four matching pairs in the kernel. [Gluing and dual isogeny example](examples/hyperelliptic-elliptic-gluing-picard.asl)

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

### Integral Satake inversion and spherical IC traces

`rvSatakeInverse(f)` recovers a spherical Hecke element when $f$ lies in the normalized Satake image over $\mathbb Z[q]$. The interpreter strips the highest dominant monomial in the Hall–Littlewood expansion on each total weight. If its coweight is $\lambda$ and its normalization exponent is $h(\lambda)$, then its leading coefficient divided by $v^{h(\lambda)}$ must have even nonnegative powers of $v$. Signed integer coefficients are allowed. The calculation converts these powers to powers of $q$ and subtracts that spherical term's image before continuing.

`rvSatSchur[[lambda0,...]]` gives the Schur character for a weakly increasing signed coweight, including the determinant shift. Its interlacing branching coefficients are all one. The spherical IC trace element `rvSphIC[[lambda0,...]]` is the inverse Satake image of $v^{h(\lambda)}s_\lambda$. Its standard-basis coefficients can be inspected with `rvSphCoeff`.

For rank three, the adjoint trace element $\mathrm{IC}_{(0,1,2)}$ is $S_{(0,1,2)}+(q+1)S_{(1,1,1)}$. Its square has the following character decomposition.
$$
\mathrm{IC}_{(0,2,4)}
+q(\mathrm{IC}_{(1,1,4)}+\mathrm{IC}_{(0,3,3)})
+2q^2\mathrm{IC}_{(1,2,3)}
+q^4\mathrm{IC}_{(2,2,2)}.
$$
The example also writes the rank-two square $S_{(0,1)}^2$ as $\mathrm{IC}_{(0,2)}+q\mathrm{IC}_{(1,1)}$. [Satake inverse and IC example](examples/rv-satake-ic.asl)

### Littlewood–Richardson hives and Horn inequalities

`RVHive([lambda],[mu],[nu])` is the finite discrete type of integral hives with three prescribed dominant weights. Enter signed weights in weakly increasing order, as for `rvSphBasis`. A hive is a triangular array $H(i,j)$ for $i,j\geq0$ and $i+j\leq n$, with $H(0,0)=0$. The directed edges from $(0,0)$ to $(n,0)$, from $(n,0)$ to $(0,n)$, and from $(0,0)$ to $(0,n)$ carry the successive differences of $\lambda$, $\mu$, and $\nu$, each read in descending order. In each unit rhombus the sum at the obtuse vertices is at least the sum at the acute vertices.

`rvHive[[lambda],[mu],[nu]]{coordinates}` supplies only the interior vertices in order of increasing $j$, then increasing $i$. The literal checks every rhombus. `rvHiveWitness[[lambda],[mu],[nu]]` constructs the lexicographically first integral point, and `rvHiveAt[i,j](h);` reads any vertex. An empty slice is an empty type. For small slices `card \h : RVHive([lambda],[mu],[nu]) => S.hi;` counts all points through normal finite-function enumeration. `rvHiveCount[[lambda],[mu],[nu]];` counts directly by interlacing successive rows and checking the third rhombus direction, including counts above the enumeration limit when the search fits its work budget.

`rvHorn[[lambda],[mu],[nu]];` reports admissibility, a total mismatch, or a violated Horn inequality. For a triple of equal-sized subsets $I,J,K$ of the descending rank-$n$ weights, it reports $\sum_I\lambda+\sum_J\mu<\sum_K\nu$ only when $c_{\tau(I),\tau(J)}^{\tau(K)}>0$ in smaller rank. Here $\tau(I)=(i_r-r,\ldots,i_1-1)$ for $I=\{i_1<\cdots<i_r\}$. The interpreter checks positivity by finding a smaller integral hive. The integral type-A hive and Horn theorems identify these admissible boundaries with the nonzero Littlewood–Richardson multiplicities.

For $\lambda=\mu=(0,1,2)$ and $\nu=(1,2,3)$, the two hives give the multiplicity behind the coefficient $2q^2$ in the spherical IC tensor-square calculation. Scaling all three weights by $k$ gives $k+1$ hives. Adding constants $a$ and $b$ to the two input weights and $a+b$ to the output shifts each vertex by $ai+(a+b)j$ and preserves the count. [Hive example](examples/rv-hive.asl)

### Level-dependent Verlinde fusion

`RVFusion(n,k)` is the integral level-$k$ $\mathrm{SU}(n)$ WZW fusion ring. Use `rvFusionBasis[k,[lambda0,...]]` for an irreducible basis element and `rvFusionOne[n,k]` for the unit. The interpreter subtracts the smallest coordinate from the weakly increasing input weight. The resulting largest coordinate must be at most $k$. This identifies determinant translates and leaves just the unit basis element at level zero.

`rvFusionMul(a,b)` multiplies elements at the same rank and level. It counts classical tensor multiplicities in integral hive slices and applies the Kac–Walton affine Weyl fold to each shifted highest weight. Weights on alcove walls contribute zero. The remaining folds contribute with the permutation sign, and basis products have nonnegative coefficients after cancellation. This level-dependent multiplication is distinct from the classical character tensor product computed by the spherical IC basis.

Elements support `zero`, `add`, `neg`, equality, signed integer multiplication by `rvFusionScale[c](a)`, and the contragredient involution `rvFusionDual(a)`. Use `rvFusionCoeff[[lambda0,...]](a);` to read an integer coefficient. The example squares the $\mathrm{SU}(3)$ adjoint class at levels two, three, and four. The level-two result is $1+8$. Level three gives $1+2\cdot 8+10+\overline{10}$, and level four adds $27$. [Verlinde fusion example](examples/rv-verlinde.asl)

### Cyclotomic Verlinde characters and conformal blocks

`RVCyclo(m)` is the integral cyclotomic ring $\mathbb Z[z]/\Phi_m(z)$. `rvCyclo[m,e]` is $z^e$ for a signed integer $e$. The ring supports `zero`, `add`, `neg`, equality, `rvCycloScale[c]`, `rvCycloMul`, and `rvCycloConj`, which sends $z$ to $z^{-1}$. Display and `rvCycloCoeff[j](a);` use the monomial coefficients of the unique remainder modulo the monic cyclotomic polynomial.

`rvVerlinde[[sigma0,...]](a)` is the exact Verlinde character of a fusion element $a$ at the integrable label $\sigma$. Read the normalized label in descending order and set $d_i=\sigma_i+n-1-i$ with $i$ zero-based. The interpreter evaluates each Schur character at $X_i=z^{nd_i-\sum_j d_j}$ in order $m=n(k+n)$. These torus coordinates have product one. `rvQDim(a)` specializes at the zero label. At rank three and level two the adjoint value $x$ satisfies $x^2=x+1$. The example checks two conjugate fundamental characters at a nonzero torus label.

`rvFusionHandle[n,k]` constructs the handle element $H=\sum_\lambda [\lambda][\lambda^*]$ over all integrable labels. `rvConformalDim[n,k,g,[[lambda0,...],...]];` returns the coefficient of the unit in $H^g$ times the product of the puncture labels. For rank three at level two, the vacuum dimensions from genus zero through five are $1,6,45,405,4050,42525$. The interpreter bounds cyclotomic orders, integrable-label enumeration, and integer arithmetic. [Cyclotomic and sewing example](examples/rv-verlinde-tqft.asl)

### Quantum Schubert calculus for Grassmannians

`RVQGr(k,n)` is the integral small quantum cohomology of $\mathrm{Gr}(k,n)$ in its Schubert basis over $\mathbb Z[q]$. Use `rvQGrBasis[k,n,[lambda0,...]]` for $\sigma_\lambda$ and `rvQGrOne[k,n]` for the unit. Enter actual weakly decreasing partitions in the $k\times(n-k)$ rectangle. You may omit trailing zeros. A uniform shift of a partition changes its Schubert class. Elements support `zero`, `add`, `neg`, equality, signed scaling `rvQGrScale[c](a)`, and multiplication by $q^d$ with `rvQGrShift[d](a)`.

`rvQGrMul(a,b)` computes classical Littlewood–Richardson coefficients by counting integral hives. For each classical output $\lambda$, form $\beta_i=\lambda_i+k-1-i$ with zero-based $i$ and let $d=\sum_i\lfloor\beta_i/n\rfloor$. Equal residues modulo $n$ give zero. Sorting the residues in decreasing order and subtracting the staircase gives a Schubert partition in the rectangle. Its contribution is multiplied by $q^d$, the sorting sign, and $(-1)^{(k-1)d}$. `rvQGrTranspose(a)` transposes the partitions and preserves $q$, giving the complementary Grassmannian ring. Products use the smaller rectangle for the hive count. Projective space uses $h^n=q$ directly.

`rvQGrCoeff[d,[lambda0,...]](a);` reads the coefficient of $q^d\sigma_\lambda$. `rvQGrClassical(a)` keeps only the $q^0$ terms, as a Schubert representative of specialization at $q=0$. `rvQGrGW[k,n,d,lambda,mu,nu];` takes the $q^d$ coefficient at the Poincaré dual of $\nu$ in $\sigma_\lambda\sigma_\mu$, yielding the three-point genus-zero Gromov–Witten invariant. For $\mathrm{Gr}(3,6)$ the example computes $\sigma_{21}\sigma_{321}=\sigma_{333}+q(2\sigma_{21}+\sigma_3+\sigma_{111})$. It also checks $\sigma_2\sigma_{11}=q$ in $\mathrm{Gr}(2,4)$ and $h^{37}=q$ in projective space and its complementary Grassmannian. [Quantum Grassmannian example](examples/rv-quantum-grassmannian.asl)

### Quantum differential operators and periods

`RVOre` is the integral ring of operators in $q$ and $\theta$ with $\theta q=q(\theta+1)$. Use `rvOre[a,j,c]` for $c\,q^a\theta^j$, where $a,j$ are nonnegative and the optional signed coefficient defaults to one. `rvOreMul` multiplies operators in that order. The ring supports `zero`, `add`, `neg`, and equality. `rvOreCoeff[a,j](L);` reads the coefficient of $q^a\theta^j$.

In the $z=1$ Grassmannian quantum differential module, `rvQGrNabla(v)` is $\nabla v=q\,dv/dq+H(q)v$ with $H(q)$ the quantum hyperplane multiplication matrix. The `rvOreApply(L,v)` action sends $q^a\theta^j$ to $q^a\nabla^j$. `rvQGrAnnihilator[k,n]` computes the least $\theta$-order operator killing the unit. It takes the first dependence of the cyclic vectors $1,\nabla 1,\nabla^2 1,\ldots$ over $\mathbb Q(q)$ by fraction-free elimination, removes their common coefficient content, and chooses the positive leading coefficient convention for the resulting integral operator.

`rvQGrPeriod[k,n,d];` gives the normalized unit-period coefficients through $q^d$. Writing $H=H_0+qH_1$, it starts with the row $\ell_0$ that takes the coefficient of the unit. It solves $\ell_j(jI-H_0)=\ell_{j-1}H_1$ and reports $\ell_j(1)$. This specifies the coefficients even when the scalar indicial polynomial has positive integral roots. The example gives $\theta^5-4q\theta-2q$ and coefficients $(2d)!/(d!)^6$ for $\mathrm{Gr}(2,4)$. For $\mathrm{Gr}(2,5)$ it gives indicial polynomial $\theta^7(\theta-1)^3$ and period starting $1,3,19/32,49/2592$. The order for $\mathrm{Gr}(3,6)$ is $14$ in Schubert rank $20$. [Quantum connection example](examples/rv-quantum-connection.asl)


#### Factorial regularization

`rvOreFactorial[h](L)` transforms an Ore operator for the coefficient substitution $c_d \mapsto (d!)^h c_d$. It replaces each term $q^aP_a(\theta)$ by $q^aP_a(\theta)\prod_{j=1}^a(\theta+j)^h$. The operation preserves multiplication in the Ore ring. For a transformed operator $T$, `rvOreEulerFactor(T)` computes the primitive polynomial $D(\theta)$ with positive leading coefficient given by the gcd of the shifted coefficient polynomials $P_a(\theta-a)$ of $T$. Then `rvOreEulerReduce(T)` returns the operator $B$ obtained by dividing each $P_a(\theta)$ by $D(\theta+a)$. The identity $T=D B$ can be checked with `rvOreMul`. On the zero operator the two operations return one and zero, respectively. Oversized intermediate shifts can be handled with modular computation and exact divisibility certificates. The requested factor or cofactor still has bounded signed 128-bit coefficients. For a single block, reduction needs only its signed content.

`rvQGrRegularized[k,n,d];` returns the coefficients $b_j=(j!)^n c_j$ of the normalized period through degree $d$. The computation rescales the vector recurrence using the Schubert degrees, so it can compute these coefficients without forming large factorials. For the transformed unit annihilator, $T G=0$ with $G=\sum_{j\ge0}b_jq^j$. From $D B G=0$, the coefficient of $q^j$ in $B G$ must vanish when $D(j)\ne 0$. `rvQGrRegularizedRhs[k,n];` finds all nonnegative integer roots of $D$ with Cauchy's bound and rational-root divisibility. It evaluates the normalized vector period at the corresponding degrees and reports the nonzero coefficients of $R=B G$ as a map from powers of $q$ to coefficients. Root search and coefficient arithmetic use the existing work and degree bounds.

For $\mathrm{Gr}(2,5)$, $B=\theta^2-q(11\theta^2+11\theta+3)-q^2(\theta+1)^2$ and $R=0$. The coefficients $1,3,19,147,1251,\ldots$ equal $\sum_{j=0}^d\binom{d}{j}^2\binom{d+j}{j}$ in degree $d$, the Apéry numbers for $\zeta(2)$. For $\mathrm{Gr}(2,4)$, $B=\theta-2q(2\theta+1)$ and the coefficients are $\binom{2d}{d}$. For projective space, the reduced equation is $(1-q)G=1$. The example also verifies Euler cancellation from order $14$ to order $4$ for $\mathrm{Gr}(3,6)$. [Quantum regularization example](examples/rv-quantum-regularization.asl)

#### Certified Wilf–Zeilberger period recurrences

`rvWZCheck[u,v,z,n0](L) {{(i,j,c),...},{(a,b,c,m),...}};` verifies an `RVOre` operator on the sum of $F(n,k)=\binom nk^u\binom{n+k}k^v z^k$ over $0\le k\le n$. Take $u\ge1$, $v\ge0$ and a nonzero integer $z$. The first list gives the numerator of $R(n,k)$ as monomials $c n^i k^j$. The second list gives denominator factors $(a n+b k+c)^m$ with positive multiplicity $m$. An empty denominator list gives denominator one.

Write $L=\sum_a q^a P_a(\theta)$, with $\theta=q\,d/dq$. The checker multiplies the known hypergeometric ratios, clears denominators, and verifies every coefficient of the bivariate polynomial identity
$\sum_a P_a(n-a)F(n-a,k)/F(n,k) = [F(n,k+1)/F(n,k)]R(n,k+1)-R(n,k)$
over the integers. Each denominator factor must be strictly positive for $n\ge n_0$ and $0\le k\le n+1$. The check uses the two affine endpoint values at $n_0$ and their slopes. The support and $R(n,0)=0$ make both telescoping boundary terms vanish. The start $n_0$ is positive and at least the largest $q$-shift. Every coefficient of $L\sum_n (\sum_k F(n,k))q^n$ below $n_0$ is checked by direct summation. Coefficients and intermediate values use checked signed 128-bit arithmetic with sparse work and storage bounds.

The examples verify the Apéry recurrence for $\zeta(2)$, the Franel cube sum, and an alternating binomial square sum. The Apéry operator is $\theta^2-q(11\theta^2+11\theta+3)-q^2(\theta+1)^2$. Its $q^0$ multiplier in degree $n$ is $n^2$, so this recurrence determines a normalized power series in every degree. [Wilf–Zeilberger period example](examples/rv-wz-period.asl)

### Laurent periods and prime-power Cartier sections

`rvperiod P[r] {([e1,...,er],c),...};` declares a Laurent polynomial $L$ in $r$ variables. Each term has signed integer exponents and a signed integer coefficient. Equal exponents are combined. `rvCTExact[P,n];` reports the constant term of $L^n$ by exact sparse multiplication. The exact integer and Laurent exponents use checked signed 128-bit arithmetic.

`rvCTMod[P,p,s,n]` returns the constant term modulo $p^s$ for a prime $p$ and positive precision $s$. The index is an unsigned 128-bit integer. `RVCTResidue(p,s)` is the residue type, and `rvCTResidue[p,s,v]` reduces a nonnegative integer to that type. Residue values support `zero`, `add` and equality, using least nonnegative representatives for display.

For the state model put $Q=1-zL$ and $M=p^{s-1}$. An `RVCTState(P,p,s)` stores the sparse numerator $V$ of $V/Q^M$ modulo $p^s$. `rvCTInitial[P,p,s]` starts with $V=Q^{M-1}$. `rvCTSection[d](v)` for $0\le d<p$ multiplies by $Q^{(p-1)M}$. It keeps monomials whose $z$ exponent is $d$ modulo $p$ and whose spatial exponents are all divisible by $p$. It subtracts $d$ from the $z$ exponent and divides all exponents by $p$. The congruence $Q(z,x)^{pM}\equiv Q(z^p,x^p)^M\pmod{p^s}$ preserves the denominator. At precision one the same section can use the $z^d$ coefficient $L^d$ of $Q^{p-1}$ directly.

Apply sections from the least significant base-$p$ digit to the most significant. `rvCTRead(v)` reads the coefficient of $z^0x^0$ in the current numerator, which equals the constant term of the corresponding rational series at $z=0$. Section results are actual state values. They support `show`, `zero`, `add` and equality, and calls such as `rvCTRead(rvCTSection[2](rvCTSection[1](v)))` preserve the intermediate numerator. The finite state type has $z$ degrees $0,\ldots,M-1$ and spatial exponent intervals given by $M$ times the coordinate bounds of $1$ and $L$. Ordinary finite enumeration applies within the interpreter's size limit. A modular query reuses repeated transitions only after matching the whole sparse state and the digit.

For $L=(1+x)(1+y)(1+x+y)/(xy)$, expanding the last factor in $L^n$ and extracting the two spatial coefficients gives $\sum_k\binom{n}{k}^2\binom{n+k}{k}$. These Apéry numbers are $(n!)^5$ times the normalized $\mathrm{Gr}(2,5)$ quantum-period coefficients. The example uses the two Cartier sections for $7$ in base $3$ to obtain $6$ modulo $9$, and evaluates the period at $3^{38}-1$ modulo $27$ and $2^{85}-1$ modulo $16$. It also shows a Laurent polynomial whose mod-two period needs a nonconstant state. [Laurent period example](examples/rv-laurent-period.asl)

### Quantum p-curvature

`RVQGrP(p,k,n)` is the polynomial Schubert module in characteristic $p$. `rvQGrMod[p](v)` reduces an integral `RVQGr(k,n)` section exactly modulo $p$, and `rvQGrPBasis[p,k,n,[lambda0,...]]` constructs a basis section. The module supports `zero`, `add`, `neg`, equality, and the `rvQGrNabla`, `rvQGrScale`, `rvQGrShift`, `rvQGrTranspose`, `rvQGrClassical`, and `rvQGrCoeff` operations.

The Euler derivation $\theta=q\,d/dq$ satisfies $\theta^{[p]}=\theta$, so the $p$-curvature is the polynomial-linear map $\Psi_p=\nabla^p-\nabla$. `rvQGrCurvatureCol[p,k,n,[lambda0,...]]` computes it on a basis section by iterating the modular connection. `rvQGrCurvature(v)` applies those columns to an arbitrary polynomial section. `rvQGrCurvatureSatake(v)` computes the same action using the signed projective-space connection with $q'=(-1)^{k-1}q$. A Schubert partition gives the increasing exterior basis indices $\lambda_{k+1-j}+j-1$ for $1\leq j\leq k$. The exterior connection's curvature acts by the sum over individual factors with their exterior signs.

`rvQGrCurvatureChar[p,k,n];` gives the coefficients of $\det(TI-\Psi_p)$ in ascending powers of $T$. Each entry lists the nonzero pairs of $q$-exponent and prime-field coefficient. `rvQGrCurvatureCoeff[p,k,n,j];` reads the coefficient of $T^j$. `rvQGrCurvatureHorizontal[p,k,n];` checks and rewrites the characteristic coefficients in $Q=q^p$. For $\mathrm{Gr}(2,4)$ the characteristic polynomial is $T^2(T^4-4q^p)$. In characteristic two it becomes $T^6$, while the unit satisfies $\Psi_2^2(1)\neq0$ and $\Psi_2^3(1)=0$. The example also gives three nonconstant coefficients for $\mathrm{Gr}(3,6)$ in characteristic five. [Quantum p-curvature example](examples/rv-quantum-pcurvature.asl)

### Genus-one WZW modular action

`RVWzw(n,k)` is the torus state module over `RVCyclo(D)`, with $L=n+k$ and $D=\operatorname{lcm}(2nL,4L)$. Its basis is indexed by integrable weights. `rvWzwBasis[n,k,[lambda0,...]]` accepts weakly increasing weights and removes the determinant translate. The module supports `zero`, `add`, `neg`, equality, cyclotomic scalar multiplication `rvWzwScale(c,v)`, and coordinate extraction `rvWzwCoord[[lambda0,...]](v)`. `rvCycloEmbed[D](c)` embeds a cyclotomic integer of order $m\mid D$ by $z_m\mapsto z_D^{D/m}$.

Write the normalized weight in descending order, set $x_{\lambda,i}=\lambda_i+n-i$ for one-based $i$, and let $X_\lambda$ be its coordinate sum. The entry `rvWzwSEntry[n,k,lambda,mu]` is the signed alternant with exponent $-(n\sum_i x_{\lambda,\sigma(i)}x_{\mu,i}-X_\lambda X_\mu)/(nL)$ in turns. The entry `rvWzwTEntry[n,k,lambda,h]` is the root of unity with exponent $h\,[Q(x_\lambda)-Q(\rho)]/(2nL)$, where $Q(x)=n\sum_i x_i^2-(\sum_i x_i)^2$ and $\rho_i=n-i$. The signed twist power defaults to one. The interpreter computes the alternant over the integral cyclic group ring and reduces it to the exact cyclotomic basis.

`rvWzwS(v)`, `rvWzwSAdjoint(v)`, and `rvWzwT[h](v)` apply these matrices and the conjugate transpose of $S$. `rvWzwCharge(v)` sends each weight to its contragredient. The raw normalization gives $S^2=aC$, $SS^*=|a|I$, and $(ST)^3=bC$, where $a=(-1)^{n(n-1)/2}nL^{n-1}$ and $b$ is the vacuum coefficient of $(ST)^3$. The example computes $b=375(z_{60}^4+z_{60}^{14})$ at rank three, level two, checks $b^2=a^3\exp(2\pi i\,k(n^2-1)/(4L))$, and recovers an $S$ entry from the vacuum amplitude and the conjugate Verlinde character. [WZW modular example](examples/rv-wzw-modular.asl)

#### Signed Galois action on the alcove

`rvWzwGalois[a](v)` applies the signed permutation $G_a$ of the integrable labels, with coefficients left fixed. The interpreter reduces the signed exponent $a$ modulo $D$ and requires a unit. For each weight it forms $a(\lambda+\rho)-\rho$ in decreasing order and uses the Kac–Walton alcove fold. The output is the folded label $\pi_a(\lambda)$ multiplied by the folding sign $\epsilon_a(\lambda)$.

Use `rvCycloGalois[a](c)` to substitute $z_D\mapsto z_D^a$ in a cyclotomic integer. `rvWzwCoeffGalois[a](v)` applies this automorphism to the coefficients of a state. Composing coefficientwise transport with $G_a$ gives the semilinear action. For the raw matrices the identities are $\sigma_a(S_{\lambda\mu})=\epsilon_a(\lambda)S_{\pi_a(\lambda),\mu}$ and $T_{\pi_a(\lambda)}=T_{\pi_a(0)}T_\lambda^{a^2}$. The signed operators satisfy $G_aG_b=G_{ab}$ and $G_{-1}=(-1)^{n(n-1)/2}C$.

The fold adds a vector in $L\mathbb{Z}^n$ of coordinate sum zero, permutes coordinates, and removes a common translate. These operations change the alternant only by the permutation sign. The quadratic form $Q$ is unchanged by a common translate. A shift by $Lz$ with $\sum_i z_i=0$ changes $Q$ by a multiple of $2nL$, since $\sum_i z_i^2$ is even. This gives the square of the Galois exponent and the vacuum phase in the twist identity.

At $\mathrm{SU}(3)_2$, the exponent $7$ sends the vacuum to the adjoint $[0,1,2]$ and sends the adjoint to the negative vacuum. Its square acts by $-1$ on the six labels. The example checks an $S$ row, the projective twist phase, and the distinct coefficient transport on a state with nonrational amplitudes. [WZW Galois example](examples/rv-wzw-galois.asl)

### Doubled WZW lens-space invariants

`rvLensTV[n,k,p,q]` computes the lens-space Turaev–Viro invariant of the unitary $\mathrm{SU}(n)_k$ theory doubled. For $p>0$, it accepts a signed $q$ coprime to $p$ and uses its residue modulo $p$. The inputs $p=0,q=\pm1$ describe $S^2\times S^1$. The normalization assigns $1$ to $S^2\times S^1$ and $1/\operatorname{Dim}\mathcal C$ to $S^3$, where $\operatorname{Dim}\mathcal C$ is the sum of the squares of the quantum dimensions.

The interpreter decomposes a gluing matrix with first column $(q,p)$ by signed Euclidean division, using $S=\left(\begin{smallmatrix}0&-1\\1&0\end{smallmatrix}\right)$ and $T=\left(\begin{smallmatrix}1&1\\0&1\end{smallmatrix}\right)$. For a word with $m$ copies of raw $S$, it returns the product of its vacuum amplitude and its complex conjugate, divided by $[n(n+k)^{n-1}]^m$. The calculation cancels common numerator and denominator factors as it proceeds. `rvWzwTV[n,k]{S,T[h],...}` uses a supplied gluing word in matrix multiplication order, acting on the vacuum from right to left. Left and right twists fix the vacuum and projective framing phases cancel in the doubled amplitude.

Both operations return `RVCycloQ(D)` with the same $D$ as `RVWzw(n,k)`. Its display gives the cyclotomic numerator coefficients and a common positive denominator, reduced by their scalar gcd. Use `rvCycloRat[d](c)` for a cyclotomic integer divided by positive $d$. Values support `zero`, `add`, `neg`, and equality. At rank two and level two, $L(8,1)$ has invariant $0$, $L(8,3)$ has invariant $1$, and $L(2,1)$ has invariant $(2-z_{16}^2+z_{16}^6)/4=(2-\sqrt2)/4$. [Lens-space invariant example](examples/rv-lens-tv.asl)

### Genus generating functions for conformal blocks

For a product of puncture labels $x$, let $D_g$ be the coefficient of the unit in $xH^g$. `rvConformalSeries[n,k,labels];` gives the reduced rational function $\sum_{g\geq 0}D_gt^g=P(t)/Q(t)$, with $Q(0)=1$. Use `rvConformalNum[n,k,labels];` or `rvConformalDen[n,k,labels];` to read one polynomial. Coefficients are listed in ascending powers of $t$. An empty numerator and denominator `[1]` give the zero series. A nonzero total central charge among the normalized labels gives this zero case.

The interpreter builds the matrix for multiplication by $H$ on the weights reachable from the unit. Cayley–Hamilton gives the initial denominator $\det(I-tH)$ on that space. Exact polynomial gcd cancellation with the numerator gives the minimal denominator for the chosen insertions. Its coefficients give the integer recurrence for the dimensions. For $\mathrm{SU}(3)_3$, inserting a fundamental and its dual reduces the denominator degree from three to two.

`rvConformalMod[n,k,g,m,labels];` calculates $D_g$ modulo any positive modulus fitting an unsigned 128-bit integer. It reduces puncture products and exponentiates the handle matrix in the residue ring. `rvConformalExact[n,k,g,labels];` uses nonnegative decimal arithmetic for the handle iterates to give dimensions beyond signed 128 bits within the work and output limits. The example includes residues at genus $10^{18}$ and the exact $\mathrm{SU}(3)_2$ vacuum dimension at genus $100$. [Conformal series example](examples/rv-conformal-series.asl)

### Temperley–Lieb braid paths and the Jones polynomial

`RVTL(k,n,h)` is the path sector for $n$ fundamental strands of $\mathrm{SU}(2)_k$ ending at height $h$. The level $k$ is positive. A basis vector `rvTLPath[k,[a0,...,an]]` has $a_0=0$, $a_n=h$, unit steps, and $0\le a_i\le k$. The coefficients are cyclotomic rationals with root $A$ of order $4(k+2)$. Use `zero`, `add`, `neg`, `rvTLScale(c,v)`, and `rvTLCoord[[a0,...,an]](v)` to work with path vectors. An impossible endpoint gives a zero sector.

Write $q=A^2$, $d_h=(-1)^h[h+1]_q$, and $\delta=d_1$. The operator `rvTLE[i](v)` gives zero on a path unless $a_{i-1}=a_{i+1}=a$. When they agree, it sums over all allowed new middle heights $b$ with coefficient $d_b/d_a$. The interpreter uses exact cyclotomic arithmetic for these ratios. The operation `rvTLBraid[[i1,...,im]](v)` uses the generators $B_i=A I + A^{-1}e_i$ and $B_i^{-1}=A^{-1}I+A e_i$. Indices satisfy $1\le |i|<n$, and the word is in matrix multiplication order, acting from right to left. `rvTLTrace[k,n,h,[i1,...,im]]` gives the trace on the indicated sector.

`rvTLJones[k,n,[i1,...,im]]` computes the Jones polynomial of the braid closure evaluated at $t=A^{-4}$. It weights each sector trace by $d_h$, divides their sum by $\delta$, and multiplies by $(-A^3)^{-w}$ for writhe $w$. This gives $1$ for the unknot and $\delta^{-1}$ for the zero-strand empty closure, following this normalization. The trefoil gives $t+t^3-t^4$ and the figure-eight gives $t^2-t+1-t^{-1}+t^{-2}$.

Use `rvCycloQMul(a,b)` for products and `rvCycloQGalois[a](c)` to substitute $A\mapsto A^a$ for $a$ coprime to the root order. At level three on three strands ending at height one, $B_1B_2^{-1}$ has determinant $1$ and trace $(3-\sqrt5)/2$. The substitution $A\mapsto A^3$ sends the trace to $(3+\sqrt5)/2>2$, so its conjugate has an eigenvalue off the unit circle. A finite-order matrix stays finite-order under a field automorphism, proving this braid matrix has infinite order. [Braid paths and Jones example](examples/rv-tl-jones.asl)


#### Jones–Wenzl projectors and colored closures

`rvTLProject[first,m](v)` applies the Jones–Wenzl projector to $m$ contiguous elementary strands of a path vector. The first index is one-based. The operator keeps the two boundary heights of the block and changes its intermediate heights. It uses $f_0=f_1=I$ and the recursion $f_{r+1}=f_r-(d_{r-1}/d_r)f_re_rf_r$, with $f_r$ embedded in the first $r$ strands of the block. At level $k$, the path quotient kills $f_{k+1}$. The computation caches exact columns on local path segments and counts both work and storage.

`rvTLColoredJones[k,n,m,[word]]` colors each of the $n$ strands by the same integrable label $0\le m\le k$. Each letter becomes a crossing of two cables of $m$ elementary strands, with $m^2$ crossings. The evaluator applies $f_m$ to every cable and takes the quantum trace with sector weights $d_h$. For colored writhe $w$, it multiplies by $\theta_m^{-w}/d_m$, where $\theta_m=(-1)^m A^{m(m+2)}$. This normalizes the unknot to one. Color zero gives the trivial invariant and color one gives `rvTLJones`. The elementary braid expansion has a limit of 4,096 crossings. The empty braid uses the trace $d_m^n$ of the identity directly.

On two colored strands, channel $h$ has crossing eigenvalue $(-1)^{(2m-h)/2}A^{(h(h+2)-2m(m+2))/2}$. The allowed $h$ are $0,2,\ldots,\min(2m,2k-2m)$. The example checks this channel formula for the color-two trefoil at level five and uses a genuine six-strand cable for the color-two figure-eight. It also checks both Markov stabilizations on eight elementary strands. [Colored Jones example](examples/rv-tl-colored.asl)

#### Component colors and framed closures

`rvTLFramed[k,n,[word],[colors],[framings]]` gives the unnormalized quantum trace of a framed braid closure as an exact `RVCycloQ(4(k+2))` value, with level $k$ and $n$ original strands. Its closure components are the cycles of the strand permutation, ordered by least original strand. Supply one integrable color and one signed integer target framing per component. Braid words act from right to left and move the cable widths along with the strands. A crossing of widths $a$ and $b$ contributes $ab$ elementary crossings, and each initial cable carries its Jones–Wenzl projector. Color zero is the empty cable. At level zero it is the only color.

For signed self-writhe $w_i$, target framing $f_i$, and color $a_i$, the evaluator multiplies the quantum trace by $\prod_i\theta_{a_i}^{f_i-w_i}$ in the signed Temperley–Lieb convention above. A zero-framed unknot of color $a$ gives $d_a$, and the empty closure gives one. `rvTLLinking[n,[word],[framings]];` reports the integer matrix with these target framings on its diagonal and half the signed intercomponent crossing sums off the diagonal. The example checks the multicolor Hopf link by its two fusion channels and verifies positive and negative Markov stabilization. [Framed closure example](examples/rv-tl-framed.asl)

#### SU(2) Kirby surgery

`rvTLSurgery[k,n,[word],[framings]]` evaluates integral surgery on a framed braid closure with nonsingular linking matrix. The result is an exact cyclotomic value with $S^3$ normalized to one. It sums the framed trace over all integrable component colors, weighted by $\prod_i d_{a_i}$, to obtain the Kirby sum $F(\Omega)$. For linking inertia $(b_+,b_-)$ it divides by $p_+^{b_+}p_-^{b_-}$, with $p_\pm=\sum_a d_a^2\theta_a^{\pm1}$ and $p_+p_-=\mathcal D^2=\sum_a d_a^2$. An integral symmetric elimination computes inertia and includes a two-by-two pivot for hyperbolic blocks.

`rvTLSurgeryTV[k,n,[word],[framings]]` also accepts singular linking matrices. For $\ell$ closure components it returns $F(\Omega)\overline{F(\Omega)}/(\mathcal D^2)^{\ell+1}$, so $S^3$ has value $1/\mathcal D^2$ and $S^1\times S^2$ has value one. The signed Temperley–Lieb twist differs from the unitary WZW twist by $(-1)^a$. The Galois substitution $A\mapsto -A$ removes that sign and fixes $d_a^2$, giving the comparison with `rvLensTV` for one-component lens surgeries. Untouched components contribute separate unknot Gauss sums, and the evaluator checks the number of remaining colorings and elementary crossing bound before expansion.

The example checks both blowups, the zero-framed Hopf presentation of $S^3$, a nonzero handle slide from framings $(2,3)$ to $(5,3)$, and the distinction between $+1$ surgery on the trefoil and the unknot. [Kirby surgery example](examples/rv-tl-surgery.asl)

#### Integral Khovanov cohomology of braid closures

Declare a braid by `khbraid K[n,[i1,...,im]];` using the same signed letters and order as `rvTLBraid`. The closure uses the orientation of the braid strands. In each positive crossing the zero resolution is vertical and the one resolution is cup–cap. At a negative crossing they are reversed. The cube order is the original order of the letters. Each resolved circle is labeled $1$ or $X$ in $\mathbb Z[X]/(X^2)$, of quantum degrees $1$ and $-1$ respectively. If the set of one resolutions has size $s$, the bigrading is $i=s-n_-$ and $j=\deg(\text{labels})+s+n_+-2n_-$.

A basis term `khTerm[K,[ones],[xs]]` lists the positions of the one resolutions and the circles bearing $X$. Both indices start at zero. Circles are ordered by their first strand segment, starting at the top and then scanning crossing outputs in the order of braid action. Give the circle indices in increasing order. An optional fourth argument sets the signed integral coefficient. The equivalent basis-coordinate constructor is `khBasis[K,i,j,g]`. Use `khD` for the differential $KhC(K,i,j)\to KhC(K,i+1,j)$. It adds the merge and split maps with sign $(-1)^{\#\{\text{earlier ones}\}}$, using $X^2=0$, $\Delta(1)=1\otimes X+X\otimes 1$ and $\Delta(X)=X\otimes X$.

`khGroup[K,i,j];` reports the chain rank, cohomology free rank, and positive torsion invariant factors greater than one. Invariant-only queries compute the ranks and incoming Smith factors directly. The integral kernel is a direct summand because the outgoing image is free. For representatives and nulls, the interpreter splits the outgoing map, expresses the incoming map in the resulting integral kernel basis, and tracks both unimodular changes of basis in its Smith form. Use `khClass` on cycles, `khRep` on classes in `KhH(K,i,j)`, and `khGen[K,i,j,g]` for a generator, with the free generators followed by the torsion generators in invariant-factor order. Zero, addition, negation, and `khScale[c]` work on chains and classes. For a torsion generator `khNull[K,i,j,g]` returns a chain $w$ with $dw=mz$, where $z$ is the generator's representative and $m$ its exact order. The braid input supports up to 128 crossings. Integer arithmetic is checked in the signed 128-bit range, and computation and sparse storage limits apply.

`khEuler[K];` lists the coefficients of the graded Euler Laurent polynomial in increasing quantum degree, including the unknot $q+q^{-1}$. The calculation sums the circle-label dimensions in each resolution. `khEulerAt[K,k]` evaluates at $q=-A^{-2}$ for the same cyclotomic parameter $A$ as the level-$k$ Temperley–Lieb calculation. The output equals $\delta$ times the repo's normalized `rvTLJones`, with $\delta=-A^2-A^{-2}$. In the trefoil example $KhH(\mathrm{Trefoil},3,7)\cong\mathbb Z/2$ and the program checks the actual integral relation $dw=2z$. It also computes the mirror torsion and compares the figure-eight and both Markov stabilizations. [Integral Khovanov example](examples/integral-khovanov.asl)

#### Filtered Lee cohomology and Rasmussen's invariant

`LeeC(K,i,f)` contains rational cochains in homological degree $i$ whose terms have quantum degree at least $f$. Use `leeLift(x)` to regard a homogeneous integral `KhC(K,i,j)` chain as an element of `LeeC(K,i,j)`. Lee cochains support `zero`, `add`, `neg`, equality, and `leeScale[num,den]` for a rational scalar with positive denominator. `leeFilter[f](x)` changes the declared filtration after checking the actual support.

`leePhi(x)` gives the integral perturbation on homogeneous Khovanov chains, in bidegree $(1,4)$. It adds $m(X,X)=1$ and the $1\otimes 1$ term in $\Delta(X)$ to the same signed cube. On Lee cochains `leeD(x)` applies $d_{\mathrm{Kh}}+\Phi$ and preserves the decreasing quantum filtration. The same circle resolutions, labels, and edge signs are used for both differentials.

`leeRank[K,i,f];` gives the dimension of the image $H^i(F^fC)\to H^i(C)$. It counts the kernel on chains of quantum degree at least $f$, subtracts the full incoming rank, and adds back the rank of the incoming image projected to degrees below $f$. `leeLevels[K,i];` lists the filtration jumps with multiplicity. Exact fraction-free elimination works over $\mathbb Q$ with checked signed 128-bit arithmetic and sparse work and storage bounds. It needs only the requested homological degree and its incident cube edges.

`leeRep[K,i,f,k]` returns the $k$th chosen cycle in the associated graded piece at the exact jump $f$, with zero-based index. The cycles form a basis of the filtered cohomology when taken across the jumps. `leeFiltration(x);` checks closure and gives the highest filtration attained by the class, or reports `boundary`. `leeIndependent(x,...);` checks independence modulo the actual boundary span for cycles of the same braid and homological degree, which may have different declared filtrations. `leeNull(x)` constructs an exact rational preimage of a boundary in degree $i-1$. Its result type shows a conservative filtration bound from the braid.

`leeS[K];` takes the mean of the two surviving jumps in homological degree zero for a knot. Links also have the full filtration queries. For the positive braid $(\sigma_1\sigma_2)^8$ closing to $T(3,8)$, the jumps are $13,15$ and $s=14$. The bound $g_4\ge |s|/2$ and the genus-seven surface from the braid give $g_4=7$. The example also verifies a nontrivial filtration lift and an exact half-integral boundary preimage on a stabilized unknot. [Lee and Rasmussen example](examples/lee-rasmussen.asl)

#### Annular Khovanov homology and the transverse cycle

A resolution circle is essential exactly when it crosses the braid closure seam an odd number of times. Its label contributes $+1$ for $1$ or $-1$ for $X$ to the annular degree $k$. Trivial circles contribute zero. `khAnnular[k](c)` selects that homogeneous component as `KhAC(K,i,j,k)`, and `khForgetK` returns its underlying integral chain. Use `khABasis[K,i,j,k,g]` for the ordered homogeneous basis. The actual signed differential has degrees $0$ and $-2$ in $k$, given by `khD0` and `khDMinus2`. The bounded `khAnnularCheck[K];` checks these two degrees and $d^2=0$ on the enhanced cube, which gives $d_0^2=\partial^2=d_0\partial+\partial d_0=0$ by degree.

`khAnnularGroup[K,i,j,k];` computes the integral homology of $d_0$ by Smith reduction. The class type is `KhAH(K,i,j,k)`. Use `khClass`, `khRep`, `khAGen` and `khANull` for classes, representatives, generators and torsion nulls. `khAnnularEuler[K];` reports the signed Euler polynomial as $((j,k),\text{coefficient})$ pairs.

`khPsi[K]` is the transverse cycle in the oriented resolution, with all circles labeled $X$. Its degrees are $(i,j,k)=(0,w-n,-n)$ for braid writhe $w$ and strand count $n$. `KhFC(K,i,j,k)` is the $\mathbb F_2$ subcomplex cut off at annular degree $k$. `khToF2[k]` reduces integer coefficients modulo two and retains the terms at or below that cutoff. It can also widen or truncate an already filtered chain. `khD` on `KhFC` uses the whole filtered differential, and on the dual type `KhFD(K,i,j,k)` it uses the transpose in the opposite homological direction. `khPair(dual,chain)` is `S.hi` when their pairing is one.

`khKappa[K];` finds the first cutoff $k$ at which $\psi$ is a boundary and reports $n+k$, or infinity if $\psi$ survives. In the finite case `khKappaNull[K]` gives its actual filtered primitive. `khKappaDual[K]` annihilates all boundaries at cutoff $k-2$ and pairs to one with $\psi$. For infinity the dual certifies nonvanishing in the full complex. The four-strand examples have $\kappa=4$ and $\kappa=2$, and the program checks both primitives and lower-cutoff duals. The invariant is defined for braid conjugacy. The illustrated positive stabilization of the second four-braid raises $\kappa$ from two to four, with both certificates checked again. [Annular example](examples/annular-khovanov.asl) · [Hubbard and Saltz](https://arxiv.org/abs/1507.06263)

#### Annular $\mathfrak{sl}_2$ and exterior current actions

`khE`, `khF` and `khH` act on both `KhAC` and `KhAH`. On each enhanced resolution, order the essential circles by their first radial seam encounter. `khE` switches $X$ to $1$, and `khF` switches $1$ to $X$, with sign $(-1)^{\mathrm{depth}}$ starting at depth zero. The respective shifts in $(j,k)$ are $(2,2)$ and $(-2,-2)$. `khH` multiplies by $k$. These operators commute with $d_0$ and satisfy $[E,F]=H$, $[H,E]=2E$ and $[H,F]=-2F$. The class operations use integral representatives and Smith projection, including torsion.

`khIrreps[K,i,jMinusK];` computes the multiplicities of the rational irreducibles as `(highest weight, multiplicity)` pairs using the actual raising maps on the free quotient. `khHighest[K,i,jMinusK,m,g]` gives the zero-based $g$th integral cocycle in highest weight $m$, killed by $E$ on chains. Their $F$ strings give a basis over $\mathbb Q$. For the identity five-braid the decomposition is $V_5\oplus4V_3\oplus5V_1$, and the longest string has $F^5z=120$ times the all-$X$ state.

The odd maps `khVMinus`, `khVZero` and `khWPlus` have degree one in $i$ and shifts $(0,-2),(2,0),(4,2)$ in $(j,k)$. They are $v_-=d_{-2}$, $v_0=[E,v_-]$ and $w_+=[E,v_0]$. `khLee0` and `khLeePlus` give the annular-degree-zero and degree-two components $\Phi_0,\Phi_+$ of the Lee perturbation on chains. In these signs $w_+=-2\Phi_+$. The two possible odd defects are $\{v_-,w_+\}=2\{d_0,\Phi_0\}$ and $v_0^2=-\{d_0,\Phi_0\}$, so the exterior current relations hold over $\mathbb Q$ on homology with $v_+=w_+/2$. For a cocycle $z$, the displayed identities give the boundary primitives $2\Phi_0 z$ and $-\Phi_0 z$. `khActionCheck[K];` checks the even and odd identities on every enhanced cube cell, subject to the work limit. The example exhibits a nonzero odd image and a nonzero chain defect with its primitive. Over $\mathbb Q$, annular link isotopy, including braid conjugation, preserves the action up to isomorphism. Positive stabilization can change the annular link. [Annular current example](examples/annular-current.asl) · [Grigsby, Licata and Wehrli](https://arxiv.org/abs/1505.04386)

#### Annular-to-ordinary Khovanov pages

`KhPC(K,p,i,j,k)` is the signed Khovanov cochain space over $\mathbb F_p$ in bidegree $(i,j)$ with annular filtration at most $k$, where the prime $p$ fits in `u32`. `khMod[p,k](c)` reduces an integral cochain and selects its terms through $k$, or changes the cutoff of a cochain over the same field. `khD` applies the whole differential. These types and the spectral page classes support addition, signed `khScale`, and bounded finite enumeration.

`KhPage(K,p,r,i,j,k)` contains $E_r^{i,j,k}$ of the increasing annular filtration. `khPageDim[K,p,r,i,j,k];` reports its dimension, and `khPageGen[K,p,r,i,j,k,g]` gives its zero-based $g$th generator. $E_0$ is the annular grading of the cochain space. $E_1$ is annular Khovanov homology over $\mathbb F_p$, so reducing integral ranks alone would miss torsion contributions. At stabilization the sum over $k$ gives ordinary mod-$p$ homology in bidegree $(i,j)$. The differential `khPageD` shifts $(i,j,k)$ by $(1,0,-2r)$.

`khPageRep(x)` returns a lifted cochain in $F_k$ whose differential lies in $F_{k-2r}$. `khPageClass[r,k](c)` checks both cutoffs on the full cochain and projects it to $E_r$. `khPageNext(x)` advances a $d_r$-cycle, and `khPagePrev(x)` chooses a representative of a later class as a cycle on the previous page. `khPageBoundary(x)` constructs a $d_r$ preimage when $x$ is a boundary. The chosen preimage has full differential equal to `khPageRep(x)` when both cochains are compared in the source cutoff. The reduction orders cells by filtration and tracks paired cochains $b$ and $db$ with distinct leading cells. A pair with filtration drop $2r$ gives an actual $d_r$ with both classes disappearing on the following page.

The positive two-braid demonstrates $d_1$ over $\mathbb F_3$ and compares it to the integral $v_-$ action. The transverse classes of the Hubbard–Saltz braids $A$ and $B$ are hit by $d_2$ and $d_1$, respectively. Positive stabilization of $B$ moves its hit to $d_2$. The program checks the actual full cochains bounding all three transverse cycles and a dual lower-cutoff obstruction for $A$. [Annular spectral sequence example](examples/annular-khovanov-pages.asl)

#### Annular Rasmussen profiles

`leeOrientation[K]` is the canonical rational Lee cycle of the braid orientation. It uses the oriented smoothing, with labels $1+X$ and $1-X$ alternating along the nested circles. For $t=p/q\in[0,2]$, give every enhanced state weight $j-tk$. The weight of a cochain is the minimum weight of its nonzero terms. `leeDt[K,p,q];` gives the highest weight of a representative of the orientation class after adding Lee boundaries.

`leeDtRep[K,p,q]` is an optimal representative and `leeDtNull[K,p,q]` is a primitive for its difference from `leeOrientation[K]`. The dual `leeDtDual[K,p,q]` has type `LeeDual(K,0)`. It pairs to one with the orientation and vanishes on all boundaries and all strictly higher-weight states. `leeDtDualD` is the transpose Lee differential, and `leeDtPair(dual,chain)` evaluates on the same braid and homological degree as an exact rational in `RVCycloQ(1)`. `leeDtWeight[p,q](value);` reports the lowest weight of a cochain or the highest weight of a dual's support.

`leeDtProfile[K];` prints pieces `((p,q),(r,s),j,k)` for the line $j-tk$ on $[p/q,r/s]$. The calculation uses every exact intersection of the finitely many relevant weight lines. Their order is constant between consecutive cuts, so fraction-free elimination identifies the active line throughout each interval. Equal adjacent pieces are merged. The figure-eight example checks both optimality witnesses at $t=1/2$, improving the raw weight $-3/2$ to $-1/2$. The four-braid $A$ has six pieces with interior breakpoints $1/2,2/3,1,4/3,3/2$. Its two witnesses at $t=2/3$ certify $d_{2/3}=5/3$, and the value at the interior point $t=3/5$ is $9/5$. For knots $d_0=s-1$, and for closed braids $d_1$ is the writhe. [Annular Rasmussen example](examples/annular-rasmussen.asl) · [Grigsby, Licata and Wehrli](https://arxiv.org/abs/1612.05953)


### Certified Laurent periods and WZ recurrences

`rvCTWZCheck[P,u,v,z,n0](L){numerator,denominator};` checks that the declared Laurent polynomial `P` is the canonical mirror of the WZ summand $\binom{n}{k}^{u}\binom{n+k}{k}^{v}z^k$. With $u-1$ variables $x_i$ and $v$ variables $y_j$, the canonical polynomial is
\[
\Lambda=\prod_i(1+x_i)\prod_j\frac{1+y_j}{y_j}
\left(1+z\prod_i x_i^{-1}\prod_j(1+y_j)\right).
\]
Choosing $k$ copies of the second term of the bracket in $\Lambda^n$ contributes $\binom{n}{k}z^k$. Each $x_i$ supplies $[x_i^k](1+x_i)^n=\binom{n}{k}$ and each $y_j$ supplies $[y_j^n](1+y_j)^{n+k}=\binom{n+k}{k}$. The resulting constant term is therefore the WZ sum in every degree.

The interpreter compares every coefficient of the expanded polynomial with the normalized declaration. It then runs the same exact WZ certificate and initial-coefficient checks as `rvWZCheck` on the Ore operator $L$. A successful check certifies the recurrence for the declared Laurent period.

The optional final argument `[[col0],...]` gives the exponent columns for the $x$ variables followed by the $y$ variables. Each column has the number of coordinates declared for `P`. The interpreter verifies full column rank over $\mathbb Q$ by fraction-free elimination. The associated map of exponent lattices is injective, so only a zero exponent can map to zero, even when the image is a proper sublattice. For $u=2,v=0,z=-1$, the single column `[[2,4]]` takes $x-x^{-1}$ to $X^2Y^4-X^{-2}Y^{-4}$ with the same period.

### A certified mirror for $\mathrm{Gr}(2,5)$

The example computes the Grassmannian unit annihilator, applies the factorial transform, and checks the Euler factorization. It uses `rvCTWZCheck` on the resulting cofactor $B=\theta^2-q(11\theta^2+11\theta+3)-q^2(\theta+1)^2$ and on the declared Laurent polynomial. This certifies both the mirror identity and the WZ recurrence. The quantum forcing term vanishes. In degree $n>0$ the coefficient of the new term is $n^2$, so the common initial value $A(0)=1$ determines the entire period.

For this mirror the general Laurent formula becomes $L=(1+x)(1+y)(1+x+y)/(xy)$. Its period is $[x^ny^n](1+x)^n(1+y)^n(1+x+y)^n$. Choosing $(1+x)^j y^{n-j}$ from the last factor gives $\binom{n}{j}[x^n](1+x)^{n+j}[y^j](1+y)^n$, hence the sum $\sum_j\binom{n}{j}^2\binom{n+j}{j}$ in every degree.

The program applies the base-three sections for $7$ to obtain $A(7)\equiv6\pmod{9}$. Modulo $27$ it checks that the whole state after the first digit $2$ is fixed by that section. All $38$ digits of $3^{38}-1=1350851717672992088$ are $2$, so reading this fixed state gives $A(1350851717672992088)\equiv19\pmod{27}$. [Grassmannian–Apéry proof example](examples/rv-grassmannian-apery-proof.asl)

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

#### De Rham cohomology of a polynomial Artin–Schreier cover

For a declared cover $y^2+y=f(x)$ with affine denominator $1$ and $\deg(f)=2g+1$, `ASDeRham(C,P)` is its $2g$-dimensional first de Rham cohomology over the binary field with modulus `P`. The Čech basis starts with $\alpha_i=x^i\,dx$ for $0\le i<g$, then uses $\eta_j$ for $1\le j\le g$. The lift $\eta_j$ has overlap function $h=y/x^j$, form $[f'x^{-j}]_{\ge0}\,dx$ on the affine chart, and that form plus $dh$ on the other chart. `rvASDeRhamGen[C,P,i]` selects a basis coordinate in this order, starting at zero. `rvASDeRhamIncl(x)` includes `ASDiff(C,P)` in this basis.

Use `zero`, `add`, and `rvASDeRhamScale[c](x)` for vector operations. Frobenius `rvASF(x)` squares scalar coefficients, and Verschiebung `rvASV(x)` takes their square roots. On the included holomorphic differentials, Verschiebung agrees with `rvASCartier`. `rvASDeRhamPair(x,y)` evaluates the alternating cup pairing in `BinaryScalar(P)`. Use `rvFieldElt[P,c]` for a field element with the same bit encoding used for differential coefficients. The basis pairing is $\langle\alpha_i,\eta_j\rangle=\delta_{i,j-1}$ and $\langle\eta_i,\eta_j\rangle=(i+j\bmod 2)[x^{i+j}]f$.

`rvASEO[C];` computes the Ekedahl–Oort final type by closing the zero and full de Rham spaces under $V$ and $F^{-1}$ and using the ranks of $V$ along this canonical flag. The final type for the polynomial covers is $\nu_i=\lfloor i/2\rfloor$ for $1\le i\le g$. The flag computation uses bounded sparse elimination. The example has a nonzero pairing between two $\eta$ basis classes, and checks semilinearity and the Cartier restriction over $\mathbb F_8$. [Artin–Schreier de Rham example](examples/rank-as-derham.asl)

#### Geometric pole orders and Ekedahl–Oort type

`rvASPoleOrders[C];` lists pairs `(order, count)` in increasing order for the geometric poles of a reduced `rvascover`. The calculation separates binary multiplicity layers in the finite denominator with polynomial gcds. The degree of each separated squarefree factor counts all its roots over the algebraic closure. A pole at infinity contributes its own order.

`rvASEO[C];` works for these covers over $\overline{\mathbb F}_2$. For geometric pole orders $m_i$, assemble one ordinary two-dimensional block for each pole after the first and a one-pole block of genus $(m_i-1)/2$ for every pole. The Frobenius and Verschiebung matrices of this direct sum give the canonical flag by taking images under $V$ and preimages under $F$. This uses the geometric decomposition of the Jacobian's 2-torsion in Elkin and Pries, Theorems 1.2–1.3. [Ekedahl–Oort strata of hyperelliptic curves in characteristic 2](https://arxiv.org/abs/1007.1226)

The genus-five curves $y^2+y=x^9+x^5+x^{-1}$ and $y^2+y=x^5+x+x^{-1}+x^{-3}+x^{-5}$ have equal zeta numerators, $p$-rank one, and $a$-number two. Their respective geometric pole orders $(1,9)$ and $(5,5)$ give the final types $[1,1,2,2,3]$ and $[1,1,1,2,3]$. [Geometric Ekedahl–Oort example](examples/rank-as-geometric-eo.asl)

`rvASWordRank[C,[F,V,...]];` computes the rank of a composition on the geometric de Rham module, with the rightmost operator acting first. The empty word acts as the identity. `rvASCommonKer[C,[[F,F],[V],...]];` gives the dimension of the intersection of the kernels, and `rvASSumImages[C,[[F,F],[V],...]];` gives the dimension of the sum of the images. An empty family has the whole space as its common kernel and zero as its image sum. The block matrices have coefficients in $\mathbb F_2$, so the twists by Frobenius and its inverse preserve these kernel and image dimensions after extending scalars.

The stable ranks of powers of either operator give the $p$-rank, and $\dim(\ker F\cap\ker V)$ is the $a$-number. For the two genus-five covers with the same zeta numerator, the ranks of $F^3$ are two and one and the dimensions of $\ker F^2\cap\ker V$ are three and four. [Geometric Frobenius and Verschiebung example](examples/rank-as-words.asl)

#### Ekedahl–Oort final types and stratum closures

`EOType(g)` is the finite type of legal final sequences of length `g`. A literal `eoType{nu1,...,nug}` starts at `nu0 = 0` and increases by zero or one at each position. `rvASEOType[C]` and `rvASOddEOType[C]` give the sequence of a declared characteristic-two or odd-prime Artin–Schreier cover as a value. `eoDim(x);` gives the dimension $\sum_i\nu_i$ of the corresponding stratum in the moduli space $\mathcal A_g$ of principally polarized abelian varieties. `eoWeyl(x);` gives its balanced Weyl permutation, with positions and images numbered from one. [Viehmann and Wedhorn, Appendix A](https://arxiv.org/abs/1011.3230)

`eoBruhat(lower,upper)` compares the sequences coordinatewise. `eoClosure(lower,upper)` tests containment of the lower stratum in the Zariski closure of the upper stratum in $\mathcal A_g$. The closure test searches for a Levi permutation $u$ such that $u^{-1}w_{\mathrm{lower}}\delta(u)\le w_{\mathrm{upper}}$ in the type-C Bruhat order, with $\delta$ conjugation by the longest Levi element. `eoClosureWitness(lower,upper);` prints that permutation or `none` when containment fails. Through genus eight it returns the lexicographically first Levi witness. In larger genus the comparison searches the downward type-C Bruhat interval from the upper permutation, stopping at the lower stratum dimension. The lower permutation is minimal in its left Levi coset, so $\ell(u^{-1}w_{\mathrm{lower}})=\ell(u)+\ell(w_{\mathrm{lower}})$ and multiplying by $\delta(u)$ can decrease length by at most $\ell(u)$. After right multiplication by the longest Levi element, the comparison reduces to matching directed cycles of signed pairs up to rotation with equal signs. This finds witnesses especially efficiently in small codimension. The computation and storage limits apply even if the resulting containment is false. [Wedhorn, Specialization of F-Zips](https://arxiv.org/abs/math/0507175)

The genus-five final types $(0,0,1,2,2)$ and $(1,1,1,1,2)$ have stratum dimensions five and six. The permutation $(1\;2)$ certifies specialization from the latter to the former even though the two types are incomparable in the coordinatewise order. [Ekedahl–Oort closure example](examples/rank-eo-closure.asl)

In genus twelve, the types $(0,0,0,0,1,2,3,4,4,5,5,5)$ and $(0,0,0,1,2,2,3,4,4,4,5,5)$ have dimensions 29 and 30 and neither stratum specializes to the other. A second genus-twelve pair exhibits a codimension-one specialization beyond the coordinatewise order. [Higher-genus EO closure example](examples/rank-eo-interval.asl)

The example with $y^3-y=x^{13}+2x^{11}+x^{10}+x$ and $y^3-y=x^{13}+x^5+x^2$ uses their actual typed final sequences as closure inputs. Their complete zeta numerators agree, and the final strata are incomparable. [Odd-prime EO closure example](examples/rank-as-odd-eo-closure.asl)

`eoFirstSlope(x);` computes Harashita’s sharp lower bound on the first Newton slope of the stratum of a positive-genus final type. The calculation extends the sequence to the full self-dual flag and finds the periodic set of the map that applies $V$ when its image is nonzero and $F^{-1}$ otherwise. The bound is the proportion of periodic indices above the middle of that flag. `eoSlopeCycles(x);` counts cycles in this periodic set. The first slope of an individual point can exceed the bound. [Harashita, Ekedahl–Oort strata and the first Newton slope strata](https://hdl.handle.net/2115/17113)

Four polynomial Artin–Schreier curves of genus eight have final type $[0,1,1,2,2,3,3,4]$. Their first slopes are $1/2$, $1/3$, $3/8$, and $1/4$. The stratum bound is $1/4$. [EO and first Newton slope example](examples/rank-eo-first-slope.asl)

#### Odd-prime Artin–Schreier characters

`rvasodd C[p] {(e0,c0),...};` declares the smooth one-pole cover $y^p-y=\sum_i c_ix^{e_i}$ over $\mathbb F_p$. The exponents are nonnegative and the coefficients are signed integers reduced modulo the odd prime `p`. Each positive exponent $pe$ reduces to $e$ by an Artin–Schreier change of variable, and terms with equal exponents combine. A constant term survives the reduction and can twist the arithmetic over $\mathbb F_p$. `rvASOddReduced[C];` prints the surviving `(exponent, coefficient)` pairs. The reduced degree `d` must be positive. `rvASOddGenus[C];` reports $(p-1)(d-1)/2$, and `rvASOddPrank[C];` reports zero.

`rvASOddCharacterL[C,j];` computes the degree-$(d-1)$ factor of the additive character $a\mapsto\zeta_p^{ja}$, for `1 <= j < p`. Each coefficient array is in the power basis of $\mathbb Z[\zeta_p]$. `rvASOddCoeff[C,j,k]` gives coefficient `k` as a `RVCyclo(p)` value. The computation uses the Euler product to express coefficient `k` as a sum over monic degree-`k` polynomials, counting roots with multiplicity. With the first $\lfloor d/2\rfloor$ reversed coefficients fixed, the rest enter the trace phase linearly. Additive character orthogonality sums those variables exactly, so the search has at most $p^{\lfloor d/2\rfloor}$ choices. `rvASOddLPolynomial[C];` takes the Galois norm of the character factor. `rvASOddPoints[C,m];` obtains the point count in positive extension degree `m` by a character power-trace recurrence. The work, storage, and signed 128-bit coefficient bounds apply, and point counts use unsigned 128-bit integers. [Lauder and Wan, Computing Zeta Functions of Artin–Schreier Curves over Finite Fields](https://doi.org/10.1112/S1461157000000681)

`rvASOddCharacterNewton[C,j];` computes the slopes of one character factor with $v_p(p)=1$ using the uniformizer $1-\zeta_p$. `rvASOddNewton[C];` computes the slopes of the full numerator. In characteristic five, the factors for $x^6$, $x^6+x^3$, and $x^6+x^4$ have first slopes $1/2$, $3/8$, and $1/4$, respectively. The example also checks signed normalization and a constant twist on a typed coefficient and on point counts. [Odd-prime Artin–Schreier example](examples/rank-as-odd.asl)

`ASOddDiff(C)` is the space of holomorphic differentials of an odd-prime cover over its ground field $\mathbb F_p$. Its basis consists of $\omega(b,r)=x^by^r\,dx$ with $0\le r\le p-2$ and $pb+dr\le 2g-2$. Use `rvASOddBasis[C];` for the `(b,r)` pairs ordered by `r` then `b`, and `rvASOddOmega[C,b,r]` for a basis differential. The space supports `zero`, `add`, and `pscale`.

`rvASOddCartier(w)` and `rvASOddCartierIter[k](w)` apply Cartier and its $k$th iterate. The column for $x^by^r\,dx$ comes from $y^r=(y^p-f)^r$ and the terms $x^{pn+p-1}\,dx$ in each coefficient of $y^{pj}$. Over $\mathbb F_p$ the map is linear. `rvASOddCartierRank[C,k];` gives the rank of the iterate by sparse reduction, and `rvASOddaNumber[C];` gives the genus minus the rank of Cartier. `rvASOddTranslate[c](w)` is the pullback under $y\mapsto y+c$, with signed `c` reduced modulo `p`. It commutes with Cartier and its `p`th iterate is the identity.

At $p=5$, the genus-ten curves with $f=x^6+x^3$ and $f=x^6+x^4$ both have $a$-number six. Their Cartier rank sequences are $[10,4,1,0]$ and $[10,4,2,1,0]$. The curve with $f=x^6$ has Cartier rank zero. [Odd-prime Cartier–Manin example](examples/rank-as-odd-cartier.asl)

`ASOddDeRham(C)` is the first de Rham cohomology of an odd-prime cover over $\mathbb F_p$. The first $g$ basis vectors are the regular differentials $\alpha_{b,r}=x^by^r\,dx$ in `rvASOddBasis` order. The next $g$ are their Čech lifts $\eta_{b,r}$ with overlap function $h=x^{-b-1}y^{p-1-r}$. The two forms of a lift satisfy $U-V=dh$. The affine form $U$ retains the nonnegative $x$-powers of $dh$. Use `rvASOddDeRhamGen[C,i]` with a zero-based index, or `rvASOddDeRhamIncl(w)` to include an `ASOddDiff` differential. Use `zero`, `add`, and `pscale` for vector arithmetic.

`rvASOddF(w)` raises the overlap function to the $p$th power and reduces the resulting Čech class. `rvASOddV(w)` applies Cartier to the two chart forms. Over $\mathbb F_p$ both maps are linear. `rvASOddDeRhamTranslate[c](w)` pulls back the three Čech entries along $y\mapsto y+c$, commutes with both operators, and preserves the cup form. `rvASOddDeRhamPair(w,z)` computes the signed cup pairing from the residue of $U_w h_z-h_w V_z$. Its value is a `PrimeScalar(p)`, whose literals `rvPrimeElt[p,c]` reduce a signed integer `c` modulo `p`. `rvASOddGram[C];` reports the pairing matrix in the same basis order.

`rvASOddEO[C];` computes the geometric final type by closing the zero and full spaces under $V$ and $F^{-1}$ over $\mathbb F_p$. The matrices and pairing use the actual reduced equation. At $p=5$ the genus-ten curves with $f=x^6+x^3$ and $f=x^6+x^4$ have equal $a$-number six and different final types. The example also shows a nonzero signed pairing of Čech lifts in characteristic three. [Odd-prime Čech de Rham example](examples/rank-as-odd-derham.asl)

Over $\mathbb F_3$, $f=x^{13}+2x^{11}+x^{10}+x$ and $f=x^{13}+x^5+x^2$ give genus-twelve curves with the same pole order, the same full primitive character polynomial, and the same $a$-number seven. Their final types differ. The example checks all thirteen primitive coefficients and reports the common zeta numerator alongside the two final types. [Equal zeta factors and distinct final types](examples/rank-as-odd-stationary.asl)

#### Tate cohomology of the wild deck action

Let $T$ be pullback under $y\mapsto y+1$ on `ASOddDeRham(C)`, with $\Delta=T-1$ and $N=\sum_{i=0}^{p-1}T^i=\Delta^{p-1}$. The two-periodic Tate complex has differential $\Delta$ in even cohomological degrees and $N$ in odd degrees. Use `rvASOddTateD[z](w)` on a de Rham vector `w` to apply that differential. The type `ASOddTate(C,z)` is its cohomology in signed degree `z`, with the usual vector operations `zero`, `add`, and `pscale`. Use `rvASOddTatedim[C,z];` for the dimension and `rvASOddTategen[C,z,k]` for a zero-based generator. Use `rvASOddTateClass[z](w)` to project a cycle and `rvASOddTateRep(x)` to obtain the representative selected by sparse elimination. Distinct signed degrees have distinct types even when their matrices coincide.

`rvASOddTateF(x)` and `rvASOddTateV(x)` use the induced Frobenius and Verschiebung. The polarization descends to a perfect pairing across the two parities, returned as `PrimeScalar(p)` by `rvASOddTatePair(x,y)` on classes of the same cover and opposite parity. The adjoint of $\Delta$ is $T^{-1}-1$, with the same kernel and image as $\Delta$, and $N$ is self-adjoint. `rvASOddTateJordan[C];` gives the descending block lengths of $\Delta$ on de Rham cohomology, including blocks of size $p$ when present. Such blocks contribute zero in both Tate parities.

For $y^3-y=x^7$, all six deck blocks have size two. Each Tate parity has dimension six, and the induced $F$ and $V$ each have rank two and square zero. The image of $F$ has dimension two while the kernel of $V$ has dimension four. The example certifies a class in $\ker V$ outside $\operatorname{im} F$ by its nonzero pairing with a class of opposite parity killed by $V$. [Wild deck Tate example](examples/rank-as-odd-tate.asl)

#### Polarized Dieudonné lattices and Newton polygons

`rvcrystal C[p,g] {(row,column,entry),...};` declares a Frobenius matrix in the standard symplectic basis $e_1,\ldots,e_g,f_1,\ldots,f_g$. Indices are zero-based, coefficients are signed integers, and omitted entries are zero. The lattice is over $W(\mathbb F_p)=\mathbb Z_p$, whose Witt Frobenius is the identity. With $J=\begin{psmallmatrix}0&I\\-I&0\end{psmallmatrix}$, the declaration verifies $F^t J F=pJ$ and computes $V=J^{-1}F^tJ$. This gives $VF=pI$, and over $\mathbb Q_p$ the left inverse $p^{-1}V$ is also the right inverse, so $FV=pI$ integrally.

`rvCrystalChar[C];` returns the integral coefficients of $\det(TI-F)$ in descending degree. `rvCrystalNewton[C];` computes slopes with multiplicities from the reciprocal polynomial using $v_p(p)=1$. `rvCrystalEOType[C]` computes the typed final sequence of the mod-$p$ module by closing its canonical flag under $V$ and $F^{-1}$. `rvCrystalMatrices[C,h];` gives the actual $F$ and $V$ matrices with least nonnegative entries modulo $p^h$. `rvCrystalCongruent[A,B,h]` compares both matrices entrywise at positive level $h$ for a common prime and genus. A true result gives a certificate in the specified symplectic bases. A symplectic change of basis can change the entrywise result. Matrix arithmetic uses checked signed 128-bit integers within the statement work and storage budgets. Congruence uses valuations of differences and works at levels for which the matrix reporting modulus is out of range.

For the cyclic signed shift $S$ sending $e_g\mapsto f_1$ and $f_g\mapsto-e_1$, set $F=S\operatorname{diag}(1_g,p_g)$. Left composition with the transvection $f_1\mapsto f_1+c e_1$ gives characteristic polynomial $T^{2g}-cT^g+p^g$. For $g=2h+1$ and $c=p^h$, the slopes are $h/g$ and $(h+1)/g$, each with multiplicity $g$. The $c=0$ lattice has all slopes $1/2$. The two lattices agree modulo $p^h$ in the named basis. The example takes $g=7$, $p=3$ and compares $c=0$ with $c=27$. Its full $F,V$ matrices agree modulo 27 and its final sequences are equal. Its Newton slopes are $1/2$ repeated 14 times in the first lattice and $3/7,4/7$ repeated seven times each in the second lattice. [Polarized Dieudonné lattice example](examples/rank-dieudonne-crystal.asl)

#### Split filtered $(\varphi,N)$-modules

`rvphodge D[p] { {(num,den),...}, {(row,column,num,den),...}, {(weight,[entries]),...} };` declares a filtered $(\varphi,N)$-module over $\mathbb Q_p$. The first list gives distinct nonzero rational eigenvalues $\lambda_i$ of $\varphi$ in the chosen eigenbasis. The length of this list determines the dimension, and each denominator is a positive integer. The second list gives rational entries of $N$ at zero-based positions, with omitted entries zero. The declaration checks $N\varphi=p\varphi N$. The third list gives integer columns forming a basis over $\mathbb Q_p$. A column annotated by weight $w$ belongs to $\mathrm{Fil}^j$ exactly for thresholds $j\le w$ in this adapted basis. The columns of weight at least $j$ span $\mathrm{Fil}^j$.

`rvPHodgeData[D];` gives the full Hodge and Newton weights, their sums $t_H,t_N$, and the degree $t_H-t_N$. Use `rvPHodgeData[D,[i0,...]];` for the induced filtration on a stable coordinate subspace. Projection of the flag columns to the complementary coordinates computes the intersection ranks over characteristic zero. This permits flags whose integer basis becomes singular modulo $p$.

The simple split eigenvalues force every $\varphi$-stable subspace to be a coordinate subspace. Nonzero entries of $N$ form disjoint chains, and a $\varphi,N$-stable subspace chooses a prefix of each chain, starting at its lowest slope. `rvPHodgeWeak[D];` checks equality on the full module and the condition $t_H\le t_N$ on every stable subspace. It reports the full equality defect or an explicit destabilizing subspace when a condition fails. `rvPHodgeHN[D];` gives the Harder–Narasimhan flag and the exact slopes of its successive quotients. At each step it selects the greatest stable subspace of maximum degree divided by dimension in the quotient. The Hodge degree is supermodular under sums and intersections, so the union of maximal-slope subspaces again has maximal slope and gives the greatest choice. All work, storage, and checked integer bounds apply to the exact computation.

The example contrasts a weakly admissible monodromy chain with the same filtration after deleting $N$ and with the same monodromy after moving its high plane. A noncoordinate flag has first stable term $\langle e_0,e_1\rangle$ and slopes $1/2,-1$. [Filtered Frobenius and monodromy example](examples/rv-filtered-phi-n.asl)

#### A chosen split quadratic embedding

`rvqpadic K[p,d,r];` chooses the embedding of $K=\mathbb Q(\sqrt d)$ into $\mathbb Q_p$ whose square root is congruent to $r$. The input requires an odd prime, a nonsquare integer $d$, and a nonzero root $r^2=d$ modulo $p$. `rvQpQuad[K,a,b,c]` is $(a+b\sqrt d)/c$ with $c>0$, as a value of type `RVQpQuad(K)`. Values support `zero`, `add`, `neg`, equality, `rvQpQuadMul`, `rvQpQuadConj`, `rvQpQuadInv`, and `rvQpQuadPow[k]` for nonnegative $k$. The embedding remains fixed, so conjugation can change the valuation. Distinct embedding declarations have distinct types.

`rvQpQuadVal(x);` takes the valuation at the chosen prime, returning `infinity` for zero. After removing a common $p$-power from the two numerator coefficients, a numerator that vanishes at the chosen root has unit conjugate. The norm then gives the remaining valuation. `rvQpQuadResidue[s](x);` gives the least nonnegative residue modulo $p^s$ for a $p$-adically integral value, including one with a written denominator divisible by $p$. It strips the powers of $p$ first. If a primitive numerator vanishes at the chosen root, its unit part is the stripped norm divided by its unit conjugate. This needs at most the requested Hensel precision, including for denominators with large powers of $p$. Stored canonical coefficients and the output modulus $p^s$ must fit in signed 128 bits. Arithmetic uses bounded wider intermediates for products, norms, and normalization, and narrows the coefficients only after cancellation.

`rvphodgequad M[K] {{eigenvalues},{(row,column,entry),...},{(weight,[column entries]),...}};` takes values of type `RVQpQuad(K)` for each eigenvalue, monodromy entry, and flag entry. The Frobenius spectrum must again be simple and split. The same `rvPHodgeData`, `rvPHodgeWeak`, and `rvPHodgeHN` queries use exact characteristic-zero flag ranks over $K$ and Newton slopes at the chosen embedding. In the example, the Frobenius roots of the Legendre curve with parameter $2$ over $\mathbb F_5$ have valuations $0$ and $1$ at $i\equiv2\bmod5$. The unit residue modulo $5^8$ is $169738$, also computed independently by `legendreDwork[5,2,8];`. Their valuations exchange at the other embedding. Weight-one lines $\langle e_\alpha\rangle$ and $\langle e_\alpha+\beta^8e_\beta\rangle$ agree modulo $5^8$. The first fails weak admissibility and the second satisfies it. These are chosen filtrations on the given Frobenius data. [Quadratic filtered Frobenius example](examples/rv-quadratic-phodge.asl)

At the root of $2$ congruent to $3$ modulo $7$, $(3-\sqrt{2})^{31}$ has valuation $31$ and angular component $6$. Its two squared norm terms exceed the coefficient range individually. Their difference is $7^{31}$, and its inverse is represented exactly. The example also computes residues of $(3-\sqrt{2})^{31}/7^{31}$ without forming a modulus of precision $s+31$. [Wide quadratic arithmetic example](examples/rv-quadratic-wide.asl)

#### Tropical Grassmannians of p-adic flags

`rvValMat[p,n]{[row0],...}` takes independent integer rows of width `n`. Any rational row space can be represented after clearing each row's denominators. For a declared split quadratic embedding `K`, `rvValMatQuad[K,n]{[row0],...}` takes rows with entries of type `RVQpQuad(K)`. Use `rvValFil[D,j]` for the actual subspace $\mathrm{Fil}^j$ of a declared rational or quadratic filtered $(\varphi,N)$-module. It uses the stored flag in its given ambient coordinates and the module's declared prime or embedding.

For $r$ rows, the result has type `RVValMat(p,n,r)`. Its coordinates are $v(I)=v_p(\det M_I)-\min_J v_p(\det M_J)$ for $r$-element subsets of the zero-based ground set. Only finite coordinates enter the minimum. Exact zero minors have valuation `infinity`. Rational minors use fraction-free elimination, and quadratic minors use exact elimination over the declared field and valuations of its pivots. If a pivot inverse or an intermediate product exceeds the narrow coefficient range, ASL clears denominators and recomputes that minor by fraction-free elimination on bounded wide integers. This permits exact cancellation before the determinant is valued and preserves exact zero support. [Wide tropical minors example](examples/rv-tropical-wide.asl) Different bases of the row space give the same normalized valuated matroid. These values can be named, compared, and used in pairs and as function results. Rank zero and full rank each have a single value at fixed $p,n$. Other ranks use the native finite-value limits for explicit support and cannot be enumerated as types.

`rvValPlucker[[I]](M);` reads one normalized coordinate. `rvValSupport(M);` lists the finite bases in lexicographic order. `rvValInitial[[w0,...]](M);` gives the minimum of $v(I)-\sum_{i\in I}w_i$ and every basis attaining it. Weights are signed 64-bit integers, and costs use signed 128-bit arithmetic. `rvValDual(M)` complements the finite bases without changing their valuations. Complementing the initial bases at $w$ gives the dual's initial bases at $-w$.

For finite $I,J$ and $i\in I\setminus J$, `rvValExchange[[I],[J],i](M);` returns the first $j\in J\setminus I$ that makes both exchanged bases finite and does not increase the sum of their valuations. The report gives $j$, the old sum and the new sum. `rvValTie[[A],[B]](M);` uses $|A|=r-1$ and $|B|=r+1$. It reports the minimum of $v(A\cup\{j\})+v(B\setminus\{j\})$ and its attaining indices $j\in B\setminus A$. A finite minimum is attained at least twice. If every term is infinite, all these indices attain the minimum.

At the embedding $i\equiv2\pmod5$, put $\beta=-1-2i$. The example uses four columns $(1,0),(0,1),(1,1),(1,1+\beta^{19})$. The last minor is zero modulo $5^{19}$ and has normalized valuation $19$, while the other five minors are units. All six are finite. At zero weight only five are initial. Assigning weight $19$ to each of the last two columns makes their pair initial and removes the first pair. The example identifies this matroid with the same flag step in a declared module. [Tropical p-adic flag example](examples/rv-tropical-flag.asl)

#### Weighted residue subspaces

`RVResSpace(p,n)` is the type of all subspaces of $\mathbb F_p^n$. Values are displayed as row spaces in canonical reduced row-echelon form. Use `rvResSpace[p,n]{[row],...}` to span signed integer vectors modulo $p$, or `zero[RVResSpace(p,n)]` for the zero space. Equality compares subspaces. `rvResDim(U);` gives the dimension, `rvResContains[[vector]](U);` tests membership, and `rvResLeq(U,V);` tests inclusion. `rvResSum(U,V)`, `rvResIntersection(U,V)`, and `rvResOrthogonal(U)` use the same ambient type, with orthogonality from the standard dot product. Small finite types can be enumerated within the usual limit. `rvResSupport(U);` lists the bases of the row-space matroid subject to the Plücker-minor limits.

`rvResSpecialize[p,n,[w...]]{[row],...}` takes an actual rational row space, with entries written as integers or pairs `(numerator,denominator)` of positive denominator. `rvResSpecializeQuad[K,n,[w...]]{[row],...}` takes entries of type `RVQpQuad(K)` and uses its selected embedding. `rvResFil[M,j,[w...]]` specializes the actual $\operatorname{Fil}^j$ of a declared filtered module. Each operation scales ambient coordinate $i$ by $p^{-w_i}$ and takes the image of the intersection with $\mathbb Z_p^n$ in $\mathbb F_p^n$. Exact valuation pivoting selects a basis of minimum $v(I)-\sum_{i\in I}w_i$. Its identity chart has integral scaled entries. Only entries with adjusted valuation zero use an angular component, so even large differences of weights require no large powers of $p$. The Plücker support of the reduction is the initial matroid. An individual line is reduced by valuation and angular component ratios. For larger charts with intermediate field coefficients beyond the scalar range, rectangular fraction-free elimination gives exact echelon rows of the weighted lattice. Each row selects a least adjusted valuation and reduces unit ratios at that pivot, including when some input rows are redundant. The fallback retains bounded 2048-bit coefficients and widens its temporary exact products and divisions. [Wide residue-echelon example](examples/rv-residue-wide-echelon.asl)

The example gives two planes with identical tropical Plücker data over $\mathbb Q_7$ and different residue planes over $\mathbb F_7$. For $\beta=-1-2i$ at $i\equiv2\bmod5$, it also exhibits a change of residue support as the weights of the last two coordinates pass the cancellation depth of $\beta^8$. [Weighted residue example](examples/rv-residue-specialization.asl)

#### Stable reduction of marked projective lines

Apply `rvStableTree(M)` to a rank-two `RVValMat` in which every pair has finite valuation, with at least three marks. The result has type `RVStableTree(n)`. Take the last mark as the anchor $*$. For other marks form $c_{ij}=v_{ij}-v_{i*}-v_{j*}$. At a cluster of height $h=\min c_{ij}$, pairs with $c_{ij}>h$ belong to the same child. Each child cluster with at least two marks gives an internal edge of length equal to the increase in height. This constructs a tree with every internal vertex at least trivalent. Row basis changes and independent rescalings of the marked columns add a constant to all the $c_{ij}$ and preserve the edge lengths. `rvStableSplits(T);` lists the sides and length of each internal edge. Both sides are sorted and compared lexicographically.

`rvStableLine[p]{[x0,y0],...}` starts from distinct ordered points $[x_i:y_i]$ with integer homogeneous coordinates. It permits rational points after clearing denominators and includes $[1:0]$. `rvStableLineQuad[K]{[x0,y0],...}` takes coordinates in the declared split quadratic field `K`. These operations return `RVStableMarks(p,n)`, which stores the tree and the marked residue configuration of each component. `rvStableTree(C)` projects to the tree. Equality of decorated reductions compares lengths, branches, and their normalized residues. The interpreter forms a projective bracket as one signed determinant before taking its valuation and angular component. A large coordinate change can cancel in the determinant, and a bracket with oversized field coefficients can still give its exact local data. [Projective gauge example](examples/rv-projective-gauge.asl)

Vertex zero is adjacent to the anchor. Number the other vertices in preorder, visiting children in order of their least leaf label. `rvStableDirections[v](T);` lists the branches at vertex `v` in order of their least mark. Choose the least mark in each of the first three branches as representatives $A,B,C$. They are sent to $\infty,0,1$ by the coordinate $[i,B][C,A]/([i,A][C,B])$, where $[i,j]=\det(P_i,P_j)$ is oriented. Reduction uses exact valuations and angular components. It verifies that all marks in a branch have the same reduction and that different branches give distinct points. `rvStableComponents(C);` lists branches and residues for every vertex. `rvStableResidue[v,i](C);` reads the residue of the branch containing the mark `i` at `v`.

The two nine-mark configurations over $\mathbb Q_7$ in the example have the same marked splits and edge lengths $7,11,19$. At one four-valent component the fourth branch has residue $2$ in the first configuration and $5$ in the second. The residue cross-ratio distinguishes them although their metric trees agree. [Stable marked reduction example](examples/rv-stable-marks.asl)

#### Hyperelliptic metric covers and component groups

`rvHyperCover(C)` takes an `RVStableMarks(p,n)` with odd $p$ and $n=2g+2\ge4$. It forms the geometric tame double cover of the stable branch tree. A direction is ramified when it contains an odd number of branch marks. A base vertex with no ramified directions has two genus-zero lifts exchanged by the deck involution. A vertex with $b>0$ ramified directions has one lift of genus $b/2-1$. An even split lifts to two edges of the original length and dilation one. An odd split lifts to a single edge of half that length and dilation two.

`rvHyperGraph(H);` lists vertices as `(base,genus,local_degree,deck_partner)`, finite edges as `(source,target,base_edge,(length_numerator,length_denominator),dilation,deck_partner)`, and branch legs as `(mark,incident_vertex,dilation)`. `rvHyperHarmonic(H);` checks the sum of the outgoing dilations in every base direction at every cover vertex. `rvHyperGenus(H);` reports the genus, toric rank and vertex genera. Equality of `RVHyperCover(p,n)` values compares the labeled metric cover with vertex genera. Its scope is the geometric reduction graph over a suitable extension, with positive-genus components recorded by genus.

`rvHyperCycles(H);` returns an integral basis of cycles with coordinates `(edge,sign)` in the listed edge orientations. `rvHyperMetric[e](H);` multiplies every length by the positive ramification multiplier $e$ and requires integral lengths, including on separating edges. `rvHyperMonodromy[e](H);` gives $Q=C\,\operatorname{diag}(\ell)\,C^{\mathsf T}$ for the cycle matrix $C$ and the resulting metric $\ell$. `rvHyperComponent[e](H);` gives the nontrivial Smith factors of $Q$, their largest factor (the exponent) and their product (the order) for the geometric component group in that integral metric. It computes the determinant exactly and uses it to bound Smith reduction in the finite quotient. The multiplier defaults to one. The monodromy pairing can change basis under relabeling, and its Smith factors are invariant.

The three pairs of branch marks in the genus-two example produce a theta graph whose paths have lengths $4,6,10$ after suppressing bivalent vertices. Its component group has factors $2,62$ and order $124$. Doubling the metric gives factors $4,124$ and order $496$. A second example has two genus-one vertices and toric rank one. Its separating odd split requires an even multiplier to give integral lengths. [Hyperelliptic metric cover example](examples/rv-hyperelliptic-graph.asl)

#### Arithmetic Frobenius and rational hyperelliptic component groups

`rvHyperArithmetic[p,c]{r0,...}` takes distinct, $p$-integral rational branch roots and a rational leading unit $c$ at odd $p$ for $y^2=c\prod_i(x-r_i)$. Rationals may be written as integers or pairs `(numerator,denominator)` with positive denominator. The affine root configuration must occupy at least three residues at depth zero. For each proper cluster with center $z$ and depth $d$, the normalized chart requires $\nu=\sum_i\min(d,v_p(z-r_i))$ even. The semistable metric must already have integral lifted lengths at multiplier one. This gives a split-branch semistable family over $\mathbb Q_p$, including nonsquare unit twists.

`rvHyperGeometric(A)` forgets the Frobenius action. For an even cluster $S$ the two edges are exchanged by arithmetic Frobenius when the unit part of $c\prod_{r_i\notin S}(z-r_i)$ is a nonsquare. At an unbranched vertex the same square class of the normalized chart exchanges its two components. The root chart uses $c$. `rvHyperFrobenius(A);` lists the Frobenius permutation of the geometric vertices and finite edges, its integral matrix $F$ in the chosen cycle basis, and its dual matrix $F^{\mathsf T}$, since $F$ is an involution. The calculation verifies the incidence map and $F^{\mathsf T}QF=Q$. Equality compares the labeled cover, its Frobenius permutations and its canonical normalized residual charts. In each chart the first two finite residue directions by mark order fix coordinates $0,1$. The leading unit is normalized to its square class after this coordinate change.

`rvHyperFixed(A);` gives the invariant factors, order and actual dual-cycle representatives of the rational component group $(\mathbb Z^b/Q\mathbb Z^b)^{F^{\mathsf T}}$. It uses the perfect pairing $x^{\mathsf T}Q^{-1}y$ to identify this group with the characters of $\mathbb Z^b/(Q\mathbb Z^b+(F^{\mathsf T}-I)\mathbb Z^b)$. The determinant of $Q$ kills this quotient, so Smith operations on the relations and character coordinates can use residues modulo that determinant. For each Smith factor $d$, it chooses a character $y/d$ with $0\le y_i<d$ and returns the integral vector $Qy/d$. The displayed order is the Tamagawa number. In the $p=7$ theta example the geometric order is $124$ and the nonsquare twist has rational order $4$. A second pair has both signs on the two cycles, and the rational orders are $12$ and $20$ even though both geometric orders are $60$. [Arithmetic hyperelliptic example](examples/rv-hyperelliptic-arithmetic.asl)

#### Semistable hyperelliptic local factors

`rvHyperLocal(A);` gives the polynomial in $T$ whose reciprocal is the local Euler factor on inertia-invariant $H^1$. The graph part is $\det(1-TF\mid H^1(\mathrm{graph}))$. The abelian part multiplies the Weil polynomials of the normalized positive-genus components. The combined degree is $2g-b$, where the graph cycle rank $b$ is also the semistable conductor exponent. The local root number is $(-1)^s$ with $s=\dim H^1(\mathrm{graph})^{F=1}$.

`rvHyperLocalFactors(A);` lists each canonical normalized chart as `(base,genus,coefficients,traces,Weil_polynomial)`, with coefficients in increasing degree and traces over $\mathbb F_{p^r}$ for $1\le r\le g_v$. A genus of `None` denotes two geometric rational components and their displayed polynomial is the constant leading unit. A connected normalized chart has a factor for each residue direction that contains an odd number of roots. It uses the same normalization at infinity as the curve equation. An odd degree gives one point at infinity. An even degree gives two or zero according to the square class of the leading unit over the counting field. Finite-field quadratic character sums compute the traces before the integral Newton identities and functional equation fill in the Weil polynomial.

`rvHyperNodalPoints[r](A);` counts points on the nodal curve given by the admissible cover. It adds points on the Frobenius-fixed normalized components and subtracts the number of fixed nodes. The report also gives these vertex and edge counts and the trace $1+p^r V_r-N_r$, an independent coefficient check on the local factor for the chosen graph.

At $p=7$ one genus-three example has two elliptic components with traces $-4,0$ and a split cycle. Its local denominator is $(1-T)(1+4T+7T^2)(1+7T^2)$. A nonsquare twist gives $(1+T)(1-4T+7T^2)(1+7T^2)$, reversing the root number. The count over $\mathbb F_{7^5}$ is $50185$ for the original and $50659$ for the twist. The mixed-sign genus-two example has the same factor $1-T^2$ for either twist, alongside different rational component groups. [Semistable local factor example](examples/rv-hyperelliptic-local.asl)




#### Two-pole odd-prime Artin–Schreier curves

`rvastwo C[p] {(signedExponent,signedCoefficient),...};` declares $y^p-y=f_+(x)+f_-(1/x)+c$ over $\mathbb F_p$ with poles at zero and infinity. Each nonzero Laurent exponent $pm$ reduces to $m$, and equal exponents combine modulo $p$. The constant survives reduction. Use `rvASTwoReduced[C];` to inspect the Laurent polynomial. If $d$ and $e$ are its positive pole orders, `rvASTwoGenus[C];` returns $(p-1)(d+e)/2$ and `rvASTwoPrank[C];` returns $p-1$.

`rvASTwoCharacterL[C,j];` gives the degree-$d+e$ character factor on $\mathbb G_m$ for $1\le j<p$. `rvASTwoCoeff[C,j,k]` returns its $k$th coefficient as a `RVCyclo(p)` value. For a monic polynomial of degree $k$ with nonzero constant $a_k$, the Euler sum uses the power sums of its roots and their inverses, counting multiplicity. The computation fixes $a_k$ and the union of the first $\lfloor d/2\rfloor$ and the last $\lfloor e/2\rfloor$ interior coefficients. The remaining coefficients enter the trace phase linearly in the polynomial and its normalized reversal. Additive character orthogonality sums them exactly. A constant twist multiplies coefficient $k$ of character $j$ by $\zeta_p^{jkc}$.

`rvASTwoLPolynomial[C];` gives the full cyclotomic norm. `rvASTwoCharacterNewton[C,j];` and `rvASTwoNewton[C];` give the character and full Newton slopes, normalized by $v_p(p)=1$. `rvASTwoPoints[C,m];` evaluates the point count in positive extension degree $m$ by a power-trace recurrence, including the two rational ramified points. All computations use the existing checked work, storage, and integer bounds.

Over $\mathbb F_5$, $f=x+x^{-1}$ gives the quadratic Kloosterman character factor $1+K_jT+5T^2$ and twelve rational points. The genus-twelve curve with $f=x^4+x^{-2}$ has primitive character slopes $0,1/4,1/2,1/2,3/4,1$ and has 16,172 points over $\mathbb F_{5^6}$. For $f=x^6+x^{-1}$, adding $x^3$ or $x^4$ preserves the pole orders, genus fourteen and $p$-rank four. The first strictly positive primitive slopes of the three curves are $1/2$, $3/8$ and $1/4$. The example also checks exchange of the poles and rescaling of $x$. [Two-pole Artin–Schreier example](examples/rank-as-two-pole.asl)

#### Cartier and the two-pole ordinary deck representation

`ASTwoDiff(C)` is the space of regular differentials of a declared two-pole cover over $\mathbb F_p$. Its basis consists of $x^iy^j\,dx$ with signed $i$, $0\le j\le p-2$, and
$pi-ej+(p-1)(e+1)\ge0$, $-pi-dj+(p-1)(d+1)-2p\ge0$.
Use `rvASTwoBasis[C];` for the ordered $(i,j)$ pairs and `rvASTwoOmega[C,i,j]` for a basis vector. The space supports `zero`, `add` and `pscale`.

`rvASTwoCartier(w)` and `rvASTwoCartierIter[k](w)` apply Cartier and its $k$th iterate via $y^j=(y^p-f)^j$. The Laurent terms selected by Cartier have exponent $m\equiv-1\pmod p$ and map to $x^{(m+1)/p-1}\,dx$. `rvASTwoTranslate[c](w)` applies the pullback $y\mapsto y+c$. The operators are $\mathbb F_p$-linear and commute. Inversion of the base coordinate sends $x^iy^j\,dx$ to $-x^{-i-2}y^j\,dx$ and exchanges the pole orders.

`rvASTwoCartierRank[C,k];` reports the rank of $C^k$ and `rvASTwoaNumber[C];` gives $g-\operatorname{rank}(C)$. `rvASTwoStableRank[C];` and `rvASTwoOrdinaryGen[C,k]` give the dimension and zero-based ambient basis vectors of the stable Cartier image. This dimension is $p-1$. The deck difference $T-1$, for $T$ the translation by one, has one block there: $(T-1)^{p-2} C^g(y^{p-2}\,dx/x)=(p-2)!\,dx/x$. Cartier fixes $dx/x$ and commutes with this block, so $C-1$ is nilpotent on the stable image. `rvASTwoOrdinaryRanks[C];` computes the rank ladders of both differences. `rvASTwoHasseWitt[C];` gives the coefficients of $\det(I-zC)$ in ascending degree over $\mathbb F_p$, using the computed nilpotence on the stable image.

Over $\mathbb F_5$, the curves with $f=x+x^{-1}$ and $f=x+x^{-1}+1$ have slopes $0$ and $1$, each with multiplicity four, and a single ordinary deck block of length four. Their ordinary $C-1$ has block lengths $(2,2)$ and $(4)$, respectively. The genus-twelve curve $f=x^4+x^{-2}$ has Cartier ranks $12,6,5,4$ before stabilization and ordinary $C-1$ ranks $4,1,0$. [Two-pole Cartier and ordinary deck example](examples/rank-as-two-cartier.asl)

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

`rvWittTowerCharacterNewton[W,j];` reports the Newton slopes of one primitive character factor, normalized by $v_2(2)=1$. The local computation uses the uniformizer $\zeta_{2^j}-1$ and keeps the actual slope multiplicities before taking the cyclotomic norm. For the height-four tower $(x,0,0,0)$, the order-16 factor has slopes $1/8,2/8,\ldots,7/8$, each once.

`rvWittTowerCharacterBaseL[W,j,m];` replaces each Frobenius reciprocal root $\alpha$ of that factor by $\alpha^m$ for `m >= 1`, returning the coefficients in the same cyclotomic power basis. `rvWittTowerFrobPoints[W,m];` recovers the point count over $\mathbb F_{2^m}$ from the traces of those powers for every conductor. These queries first compute the character factors over $\mathbb F_2$. They use the existing signed 128-bit coefficient, storage, and work bounds, and the trace point count allows extension degrees up to 127. The example uses the trace recurrence for point counts through extension degree 96. [Character Newton polygons and Frobenius example](examples/rank-witt-frobenius.asl)

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

#### Untruncated cyclic homology and Connes SBI

`RVHNChain(J,b,[mu])`, `RVHPChain(J,b,[mu])` and `RVHCChain(J,b,[mu])` are the negative, periodic and ordinary cyclic chain groups. They use $C[[u]]$, $C((u))$ and $C((u))/uC[[u]]$, respectively, for the normalized Hochschild complex $C$. Their allowed powers of $u$ are nonnegative, arbitrary integers and nonpositive. The differential `rvCyclicD` is $d+uB$. In the ordinary cyclic quotient it discards the $uB$ term of the constant coefficient.

For fixed total degree $b$ and weight $\mu$, only powers $j$ with $0\le b+2j\le2|\mu|$ can occur. The interpreter intersects this range with the allowed exponents before building the degree slice. Use `rvHNTerm[j]`, `rvHPTerm[j]` or `rvHCTerm[j]` to insert a Hochschild chain in coefficient $j$, and `rvCyclicPart[j]` to extract it.

The corresponding homology types are `RVHN`, `RVHP` and `RVHC`, with the same arguments. `rvCyclicClass` reduces a cycle and `rvCyclicRep` gives a chosen cycle representative. The dimension queries are `rvHNdim`, `rvHPdim` and `rvHCdim`. The matching generator queries end in `gen` and append a zero-based basis index. Multiplication by $u^k$ on periodic chains and classes is `rvCyclicU[k]` for signed $k$. Its inverse is `rvCyclicU[-k]`.

The Connes sequence is
$\cdots\longrightarrow HH_b\overset{I}{\longrightarrow}HC_b\overset{S}{\longrightarrow}HC_{b-2}\overset{B}{\longrightarrow}HH_{b-1}\longrightarrow\cdots$.
Use `rvCyclicI` to include the constant coefficient and `rvCyclicS` for multiplication by $u$ in the ordinary cyclic quotient. The connecting map `rvCyclicB` applies the Hochschild Connes operator to the constant coefficient of an ordinary cyclic chain and takes the resulting Hochschild class. All three maps also act on the underlying chains.

For one exterior generator in odd positive weight $n$, the periodic groups vanish. The negative group is generated by $1[e|\cdots|e]$ in degree $2n$ and the ordinary cyclic group by $e[e|\cdots|e]$ in degree $2n-1$. The connecting map sends the latter to the upper Hochschild class. In even positive weight, each periodic degree has rank one. In weight zero the ordinary, negative and periodic groups occur in even nonnegative degrees, even nonpositive degrees and all even degrees, respectively. [Connes SBI example](examples/rank-cyclic.asl)

### Declared mixed complexes and their u-adic pages

`mixedcomplex C[p] { ... } { ... } { ... };` declares a finite graded vector space over $\mathbb F_p$ and two operators $b$ and $B$ of homological degrees $-1$ and $+1$. List generators and signed degrees first, with entries such as `x:0, y:1`. The second and third blocks list the columns of $b$ and $B$. A column such as `x => [(y,-1)]` sends $x$ to $-y$. Omitted columns vanish, and coefficients are reduced modulo $p$. The interpreter checks homogeneity and the three identities $b^2=B^2=bB+Bb=0$ on every generator.

`MixedChain(C,b)` is the underlying degree-$b$ chain space. Construct its elements with `mixchain[C,b][(generator,coefficient),...]` and apply the two operators with `mixb` and `mixB`. The total degree-$b$ space of $C[[u]]$ is `MixedU(C,b)`, with $|u|=-2$. The operation `mixU[f](x)` places an underlying chain at $u^f$ or multiplies a total chain by $u^f$. Use `mixD` for its differential $b+uB$. Each underlying basis element can appear in at most one power of $u$ in a fixed total degree.

`MixedE(C,r,b,f)` is page $r$ of the decreasing filtration by powers of $u$. A page representative in $F^f$ has differential in $F^{f+r}$. `mixPageClass[r,f]` checks that condition and reduces by actual page boundaries, and `mixPageRep` returns a corrected representative. `mixD` on page classes takes $(b,f)$ to $(b-1,f+r)$. `mixPrev` maps a chosen lift to the preceding page, and `mixNext` passes a $d_r$-cycle to page $r+1$. A class hit by $d_r$ maps to zero. Use `mixdim[C,r,b,f];` for the dimension and `mixgen[C,r,b,f,k]` to select a computed generator. Chain and page operations are linear over the declared prime field, including `add` and `pscale`.

This is the $u$-adic spectral sequence of a mixed complex. For a Hochschild mixed complex it is the Hodge-to-de Rham spectral sequence. In the staircase example, $B x_{i-1}=y_i$ and $b x_i=-y_i$. The corrected lift $\sum_{i=0}^{r-1}u^i x_i$ has total differential $u^r y_r$. [Small worked example](examples/mixed-pages.asl) · [Staircase generator](examples/mixed-page-staircase.awk)

```sh
awk -v r=720 -v p=3 -f examples/mixed-page-staircase.awk | nix run path:. -- -
```

The program checks every preceding differential, the nonzero $d_{720}$, and the vanishing of its source and target on page $721$. Its final line is `MixedE(Stair, 720, -1, 720){g0}`. With `-v control=1`, the last $B$ column is omitted and the corrected class survives onto page $721$.

#### Importing a Hochschild weight

`rvHHMixed C[J,[mu]];` materializes the normalized Hochschild complex of $J$ at weight $\mu$ as a finite mixed complex over $\mathbb F_2$. Its generators `w0`, `w1`, and so on follow the ordered Hochschild word basis at that weight. The interpreter computes the columns of $b$ and $B$ on those words and checks the mixed-complex identities before registering the declaration. The computation is bounded by the finite basis and storage limits.

`rvHHToMixed[C](x)` takes an `RVHHChain` in the imported ideal and weight to `MixedChain(C,b)`. The inverse is `rvMixedToHH`. For a Hochschild homology class of degree $b$, `rvHHToMixed[C,f]` gives the corresponding class of `MixedE(C,1,b-2f,f)`. The filtration $f$ defaults to zero. The inverse on page one restores ordinary degree $b$. Under these identifications $d_1$ from filtration $f$ to $f+1$ is the induced Connes operator.

In odd exterior weight the example gives a nonzero $d_1$ at filtration four. In even weight it gives a class that survives to page two. In the Flag weight $(1,\ldots,1)$ the imported differential keeps the actual $m_3$ term. The only nonzero Hochschild groups are in degrees four, five and six with dimensions $1,3,2$, and the Connes ranks from four to five and from five to six are $1,2$. Thus the $B$ complex on these groups is exact and page two is supported at filtration zero in this weight. Every later differential there vanishes. [Hochschild import example](examples/rank-hoch-mixed.asl)

#### Mixed-complex homology as a graded module over the power variable

`mixbars[C];` decomposes the negative cyclic homology of a finite declared mixed complex over $\mathbb F_p[u]$ with $|u|=-2$. A free generator of birth degree $b$ contributes in degrees $b-2k$ for all $k\ge0$. A torsion bar $(b,r)$ is a generator of degree $b$ killed by $u^r$, and contributes for $0\le k<r$. The grading makes each degree calculation finite and the same calculation applies to the $u$-adically completed complex.

`mixBarRep[C,i]` gives the cycle representing the $i$th bar, taking the free list first. For a torsion bar, `mixBarNull[C,i]` returns a chain of degree $b-2r+1$ with differential equal to `mixU[r]` of the representative. The reduction orders generators by decreasing raw degree and tracks homogeneous column operations on $b+B$ at $u=1$. Undoing that specialization restores the unique power of $u$ for each coefficient. The distinct pivots give a triangular kernel basis and the incoming image is generated by the indicated powers of the torsion basis cycles. This gives the exact order of each torsion generator.

`mixbardim[C,b];` gives the negative cyclic dimension in signed total degree $b$. The query `mixbarrank[C,b,c];` gives the rank of multiplication by $u^{(b-c)/2}$ for $c\le b$ of the same parity. `mixbarperiodicdim[C,b];` counts the free bars of parity $b$ after inverting $u$. A torsion bar of length $r>0$ gives a primitive $d_r$ from total degree $b-2r+1$ at filtration zero to the representative at filtration $r$. The unit pivots account for $d_0$ pairs and contribute no bars. [Signed mixed-complex example](examples/mixed-bars.asl)

Use `-v bars=1` with the staircase generator to check the bar certificate alongside the page calculation. At $r=720$ the final line is `mixbars[Stair] = free [], torsion [(1439, 720)]`. With `-v control=1` it is `mixbars[Stair] = free [0, 1439], torsion []`. The `direct=1` mode checks the immediately preceding page and the selected differential and its target at the next page. It verifies $d_{4096}$ with the page calculation alone.

#### Mixed Hom complexes and coherent maps

`mixhom H[C,D];` constructs the finite graded Hom complex between two declared mixed complexes over the same prime field. The matrix unit from a source generator $x$ to a target generator $y$ has degree $|y|-|x|$. Construct that unit with `mixunit[H,x,y]`. For a homogeneous map $f$ the two mixed operators are
$b_H(f)=b_Df-(-1)^{|f|}fb_C$ and $B_H(f)=B_Df-(-1)^{|f|}fB_C$.

Use `mixEval(f,x)` to evaluate a map and `mixComp[K](g,f)` to compose $g$ after $f$, with $K$ a declared Hom complex for the resulting endpoints. `mixId[H]` is the degree-zero identity for a declared endomorphism complex. Evaluation and composition work on both ordinary chains and total chains, with powers of $u$ adding under composition. The signed Leibniz identity holds for both operators and for the total differential. [Mixed Hom example](examples/mixed-hom.asl)

The pages and bar queries on $H$ compute the obstruction to lifting a leading map to a homogeneous $\mathbb F_p[u]$-linear chain map, and the negative cyclic homology class of a closed lift. A degree-zero class represents a chain map modulo the boundaries from degree-one maps. The Hom staircase sends the single source generator to the corrected sum $\sum_{i<r}u^i x_i$. Its differential evaluates on that generator to $u^r y_r$. With the split control the same corrected map is closed. [Staircase generator](examples/mixed-page-staircase.awk)

```sh
awk -v r=720 -v p=3 -v hom=1 -v bars=1 -f examples/mixed-page-staircase.awk | nix run path:. -- -
```

This verifies composition with the source identity, evaluation of the actual obstruction, and the nonzero $d_{720}$ in the Hom complex. The matching torsion bar is $(1439,720)$. Add `-v control=1` for the closed map and two free bars.

#### Minimal models of mixed complexes

`mixminimal M[C];` chooses a contraction of $(C,b)$ onto its homology. The ordered boundary columns choose linear preimages, and reduction modulo cycles and boundaries gives the maps $i$, $p$ and $h$. On the underlying chains these satisfy $bh+hb=1-ip$, $pi=1$ and $h^2=ph=hi=0$.

`MixedMin(M,b)` is the degree-$b$ homology vector space. The query `mixmindim[M,b];` gives the dimension and `mixMinGen[M,b,k]` selects a basis class. `mixDelta[r]` acts from degree $b$ to $b+2r-1$. The coefficients of the transferred differential are
$\delta_r=pB(-hB)^{r-1}i$ for $r\ge1$ and $\delta_0=0$.
Together they give a higher $S^1$ action on the minimal model. The higher operators can be nonzero beyond $r=1$, whereas a strict mixed complex only has $b$ and $B$.

`MixedMinU(M,b)` is the total space with $|u|=-2$. The operation `mixU[f]` inserts a class at $u^f$ or multiplies a total chain by $u^f$. On these chains, `mixD` is $D_{\min}=\sum_{r\ge1}u^r\delta_r$ and satisfies $D_{\min}^2=0$. On underlying chains, `mixMinI`, `mixMinP[M]` and `mixMinH[M]` are $i,p,h$. On total chains they are the corrected maps
$I=(1+u hB)^{-1}i$, $P=p(1+u Bh)^{-1}$ and $H=h(1+u Bh)^{-1}$.
The inverse series terminate because their iterates raise raw degree. The maps $I$ and $P$ commute with the corresponding differentials and satisfy $PI=1$ and $DH+HD=1-IP$. [Minimal mixed-model example](examples/mixed-minimal.asl)

```sh
awk -v r=720 -v p=3 -v minimal=1 -f examples/mixed-page-staircase.awk | nix run path:. -- -
```

The program checks the corrected 720-term inclusion and prints `MixedMin(Minimal, 1439){g1}` for the sole positive higher operator on the leading class. The split control gives zero. With `-v direct=1 -v r=2048`, the program verifies the same identities for the $2048$th operator.

#### Tensor products and Künneth certificates for mixed complexes

`mixtensor T[C,D];` takes two finite mixed complexes over the same prime field. In the finite tensor complex, both operators use the signed rule
$F(c\otimes e)=F_C(c)\otimes e+(-1)^{|c|}c\otimes F_D(e)$ for $F=b,B$.
Apply `mixTensor[T](x,y)` to the two underlying chains or the two total chains. Powers of $u$ add in a total tensor product. `mixSwap[S]` sends a tensor to the declared complex with reversed factors using the sign $(-1)^{|c||d|}$.

The chosen bar cycles give explicit Künneth classes over $\mathbb F_p[u]$. If their birth degrees are $a,b$, `mixTensorCycle[T,i,j]` gives their product $A\otimes B$ in degree $a+b$. Two free factors give a free class. A free factor with a length-$m$ torsion factor gives a class of order $u^m$. For two torsion factors of lengths $l,m$, let $r=\min(l,m)$. The product has order $u^r$. The secondary class also has order $u^r$ and is born in degree $a+b-2\max(l,m)+1$.

`mixTensorNull[T,i,j]` gives a chain whose differential is the first killing power of the product. For two torsion factors, write $D X=u^l A$ and $D Y=u^m B$ for their chosen nulls. The operation `mixTensorTor[T,i,j]` constructs the secondary cycle
$W=u^{l-r}A\otimes Y-(-1)^a u^{m-r}X\otimes B$.
The witness `mixTensorTorNull[T,i,j]` is $X\otimes Y$, with differential $u^r W$. The secondary generator uses the chosen cycles and nulls. Changing those choices can change the representative within the Künneth extension.

The length-two and length-three staircase tensor has two length-two bars born in degrees five and ten. The example checks the signed secondary cycle, both nulls, and the surviving first multiples on page three. The test suite also tensors two actual imported exterior Hochschild weight complexes. [Tensor and Künneth example](examples/mixed-tensor.asl)

#### Higher complexes and mapping cones

`mixedhigher K[p] {x:a,y:b,...} {x => [(r,y,c),...],...};` declares a finite free graded $\mathbb F_p[u]$ complex with $|u|=-2$. Each entry gives the term $c u^r y$ of $D(x)$ and must have $|y|=|x|+2r-1$. Missing columns are zero. The interpreter checks $D^2=0$ coefficientwise. Indeed, the source and target degrees determine the only possible exponent in each matrix entry, so checking at $u=1$ detects every coefficient of the homogeneous relation. Large exponents occupy a single sparse entry.

Use `MixedChain` and `MixedU` for underlying and total chains, and `mixDelta[r]` for the raw operator of degree $2r-1$. On total chains `mixD` computes $D=\sum_{r\ge0}u^r\delta_r$. The existing page and bar operations use this full differential, and `mixhom` uses the corresponding signed total Hom differential. The individual mixed operators `mixb` and `mixB`, strict mixed tensor, and strict mixed minimal-model declarations require a complex with no higher components.

`mixcone K(f);` takes a closed degree-zero total chain $f$ in a declared Hom complex from $C$ to $E$. Its cone has underlying generators from $E$ followed by the degree-one suspension of those from $C$. It has differential
$D(y,sx)=(D_Ey+f(x),-sD_Cx)$.
The maps `mixConeI[K]`, `mixConeLift[K]` and `mixConeP[K]` send $y$ to $(y,0)$, $x$ to $(0,sx)$ and $(y,sx)$ to $x$, respectively. Projection lowers degree by one and anticommutes with $D$. For $g=f+D_Hh$ with $|h|=1$, `mixConeShear[L](h,z)` gives the isomorphism from the cone of $f$ into the cone $L$ of $g$ by $(y,sx)\mapsto(y-h(x),sx)$.

The cone of $x\mapsto u^{1000000000}y$ on two generators verifies a nonzero $d_{1000000000}$ and gives the torsion bar $(2000000000,1000000000)$ over $\mathbb F_3$. The zero map gives two free bars in degrees $1$ and $2000000000$. A second test uses a three-term coherent map whose cone has $\delta_1^2\ne0$, with its contribution to $D^2$ canceled by $\delta_0\delta_2$. [Two-generator cone example](examples/mixed-cone.asl)

#### Graded negative cyclic homology and the cone triangle

`MixedH(C,b)` is the actual homology of `MixedU(C,b)` with differential `mixD`. It applies to strict mixed complexes and sparse higher complexes. In a fixed signed degree the polynomial and completed complexes have the same finite basis. `mixHdim[C,b];` gives the dimension. The coordinate `mixHgen[C,b,i]` uses the subsequence of bars alive in degree $b$, in the order of `mixbars`. Each bar with birth $a$ and representative $z$ contributes $[u^{(a-b)/2}z]$ when the exponent exists and lies below its torsion length. A free bar has no upper cutoff.

`mixHClass` takes a total cycle to its coordinates and rejects a noncycle. `mixHRep` returns the corresponding sum of shifted bar representatives. `mixHPrimitive` gives a deterministic preimage of an exact total chain under $D$ and rejects chains outside the image. It reduces the requested degree against the incoming columns in raw generator order. The operations `add` and `pscale` are linear on these coordinates. `mixU[k]` sends a living coordinate to its further shift if the bar still lives, with zero for a killed coordinate.

`mixInduce(f,x)` takes a homogeneous closed total Hom map and a source homology class. It evaluates the map on the representative and reduces the image into target coordinates. The output degree is the sum of the map and class degrees. On homology `mixConeI[K]` and `mixConeP[K]` induce inclusion and projection, and `mixConeConnect[K]` is the connecting map given by the map of the cone. They give
$H_b(X)\xrightarrow{f_*}H_b(Y)\xrightarrow{I_*}H_b(K)\xrightarrow{P_*}H_{b-1}(X)\xrightarrow{f_*}H_{b-1}(Y)$.
Projection has degree $-1$. The two-generator and rank-one examples test the maps on actual cycles and boundary witnesses.

A map can induce zero on homology while its cone records a nonsplit extension of the two homology modules. In the four-generator example $f$ sends $s_0$ to $u^j y$ and kills $x_0$. Its class in degree zero of the Hom complex has exact $u$-order $\min(l,m)-j$. The cone relations are $u^j y+u^l x$ and $u^m y$, with $x=s x_0$. For $0\le j<\min(l,m)$ the cone has a short bar of length $j$ at birth $2j$, when $j>0$, and a long bar of length $l+m-j$ at birth $2l$. With $l=m=100000000$ and $j=17$, the example verifies the primitive for the long bar and the survival of its preceding power. The zero map gives separate bars of length $100000000$. [Cone triangle example](examples/mixed-triangle.asl)

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
