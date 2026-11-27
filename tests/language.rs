mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn cli_accepts_the_full_syntax_depth_and_reports_overflow_as_a_diagnostic() {
    let nested_add = |depth| {
        let mut term = "v".to_owned();
        for _ in 0..depth {
            term = format!("add({term}, v)");
        }
        format!("let v = zero[One, 0]; show {term};")
    };
    let grouped = |depth| format!("show {}One.pt{};", "(".repeat(depth), ")".repeat(depth));

    assert_eq!(stdout(interpret(&nested_add(64))), "H(One, 0){0}\n");
    assert_eq!(stdout(interpret(&grouped(64))), "One.pt\n");
    for source in [nested_add(65), grouped(65)] {
        let result = interpret(&source);
        assert_eq!(result.status.code(), Some(1));
        let error = diagnostic(result);
        assert!(
            error.contains("syntax nesting exceeds 64 levels"),
            "{error}"
        );
    }
}

#[test]
fn named_unary_parsing_respects_the_full_syntax_depth() {
    let frobenius = |depth| {
        format!(
            "rvunram F[3,1]{{0,1}}; show {}rvUnramElem[F]{{1}}{};",
            "rvUnramFrob[1](".repeat(depth),
            ")".repeat(depth)
        )
    };

    assert_eq!(stdout(interpret(&frobenius(64))), "RVUnram(F)[1]\n");
    let result = interpret(&frobenius(65));
    assert_eq!(result.status.code(), Some(1));
    let error = diagnostic(result);
    assert!(
        error.contains("syntax nesting exceeds 64 levels"),
        "{error}"
    );
}

#[test]
fn rotation_descends_and_its_finite_homotopy_class_is_rigid() {
    assert_eq!(
        stdout(example("descent")),
        "Circle.bc\nbetti(Circle) = [1, 1]\nbetti(Circle * Circle) = [1, 2, 1]\nhomotopy = no\n"
    );
}

#[test]
fn a_point_evaluates_opens_and_meets_compute_intersections() {
    assert_eq!(
        stdout(example("observables")),
        "S.hi\nS.lo\nhomotopy = yes (2 steps)\n"
    );
}

#[test]
fn torus_projections_generate_the_top_class_and_pullbacks_preserve_products() {
    assert_eq!(stdout(example("torus-cup")), "H(Circle * Circle, 2){g0}\n");
}

#[test]
fn mixed_degree_cup_products_associate_on_the_three_torus() {
    let program = r"
        space C { p < a, q < a, p < b, q < b };
        let u = gen[C, 1, 0];
        let a = pull((\p : C * C * C => fst(fst(p))), u);
        let b = pull((\p : C * C * C => snd(fst(p))), u);
        let c = pull((\p : C * C * C => snd(p)), u);
        let volume = cup(cup(a, b), c);
        assert volume == cup(a, cup(b, c));
        assert volume == gen[C * C * C, 3, 0];
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn cocycle_representatives_are_equal_modulo_coboundaries() {
    let program = r"
        space T { top, bottom < middle < top };
        assert cocycle[T, 1] { bottom < middle, bottom < top } == zero[T, 1];
        assert cocycle[T, 2] { bottom < middle < top } == zero[T, 2];
        assert cocycle[T, 0] { top, top } == zero[T, 0];
        let unit = cocycle[T, 0] { bottom, middle, top };
        assert unit == gen[T, 0, 0];
        let inclusion = res[middle](\x : T => x);
        assert pull(inclusion, unit) == cocycle[Star(T, middle), 0] { middle, top };
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn degree_zero_ring_and_discrete_class_spaces_are_first_class() {
    let program = r"
        space D { a, b };
        let a = cocycle[D, 0] { a };
        let b = cocycle[D, 0] { b };
        let one = cocycle[D, 0] { a, b };
        let plus = \x : H(D, 0) => \y : H(D, 0) => add(x, y);
        assert plus(a)(b) == one;
        assert plus(one)(a) == b;
        assert cup(one, a) == a;
        assert cup(a, b) == zero[D, 0];
        assert (\x : H(D, 0) => cup(x, x)) == (\x : H(D, 0) => x);
        let swap = tab D -> D { a => D.b, b => D.a };
        assert pull(swap, a) == b;
        let translate = \p : H(D, 0) * H(D, 0) => add(fst(p), snd(p));
        assert translate((a, b)) == one;
        show add(gen[D, 0, 0], gen[D, 0, 1]);
        betti H(D, 0);
        homotopy (\x : One => a) ~ (\x : One => b);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "H(D, 0){g0 + g1}\nbetti(H(D, 0)) = [4]\nhomotopy = no\n"
    );
}

#[test]
fn pull_is_continuous_in_the_function_argument_and_handles_degenerate_images() {
    let program = r"
        space C { p < a, q < a, p < b, q < b };
        let u = cocycle[C, 1] { p < a };
        let degree = \f : C -> C => pull(f, u);
        let id = \x : C => x;
        let collapse = \x : C => C.p;
        assert degree(id) == u;
        assert degree(collapse) == zero[C, 1];
        let turn = tab C -> C { p => C.q, q => C.p, a => C.b, b => C.a };
        assert degree(turn) == u;
        assert pull((\x : C => turn(collapse(x))), u) == pull(collapse, pull(turn, u));
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn empty_spaces_high_degrees_and_iterated_cohomology_have_the_right_zero() {
    let program = r"
        space Void {};
        betti H(Void, 0);
        betti H(H(Void, 0), 0);
        betti H(S -> S, 0);
        betti H(One, 1000000000);
        assert zero[Void, 0] == cocycle[Void, 0] {};
        assert pull((\x : Void => One.pt), gen[One, 0, 0]) == zero[Void, 0];
        assert cup(zero[One, 1000000000], zero[One, 1000000000]) == zero[One, 2000000000];
    ";
    let out = stdout(interpret(program));
    let ranks: Vec<_> = out
        .lines()
        .map(|s| s.split_once(" = ").unwrap().1)
        .collect();
    assert_eq!(ranks, ["[1]", "[2]", "[2]", "[1]"]);
}

#[test]
fn a_nonclosed_literal_names_the_simplex_that_detects_its_coboundary() {
    let declaration = "space T { bottom < middle < top };";
    let bad = "cocycle[T, 1] { bottom < middle }";
    for expr in [bad.to_string(), format!(r"\x : Void => {bad}")] {
        let error = diagnostic(interpret(&format!(
            "{declaration} space Void {{}};\nshow {expr};"
        )));
        for witness in [":2:", "bottom", "middle", "top"] {
            assert!(error.contains(witness), "{error}");
        }
    }
}

#[test]
fn cohomology_requires_matching_spaces_degrees_and_strict_chains() {
    for expr in [
        "add(zero[S, 0], zero[S, 1])",
        "cup(zero[T, 0], zero[S, 1])",
        r"pull((\x : One => T.a), zero[S, 1])",
        "pull(S.lo, zero[S, 0])",
        "cup(S.lo, S.hi)",
        "gen[One, 2, 0]",
        "gen[S, 0, 1]",
        "cocycle[T, 1] { a < b }",
        "cocycle[T, 1] { a < a }",
        "cocycle[T, 1] { c < a }",
        "cocycle[T, 1] { a }",
        "cocycle[T, 1] { a < absent }",
        "cocycle[T * T, 1] {}",
        "tab S -> H(One, 0) { lo => zero[One, 0], hi => gen[One, 0, 0] }",
    ] {
        let error = diagnostic(interpret(&format!(
            "space T {{ a < c, b < c }};\nshow {expr};"
        )));
        assert!(error.contains(":2:"), "{expr}\n{error}");
    }
}

#[test]
fn individual_classes_work_beyond_the_class_space_enumeration_limit() {
    let large = "space D { a, b, c, d, e, f, g, h, i, j, k, l, m };";
    let source = format!("{large} let u = gen[D, 0, 12]; assert cup(u, u) == u; show add(u, u);");
    assert_eq!(stdout(interpret(&source)), "H(D, 0){0}\n");
    let error = diagnostic(interpret(&format!("{large}\nbetti H(D, 0);")));
    assert!(error.contains(":2:") && error.contains("limit"), "{error}");
}

#[test]
fn a_saturated_height_estimate_cannot_prove_a_group_vanishes() {
    // For D discrete, D -> C is a product of one circle per point of D.
    // Here D has more points than the largest representable degree.
    let labels = (0..usize::BITS)
        .map(|i| format!("p{i}"))
        .collect::<Vec<_>>()
        .join(", ");
    let source = format!(
        "space Bits {{ {labels} }}; space C {{ p < a, q < a, p < b, q < b }};\n\
         betti H(H(Bits, 0) -> C, {});",
        usize::MAX
    );
    let error = diagnostic(interpret(&source));
    assert!(error.contains(":2:") && error.contains("limit"), "{error}");
}

#[test]
fn earlier_names_remain_usable_around_cohomology_syntax() {
    let program = r"
        space H { p };
        let zero = H.p;
        let gen = \x : H => x;
        let cocycle = gen;
        let cup = cocycle;
        let add = cup;
        let pull = add;
        assert pull(add(cup(cocycle(gen(zero))))) == H.p;
        show gen[H, 0, 0];
    ";
    assert_eq!(stdout(interpret(program)), "H(H, 0){g0}\n");
}

#[test]
fn double_cover_sections_disagree_on_the_last_overlap() {
    let error = diagnostic(example("rejected/double-cover"));
    for witness in ["double-cover.asl:", "ca0", "ca1"] {
        assert!(error.contains(witness), "{error}");
    }
}

#[test]
fn currying_and_uncurrying_are_inverse_on_the_whole_exponential() {
    let program = r"
        let curry = \k : (S * S) -> S => \x : S => \y : S => k((x, y));
        let uncurry = \h : S -> S -> S => \p : S * S => h(fst(p))(snd(p));
        assert (\k : (S * S) -> S => uncurry(curry(k))) == (\k : (S * S) -> S => k);
        assert (\h : S -> S -> S => curry(uncurry(h))) == (\h : S -> S -> S => h);
        let eval = \p : (S -> S) * S => fst(p)(snd(p));
        show eval(((\x : S => x), S.hi));
    ";
    assert_eq!(stdout(interpret(program)), "S.hi\n");
}

#[test]
fn restriction_and_gluing_compose_on_nested_stars() {
    let program = r"
        space Chain { a < b < c };
        let id = \x : Chain => x;
        let local = res[b](id);
        let germ = res[c](local);
        assert germ(Star(Star(Chain, b), c).c) == Chain.c;
        assert (sew Star(Chain, b) -> Chain { b => res[b](local) }) == local;
        show local(Star(Chain, b).c);
    ";
    assert_eq!(stdout(interpret(program)), "Chain.c\n");
}

#[test]
fn overlap_agreement_is_extensional_for_function_valued_charts() {
    let program = r"
        space Joint { w < j, e < j };
        let id = \x : S => x;
        let eta = \y : S => id(y);
        let bottom = \x : S => S.lo;
        let west = tab Star(Joint, w) -> (S -> S) { w => bottom, j => eta };
        let east = \x : Star(Joint, e) => id;
        let glued = sew Joint -> (S -> S) { e => east, w => west };
        assert res[w](glued) == west;
        show glued(Joint.j)(S.hi);
    ";
    assert_eq!(stdout(interpret(program)), "S.hi\n");
}

#[test]
fn empty_domains_and_discrete_targets_have_their_categorical_counts() {
    let program = r"
        space Void {};
        space Two { a, b };
        let empty = tab Void -> Void {};
        assert empty == (sew Void -> Void {});
        assert empty == (\x : Void => x);
        betti Void;
        betti Void * S;
        betti Void -> Void;
        betti S -> Void;
        betti Two -> Two;
        betti S -> Two;
    ";
    let output = stdout(interpret(program));
    let ranks: Vec<_> = output
        .lines()
        .map(|line| line.split_once(" = ").unwrap().1)
        .collect();
    assert_eq!(ranks, ["[]", "[]", "[1]", "[]", "[4]", "[2]"]);
}

#[test]
fn homotopy_measures_shortest_comparability_zigzags() {
    let program = r"
        space Fence { a < b, c < b, c < d, e < d, isolated };
        let at = \x : Fence => \u : One => x;
        homotopy at(Fence.a) ~ at(Fence.a);
        homotopy at(Fence.a) ~ at(Fence.b);
        homotopy at(Fence.a) ~ at(Fence.e);
        homotopy at(Fence.b) ~ at(Fence.isolated);
        space D { x, y };
        let left = tab D -> S { x => S.hi, y => S.lo };
        let right = tab D -> S { x => S.lo, y => S.hi };
        homotopy left ~ right;
    ";
    assert_eq!(
        stdout(interpret(program)),
        "homotopy = yes (0 steps)\nhomotopy = yes (1 step)\nhomotopy = yes (4 steps)\nhomotopy = no\nhomotopy = yes (2 steps)\n"
    );
}

#[test]
fn individual_maps_work_beyond_the_function_space_enumeration_limit() {
    let domain = "space D { a, b, c, d, e, f, g, h, i, j, k, l, m };";
    let program = format!(
        r"{domain}
        let low = \x : D => S.lo;
        let high = \x : D => S.hi;
        show high(D.m);
        homotopy low ~ high;"
    );
    assert_eq!(
        stdout(interpret(&program)),
        "S.hi\nhomotopy = yes (1 step)\n"
    );

    let error = diagnostic(interpret(&format!("{domain}\nbetti D -> S;")));
    assert!(error.contains(":2:") && error.contains("limit"), "{error}");
}

#[test]
fn discontinuity_points_to_the_reversed_comparison() {
    let error = diagnostic(interpret(
        "\nlet swap = tab S -> S { lo => S.hi, hi => S.lo };",
    ));
    for witness in [":2:", "S.lo", "S.hi"] {
        assert!(error.contains(witness), "{error}");
    }
}

#[test]
fn space_and_star_names_survive_type_checking() {
    for program in [
        r"space Copy { lo < hi }; show (\x : S => x)(Copy.lo);",
        r"let f = res[lo](\x : S => x); show f(S.hi);",
        r"let f = \x : S => x; let bad = sew S -> S { lo => f };",
    ] {
        diagnostic(interpret(program));
    }
}

#[test]
fn empty_binders_still_require_well_typed_bodies() {
    for body in [
        "ghost",
        "S.absent",
        "fst(S.lo)",
        "tab S -> S { hi => S.hi }",
    ] {
        let program = format!(r"space Void {{}}; let f = \x : Void => {body};");
        diagnostic(interpret(&program));
    }
}

#[test]
fn invalid_orders_incomplete_atlases_and_false_assertions_fail() {
    for program in [
        "space Cycle { a < b < c < a };",
        "space D { a, b }; let bad = tab D -> S { a => S.hi, a => S.lo };",
        r"space D { a, b }; let f = \x : D => x; let bad = sew D -> D { a => res[a](f) };",
        "assert S.lo == S.hi;",
        r"homotopy (\x : S => x) ~ (\x : One => S.lo);",
    ] {
        diagnostic(interpret(program));
    }
}
