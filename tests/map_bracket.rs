mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn the_sphere_triple_is_a_nonzero_singleton_and_a_nonexact_triple_is_empty() {
    assert_eq!(stdout(example("map-massey")), "card = 1 / 2\n");
    let program = "
        let h0 = mapextgen[One, One, 1, 1, 0];
        let h1 = mapextgen[One, One, 1, 2, 0];
        let square = ecomp(h1, h1);
        assert eInBracket(h0, h1, h0, square) == S.hi;
        let empty = eBracket(h0, h0, h0);
        assert empty(zero[MapExt(One, One, 2, 3)]) == S.lo;
        assert eInBracket(h0, h0, h0, zero[MapExt(One, One, 2, 3)]) == S.lo;
        card empty;
    ";
    assert_eq!(stdout(interpret(program)), "card = 0 / 1\n");
}

#[test]
fn the_projective_plane_coaction_gives_a_singleton_sphere_identity() {
    let calculation = "
        let f = mapextgen[One, RP2, 0, 1, 0];
        let g = mapextgen[RP2, RP2, 0, -2, 0];
        let h = mapextgen[RP2, One, 1, 1, 0];
        let bottom = mword[One, 0, 0, RP2, 1, 0](cword{});
        let top = mword[RP2, 2, 0, RP2, 0, 0](cword{});
        let scalar = mword[RP2, 0, 0, One, 0, 0](cword{[1]});
        let primitive = mword[RP2, 1, 0, One, 0, 0](cword{});
        assert eRep(f) == bottom;
        assert eRep(g) == top;
        assert eRep(h) == scalar;
        assert ecomp(top, bottom) == zero[MapCobar(One, RP2, 0, -1)];
        assert ecomp(scalar, top) == eD(primitive);
        assert ecomp(primitive, bottom) == eRep(eid[One]);

        let triple = eBracket(h, g, f);
        assert triple(eid[One]) == S.hi;
        assert eInBracket(h, g, f, eid[One]) == S.hi;
        assert triple(zero[MapExt(One, One, 0, 0)]) == S.lo;
        card triple;
    ";
    let source = format!(
        "{}\n{calculation}",
        include_str!("../examples/map-cobar.asl")
    );
    assert!(stdout(interpret(&source)).ends_with("card = 1 / 2\n"));
}

#[test]
fn mixed_endpoints_and_rank_one_outer_maps_give_full_and_proper_indeterminacy() {
    let program = r"
        space Two { a, b };
        let middleMixed = zero[MapExt(One, Two, 1, 0)];
        let full = eBracket(eid[Two], middleMixed, eid[One]);
        let inclusion = emap(\x : One => Two.a);
        assert full(inclusion) == S.hi;
        assert eInBracket(eid[Two], middleMixed, eid[One], inclusion) == S.hi;
        card full;

        let projector = emap(\x : Two => Two.a);
        let middle = zero[MapExt(Two, Two, 1, 0)];
        let proper = eBracket(projector, middle, projector);
        assert proper(projector) == S.hi;
        assert proper(eid[Two]) == S.lo;
        let direct = \x : MapExt(Two, Two, 0, 0) => eInBracket(projector, middle, projector, x);
        assert direct == proper;
        card proper;
    ";
    assert_eq!(stdout(interpret(program)), "card = 4 / 4\ncard = 8 / 16\n");
}

#[test]
fn a_negative_primitive_degree_allows_a_zero_composition_only() {
    let program = "
        let id = eid[One];
        let positive = zero[MapExt(One, One, 1, 0)];
        let zeroMap = zero[MapExt(One, One, 0, 0)];
        assert eBracket(positive, id, id)(zeroMap) == S.lo;
        assert eInBracket(positive, id, id, zeroMap) == S.lo;
        assert eBracket(positive, zeroMap, id)(id) == S.hi;
        assert eInBracket(positive, zeroMap, id, id) == S.hi;
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn direct_membership_handles_a_target_too_large_to_tabulate() {
    let setup = "
        space Five { a, b, c, d, e };
        let id = eid[Five];
        let middle = zero[MapExt(Five, Five, 1, 0)];
        let zeroMap = zero[MapExt(Five, Five, 0, 0)];
    ";
    let direct = format!(
        "{setup}
        assert eInBracket(id, middle, id, id) == S.hi;
        assert eInBracket(zeroMap, middle, zeroMap, id) == S.lo;"
    );
    assert!(stdout(interpret(&direct)).is_empty());
    let error = diagnostic(interpret(&format!("{setup}show eBracket(id, middle, id);")));
    assert!(error.contains("limit"), "{error}");
}

#[test]
fn triple_degree_checks_widen_signed_sums_and_validate_the_result() {
    let max_t = i64::MAX;
    let min_t = i64::MIN;
    let max_s = usize::MAX;
    for (h, g, f, target) in [
        ((1, -1), (0, 1), (0, max_t), (0, max_t)),
        ((1, 1), (0, -1), (0, min_t), (0, min_t)),
        ((0, 0), (1, 0), (max_s, 0), (max_s, 0)),
    ] {
        let zero = |(s, t): (usize, i64)| format!("zero[MapExt(One, One, {s}, {t})]");
        let program = format!(
            "assert eInBracket({}, {}, {}, {}) == S.hi;",
            zero(h),
            zero(g),
            zero(f),
            zero(target)
        );
        assert!(stdout(interpret(&program)).is_empty());
    }
    let invalid = format!(
        "show eInBracket(zero[MapExt(One, One, 1, {max_t})], \
         zero[MapExt(One, One, 0, 1)], eid[One], eid[One]);"
    );
    assert!(diagnostic(interpret(&invalid)).contains("internal degree is out of range"));
}

#[test]
fn triple_operations_check_arity_endpoints_and_query_type() {
    let h = "zero[MapExt(One, One, 1, 0)]";
    let expressions = [
        (format!("eBracket({h}, eid[One])"), "three"),
        (format!("eInBracket({h}, eid[One], eid[One])"), "four"),
        ("eBracket(eid[One], eid[One], eid[One])".into(), "negative"),
        (format!("eBracket({h}, eid[S], eid[One])"), "composable"),
        (format!("eBracket({h}, S.lo, eid[One])"), "MapExt"),
        (
            format!("eInBracket({h}, eid[One], eid[One], eid[S])"),
            "MapExt",
        ),
    ];
    for (expression, expected) in expressions {
        let error = diagnostic(interpret(&format!("show {expression};")));
        assert!(error.contains(expected), "{error}");
    }
}
