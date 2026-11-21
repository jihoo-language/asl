mod common;

use common::{diagnostic, example, interpret, stdout};

const SPHERE: &str = "
    let h0 = mapextgen[One, One, 1, 1, 0];
    let h1 = mapextgen[One, One, 1, 2, 0];
    let h2 = mapextgen[One, One, 1, 4, 0];
";

#[test]
fn the_multi_object_transfer_uses_composable_matrix_paths() {
    assert_eq!(
        stdout(example("map-transfer")),
        "MapExt(One, One, 2, 8){g0}\n"
    );
    let source = r"
        space Two { a, b };
        let swap = emap(tab Two -> Two { a => Two.b, b => Two.a });
        let first = emap(\x : Two => Two.a);
        let second = emap(\x : Two => Two.b);
        assert eHigher(swap, first) == second;
        assert eHigher(first, swap) == first;
        assert eHigher(swap, swap) == eid[Two];
        let h0 = add(mapextgen[Two, Two, 1, 1, 0], mapextgen[Two, Two, 1, 1, 3]);
        let h1 = add(mapextgen[Two, Two, 1, 2, 0], mapextgen[Two, Two, 1, 2, 3]);
        assert eHigher(h0, h1) == ecomp(h0, h1);
        assert eHigher(h0, h1, h0) == ecomp(h1, h1);
        assert eHigher(eid[Two], h0) == h0;
        assert eHigher(h0, eid[Two]) == h0;
        assert eLift(h0, eid[Two]) == zero[MapCobar(Two, Two, 0, 1)];
        assert eHigher(h0, eid[Two], h1) == zero[MapExt(Two, Two, 1, 3)];
        assert eLift(h0, eid[Two], h1) == zero[MapCobar(Two, Two, 0, 3)];
        assert eHigher(h0, swap, first) == zero[MapExt(Two, Two, 0, 1)];
    ";
    assert!(stdout(interpret(source)).is_empty());
}

#[test]
fn strict_units_do_not_discard_nonidentity_maps_of_filtration_zero() {
    let source = format!(
        r"{}
        let outA = mapextgen[RP2, One, 1, 1, 0];
        let crossA = mapextgen[RP2, RP2, 0, -2, 0];
        let inA = mapextgen[One, RP2, 0, 1, 0];
        assert eHigher(outA, crossA, inA) == eid[One];
        let outB = mapextgen[RP2, One, 1, 0, 0];
        let crossB = mapextgen[RP2, RP2, 0, 1, 0];
        let inB = mapextgen[One, RP2, 1, 1, 0];
        assert eHigher(outB, crossB, inB) == mapextgen[One, One, 1, 2, 0];
        assert eHigher(outB, eid[RP2], inB) == zero[MapExt(One, One, 1, 1)];
        assert eLift(outB, eid[RP2], inB) == zero[MapCobar(One, One, 0, 1)];
    ",
        include_str!("../examples/projective-plane.asl")
    );
    stdout(interpret(&source));
}

#[test]
fn lifts_satisfy_the_morphism_identity_on_a_nonclosed_arity_four_input() {
    let source = format!(
        "{SPHERE}
        assert eHigher(h2, h0) == ecomp(h2, h0);
        assert eHigher(h0, h2) == ecomp(h0, h2);
        assert eD(eLift(h0, h2)) == add(ecomp(eRep(h0), eRep(h2)), eRep(eHigher(h0, h2)));
        assert eLift(h1, h0) == mword[One, 0, 0, One, 0, 0](add(cword{{[0, 1]}}, cword{{[3]}}));

        let fourth = eHigher(h0, h2, h1, h0);
        assert fourth == ecomp(h2, h2);
        let f4 = eLift(h0, h2, h1, h0);
        let expected = add(cword{{[1, 0, 1]}}, add(cword{{[2, 2]}}, cword{{[5, 1]}}));
        assert f4 == mword[One, 0, 0, One, 0, 0](expected);
        let u4a = ecomp(eRep(h0), eLift(h2, h1, h0));
        let u4b = ecomp(eLift(h0, h2), eLift(h1, h0));
        let u4c = ecomp(eLift(h0, h2, h1), eRep(h0));
        let u4 = add(u4a, add(u4b, u4c));
        let defectA = mword[One, 0, 0, One, 0, 0](cword{{[1], [4]}});
        let defectB = mword[One, 0, 0, One, 0, 0](cword{{[1], [0, 1], [4]}});
        assert eD(u4) == add(ecomp(defectA, eLift(h1, h0)), defectB);

        let a = eLift(h0, h2, eHigher(h1, h0));
        let b = eLift(h0, eHigher(h2, h1), h0);
        let c = eLift(eHigher(h0, h2), h1, h0);
        let d = eLift(h0, eHigher(h2, h1, h0));
        let e = eLift(eHigher(h0, h2, h1), h0);
        assert eD(f4) == add(u4, add(eRep(fourth), add(add(a, b), add(c, add(d, e)))));
    "
    );
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn higher_operations_satisfy_a_nontrivial_stasheff_relation() {
    let source = format!(
        "{SPHERE}
        let left = eHigher(h0, h1, h0);
        let right = eHigher(h1, h0, h1);
        assert left == ecomp(h1, h1);
        assert right == ecomp(h2, h0);
        assert eHigher(h1, left) == mapextgen[One, One, 3, 6, 0];
        assert eHigher(h1, left) == eHigher(right, h0);
        let a = eHigher(h0, h2, eHigher(h1, h0));
        let b = eHigher(h0, eHigher(h2, h1), h0);
        let c = eHigher(eHigher(h0, h2), h1, h0);
        let d = eHigher(h0, eHigher(h2, h1, h0));
        let e = eHigher(eHigher(h0, h2, h1), h0);
        assert add(add(a, b), add(c, add(d, e))) == zero[MapExt(One, One, 3, 8)];
    "
    );
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn transferred_grades_use_wide_signed_sums_and_diagnose_invalid_results() {
    let zero = |(s, t): (usize, i64)| format!("zero[MapExt(One, One, {s}, {t})]");
    for (left, middle, right, t) in [
        ((1, i64::MAX), (1, 1), (1, -1), i64::MAX),
        ((1, i64::MIN), (1, -1), (1, 1), i64::MIN),
    ] {
        let inputs = format!("{}, {}, {}", zero(left), zero(middle), zero(right));
        let source = format!(
            "assert eHigher({inputs}) == {}; assert eLift({inputs}) == zero[MapCobar(One, One, 1, {t})];",
            zero((2, t)),
        );
        assert!(stdout(interpret(&source)).is_empty());
    }
    let s = usize::MAX;
    let inputs = format!("{}, {}, {}", zero((s, 0)), zero((1, 0)), zero((1, 0)));
    assert!(stdout(interpret(&format!(
        "assert eLift({inputs}) == zero[MapCobar(One, One, {s}, 0)];"
    )))
    .is_empty());
    for (source, message) in [
        ("show eHigher();".into(), "at least two"),
        ("show eLift(eid[One]);".into(), "at least two"),
        ("show eHigher(eid[S], eid[One]);".into(), "composable"),
        ("show eHigher(S.lo, eid[One]);".into(), "MapExt"),
        ("show eLift(eid[One], eid[One]);".into(), "negative"),
        (
            "show eHigher(eid[One], eid[One], eid[One]);".into(),
            "negative",
        ),
        (
            format!("show eHigher({inputs});"),
            "cohomological degree is out of range",
        ),
        (
            format!("show eHigher({}, {});", zero((0, i64::MAX)), zero((0, 1))),
            "internal degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&source));
        assert!(error.contains(message), "{error}");
    }
}

#[test]
fn unary_calls_can_use_bound_names_for_higher_operations() {
    let source = r"
        let eHigher = \x : One => x;
        let eLift = \x : One => eHigher(x);
        assert eLift(One.pt) == One.pt;
    ";
    assert!(stdout(interpret(source)).is_empty());
}
