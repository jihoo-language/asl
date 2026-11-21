mod common;

use common::{diagnostic, example, interpret, stdout};

const SPHERE: &str = "
    let h0 = mapextgen[One, One, 1, 1, 0];
    let h1 = mapextgen[One, One, 1, 2, 0];
    let a = mword[One, 0, 0, One, 0, 0](cword{[3]});
    let b = mword[One, 0, 0, One, 0, 0](cword{[0, 1]});
    let u = add(a, b);
";

#[test]
fn explicit_sphere_witnesses_give_the_closed_raw_top_and_its_class() {
    assert_eq!(
        stdout(example("map-system")),
        "MapExt(One, One, 2, 4){g0}\n"
    );
    let source = format!(
        "{SPHERE}
        let top = eCurvature(h0, h1, eRep(h0)){{1 < 3 => b, 0 < 2 => u}};
        let one = cword{{[3], [1]}};
        let two = cword{{[0, 1], [1]}};
        let three = cword{{[1], [0, 1]}};
        let expected = mword[One, 0, 0, One, 0, 0](add(one, add(two, three)));
        assert top == expected;
        assert eD(top) == zero[MapCobar(One, One, 3, 4)];
        assert eCurvature(h1, h0){{}} == ecomp(eRep(h1), eRep(h0));
        assert eMassey(h0, h0){{}} == ecomp(h0, h0);
        let exact = eD(u);
        assert eCurvature(eid[One], exact){{}} == exact;
        assert eMassey(eid[One], exact){{}} == zero[MapExt(One, One, 2, 3)];
    "
    );
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn a_projective_plane_defining_system_uses_a_single_signed_degree_witness() {
    let source = format!(
        "{}
        let f = mapextgen[One, RP2, 0, 1, 0];
        let g = mapextgen[RP2, RP2, 0, -2, 0];
        let h = mapextgen[RP2, One, 1, 1, 0];
        let v = mword[RP2, 1, 0, One, 0, 0](cword{{}});
        assert eCurvature(h, g, f){{1 < 3 => v}} == eRep(eid[One]);
        assert eMassey(h, g, f){{1 < 3 => v}} == eid[One];
    ",
        include_str!("../examples/projective-plane.asl")
    );
    stdout(interpret(&source));
}

#[test]
fn a_fourfold_system_can_choose_a_matrix_product_with_zero_adjacent_classes() {
    let source = "
        space Two { x, y };
        let f = zero[MapExt(Two, One, 1, 0)];
        let g = zero[MapExt(One, One, 0, 0)];
        let h = zero[MapExt(One, One, 1, 0)];
        let k = zero[MapExt(One, Two, 0, 0)];
        let row = mword[Two, 0, 0, One, 0, 0](cword{});
        let column = mword[One, 0, 0, Two, 0, 1](cword{});
        let rankOne = mword[Two, 0, 0, Two, 0, 1](cword{});
        assert eCurvature(k, h, g, f){} == zero[MapCobar(Two, Two, 0, 0)];
        assert eCurvature(k, h, g, f){0 < 2 => row, 2 < 4 => column} == rankOne;
        let chosen = eMassey(k, h, g, f){0 < 2 => row, 2 < 4 => column};
        assert chosen == eClass(rankOne);
        assert eInBracket(k, h, g, f, chosen) == S.hi;

        let z = zero[MapExt(One, One, 1, 0)];
        let w = zero[MapExt(One, One, 0, 0)];
        let id = eRep(eid[One]);
        assert eCurvature(z, w, z, w, z){0 < 2 => id, 2 < 5 => id} == id;
        assert eCurvature(z, w, z, w, z){0 < 2 => id, 2 < 5 => id,
                                        0 < 3 => id, 3 < 5 => id}
            == zero[MapCobar(One, One, 0, 0)];
    ";
    assert!(stdout(interpret(source)).is_empty());
}

#[test]
fn invalid_intervals_types_and_equations_report_their_source_spans() {
    let cases = [
        ("show eCurvature(h0){};", "at least two", 6),
        (
            "show eCurvature(h0, h1, h0){\n  0 < 4 => b};",
            "invalid interval [0,4]",
            3,
        ),
        (
            "show eCurvature(h0, h1, h0){\n  0 < 1 => b};",
            "adjacent interval [0,1]",
            3,
        ),
        (
            "show eMassey(h0, h1, h0){\n  0 < 3 => b};",
            "top interval [0,3]",
            3,
        ),
        (
            "show eMassey(h0, h1, h0){0 < 2 => u,\n  0 < 2 => b};",
            "duplicate witness for interval [0,2]",
            3,
        ),
        (
            "show eMassey(h0, h1, h0){\n  0 < 2 => zero[MapExt(One, One, 1, 3)]};",
            "expected type MapCobar(One, One, 1, 3)",
            12,
        ),
        (
            "show eMassey(h0, h1, h0){};",
            "equation fails on interval [0,2]",
            6,
        ),
        (
            "show eCurvature(h0, h1, h0){0 < 2 => u,\n  1 < 3 => a};",
            "equation fails on interval [1,3]",
            3,
        ),
        (
            "show eCurvature(eid[One],\n  a){};",
            "nonzero differential on interval [0,1]",
            3,
        ),
        (
            "show eCurvature(eid[S], h0){};",
            "requires composable map endpoints",
            17,
        ),
        (
            "show eCurvature(S.hi, h0){};",
            "requires MapExt or MapCobar inputs",
            17,
        ),
    ];
    for (program, message, column) in cases {
        let source = format!("{SPHERE}\n{program}");
        let error = diagnostic(interpret(&source));
        let line = source.lines().count();
        let location = format!(":{line}:{column}:");
        assert!(
            error.contains(message) && error.contains(&location),
            "{program}\n{error}"
        );
    }
}

#[test]
fn omitted_negative_intervals_and_wide_signed_sums_are_valid() {
    let source = format!(
        "let a = zero[MapExt(One, One, 1, {})];
        let b = zero[MapCobar(One, One, 0, 1)];
        let c = zero[MapCobar(One, One, 0, -1)];
        assert eCurvature(c, b, a){{}} == zero[MapCobar(One, One, 0, {})];
        assert eMassey(c, b, a){{}} == zero[MapExt(One, One, 0, {})];",
        i64::MAX,
        i64::MAX,
        i64::MAX,
    );
    assert!(stdout(interpret(&source)).is_empty());
    let negative =
        format!("{source}\nshow eCurvature(c, b, a){{1 < 3 => zero[MapCobar(One, One, 0, 0)]}};");
    assert!(diagnostic(interpret(&negative)).contains("negative cobar length"));
    let wide_interval =
        format!("{source}\nshow eCurvature(c, b, a){{0 < 2 => zero[MapCobar(One, One, 0, 0)]}};");
    assert!(diagnostic(interpret(&wide_interval)).contains("witness internal degree"));
    let invalid = format!("{source}\nshow eCurvature(b, a){{}};");
    assert!(diagnostic(interpret(&invalid)).contains("internal degree is out of range"));
    let negative = format!("{source}\nshow eMassey(c, b, b){{}};");
    assert!(diagnostic(interpret(&negative)).contains("negative cohomological degree"));
    let missing = "show eCurvature(zero[MapExt(One, One, 1, 0)], eid[One], eid[One]){};";
    assert!(diagnostic(interpret(missing)).contains("equation fails on interval [0,2]"));
}

#[test]
fn sparse_high_arity_curvature_uses_only_nonzero_entries_and_needs_no_final_basis() {
    let arity = 8192;
    let half = arity / 2;
    let inputs = vec!["z"; arity].join(",");
    let source = format!(
        "let z = zero[MapCobar(One, One, 1, 1)];
        let q = mword[One, 0, 0, One, 0, 0](cword{{[{half}]}});
        assert eCurvature({inputs}){{}} == zero[MapCobar(One, One, 2, {arity})];
        assert eCurvature({inputs}){{0 < {half} => q, {half} < {arity} => q}} == ecomp(q, q);"
    );
    assert!(stdout(interpret(&source)).is_empty());
}
