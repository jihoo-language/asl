mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn projective_four_space_detects_the_localized_torsion_relation() {
    assert_eq!(
        stdout(example("v-koszul")),
        concat!(
            "vextdim[One, RP4, [0, 1], 0, 0] = 1\n",
            "vextdim[One, RP4, [0, 1], 0, 1] = 1\n",
            "vextdim[One, RP4, [0, 1], 0, 2] = 0\n",
            "vextdim[One, RP4, [0, 1], 0, 3] = 1\n",
            "vextdim[One, RP4, [0, 1], 0, 4] = 0\n"
        )
    );
    let program = format!(
        "{}\n{}",
        include_str!("../examples/v-koszul.asl"),
        "vextdim[One, RP4, [1], 0, 1];
         vextdim[One, RP4, [1], 0, 2];
         vextdim[One, RP4, [1], -1, -1];
         vextdim[One, RP4, [1], -1, 0];"
    );
    let output = stdout(interpret(&program));
    assert!(
        output.ends_with(concat!(
            "vextdim[One, RP4, [1], 0, 1] = 0\n",
            "vextdim[One, RP4, [1], 0, 2] = 1\n",
            "vextdim[One, RP4, [1], -1, -1] = 1\n",
            "vextdim[One, RP4, [1], -1, 0] = 1\n"
        )),
        "{output}"
    );
}

#[test]
fn localization_is_functorial_and_top_powers_are_invertible() {
    let program = r"
        space Two { a, b };
        let left = qRes[1](emap(\x : Two => One.pt));
        let right = qRes[1](emap(\x : One => Two.a));
        assert ecomp(vLocal(left), vLocal(right)) == vLocal(ecomp(left, right));
        assert ecomp(vLocal(left), vLocal(right)) == vLocal(qRes[1](eid[One]));

        let u = vLocal(qRes[1](eid[One]));
        let lower = vShift[-1](u);
        let higher = vShift[1](u);
        assert ecomp(higher, lower) == u;
        assert ecomp(lower, higher) == u;
        assert lower == vextgen[One, One, [1], -1, -3, 0];
        assert qClass(qRep(lower)) == lower;
        assert vShift[-1](qRep(u)) == qRep(lower);
        assert ecomp(qRep(lower), qRep(higher)) == qRep(u);
        assert qD(qRep(lower)) == zero[VKoszul(One, One, [1], 0, -3)];
        let translate = \x : VKoszul(One, One, [1], -1, -3) => add(x, qRep(lower));
        assert translate(translate(qRep(lower))) == qRep(lower);
        assert vShift[0](u) == u;
        assert vShift[2](ecomp(vLocal(left), vLocal(right))) ==
            ecomp(vShift[1](vLocal(left)), vShift[1](vLocal(right)));
        show qRep(lower);
        vextdim[One, One, [1], -1, -2];
        vextdim[One, One, [0, 1, 2], 0, -12];
        betti VExt(One, One, [0, 1, 2], 0, -12);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "VKoszul(One, One, [1], -1, -3){[0, 0 -> 0, 0] v1^-1}\n",
            "vextdim[One, One, [1], -1, -2] = 0\n",
            "vextdim[One, One, [0, 1, 2], 0, -12] = 2\n",
            "betti(VExt(One, One, [0, 1, 2], 0, -12)) = [4]\n"
        )
    );
}

#[test]
fn signed_extremes_are_sparse_and_checked_after_widening() {
    let min = i64::MIN;
    let max = i64::MAX;
    let program = format!(
        "
        vextdim[One, One, [0], {min}, {min}];
        vextdim[One, One, [0, 1], {min}, {min}];
        show qRep(vextgen[One, One, [0, 1], {max}, -{max}, 0]);
        let base = vextgen[One, One, [1], -1, -3, 0];
        assert vShift[3074457345618258603](base) ==
            vextgen[One, One, [1], 3074457345618258602, 9223372036854775806, 0];
        "
    );
    assert_eq!(
        stdout(interpret(&program)),
        format!(
            "vextdim[One, One, [0], {min}, {min}] = 1\n\
             vextdim[One, One, [0, 1], {min}, {min}] = 0\n\
             VKoszul(One, One, [0, 1], {max}, -{max}){{[0, 0 -> 0, 0] v0^18446744073709551614 v1^-{max}}}\n"
        )
    );
}

#[test]
fn localized_families_check_indices_endpoints_and_indices_before_reduction() {
    let max = i64::MAX;
    let min = i64::MIN;
    let programs = [
        ("show vLocal(qRes[](cword{}));".into(), "nonempty"),
        ("vextdim[One, One, [], 0, 0];".into(), "nonempty"),
        ("vextdim[One, One, [1, 0], 0, 0];".into(), "increasing"),
        ("show vShift[-1](qRes[0](cword{}));".into(), "VKoszul"),
        (
            format!("show vShift[-1](zero[VExt(One, One, [0], {min}, 0)]);"),
            "out of range",
        ),
        (
            format!("show vShift[1](zero[VExt(One, One, [1], 0, {max})]);"),
            "out of range",
        ),
        (
            format!("show qD(zero[VKoszul(One, One, [0], {max}, 0)]);"),
            "out of range",
        ),
        (
            "show ecomp(vLocal(qRes[0](eid[One])), vLocal(qRes[0](eid[S])));".into(),
            "endpoints",
        ),
        (
            "show ecomp(vLocal(qRes[0](cword{})), qRes[0](cword{}));".into(),
            "family",
        ),
        (
            "show ecomp(vLocal(qRes[0](cword{})), vLocal(qRes[1](cword{})));".into(),
            "same Milnor primitives",
        ),
    ];
    for (program, expected) in programs {
        let error = diagnostic(interpret(&program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
