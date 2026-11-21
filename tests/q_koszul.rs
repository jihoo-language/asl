mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn sphere_koszul_products_and_finite_steenrod_restrictions() {
    assert_eq!(
        stdout(example("q-koszul")),
        concat!(
            "qextdim[One, One, [0, 1, 2], 3, 9] = 2\n",
            "QExt(One, One, [1], 4, 12){g0}\n",
            "qextdim[Two, Two, [], 0, 0] = 4\n",
        )
    );
}

#[test]
fn projective_four_space_detects_both_milnor_primitives_at_both_endpoints() {
    // Identify antipodal faces of the five-dimensional cross-polytope.
    fn face(axes: u32, negative: u32) -> String {
        let mut negative = negative & axes;
        if negative & (1 << axes.trailing_zeros()) != 0 {
            negative ^= axes;
        }
        format!("f{axes}_{negative}")
    }
    let mut relations = Vec::new();
    for axes in 1u32..32 {
        let first = 1 << axes.trailing_zeros();
        for negative in 0..32 {
            if negative & !axes != 0 || negative & first != 0 {
                continue;
            }
            for axis in (0..5)
                .map(|i| 1 << i)
                .filter(|bit| axes & bit != 0 && axes != *bit)
            {
                relations.push(format!(
                    "{} < {}",
                    face(axes ^ axis, negative),
                    face(axes, negative)
                ));
            }
        }
    }
    let program = format!(
        "space P {{ {} }};\n{}",
        relations.join(", "),
        r"
        let topRaw = mword[One, 0, 0, P, 4, 0](cword{});
        let top = qRes[0, 1](topRaw);
        let topD = qRes[0, 1](add(mword[One, 0, 0, P, 3, 0](cword{[1]}),
                                  mword[One, 0, 0, P, 1, 0](cword{[0, 1]})));
        assert qD(top) == topD;
        assert topD == qRes[0, 1](eD(topRaw));
        assert qD(qRes[1](top)) == qRes[1](topD);
        assert qD(qRes[0](top)) == qRes[0](topD);

        let downRaw = mword[P, 1, 0, One, 0, 0](cword{});
        let down = qRes[0, 1](downRaw);
        let downD = qRes[0, 1](add(mword[P, 2, 0, One, 0, 0](cword{[1]}),
                                   mword[P, 4, 0, One, 0, 0](cword{[0, 1]})));
        assert qD(down) == downD;
        assert downD == qRes[0, 1](eD(downRaw));
        assert ecomp(downD, top) == qRes[0, 1](cword{[0, 1]});
        assert ecomp(down, topD) == qRes[0, 1](cword{[0, 1]});
        assert qD(ecomp(down, top)) == add(ecomp(downD, top), ecomp(down, topD));
        let both = qRes[0, 1](mword[P, 1, 0, P, 4, 0](cword{}));
        assert qD(qD(both)) == zero[QKoszul(P, P, [0, 1], 2, 3)];
        qextdim[One, P, [0, 1], 0, 4];
        qextdim[P, One, [0, 1], 0, -1];
        qextdim[P, One, [], 0, -1];
    "
    );
    assert_eq!(
        stdout(interpret(&program)),
        concat!(
            "qextdim[One, P, [0, 1], 0, 4] = 0\n",
            "qextdim[P, One, [0, 1], 0, -1] = 0\n",
            "qextdim[P, One, [], 0, -1] = 1\n",
        )
    );
}

#[test]
fn raw_and_class_spaces_use_finite_functions_and_signed_grades() {
    let program = r"
        let u = qRes[0, 1](cword{[1], [0, 1]});
        let translate = \x : QKoszul(One, One, [0, 1], 2, 4) => add(x, u);
        assert translate(u) == zero[QKoszul(One, One, [0, 1], 2, 4)];
        let identity = \x : QExt(One, One, [0, 1, 2], 3, 9) => qClass(qRep(x));
        assert identity == (\x : QExt(One, One, [0, 1, 2], 3, 9) => x);
        betti QKoszul(One, One, [0, 1], 2, 4);
        betti QExt(One, One, [0, 1, 2], 3, 9);
    ";
    assert_eq!(stdout(interpret(program)), "betti(QKoszul(One, One, [0, 1], 2, 4)) = [2]\nbetti(QExt(One, One, [0, 1, 2], 3, 9)) = [4]\n");
    let large = (i64::MAX as u128).min(usize::MAX as u128);
    let small = i64::MIN;
    assert_eq!(stdout(interpret(&format!(
        "qextdim[One, One, [0], {large}, {large}];\nqextdim[One, One, [0], 3, {small}];\nqextdim[One, One, [0, 2, 4], 10000000, 150000002];"
    ))), format!("qextdim[One, One, [0], {large}, {large}] = 1\nqextdim[One, One, [0], 3, {small}] = 0\nqextdim[One, One, [0, 2, 4], 10000000, 150000002] = 0\n"));
}

#[test]
fn incompatible_primitive_selections_and_koszul_bidegrees_are_rejected() {
    for (program, expected) in [
        ("show qRes[1, 0](cword{});".into(), "strictly increasing"),
        (
            "qextdim[One, One, [0, 0], 0, 0];".into(),
            "strictly increasing",
        ),
        (
            format!("qextdim[One, One, [{}], 0, 0];", usize::BITS),
            "degree is out of range",
        ),
        (
            "show qRes[1](resA[0](cword{}));".into(),
            "source Steenrod subalgebra",
        ),
        ("show qRes[0, 1](qRes[1](cword{}));".into(), "subset"),
        (
            "show ecomp(qRes[0](cword{}), qRes[1](cword{}));".into(),
            "same Milnor primitives",
        ),
        (
            "show ecomp(qRes[0](eid[One]), qRes[0](eid[S]));".into(),
            "composable Koszul endpoints",
        ),
        (
            format!("show qD(zero[QKoszul(One, One, [0], {}, 0)]);", usize::MAX),
            "qD degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
