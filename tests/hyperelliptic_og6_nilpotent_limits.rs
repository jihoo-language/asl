mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn setup() -> &'static str {
    include_str!("../examples/hyperelliptic-og6-nilpotent-limits.asl")
        .split_once("// A central")
        .unwrap()
        .0
}

const DECLARATIONS: &str = r"
    rvog6nilpotentii A[P,3];
    rvog6nilpotentiii B3[P,0](u);
";

#[test]
fn runnable_example_computes_both_types_and_their_image_indices() {
    assert_eq!(
        stdout(example("hyperelliptic-og6-nilpotent-limits")),
        concat!(
            "rvOG6NilpotentData[A] = type=II, H_rank=7, graded=[0, 2, 3, 2, 0], imN_smith=[3, 3], imN2_smith=[], imN_index=9, imN2_index=1\n",
            "rvOG6NilpotentData[B3] = type=III, H_rank=7, graded=[1, 0, 5, 0, 1], imN_smith=[1, 1], imN2_smith=[2], imN_index=1, imN2_index=2\n",
            "rvOG6NilpotentInWeight[A, 1] = true\n",
            "rvOG6NilpotentInWeight[B3, 1] = false\n",
            "rvOG6NilpotentInWeight[B3, 2] = true\n",
        )
    );
}

#[test]
fn smith_uses_the_polarized_lattice_when_an_axis_has_divisibility_two() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{}
            let h2=add(e,add(h,h));
            rvog6cusp Q(h2,e,f);
            rvog6nilpotentii C2[Q,2];
            rvog6nilpotentiii C3[Q,0](u);
            rvog6nilpotentiii Triple[Q,0](add(u,add(u,u)));
            let x=add(add(eDual,eDual),neg(hF));
            rvOG6Pair(x,h2);
            rvOG6Pair(x,e);
            assert rvOG6NilpotentPower[C2,1](x)==add(add(f,f),add(f,f));
            assert rvOG6NilpotentPower[C2,1](fDual)==neg(add(e,e));
            assert rvOG6NilpotentPower[C3,1](fDual)==neg(e);
            assert rvOG6NilpotentPower[C3,2](x)==neg(add(add(e,e),add(e,e)));
            rvOG6NilpotentData[C2];
            rvOG6NilpotentData[C3];
            rvOG6NilpotentData[Triple];
            ",
            setup()
        ))),
        concat!(
            "rvOG6Pair = 0\n",
            "rvOG6Pair = 2\n",
            "rvOG6NilpotentData[C2] = type=II, H_rank=7, graded=[0, 2, 3, 2, 0], imN_smith=[2, 4], imN2_smith=[], imN_index=8, imN2_index=1\n",
            "rvOG6NilpotentData[C3] = type=III, H_rank=7, graded=[1, 0, 5, 0, 1], imN_smith=[1, 2], imN2_smith=[4], imN_index=2, imN2_index=4\n",
            "rvOG6NilpotentData[Triple] = type=III, H_rank=7, graded=[1, 0, 5, 0, 1], imN_smith=[3, 6], imN2_smith=[36], imN_index=18, imN2_index=36\n",
        )
    );
}

#[test]
fn saturated_weights_separate_mixed_classes_and_ignore_image_index() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{}{DECLARATIONS}
            let mix=add(e,delta);
            let transverse=add(mix,eDual);
            rvOG6NilpotentInWeight[A,0](e);
            rvOG6NilpotentInWeight[A,1](e);
            rvOG6NilpotentInWeight[A,1](mix);
            rvOG6NilpotentInWeight[A,2](mix);
            rvOG6NilpotentInWeight[A,2](transverse);
            rvOG6NilpotentInWeight[A,3](transverse);
            rvOG6NilpotentInWeight[B3,0](e);
            rvOG6NilpotentInWeight[B3,1](mix);
            rvOG6NilpotentInWeight[B3,2](mix);
            rvOG6NilpotentInWeight[B3,3](transverse);
            rvOG6NilpotentInWeight[B3,4](transverse);
            rvOG6NilpotentInWeight[B3,-1](e);
            rvOG6NilpotentInWeight[B3,-4](zero[RVOG6(M)]);
            rvOG6NilpotentInWeight[B3,19](transverse);
            ",
            setup()
        ))),
        concat!(
            "rvOG6NilpotentInWeight[A, 0] = false\n",
            "rvOG6NilpotentInWeight[A, 1] = true\n",
            "rvOG6NilpotentInWeight[A, 1] = false\n",
            "rvOG6NilpotentInWeight[A, 2] = true\n",
            "rvOG6NilpotentInWeight[A, 2] = false\n",
            "rvOG6NilpotentInWeight[A, 3] = true\n",
            "rvOG6NilpotentInWeight[B3, 0] = true\n",
            "rvOG6NilpotentInWeight[B3, 1] = false\n",
            "rvOG6NilpotentInWeight[B3, 2] = true\n",
            "rvOG6NilpotentInWeight[B3, 3] = false\n",
            "rvOG6NilpotentInWeight[B3, 4] = true\n",
            "rvOG6NilpotentInWeight[B3, -1] = false\n",
            "rvOG6NilpotentInWeight[B3, -4] = true\n",
            "rvOG6NilpotentInWeight[B3, 19] = true\n",
        )
    );
}

#[test]
fn nilpotent_powers_and_exponentials_match_center_and_eichler_laws() {
    let mut source = format!(
        r"{}{DECLARATIONS}
        rvog6nilpotentiii Shifted[P,0](add(u,e));
        rvog6nilpotentiii OnF[P,1](add(e,eDual));
        assert rvOG6NilpotentPower[OnF,2](fDual)==neg(add(f,f));
        assert rvOG6NilpotentPower[A,1](h)==zero[RVOG6(M)];
        assert rvOG6NilpotentExp[B3,5](h)==h;
        assert rvOG6NilpotentPower[B3,1](u)==neg(add(e,e));
        let a=rvOG6CuspTransvection[P,1](neg(delta),eDual);
        let b=rvOG6CuspTransvection[P,0](neg(delta),a);
        let c=rvOG6CuspTransvection[P,1](delta,b);
        let comm=rvOG6CuspTransvection[P,0](delta,c);
        rvog6nilpotentii Commutator[P,-2];
        assert rvOG6NilpotentExp[Commutator,1](eDual)==comm;
        rvOG6Pair(rvOG6NilpotentPower[B3,1](eDual),u);
        rvOG6Pair(eDual,rvOG6NilpotentPower[B3,1](u));
        rvOG6Pair(rvOG6NilpotentExp[B3,-3](eDual),rvOG6NilpotentExp[B3,-3](u));
        rvOG6Pair(eDual,u);
        ",
        setup(),
    );
    for index in 0..8 {
        writeln!(source, r"
            let x{index}=rvOG6Basis[M,{index}];
            assert rvOG6NilpotentPower[B3,0](x{index})==x{index};
            assert rvOG6NilpotentPower[A,2](x{index})==zero[RVOG6(M)];
            assert rvOG6NilpotentPower[B3,2](x{index})==rvOG6NilpotentPower[B3,1](rvOG6NilpotentPower[B3,1](x{index}));
            assert rvOG6NilpotentPower[B3,4](x{index})==zero[RVOG6(M)];
            assert rvOG6NilpotentPower[Shifted,1](x{index})==rvOG6NilpotentPower[B3,1](x{index});
            assert rvOG6NilpotentExp[Shifted,-2](x{index})==rvOG6NilpotentExp[B3,-2](x{index});
            assert rvOG6NilpotentExp[B3,2](x{index})==rvOG6NilpotentExp[B3,-3](rvOG6NilpotentExp[B3,5](x{index}));
            assert rvOG6NilpotentExp[A,4](x{index})==rvOG6CuspCentral[P,12](x{index});
        ").unwrap();
    }
    assert_eq!(
        stdout(interpret(&source)),
        "rvOG6Pair = 2\nrvOG6Pair = -2\nrvOG6Pair = 0\nrvOG6Pair = 0\n"
    );
}

#[test]
fn huge_smith_indices_and_shifted_logarithms_remain_exact() {
    let big = i128::MAX;
    let min = i128::MIN;
    let neg_big = -big;
    assert_eq!(
        stdout(interpret(&format!(
            r"{}{DECLARATIONS}
            rvog6nilpotentii Huge[P,{min}];
            rvOG6NilpotentData[Huge];
            let multiple=rvOG6[M,0](rvAbMukai[B,0,0]{{{big},0,0,0,0,0}});
            let dualMultiple=rvOG6[M,0](rvAbMukai[B,0,0]{{0,0,0,0,0,{neg_big}}});
            rvog6nilpotentiii Shift[P,0](add(u,multiple));
            assert rvOG6NilpotentPower[Shift,1](dualMultiple)==rvOG6NilpotentPower[B3,1](dualMultiple);
            assert rvOG6NilpotentExp[B3,{min}](delta)==delta;
            assert rvOG6NilpotentExp[Huge,0](eDual)==eDual;
            ",
            setup()
        ))),
        concat!(
            "rvOG6NilpotentData[Huge] = type=II, H_rank=7, graded=[0, 2, 3, 2, 0], ",
            "imN_smith=[170141183460469231731687303715884105728, 170141183460469231731687303715884105728], ",
            "imN2_smith=[], imN_index=28948022309329048855892746252171976963317496166410141009864396001978282409984, imN2_index=1\n",
        )
    );
}

#[test]
fn rejects_invalid_directions_and_unpolarized_weights() {
    let source = format!("{}{DECLARATIONS}", setup());
    for (input, fragment) in [
        ("rvog6nilpotentii Bad[P,0];", "requires a nonzero central coefficient"),
        ("rvog6nilpotentiii Bad[P,2](u);", "axis must be 0 or 1"),
        ("rvog6nilpotentiii Bad[P,0](eDual);", "perpendicular to the polarization and the selected axis"),
        ("rvog6nilpotentiii Bad[P,0](h);", "perpendicular to the polarization and the selected axis"),
        ("rvog6nilpotentiii Bad[P,0](delta);", "requires a direction of positive square"),
        ("rvog6nilpotentiii Bad[P,0](f);", "requires a direction of positive square"),
        ("rvog6nilpotentii A[P,1];", "is already defined"),
        ("rvog6nilpotentii Bad[Missing,1];", "unknown OG6 cusp"),
        ("rvOG6NilpotentData[Missing];", "unknown OG6 nilpotent"),
        ("rvOG6NilpotentInWeight[B3,4](h);", "requires a class perpendicular to the polarization"),
        ("rvog6 M2(rvAbMukai[B,1,-1]{0,0,0,0,0,0}); show rvOG6NilpotentPower[A,1](zero[RVOG6(M2)]);", "expected type RVOG6(M)"),
        ("rvog6 M2(rvAbMukai[B,1,-1]{0,0,0,0,0,0}); rvog6nilpotentiii Bad[P,0](zero[RVOG6(M2)]);", "expected type RVOG6(M)"),
        ("show rvOG6NilpotentExp[B3,170141183460469231731687303715884105727](eDual);", "abelian Mukai coordinate exceeds"),
    ] {
        let error = diagnostic(interpret(&format!("{source}{input}")));
        assert!(error.contains(fragment), "{error}");
    }
}
