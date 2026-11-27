mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn signed_hom_differentials_composition_and_higher_coefficients() {
    assert_eq!(stdout(example("mixed-hom")), "MixedU(Z, 0){u^4*q}\n",);
}

#[test]
fn mixed_hom_coordinates_and_degree_bounds_are_checked() {
    let fixture = include_str!("../examples/mixed-hom.asl");
    for (extra, expected) in [
        (
            "show mixunit[End,a,unknown];",
            "unknown mixed complex generator",
        ),
        ("show mixunit[V,a,b];", "defined by mixhom"),
        (
            "show mixEval(f,mixchain[X,0][(x,1)]);",
            "wrong source complex",
        ),
        ("show mixEval(f,mixU[0](a));", "same form"),
        ("show mixComp[End](mixU[0](f),f);", "same form"),
        (
            "show mixComp[XX](post,map);",
            "endpoint complexes do not match",
        ),
        ("show mixId[XY];", "requires a mixed endomorphism complex"),
        (
            "mixedcomplex W[7] {z:0}{}{}; mixhom ZW[Z,W];",
            "same prime field",
        ),
        (
            "mixedcomplex L[5] {z:-9223372036854775808}{}{}; mixhom LZ[L,Z];",
            "basis degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{fixture} {extra}")));
        assert!(error.contains(expected), "{extra}: {error}");
    }
    let gens = (0..300)
        .map(|i| format!("x{i}:0"))
        .collect::<Vec<_>>()
        .join(",");
    let oversized = format!("mixedcomplex L[2] {{{gens}}} {{}} {{}}; mixhom H[L,L];");
    assert!(
        diagnostic(interpret(&oversized)).contains("mixed Hom complex exceeds the storage limit")
    );
}

#[test]
fn an_empty_hom_complex_still_has_the_correct_actions_and_identity() {
    let source = r"
        mixedcomplex A[2] {} {} {};
        mixedcomplex B[2] {b:0} {} {};
        mixhom AB[A,B];
        mixhom BA[B,A];
        mixhom EE[A,A];
        assert mixId[EE] == zero[MixedChain(EE,0)];
        assert mixEval(zero[MixedChain(BA,0)],mixchain[B,0][(b,1)]) == zero[MixedChain(A,0)];
        assert mixComp[EE](zero[MixedChain(BA,0)],zero[MixedChain(AB,0)]) == mixId[EE];
        mixbars[AB];
        ";
    assert_eq!(
        stdout(interpret(source)),
        "mixbars[AB] = free [], torsion []\n"
    );
}

#[test]
fn corrected_map_has_the_actual_staircase_obstruction_or_a_closed_lift() {
    for (last_column, twisted) in [("", false), (",x2=>[(y3,1)]", true)] {
        let mut source = format!(
            r"mixedcomplex O[3] {{ a:0 }} {{}} {{}};
            mixedcomplex S3[3] {{ x0:0,x1:2,x2:4,y1:1,y2:3,y3:5 }}
                {{ x1=>[(y1,-1)],x2=>[(y2,-1)] }}
                {{ x0=>[(y1,1)],x1=>[(y2,1)]{last_column} }};
            mixhom H[O,S3];
            let lead = mixgen[H,3,0,0,0];
            let f = mixPageRep(lead);
            let a = mixU[0](mixchain[O,0][(a,1)]);
            let y = mixU[3](mixchain[S3,5][(y3,1)]);
            assert f == add(add(mixU[0](mixunit[H,a,x0]),mixU[1](mixunit[H,a,x1])),mixU[2](mixunit[H,a,x2]));
            assert mixD(mixEval(f,a)) == mixEval(mixD(f),a);
            ",
        );
        if twisted {
            source.push_str(
                "assert mixEval(mixD(f),a) == y;
                 assert mixD(lead) == mixgen[H,3,-1,3,0];
                 show mixD(lead);",
            );
        } else {
            source.push_str(
                "assert mixD(f) == zero[MixedU(H,-1)];
                 assert mixNext(lead) == mixgen[H,4,0,0,0];
                 show mixNext(lead);",
            );
        }
        assert_eq!(
            stdout(interpret(&source)),
            if twisted {
                "MixedE(H, 3, -1, 3){g0}\n"
            } else {
                "MixedE(H, 4, 0, 0){g0}\n"
            }
        );
    }
}
