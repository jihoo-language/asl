mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn quadratic_real_multiplication_has_the_actual_relations_and_pairing() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-verbitsky-hodge-ring")),
        concat!(
            "rvVerbitskyRank[E, 2, 0] = 1\n",
            "rvVerbitskyRank[E, 2, 2] = 1\n",
            "rvVerbitskyRank[E, 2, 4] = 3\n",
            "rvVerbitskyRank[E, 2, 6] = 1\n",
            "rvVerbitskyRank[E, 2, 8] = 1\n",
            "rvVerbitskyRelations[E, 2, 6] = [{-1/2*D[0]*Q[0], 1*D[0]*Q[1]}, {-3/11*D[0]*Q[0], 1*D[0]^3}]\n",
            "rvVerbitskyPair = 220\n",
            "rvVerbitskySignature[4] = full=[3, 0, 0] primitive=[2, 0, 0]\n",
        )
    );
}

#[test]
fn canonical_classes_reduce_noncoordinate_products_and_cohomology_equality() {
    let formal = diagnostic(interpret(
        r"
        rvhodgewick E[11]{-1,-1,1}{[2]};
        let d=rvHodgeD[E,2,0];
        let center=rvHodgeQ[E,2]{-1,2};
        assert rvHodgeCup(d,center)==zero[RVHodgeTensor(E,2,6)];
    ",
    ));
    assert!(formal.contains("assert"), "{formal}");
    let result = stdout(interpret(
        r"
        rvhodgewick E[11]{-1,-1,1}{[2]};
        let d=rvHodgeD[E,2,0];
        let center=rvHodgeQ[E,2]{-1,2};
        let relation=rvHodgeCup(d,center);
        assert rvVerbitskyClass(relation)==zero[RVVerbitsky(E,2,6)];
        let D=rvVerbitskyClass(d);
        let p=rvVerbitskyPrimitive[4,1](D);
        assert rvVerbitskyClass(rvVerbitskyRepresentative(p))==p;
        assert rvVerbitskyCup(D,p)==zero[RVVerbitsky(E,2,6)];
        let D3=rvVerbitskyCup(rvVerbitskyCup(D,D),D);
        let Q=rvVerbitskyBasis[E,2,6,0];
        assert rvVerbitskyScale[11,3](D3)==Q;
        let P=rvVerbitskyClass(center);
        let u=rvVerbitskyBasis[E,2,0,0];
        assert rvVerbitskyCup(u,P)==P;
        let theta=rvVerbitskyClass(rvHodgeQ[E,2]{0,1});
        assert rvVerbitskyCup(rvVerbitskyCup(D,theta),D)==rvVerbitskyCup(D,rvVerbitskyCup(theta,D));
        let top=rvVerbitskyCup(P,P);
        assert rvVerbitskyCup(top,D)==zero[RVVerbitsky(E,2,10)];
        rvVerbitskyIntegral(top);
        rvHodgeIntegral(rvVerbitskyRepresentative(top));
        ",
    ));
    assert_eq!(result, "rvVerbitskyIntegral = 220\nrvHodgeIntegral = 220\n");
}

#[test]
fn primitive_form_accepts_a_noncoordinate_marked_lefschetz_divisor() {
    let text = stdout(interpret(
        r"
        rvhodgewick B[21]{0,1}{[2,0],[0,-2]};
        let D=rvVerbitskyClass(rvHodgeD[B,2,0]);
        let E=rvVerbitskyClass(rvHodgeD[B,2,1]);
        let h=add(D,rvVerbitskyScale[1,2](E));
        rvVerbitskySignature[2](D);
        rvVerbitskySignature[2](h);
        rvVerbitskySignature[4](h);
        let p=rvVerbitskyPrimitive[2,0](h);
        assert p==add(E,rvVerbitskyScale[1,2](D));
        let h2=rvVerbitskyCup(h,h);
        let h3=rvVerbitskyCup(h2,h);
        assert rvVerbitskyCup(h3,p)==zero[RVVerbitsky(B,2,8)];
        rvVerbitskyPair(rvVerbitskyCup(h2,p),p);

        rvhodgewick S[21]{0,1}{[2]};
        let h1=rvVerbitskyClass(rvHodgeD[S,1,0]);
        assert rvVerbitskyClass(rvHodgeQ[S,1]{1})==rvVerbitskyScale[21,2](rvVerbitskyCup(h1,h1));
        rvVerbitskySignature[2](h1);
        ",
    ));
    assert_eq!(
        text,
        concat!(
            "rvVerbitskySignature[2] = full=[1, 1, 0] primitive=[0, 1, 0]\n",
            "rvVerbitskySignature[2] = full=[1, 1, 0] primitive=[0, 1, 0]\n",
            "rvVerbitskySignature[4] = full=[3, 1, 0] primitive=[2, 0, 0]\n",
            "rvVerbitskyPair = -9/4\n",
            "rvVerbitskySignature[2] = full=[1, 0, 0] primitive=[0, 0, 0]\n",
        )
    );
}

#[test]
fn graded_ranks_match_the_verbitsky_irrep_invariant_dimensions() {
    let mut program = "rvhodgewick E[11]{-1,-1,1}{[2]};\n".to_owned();
    for n in 2..=7 {
        for k in 0..=2 * n {
            let d = 2 * k;
            writeln!(
                program,
                "rvVerbitskyRank[E,{n},{d}]; rvHilbLLVHodgeInvariant[{n},{d},1,real,2]{{{n}}};"
            )
            .unwrap();
        }
    }
    let text = stdout(interpret(&program));
    let lines = text.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 2 * (2..=7).map(|n| 2 * n + 1).sum::<usize>());
    for pair in lines.chunks_exact(2) {
        let rank = pair[0].split_once(" = ").unwrap().1;
        let contribution = pair[1].rsplit_once(" contribution=").unwrap().1;
        assert_eq!(rank, contribution);
    }
}

#[test]
fn malformed_marks_grades_and_isotropic_polarizations_are_diagnosed() {
    let tests = [
        ("rvVerbitskyRank[N,2,4];", "unknown Hodge Wick marking"),
        ("rvhodgewick B[21]{0,1}{[2,0],[0,-2]}; rvVerbitskyRank[B,1,2];", "requires total Beauville rank 22"),
        ("rvhodgewick B[22]{0,1}{[2]}; rvVerbitskyRank[B,2,1];", "nonnegative even"),
        ("rvhodgewick B[22]{0,1}{[2]}; rvVerbitskyPair(rvVerbitskyBasis[B,2,2,0],rvVerbitskyBasis[B,2,4,0]);", "complementary degrees"),
        ("rvhodgewick B[22]{0,1}{[2]}; rvVerbitskySignature[6](rvVerbitskyBasis[B,2,2,0]);", "at most the middle"),
        ("rvhodgewick B[21]{0,1}{[2,0],[0,-2]}; rvVerbitskySignature[4](add(rvVerbitskyBasis[B,2,2,0],rvVerbitskyBasis[B,2,2,1]));", "nonisotropic"),
        ("rvhodgewick B[22]{0,1}{[2]}; rvVerbitskySignature[4](rvVerbitskyBasis[B,2,4,0]);", "needs a degree-two divisor"),
        ("rvhodgewick B[22]{0,1}{[2]}; rvVerbitskyIntegral(rvHodgeInverseBB[B,2]);", "RVVerbitsky class"),
        ("rvhodgewick B[22]{0,1}{[2]}; rvVerbitskyRank[B,1000000000000000,2000000000000000];", "storage limit"),
    ];
    for (program, reason) in tests {
        assert!(diagnostic(interpret(program)).contains(reason), "{program}");
    }
    let valid=stdout(interpret("rvhodgewick B[22]{0,1}{[2]}; rvVerbitskyRank[B,1000000000000000,4000000000000000]; rvVerbitskyRank[B,2,10];"));
    assert_eq!(valid,"rvVerbitskyRank[B, 1000000000000000, 4000000000000000] = 1\nrvVerbitskyRank[B, 2, 10] = 0\n");
}

#[test]
fn higher_lefschetz_power_and_one_dimensional_etale_blocks() {
    let mut program = r"
        rvhodgewick E[11]{-1,-1,1}{[2]};
        let h=rvVerbitskyClass(rvHodgeD[E,3,0]);
        let p=rvVerbitskyClass(rvHodgeQ[E,3]{-1,2});
        let h2=rvVerbitskyCup(h,h);
        let h3=rvVerbitskyCup(h2,h);
        assert rvVerbitskyCup(h3,p)==zero[RVVerbitsky(E,3,10)];
        rvVerbitskyPair(rvVerbitskyCup(h2,p),p);
        rvVerbitskySignature[4](h);
    "
    .to_owned();
    let rows = (0..21)
        .map(|i| {
            let row = (0..21)
                .map(|j| {
                    if i == j {
                        if i == 0 {
                            "2"
                        } else {
                            "-2"
                        }
                    } else {
                        "0"
                    }
                })
                .collect::<Vec<_>>()
                .join(",");
            format!("[{row}]")
        })
        .collect::<Vec<_>>()
        .join(",");
    writeln!(program, "rvhodgewick S[1]{{-1,0,1}}{{{rows}}};").unwrap();
    program.push_str(
        r"
        rvVerbitskyRank[S,2,4];
        let r=rvVerbitskyClass(rvHodgeQ[S,2]{(1,2),(1,2)});
        let s=rvVerbitskyClass(rvHodgeQ[S,2]{(1,2),(-1,2)});
        rvVerbitskyPair(r,r);
        rvVerbitskyPair(r,s);
    ",
    );
    assert_eq!(
        stdout(interpret(&program)),
        concat!(
            "rvVerbitskyPair = 440\n",
            "rvVerbitskySignature[4] = full=[3, 0, 0] primitive=[2, 0, 0]\n",
            "rvVerbitskyRank[S, 2, 4] = 233\n",
            "rvVerbitskyPair = 3\n",
            "rvVerbitskyPair = 1\n",
        )
    );
}
