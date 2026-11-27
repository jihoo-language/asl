mod common;

use common::{diagnostic, example, interpret, stdout};

const F3: &str = r"
    rvunram F[3,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    rvhyperelliptic C[F]{z,t,z,t,o,o};
";

#[test]
fn factored_zeta_and_supersingular_resolution_match_exact_values() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-k3")),
        concat!(
            "rvKummerK3Zeta[C] = 1/((1-T)(1-9*T)(1-3*T)^2(1+8*T+33*T^2+72*T^3+81*T^4)(1-3*T)^4(1-(3*T)^3)^4), node orbit lengths=[1, 1, 1, 1, 3, 3, 3, 3]\n",
            "rvKummerK3Points[C, 1] = singular=8, nodes=4, resolved=20\n",
            "rvKummerK3Points[C, 3] = singular=848, nodes=16, resolved=1280\n",
            "rvKummerK3Picard[C, 1] = rank=10, geometric rank=20\n",
            "rvKummerK3Picard[C, 3] = rank=20, geometric rank=20\n",
            "rvKummerK3Zeta[S] = 1/((1-T)(1-625*T)(1-25*T)^2(1-100*T+3750*T^2-62500*T^3+390625*T^4)(1-25*T)^16), node orbit lengths=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]\n",
            "rvKummerK3Points[S, 1] = singular=776, nodes=16, resolved=1176\n",
            "rvKummerK3Picard[S, 1] = rank=22, geometric rank=22\n",
        )
    );
}

#[test]
fn point_and_node_counts_agree_with_small_quartic_and_jacobian_enumerations() {
    let output = stdout(interpret(&format!(
        r"{F3}
        rvunram W[3,1]{{1,0,1}};
        rvjacbase L[C,W];
        rvKummerK3Points[C,1];
        rvKummerK3Points[C,2];
        rvKummerK3Points[L,1];
        show (\p:RVKummer(C)=>p);
        show (\p:RVKummer(L)=>p);
        show (\d:RVJac(C)=>add(d,d));
        show (\d:RVJac(L)=>add(d,d));
        "
    )));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(
        lines[0],
        "rvKummerK3Points[C, 1] = singular=8, nodes=4, resolved=20"
    );
    assert_eq!(
        lines[1],
        "rvKummerK3Points[C, 2] = singular=98, nodes=4, resolved=134"
    );
    assert_eq!(
        lines[2],
        "rvKummerK3Points[L, 1] = singular=98, nodes=4, resolved=134"
    );
    assert_eq!(lines[3].matches("RVKummer(C)").count(), 2 * 8);
    assert_eq!(lines[4].matches("RVKummer(L)").count(), 2 * 98);
    assert_eq!(lines[5].matches("=> RVJac(C){u=[[1]], v=[]}").count(), 4);
    assert_eq!(lines[6].matches("=> RVJac(L){u=[[1, 0]], v=[]}").count(), 4);
}

#[test]
fn picard_rank_distinguishes_cycles_from_nodes_and_handles_huge_degrees() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F3}
            rvhyperelliptic G[F]{{z,o,z,z,o,o}};
            rvKummerK3Picard[G,1];
            rvKummerK3Picard[G,3];
            rvKummerK3Picard[C,2];
            rvKummerK3Picard[C,340282366920938463463374607431768211454];
            rvKummerK3Picard[C,340282366920938463463374607431768211455];
            "
        ))),
        concat!(
            "rvKummerK3Picard[G, 1] = rank=10, geometric rank=18\n",
            "rvKummerK3Picard[G, 3] = rank=18, geometric rank=18\n",
            "rvKummerK3Picard[C, 2] = rank=10, geometric rank=20\n",
            "rvKummerK3Picard[C, 340282366920938463463374607431768211454] = rank=10, geometric rank=20\n",
            "rvKummerK3Picard[C, 340282366920938463463374607431768211455] = rank=20, geometric rank=20\n",
        )
    );
}

#[test]
fn point_counts_use_the_unsigned_range_and_diagnose_overflow() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[19,1]{0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            rvhyperelliptic C[F]{z,o,z,z,z,o};
            rvKummerK3Points[C,15];
            "
        )),
        "rvKummerK3Points[C, 15] = singular=230466617897195215069474308673820901804, nodes=4, resolved=230466617897195215130198816793320095000\n",
    );
    assert!(
        diagnostic(interpret(&format!("{F3} rvKummerK3Points[C,41];")))
            .contains("point count exceeds unsigned 128-bit range")
    );
}

#[test]
fn rejects_zero_degrees_and_unsupported_models() {
    for (source, expected) in [
        (
            "rvKummerK3Picard[C,0];",
            "extension degree must be positive",
        ),
        (
            "rvKummerK3Points[C,0];",
            "extension degree must be positive",
        ),
        (
            "rvhyperelliptic E[F]{o,o,z,o}; rvKummerK3Zeta[E];",
            "smooth monic quintic",
        ),
        (
            "rvKummerK3Zeta[Missing];",
            "unknown residue hyperelliptic curve",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{F3}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
