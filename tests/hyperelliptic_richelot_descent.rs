mod common;

use common::{diagnostic, example, interpret, stdout};

const TRANSPOSED: &str = r"
    rvunram F[3,1]{0,1};
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    rvhyperelliptic C[F]{z,t,t,z,t,o};
    rvjacbase CL[C,W];
    rvrichelot DL[CL]{[1,1],[1],[1]}{[1,2],[1],[1]}{[],[1]};
";

const CUBIC: &str = r"
    rvunram F[5,1]{0,1};
    rvunram W[5,1]{1,1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    let four=rvUnramElem[F]{4};
    rvhyperelliptic C[F]{o,four,o,t,four,z,o};
    let zw=rvUnramElem[W]{};
    let ow=rvUnramElem[W]{1};
    let tw=rvUnramElem[W]{2};
    let fourw=rvUnramElem[W]{4};
    rvhyperelliptic CL[W]{ow,fourw,ow,tw,fourw,zw,ow};
    // The three quadratics are x²+a x+a² and its successive Frobenius conjugates.
    rvrichelot DL[CL]{[0,0,1],[0,1],[1]}{[3,3,3],[1,1,4],[1]}{[0,2,1],[4,3,1],[1]};
    rvrichelotdescend D[C,DL];
";

#[test]
fn a_transposition_has_rational_kernel_order_two_and_image_order_ten() {
    assert_eq!(
        stdout(example("hyperelliptic-richelot-descent")),
        concat!(
            "rvRichelotDescent[D] = relative degree=2, permutation=[1, 0, 2], rational kernel order=2\n",
            "rvHWWWeil[C] = [1, 2, 2, 6, 9]\n",
            "rvHWWWeil[D] = [1, 2, 2, 6, 9]\n",
            "rvJacCount[C] = 20\n",
            "rvJacCount[D] = 20\n",
            "RVJac(D){u=[[2], [0], [1]], v=[[2], [1]]}\n",
            "rvJacOrder = 10\n",
        )
    );
}

#[test]
fn descended_push_is_additive_and_factor_scaling_preserves_the_equation() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{TRANSPOSED}
            rvrichelotdescend D[C,DL];
            // The factors are x, (1+a)^(-1) G2, and (1+a) G1.
            rvrichelot RL[CL]{{[],[1]}}{{[0,2],[2,1],[2,1]}}{{[0,2],[1,1],[1,1]}};
            rvrichelotdescend R[C,RL];
            rvRichelotDescent[R];
            rvhyperelliptic Known[F]{{o,t,t,o,t,o}};
            assert rvJacExtend[Known](zero[RVJac(D)])==zero[RVJac(Known)];
            let push=\p:RVJac(C)=>rvRichelotDescPush[D](p);
            assert push==(\p:RVJac(C)=>rvJacExtend[D](rvRichelotDescPush[R](p)));
            assert (\p:RVJac(C)=>\q:RVJac(C)=>push(add(p,q)))
                == (\p:RVJac(C)=>\q:RVJac(C)=>add(push(p),push(q)));
            "
        ))),
        "rvRichelotDescent[R] = relative degree=2, permutation=[0, 2, 1], rational kernel order=2\n"
    );
}

#[test]
fn a_cubic_frobenius_orbit_has_trivial_rational_kernel() {
    assert_eq!(
        stdout(interpret(&format!(
            "{CUBIC} rvRichelotDescent[D]; rvHWWWeil[C]; rvHWWWeil[D];"
        ))),
        concat!(
            "rvRichelotDescent[D] = relative degree=3, permutation=[1, 2, 0], rational kernel order=1\n",
            "rvHWWWeil[C] = [1, 2, 7, 10, 25]\n",
            "rvHWWWeil[D] = [1, 2, 7, 10, 25]\n",
        )
    );
    assert!(diagnostic(interpret(&format!(
        "{CUBIC} show rvRichelotDescPush[D](o);"
    )))
    .contains("requires a monic quintic source"));
}

#[test]
fn rejects_a_splitting_with_genuinely_noninvariant_branch_pairs() {
    let prefix = r"
        rvunram F[5,1]{0,1};
        rvunram W[5,1]{-2,0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let t=rvUnramElem[F]{2};
        rvhyperelliptic C[F]{z,o,o,z,t,o};
        rvjacbase CL[C,W];
        // G1 and G2 pair 1 and 2 with opposite square roots of 2.
        rvrichelot DL[CL]{[0,1],[4,4],[1]}{[0,3],[3,1],[1]}{[],[1]};
    ";
    assert!(stdout(interpret(&format!("{prefix} rvRichelotData[DL];")))
        .starts_with("rvRichelotData[DL] = source=CL"));
    let problem = diagnostic(interpret(&format!("{prefix} rvrichelotdescend D[C,DL];")));
    assert!(
        problem.contains("Frobenius does not preserve the Richelot splitting"),
        "{problem}"
    );
}

#[test]
fn rejects_wrong_base_equations_and_nominal_endpoints() {
    for (code, expected) in [
        ("rvhyperelliptic B[F]{z,o,z,o,t,o}; rvrichelotdescend D[B,DL];", "not the declared base extension"),
        ("rvrichelotdescend D[C,CL];", "has no stored Richelot splitting"),
        ("rvrichelotdescend D[C,DL]; rvrichelotdescend D[C,DL];", "already defined"),
        ("rvRichelotDescent[DL];", "has no stored Richelot descent"),
        ("rvrichelotdescend D[C,DL]; rvhyperelliptic B[F]{z,t,t,z,t,o}; show rvRichelotDescPush[D](zero[RVJac(B)]);", "expected type RVJac(C)"),
    ] {
        let problem = diagnostic(interpret(&format!("{TRANSPOSED}{code}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
