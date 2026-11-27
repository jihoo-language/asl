mod common;

use common::{diagnostic, example, interpret, stdout};

const F7: &str = r"
    rvunram F[7,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    let h=rvUnramElem[F]{4};
    let e=rvUnramElem[F]{5};
    rvhyperelliptic C[F]{z,h,t,e,t,o};
    rvrichelot D[C]{[4],[2],[1]}{[],[1]}{[1],[],[1]};
    rvjacchart N[D](t);
";

#[test]
fn split_fiber_dual_composition_and_weil_polynomial() {
    assert_eq!(
        stdout(example("hyperelliptic-richelot-point")),
        concat!(
            "RVJac(D){u=[[0, 2], [0, 0], [1, 0]], v=[[0, 2], [2, 2]]}\n",
            "rvHWWWeil[C] = [1, -2, 2, -18, 81]\n",
            "rvHWWWeil[D] = [1, -2, 2, -18, 81]\n",
            "RVJac(N){u=[[0], [0], [1]], v=[[1], [5]]}\n",
        )
    );
}

#[test]
fn points_on_kernel_and_branch_fibers_under_a_change_of_chart() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F7}
            let kernel=rvRichelotPoint[N](z,z);
            assert kernel==zero[RVJac(N)];
            let a=rvRichelotPoint[N](o,z);
            let b=rvRichelotPoint[N](h,z);
            assert a==rvjac[N]{{[3],[1]}}{{}};
            assert a==b;
            assert add(a,b)==zero[RVJac(N)];
            assert add(rvJacPoint[C](o,z),rvJacPoint[C](h,z))==rvjac[C]{{[4],[2],[1]}}{{}};

            rvrichelot R[C]{{[],[1]}}{{[1],[],[1]}}{{[4],[2],[1]}};
            rvjacchart M[R](t);
            let oldInfinity=rvRichelotPoint[N](t,o);
            assert rvJacExtend[M](oldInfinity)==rvRichelotPoint[M](t,o);
            assert rvRichelotPoint[N](t,neg(o))==neg(oldInfinity);
            assert rvRichelotPoint[N](e,rvUnramElem[F]{{3}})==rvjac[N]{{[2],[3],[1]}}{{[],[4]}};
            "
        ))),
        ""
    );
}

#[test]
fn nonbranch_fiber_can_meet_new_infinity() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[7,1]{0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            let t=rvUnramElem[F]{2};
            let a=rvUnramElem[F]{3};
            let h=rvUnramElem[F]{4};
            rvhyperelliptic C[F]{z,a,t,h,t,o};
            rvrichelot D[C]{[1],[],[1]}{[3],[2],[1]}{[],[1]};
            rvjacchart N[D](a);
            let image=rvRichelotPoint[N](a,o);
            assert image==rvJacPoint[N](o,a);
            assert rvRichelotPoint[N](a,neg(o))==neg(image);
            "
        )),
        ""
    );
}

#[test]
fn a_chart_on_a_monic_quintic_target_agrees_with_jacobian_transport() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram U[3,1]{1,0,1};
            let z=rvUnramElem[U]{};
            let o=rvUnramElem[U]{1};
            let a=rvUnramElem[U]{0,1};
            let b=rvUnramElem[U]{1,1};
            let c=rvUnramElem[U]{2,1};
            rvhyperelliptic C[U]{z,b,z,c,z,o};
            rvrichelot D[C]{[1],[],[1]}{[1,1],[],[1]}{[],[1]};
            rvjacchart N[D](o);
            assert rvRichelotPoint[N](b,a)==rvJacRechart[N](rvRichelotPoint[D](b,a));
            assert rvRichelotPoint[N](a,z)==rvJacRechart[N](rvRichelotPoint[D](a,z));
            assert rvRichelotPoint[N](neg(a),z)==rvRichelotPoint[N](a,z);
            "
        )),
        ""
    );
}

#[test]
fn enforces_membership_monic_source_and_nominal_endpoints() {
    for (input, expected) in [
        (
            "show rvRichelotPoint[N](z,o);",
            "does not lie on the Richelot source curve",
        ),
        ("show rvRichelotPoint[D](t,o);", "monic odd-degree"),
        (
            "show rvRichelotPoint[C](t,o);",
            "has no stored Richelot splitting",
        ),
        (
            "rvrichelotdual E[D]; rvjacchart M[E](z); show rvRichelotPoint[M](z,z);",
            "requires a monic quintic source",
        ),
        (
            "rvunram Q[7,1]{0,1}; show rvRichelotPoint[N](t,rvUnramElem[Q]{1});",
            "expected type RVUnram(F)",
        ),
        (
            "show add(rvRichelotPoint[N](t,o),zero[RVJac(C)]);",
            "expected type RVJac(N)",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{F7}{input}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
