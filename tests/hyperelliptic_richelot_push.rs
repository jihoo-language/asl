mod common;

use common::{diagnostic, example, interpret, stdout};

const F7: &str = r"
    rvunram F[7,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    let a=rvUnramElem[F]{3};
    let h=rvUnramElem[F]{4};
    let e=rvUnramElem[F]{5};
    let s=rvUnramElem[F]{6};
";

#[test]
fn nonsplit_and_tangent_support_have_exact_images() {
    assert_eq!(
        stdout(example("hyperelliptic-richelot-push")),
        concat!(
            "RVJac(N){u=[[2], [3], [1]], v=[[0], [3]]}\n",
            "RVJac(N){u=[[4], [4], [1]], v=[[6]]}\n",
        )
    );
}

#[test]
fn push_is_additive_on_every_pair_in_a_charted_jacobian() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F7}
            rvhyperelliptic C[F]{{z,s,s,z,s,o}};
            rvrichelot D[C]{{[],[1]}}{{[1],[],[1]}}{{[6],[6],[1]}};
            rvjacchart N[D](o);
            rvJacCount[C];
            rvJacCount[N];
            let push=\p:RVJac(C)=>rvRichelotPush[N](p);
            assert (\p:RVJac(C)=>\q:RVJac(C)=>push(add(p,q)))
                ==(\p:RVJac(C)=>\q:RVJac(C)=>add(push(p),push(q)));
            "
        ))),
        "rvJacCount[C] = 20\nrvJacCount[N] = 20\n"
    );
}

#[test]
fn dual_composition_and_chart_transport_on_all_sixty_four_classes() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram U[3,1]{1,0,1};
            let z=rvUnramElem[U]{};
            let o=rvUnramElem[U]{1};
            let b=rvUnramElem[U]{1,1};
            let c=rvUnramElem[U]{2,1};
            rvhyperelliptic C[U]{z,b,z,c,z,o};
            rvrichelot D[C]{[1],[],[1]}{[1,1],[],[1]}{[],[1]};
            rvrichelotdual E[D];
            rvJacCount[C];
            let forward=\p:RVJac(C)=>rvRichelotPush[D](p);
            let backward=\p:RVJac(D)=>rvRichelotPush[E](p);
            // The chart to E multiplies y by 2, which is negation in F9.
            assert (\p:RVJac(C)=>backward(forward(p)))
                ==(\p:RVJac(C)=>rvJacExtend[E](rvJacScale[-2](p)));

            rvjacchart N[D](o);
            assert (\p:RVJac(C)=>rvRichelotPush[N](p))
                ==(\p:RVJac(C)=>rvJacRechart[N](forward(p)));
            "
        )),
        "rvJacCount[C] = 64\n"
    );
}

#[test]
fn split_support_handles_a_branch_point_and_a_degree_drop_in_separate_cases() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F7}
            rvhyperelliptic C[F]{{z,h,t,e,t,o}};
            rvrichelot D[C]{{[4],[2],[1]}}{{[],[1]}}{{[1],[],[1]}};
            rvjacchart N[D](t);
            let mixed=add(rvJacPoint[C](t,o),rvJacPoint[C](o,z));
            assert rvRichelotPush[N](mixed)==add(rvRichelotPoint[N](t,o),rvRichelotPoint[N](o,z));

            rvhyperelliptic B[F]{{z,a,t,h,t,o}};
            rvrichelot G[B]{{[1],[],[1]}}{{[3],[2],[1]}}{{[],[1]}};
            rvjacchart M[G](a);
            let two=add(rvJacPoint[B](a,o),rvJacPoint[B](h,a));
            let image=rvRichelotPush[M](two);
            assert image==add(rvRichelotPoint[M](a,o),rvRichelotPoint[M](h,a));
            assert image==rvjac[M]{{[1],[5],[1]}}{{[2],[2]}};
            "
        ))),
        ""
    );
}

#[test]
fn push_has_the_three_split_branch_pair_classes_in_its_kernel() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F7}
            rvhyperelliptic C[F]{{z,a,s,z,h,o}};
            rvrichelot D[C]{{[],[1]}}{{[2],[4],[1]}}{{[5],[],[1]}};
            rvjacchart N[D](a);
            let first=rvJacPoint[C](z,z);
            let second=add(rvJacPoint[C](o,z),rvJacPoint[C](t,z));
            let third=add(rvJacPoint[C](a,z),rvJacPoint[C](h,z));
            assert add(first,second)==third;
            assert rvRichelotPush[N](zero[RVJac(C)])==zero[RVJac(N)];
            assert rvRichelotPush[N](first)==zero[RVJac(N)];
            assert rvRichelotPush[N](second)==zero[RVJac(N)];
            assert rvRichelotPush[N](third)==zero[RVJac(N)];
            let mixed=add(rvJacPoint[C](o,z),rvJacPoint[C](a,z));
            let image=rvRichelotPush[N](mixed);
            assert image==add(rvRichelotPoint[N](o,z),rvRichelotPoint[N](a,z));
            assert image==rvjac[N]{{[6],[4],[1]}}{{}};
            "
        ))),
        ""
    );
}

#[test]
fn push_checks_nominal_source_and_target_models() {
    for (input, expected) in [
        (
            "rvhyperelliptic Twin[F]{z,h,t,e,t,o}; show rvRichelotPush[N](zero[RVJac(Twin)]);",
            "expected type RVJac(C)",
        ),
        (
            "show rvRichelotPush[D](zero[RVJac(C)]);",
            "monic odd-degree",
        ),
        (
            "show rvRichelotPush[C](zero[RVJac(C)]);",
            "has no stored Richelot splitting",
        ),
        (
            "show add(rvRichelotPush[N](zero[RVJac(C)]),zero[RVJac(C)]);",
            "expected type RVJac(N)",
        ),
    ] {
        let problem = diagnostic(interpret(&format!(
            r"{F7}
            rvhyperelliptic C[F]{{z,h,t,e,t,o}};
            rvrichelot D[C]{{[4],[2],[1]}}{{[],[1]}}{{[1],[],[1]}};
            rvjacchart N[D](t);
            {input}"
        )));
        assert!(problem.contains(expected), "{problem}");
    }
}
