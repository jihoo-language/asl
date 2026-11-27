mod common;

use common::{diagnostic, example, interpret, stdout};

const F5: &str = r"
    rvunram F[5,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    let three=rvUnramElem[F]{3};
    let four=rvUnramElem[F]{4};
    rvhyperelliptic C[F]{o,o,z,three,o,o};
";

#[test]
fn even_and_nonmonic_models_gain_an_odd_jacobian_chart() {
    assert_eq!(
        stdout(example("hyperelliptic-chart")),
        concat!(
            "rvHWWWeil[E] = [1, 3, 8, 15, 25]\n",
            "rvHWWWeil[M] = [1, 3, 8, 15, 25]\n",
            "rvJacCount[M] = 52\n",
        )
    );
    assert_eq!(
        stdout(interpret(&format!(
            r"{F5}
            rvhyperelliptic E[F]{{four,o,z,three,z,o,t}};
            rvjacchart M[E](t);
            rvhyperelliptic ExpectedEven[F]{{t,z,z,three,three,o}};
            let even=rvJacPoint[M](o,t);
            assert rvJacExtend[ExpectedEven](even)==rvJacPoint[ExpectedEven](o,t);
            assert add(even,neg(even))==zero[RVJac(M)];

            rvhyperelliptic O[F]{{z,o,z,three,z,t}};
            rvjacchart N[O](t);
            rvhyperelliptic ExpectedOdd[F]{{z,o,z,o,three,o}};
            assert rvJacExtend[ExpectedOdd](zero[RVJac(N)])==zero[RVJac(ExpectedOdd)];
        "
        ))),
        ""
    );
}

#[test]
fn transport_handles_poles_repeated_points_and_nonsplit_support() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F5}
            rvjacchart N[C](t);
            rvhyperelliptic Expected[F]{{z,four,o,four,t,o}};
            assert rvJacExtend[Expected](zero[RVJac(N)])==zero[RVJac(Expected)];
            let branch=rvJacPoint[C](t,z);
            let infinity=rvJacPoint[N](z,z);
            let p=rvJacPoint[C](z,o);
            let q=rvJacPoint[C](three,t);
            let mapped=rvJacPoint[N](three,three);
            assert rvJacRechart[N](branch)==infinity;
            assert rvJacUnchart[N](infinity)==branch;
            assert rvJacRechart[N](p)==add(mapped,infinity);
            assert rvJacRechart[N](add(branch,p))==mapped;
            assert rvJacUnchart[N](mapped)==add(branch,p);
            assert rvJacRechart[N](add(p,q))==add(rvJacRechart[N](p),rvJacRechart[N](q));
            let repeated=rvjac[C]{{[],[],[1]}}{{[1],[3]}};
            let repeatedTarget=rvjac[N]{{[4],[4],[1]}}{{[1],[4]}};
            assert repeated==add(p,p);
            assert rvJacRechart[N](repeated)==repeatedTarget;
            assert rvJacUnchart[N](repeatedTarget)==repeated;
            let nonsplit=rvjac[C]{{[2],[1],[1]}}{{[1],[2]}};
            let nonsplitTarget=rvjac[N]{{[2],[],[1]}}{{[1]}};
            assert rvJacRechart[N](nonsplit)==nonsplitTarget;
            assert rvJacUnchart[N](nonsplitTarget)==nonsplit;
        "
        ))),
        ""
    );
}

#[test]
fn the_two_transports_are_inverse_on_the_full_small_group() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[3,1]{0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            let t=rvUnramElem[F]{2};
            rvhyperelliptic C[F]{z,o,o,o,t,o};
            rvjacchart N[C](t);
            let forward=\d:RVJac(C)=>rvJacRechart[N](d);
            let backward=\d:RVJac(N)=>rvJacUnchart[N](d);
            assert (\d:RVJac(C)=>backward(forward(d)))==(\d:RVJac(C)=>d);
            assert (\d:RVJac(N)=>forward(backward(d)))==(\d:RVJac(N)=>d);
            assert (\a:RVJac(C)=>\b:RVJac(C)=>forward(add(a,b)))
                ==(\a:RVJac(C)=>\b:RVJac(C)=>add(forward(a),forward(b)));
            rvjacchart M[N](z);
            let multicharted=\d:RVJac(C)=>rvJacRechart[M](forward(d));
            assert (\d:RVJac(C)=>backward(rvJacUnchart[M](multicharted(d))))==(\d:RVJac(C)=>d);
        "
        )),
        ""
    );
}

#[test]
fn branch_and_derivative_can_lie_in_an_extension_field() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[3,1]{1,0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            let a=rvUnramElem[F]{0,1};
            let d=rvUnramElem[F]{0,2};
            let b=rvUnramElem[F]{1,1};
            rvhyperelliptic C[F]{b,z,o,o};
            rvjacchart N[C](a);
            rvhyperelliptic Expected[F]{z,d,o,o};
            assert rvJacExtend[Expected](zero[RVJac(N)])==zero[RVJac(Expected)];
            let p=rvJacPoint[C](o,rvUnramElem[F]{1,2});
            let im=rvUnramElem[F]{2,1};
            let q=add(rvJacPoint[N](im,im),rvJacPoint[N](z,z));
            assert rvJacRechart[N](p)==q;
            assert rvJacUnchart[N](q)==p;
        "
        )),
        ""
    );
}

#[test]
fn charts_require_smooth_curves_branch_values_and_nominal_endpoints() {
    for (source, expected) in [
        ("rvjacchart N[C](o);", "finite branch point"),
        ("rvhyperelliptic Bad[F]{z,z,o}; rvjacchart N[Bad](z);", "squarefree"),
        ("rvunram W[5,2]{0,1}; let w=rvUnramElem[W]{1}; rvhyperelliptic D[W]{w,w}; rvjacchart N[D](neg(w));", "precision one"),
        ("rvunram G[5,1]{0,1}; rvjacchart N[C](rvUnramElem[G]{2});", "expected type RVUnram(F)"),
        ("rvjacchart N[C](t); rvhyperelliptic Twin[F]{o,o,z,three,o,o}; show rvJacRechart[N](zero[RVJac(Twin)]);", "expected type RVJac(C)"),
        ("rvjacchart N[C](t); rvjacchart Twin[C](t); show rvJacUnchart[N](zero[RVJac(Twin)]);", "expected type RVJac(N)"),
        ("show rvJacRechart[C](zero[RVJac(C)]);", "unknown Jacobian chart"),
        ("rvjacchart N[C](t); rvjacchart N[C](t);", "already defined"),
        ("rvhyperelliptic E[F]{four,o,z,three,z,o,t}; rvjacchart M[E](t); show rvJacUnchart[M](zero[RVJac(M)]);", "monic odd-degree"),
        ("rvhyperelliptic E[F]{four,o,z,three,z,o,t}; rvjacchart M[E](t); show rvJacRechart[M](zero[RVJac(C)]);", "monic odd-degree"),
    ] {
        let problem = diagnostic(interpret(&format!("{F5}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
