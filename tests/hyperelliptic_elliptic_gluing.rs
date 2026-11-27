mod common;

use common::{diagnostic, example, interpret, stdout};

const F5: &str = r"
    rvunram K[5,1]{0,1};
    let z=rvUnramElem[K]{};
    let o=rvUnramElem[K]{1};
    let two=rvUnramElem[K]{2};
    let three=rvUnramElem[K]{3};
    let four=rvUnramElem[K]{4};
    rvhyperelliptic E[K]{o,o,z,o};
    rvhyperelliptic F[K]{o,z,o,o};
    rvellglue C[E,F](z,o,o,z);
";

#[test]
fn nonsplit_cubics_and_independent_frobenius_polynomials() {
    assert_eq!(
        stdout(example("hyperelliptic-elliptic-gluing")),
        concat!(
            "RVJac(E){u=[[1], [1]], v=[[2]]}\n",
            "RVJac(F){u=[[1], [1]], v=[[4]]}\n",
            "rvHWWWeil[C] = [1, 2, 7, 10, 25]\n",
            "rvHWWWeil[E] = [1, 3, 5]\n",
            "rvHWWWeil[F] = [1, -1, 5]\n",
        ),
    );
}

#[test]
fn both_shifts_and_a_nonsquare_scale_survive_projective_rescaling() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F5}
            rvhyperelliptic G[K]{{four,four,two,o}};
            rvhyperelliptic H[K]{{three,four,o,o}};
            // The parameters are A=1, B=2, m=2, k=3.
            rvellglue D[G,H](two,z,o,four);
            rvellglue Scaled[G,H](o,z,three,two);
            let left=rvEllGlueLeft[D](o,o);
            let right=rvEllGlueRight[D](o,o);
            assert left==rvJacPoint[G](four,o);
            assert right==rvJacPoint[H](o,two);
            assert rvEllGlueLeft[Scaled](o,o)==left;
            assert rvEllGlueRight[Scaled](o,o)==right;
            assert rvEllGlueLeft[D](four,o)==left;
            assert rvEllGlueRight[D](four,o)==neg(right);
            assert rvEllGlueLeft[D](o,four)==neg(left);
            assert rvEllGlueRight[D](o,four)==neg(right);
            assert rvEllGlueRight[D](four,four)==right;
            assert rvEllGlueLeft[D](z,o)==rvJacPoint[G](o,o);
            assert rvEllGlueRight[D](z,o)==zero[RVJac(H)];
            rvHWWWeil[D];
            rvHWWWeil[Scaled];
            rvHWWWeil[G];
            rvHWWWeil[H];
        "
        ))),
        concat!(
            "rvHWWWeil[D] = [1, 4, 13, 20, 25]\n",
            "rvHWWWeil[Scaled] = [1, 4, 13, 20, 25]\n",
            "rvHWWWeil[G] = [1, 3, 5]\n",
            "rvHWWWeil[H] = [1, 1, 5]\n",
        ),
    );
}

#[test]
fn supersingular_factor_and_branch_chart() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram K[7,1]{0,1};
            let z=rvUnramElem[K]{};
            let o=rvUnramElem[K]{1};
            let three=rvUnramElem[K]{3};
            let four=rvUnramElem[K]{4};
            let m=neg(o);
            rvhyperelliptic E[K]{z,o,z,o};
            rvhyperelliptic F[K]{three,three,m,o};
            rvellglue C[E,F](z,three,o,m);
            assert rvEllGlueLeft[C](three,z)==rvJacPoint[E](z,z);
            assert rvEllGlueRight[C](three,z)==rvJacPoint[F](four,z);
            assert rvEllGlueRight[C](z,three)==zero[RVJac(F)];
            rvjacchart N[C](three);
            rvHWWWeil[C];
            rvHWWWeil[N];
            rvHWWWeil[E];
            rvHWWWeil[F];
        ",
        )),
        concat!(
            "rvHWWWeil[C] = [1, -4, 14, -28, 49]\n",
            "rvHWWWeil[N] = [1, -4, 14, -28, 49]\n",
            "rvHWWWeil[E] = [1, 0, 7]\n",
            "rvHWWWeil[F] = [1, -4, 7]\n",
        ),
    );
}

#[test]
fn extension_field_parameters_and_images_have_nonprime_coordinates() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram K[5,1]{3,0,1};
            let o=rvUnramElem[K]{1};
            let two=rvUnramElem[K]{2};
            let alpha=rvUnramElem[K]{0,1};
            let beta=rvUnramElem[K]{1,1};
            let minusAlpha=neg(alpha);
            let twiceAlpha=add(alpha,alpha);
            rvhyperelliptic E[K]{rvUnramElem[K]{1,2},two,twiceAlpha,o};
            rvhyperelliptic F[K]{rvUnramElem[K]{2,4},twiceAlpha,rvUnramElem[K]{4,2},o};
            rvellglue C[E,F](beta,minusAlpha,o,minusAlpha);
            assert rvEllGlueLeft[C](o,o)==rvJacPoint[E](rvUnramElem[K]{3,1},o);
            assert rvEllGlueRight[C](o,o)==rvJacPoint[F](alpha,two);
        ",
        )),
        "",
    );
}

#[test]
fn rejects_uncertified_matches_noncurves_and_nominal_mismatches() {
    for (source, expected) in [
        ("rvellglue D[E,F](o,o,o,o);", "determinant must be nonzero"),
        ("rvellglue D[E,E](o,z,z,o);", "non-affine branch matching"),
        ("rvellglue D[E,F](o,o,o,z);", "does not match the finite branch roots"),
        ("rvhyperelliptic R[K]{z,o,z,o}; rvellglue D[R,F](z,o,o,z);", "f(A) and g(B) to be nonzero"),
        ("rvhyperelliptic R[K]{z,o,z,o}; rvellglue D[E,R](z,o,o,z);", "f(A) and g(B) to be nonzero"),
        ("rvellglue E[E,F](z,o,o,z);", "curve is already defined"),
        ("rvellglue D[C,F](z,o,o,z);", "must be a smooth monic cubic"),
        ("rvunram L[5,1]{0,1}; let w=rvUnramElem[L]{1}; let v=rvUnramElem[L]{}; rvhyperelliptic G[L]{w,v,w,w}; rvellglue D[E,G](z,o,o,z);", "same named field"),
        ("rvunram L[5,1]{0,1}; rvellglue D[E,F](z,rvUnramElem[L]{1},o,z);", "rvellglue matching coefficient"),
        ("show rvEllGlueRight[C](z,z);", "point does not lie on the glued genus-two curve"),
        ("show rvEllGlueLeft[E](z,o);", "has no stored elliptic gluing"),
        ("show add(rvEllGlueLeft[C](two,two),rvEllGlueRight[C](two,two));", "expected type RVJac(E)"),
    ] {
        let problem = diagnostic(interpret(&format!("{F5}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
