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
    rvhyperelliptic T[F]{three,o,z,t,three,o};
    let back=\k:RVKummer(T)=>rvKummer[C](rvKummerCoord[0](k),rvUnramMul(t,rvKummerCoord[1](k)),rvUnramMul(four,rvKummerCoord[2](k)),rvUnramMul(three,rvKummerCoord[3](k)));
    let p=rvJacPoint[T](three,four);
    let q=rvJacPoint[T](t,four);
    let a=back(rvJacKummerPoint(p));
    let b=back(rvJacKummerPoint(q));
    let id=zero[RVKummer(C)];
";

#[test]
fn explicit_twist_gives_both_differential_addition_candidates() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-arithmetic")),
        "RVKummer(C)[[1], [0], [4], [0]]\nRVKummer(C)[[1], [0], [4], [2]]\n"
    );
}

#[test]
fn interior_twist_points_and_signed_scalar_recurrences() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F5}
            let sum=back(rvJacKummerPoint(add(p,q)));
            let diff=back(rvJacKummerPoint(add(p,neg(q))));
            assert rvKummerLiftable(sum)==S.lo;
            assert rvKummerLiftable(diff)==S.lo;
            assert rvKummerScale[0](sum)==id;
            assert rvKummerScale[1](diff)==diff;
            assert rvKummerScale[-1](diff)==diff;
            assert rvKummerScale[7](sum)==back(rvJacKummerPoint(rvJacScale[7](add(p,q))));
            assert rvKummerScale[-7](diff)==back(rvJacKummerPoint(rvJacScale[7](add(p,neg(q)))));
            assert rvKummerScale[6](diff)==rvKummerScale[3](rvKummerScale[2](diff));
            let twice=rvKummerScale[2](sum);
            let thrice=rvKummerScale[3](sum);
            assert rvKummerDiffAdd(thrice,twice,sum)==rvKummerScale[5](sum);
            assert rvKummerDiffAdd(thrice,twice,rvKummerScale[5](sum))==sum;
        "
        ))),
        ""
    );
}

#[test]
fn addition_handles_nodes_identity_and_repeated_inputs() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F5}
            let node=rvJacKummerPoint(rvJacPoint[C](t,z));
            let tn=rvJacPoint[T](o,z);
            let shifted=back(rvJacKummerPoint(add(p,tn)));
            assert rvKummerScale[3](node)==node;
            assert rvKummerScale[2](node)==id;
            assert rvKummerScale[-1](id)==id;
            assert rvKummerDiffAdd(a,id,a)==a;
            assert rvKummerDiffAdd(id,a,a)==a;
            assert rvKummerDiffAdd(node,node,id)==id;
            assert rvKummerDiffAdd(a,node,shifted)==shifted;
            assert rvKummerDiffAdd(node,a,shifted)==shifted;
            let twice=rvKummerScale[2](a);
            assert rvKummerCoord[0](twice)==o;
            assert rvKummerDiffAdd(a,a,id)==twice;
            assert rvKummerDiffAdd(a,a,twice)==id;
        "
        ))),
        ""
    );
}

#[test]
fn rational_lifts_match_cantor_scaling_on_the_full_small_group() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[3,1]{0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            let t=rvUnramElem[F]{2};
            rvhyperelliptic C[F]{z,o,o,o,z,o};
            let viaKummer=\d:RVJac(C)=>rvKummerScale[-7](rvJacKummerPoint(d));
            let viaJacobian=\d:RVJac(C)=>rvJacKummerPoint(rvJacScale[7](d));
            assert viaKummer==viaJacobian;
            let signIndependent=\k:RVKummer(C)=>rvKummerScale[6](k);
            assert signIndependent==(\k:RVKummer(C)=>rvKummerScale[-6](k));
            let p=rvJacPoint[C](t,o);
            assert rvKummerScale[-170141183460469231731687303715884105728](rvJacKummerPoint(p))
              ==rvJacKummerPoint(rvJacScale[-170141183460469231731687303715884105728](p));
        "
        )),
        ""
    );
}

#[test]
fn differential_addition_reconciles_distinct_radicands_over_a_nonprime_field() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[3,1]{1,0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            let t=rvUnramElem[F]{2};
            let alpha=rvUnramElem[F]{0,1};
            let minusAlpha=rvUnramElem[F]{0,2};
            let k=rvUnramElem[F]{2,1};
            let k3=rvUnramElem[F]{2,2};
            let b=rvUnramElem[F]{1,1};
            rvhyperelliptic C[F]{t,z,z,z,z,o};
            rvhyperelliptic T[F]{b,z,z,z,z,o};
            let back=\v:RVKummer(T)=>rvKummer[C](rvKummerCoord[0](v),rvUnramMul(k,rvKummerCoord[1](v)),rvUnramMul(alpha,rvKummerCoord[2](v)),rvUnramMul(k3,rvKummerCoord[3](v)));
            let d=rvJacPoint[T](k,minusAlpha);
            let e=rvJacPoint[T](rvUnramElem[F]{1,2},b);
            let p=rvKummer[C](z,o,alpha,t);
            let q=rvKummer[C](z,o,minusAlpha,t);
            assert rvKummerLiftable(p)==S.lo;
            assert rvKummerLiftable(q)==S.lo;
            assert p==back(rvJacKummerPoint(d));
            assert q==back(rvJacKummerPoint(e));
            let sum=back(rvJacKummerPoint(add(d,e)));
            let diff=back(rvJacKummerPoint(add(d,neg(e))));
            assert rvKummerDiffAdd(p,q,diff)==sum;
            assert rvKummerDiffAdd(q,p,sum)==diff;
        "
        )),
        ""
    );
}

#[test]
fn rejects_nominal_mismatches_mixed_twists_and_incompatible_differences() {
    for (source, expected) in [
        (
            "show rvKummerScale[2](zero[RVJac(C)]);",
            "an RVKummer point",
        ),
        (
            "show rvKummerDiffAdd(zero[RVJac(C)],a,b);",
            "an RVKummer point",
        ),
        (
            "show rvKummerDiffAdd(a,zero[RVKummer(T)],id);",
            "expected type RVKummer(C)",
        ),
        (
            "show rvKummerDiffAdd(a,a,zero[RVKummer(T)]);",
            "expected type RVKummer(C)",
        ),
        ("show rvKummerDiffAdd(a,a,a);", "difference is incompatible"),
        (
            "let rational=rvJacKummerPoint(rvJacPoint[C](z,o)); show rvKummerDiffAdd(a,rational,id);",
            "different quadratic twist classes",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{F5}{source}")));
        assert!(error.contains(expected), "{error}");
    }
}
