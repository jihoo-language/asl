mod common;

use common::{diagnostic, example, interpret, stdout};

const SPLIT: &str = include_str!("../examples/hyperelliptic-elliptic-gluing-picard.asl");

const NONSPLIT: &str = r"
    rvunram K[7,1]{0,1};
    let z=rvUnramElem[K]{};
    let o=rvUnramElem[K]{1};
    let three=rvUnramElem[K]{3};
    let four=rvUnramElem[K]{4};
    rvhyperelliptic E[K]{z,o,z,o};
    rvhyperelliptic F[K]{three,three,neg(o),o};
    rvellglue C[E,F](z,three,o,neg(o));
    rvjacchart S[C](three);
";

#[test]
fn both_compositions_double_on_the_full_split_groups() {
    assert_eq!(
        stdout(example("hyperelliptic-elliptic-gluing-picard")),
        "rvJacCount[E] = 8\nrvJacCount[F] = 4\nrvJacCount[S] = 32\n"
    );
}

#[test]
fn the_split_kernel_is_exactly_the_four_matching_pairs_in_either_chart() {
    let text = stdout(interpret(&format!(
        r"{SPLIT}
        show glue;
        // Another base point gives the same graph kernel.
        rvjacchart T[C](two);
        let other=\pq:RVJac(E)*RVJac(F)=>rvEllGlueMap[T](fst(pq),snd(pq));
        let id=zero[RVJac(T)];
        assert other((zero[RVJac(E)],zero[RVJac(F)]))==id;
        assert other((rvJacPoint[E](z,z),rvJacPoint[F](two,z)))==id;
        assert other((rvJacPoint[E](o,z),rvJacPoint[F](z,z)))==id;
        assert other((rvJacPoint[E](two,z),rvJacPoint[F](four,z)))==id;
        show other;
        "
    )));
    assert_eq!(text.matches("=> RVJac(S){u=[[1]], v=[]}").count(), 4);
    assert_eq!(text.matches("=> RVJac(T){u=[[1]], v=[]}").count(), 4);
}

#[test]
fn a_nonsplit_matching_has_two_rational_kernel_pairs_and_four_over_f49() {
    let text = stdout(interpret(&format!(
        r"{NONSPLIT}
        rvJacCount[E];
        rvJacCount[F];
        rvJacCount[S];
        let glue=\pq:RVJac(E)*RVJac(F)=>rvEllGlueMap[S](fst(pq),snd(pq));
        let dual=\d:RVJac(S)=>rvEllGlueDual[S](d);
        assert (\pq:RVJac(E)*RVJac(F)=>dual(glue(pq)))
            ==(\pq:RVJac(E)*RVJac(F)=>(rvJacScale[2](fst(pq)),rvJacScale[2](snd(pq))));
        assert (\d:RVJac(S)=>glue(dual(d)))==(\d:RVJac(S)=>rvJacScale[2](d));
        assert glue((zero[RVJac(E)],zero[RVJac(F)]))==zero[RVJac(S)];
        assert glue((rvJacPoint[E](z,z),rvJacPoint[F](four,z)))==zero[RVJac(S)];
        show glue;

        rvunram W[7,1]{{1,0,1}};
        let zw=rvUnramElem[W]{{}};
        let ow=rvUnramElem[W]{{1}};
        let tw=rvUnramElem[W]{{3}};
        let fw=rvUnramElem[W]{{4}};
        let i=rvUnramElem[W]{{0,1}};
        let a=rvUnramElem[W]{{2,2}};
        let b=rvUnramElem[W]{{2,5}};
        rvjacbase EL[E,W];
        rvjacbase FL[F,W];
        rvellglue L[EL,FL](zw,tw,ow,neg(ow));
        // This branch lies over x=i and requires the extension field.
        rvjacchart T[L](rvUnramElem[W]{{2,3}});
        let identity=zero[RVJac(T)];
        assert rvEllGlueMap[T](zero[RVJac(EL)],zero[RVJac(FL)])==identity;
        assert rvEllGlueMap[T](rvJacPoint[EL](zw,zw),rvJacPoint[FL](fw,zw))==identity;
        assert rvEllGlueMap[T](rvJacPoint[EL](i,zw),rvJacPoint[FL](a,zw))==identity;
        assert rvEllGlueMap[T](rvJacPoint[EL](neg(i),zw),rvJacPoint[FL](b,zw))==identity;
        let p=rvJacPoint[EL](ow,tw);
        let image=rvEllGlueMap[T](p,zero[RVJac(FL)]);
        assert rvEllGlueDual[T](image)==(rvJacScale[2](p),zero[RVJac(FL)]);

        // The chart derivative is 4+6i. These are the images of the
        // raw points (1,3i) and (-1,3i) in the chart at 2+3i.
        let y=rvUnramElem[W]{{0,3}};
        let first=rvJacPoint[T](a,rvUnramElem[W]{{1,4}});
        let second=rvJacPoint[T](rvUnramElem[W]{{3,2}},rvUnramElem[W]{{0,2}});
        let left=rvEllGlueLeft[L](ow,y);
        let right=rvEllGlueRight[L](ow,y);
        assert rvEllGlueMap[T](left,zero[RVJac(FL)])==add(first,second);
        assert rvEllGlueMap[T](zero[RVJac(EL)],right)==add(first,neg(second));
        assert rvEllGlueDual[T](first)
            ==(add(left,rvJacPoint[EL](i,zw)),add(right,rvJacPoint[FL](a,zw)));
        "
    )));
    assert!(text.starts_with("rvJacCount[E] = 8\nrvJacCount[F] = 4\nrvJacCount[S] = 32\n"));
    assert_eq!(text.matches("=> RVJac(S){u=[[1]], v=[]}").count(), 2);
}

#[test]
fn the_same_nominal_elliptic_curve_can_fill_both_slots_and_be_reused() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{SPLIT}
            // This self-match interchanges the roots 0 and 1 and fixes 2.
            rvellglue R[E,E](three,two,o,two);
            rvjacchart T[R](two);
            rvJacCount[T];
            let selfGlue=\pq:RVJac(E)*RVJac(E)=>rvEllGlueMap[T](fst(pq),snd(pq));
            let selfDual=\d:RVJac(T)=>rvEllGlueDual[T](d);
            assert (\pq:RVJac(E)*RVJac(E)=>selfDual(selfGlue(pq)))
                ==(\pq:RVJac(E)*RVJac(E)=>(rvJacScale[2](fst(pq)),rvJacScale[2](snd(pq))));
            assert (\d:RVJac(T)=>selfGlue(selfDual(d)))==(\d:RVJac(T)=>rvJacScale[2](d));
            assert selfGlue((rvJacPoint[E](z,z),rvJacPoint[E](o,z)))==zero[RVJac(T)];
            assert rvEllGlueMap[S](rvJacPoint[E](z,z),rvJacPoint[F](two,z))==identity;
            "
        ))),
        "rvJacCount[E] = 8\nrvJacCount[F] = 4\nrvJacCount[S] = 32\nrvJacCount[T] = 64\n"
    );
}

#[test]
fn both_maps_require_the_exact_record_of_the_named_direct_chart() {
    for (input, expected) in [
        (
            "rvhyperelliptic Twin[K]{z,two,two,o}; show rvEllGlueMap[S](zero[RVJac(Twin)],zero[RVJac(F)]);",
            "expected type RVJac(E), found RVJac(Twin)",
        ),
        (
            "show rvEllGlueMap[S](zero[RVJac(E)],zero[RVJac(E)]);",
            "expected type RVJac(F), found RVJac(E)",
        ),
        (
            "rvjacchart T[C](o); show rvEllGlueDual[S](zero[RVJac(T)]);",
            "expected type RVJac(S), found RVJac(T)",
        ),
        (
            "rvhyperelliptic Twin[K]{o,z,three,z,four,z,two}; rvjacchart T[Twin](o); show rvEllGlueMap[T](zero[RVJac(E)],zero[RVJac(F)]);",
            "curve 'Twin' has no stored elliptic gluing",
        ),
        (
            "rvjacchart T[S](four); show rvEllGlueDual[T](zero[RVJac(T)]);",
            "curve 'S' has no stored elliptic gluing",
        ),
        (
            "show rvEllGlueDual[C](zero[RVJac(S)]);",
            "require a direct branch chart",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{SPLIT}{input}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
