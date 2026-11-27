mod common;

use std::collections::BTreeSet;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const SPLIT_F11: &str = r"
    rvunram F[11,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic C[F]{o,z,z,z,z,o};
";

const F3: &str = r"
    rvunram F[3,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    rvhyperelliptic C[F]{z,o,o,o,z,o};
    rvhyperelliptic D[F]{o,z,z,o,t,o};
";

#[test]
fn split_genus_two_fiber_has_sixteen_distinct_checked_halves() {
    assert_eq!(
        stdout(example("hyperelliptic-half-fiber")),
        concat!(
            "rvJacHalfFiber[C] = F2 dimension 4, count 16\n",
            "RVJac(C){u=[[10], [9], [1]], v=[[8], [2]]}\n",
            "RVJac(C){u=[[10], [3], [1]], v=[[3], [5]]}\n",
        ),
    );
    let mut source = format!("{SPLIT_F11} let p=rvJacPoint[C](z,o);");
    for index in 0..16 {
        writeln!(
            source,
            "let h{index}=rvJacHalfAt[C,{index}](z,o); assert add(h{index},h{index})==p; show h{index};",
        )
        .unwrap();
    }
    let output = stdout(interpret(&source));
    assert_eq!(output.lines().collect::<BTreeSet<_>>().len(), 16);
}

#[test]
fn nonsplit_and_ramified_fibers_agree_with_witnesses_and_descent() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F3}
            let p=rvJacPoint[C](o,o);
            let r=rvJacPoint[C](z,z);
            assert rvJacDivisible2(p)==S.hi;
            assert rvJacDivisible2(r)==S.hi;
            rvJacHalfFiber[C](o,o);
            rvJacHalfFiber[C](z,z);
            let h=rvJacHalfAt[C,0](o,o);
            let h1=rvJacHalfAt[C,1](o,o);
            assert h==rvJacHalfPoint[C](o,o){{[2],[2],[1],[2],[1]}};
            assert add(h,r)==h1;
            assert add(h,h)==p;
            assert add(h1,h1)==p;
            let q=rvJacHalfAt[C,0](z,z);
            let q1=rvJacHalfAt[C,1](z,z);
            assert q1==rvJacHalfPoint[C](z,z){{[],[1],[2],[],[1]}};
            assert add(q,q)==r;
            assert q1==neg(q);
            show h;
            show h1;
            show q;
            show q1;
        "
        ))),
        concat!(
            "rvJacHalfFiber[C] = F2 dimension 1, count 2\n",
            "rvJacHalfFiber[C] = F2 dimension 1, count 2\n",
            "RVJac(C){u=[[2], [2], [1]], v=[[1], [1]]}\n",
            "RVJac(C){u=[[2], [0], [1]], v=[[2]]}\n",
            "RVJac(C){u=[[2], [0], [1]], v=[[0], [1]]}\n",
            "RVJac(C){u=[[2], [0], [1]], v=[[0], [2]]}\n",
        ),
    );
}

#[test]
fn self_paired_factors_give_empty_fibers_on_a_curve_and_at_a_branch_point() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F3}
            rvhyperelliptic E[F]{{z,o,o,o,t,o}};
            assert rvJacDivisible2(rvJacPoint[D](t,o))==S.lo;
            assert rvJacDivisible2(rvJacPoint[E](t,z))==S.lo;
            rvJacHalfFiber[D](t,o);
            rvJacHalfFiber[E](t,z);
        "
        ))),
        "rvJacHalfFiber[D] = empty\nrvJacHalfFiber[E] = empty\n",
    );
}

#[test]
fn extension_field_coefficients_and_genus_zero() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[3,1]{1,0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            let t=rvUnramElem[F]{2};
            let a=rvUnramElem[F]{0,1};
            let b=add(a,neg(o));
            let c=neg(b);
            rvhyperelliptic C[F]{b,neg(a),c,z,neg(a),o};
            let p=rvJacPoint[C](c,o);
            let h=rvJacHalfAt[C,3](c,o);
            assert add(h,h)==p;
            rvJacHalfFiber[C](c,o);
            rvhyperelliptic L[F]{o,o};
            rvJacHalfFiber[L](z,t);
            rvJacHalfFiber[L](t,z);
            assert rvJacHalfAt[L,0](z,t)==zero[RVJac(L)];
            assert rvJacHalfAt[L,00](t,z)==zero[RVJac(L)];
        ",
        )),
        concat!(
            "rvJacHalfFiber[C] = F2 dimension 2, count 4\n",
            "rvJacHalfFiber[L] = F2 dimension 0, count 1\n",
            "rvJacHalfFiber[L] = F2 dimension 0, count 1\n",
        ),
    );
}

#[test]
fn factorization_when_base_or_factor_field_cardinality_exceeds_u64() {
    for modulus in ["4,1,0,1", "3,1,0,0,0,1"] {
        assert_eq!(
            stdout(interpret(&format!(
                r"
                rvunram F[65537,1]{{{modulus}}};
                let z=rvUnramElem[F]{{}};
                let o=rvUnramElem[F]{{1}};
                rvhyperelliptic C[F]{{o,z,z,o}};
                let p=rvJacPoint[C](z,o);
                let h0=rvJacHalfAt[C,0](z,o);
                let h1=rvJacHalfAt[C,1](z,o);
                assert add(h0,h0)==p;
                assert add(h1,h1)==p;
                assert add(h0,rvJacPoint[C](neg(o),z))==h1;
                rvJacHalfFiber[C](z,o);
            "
            ))),
            "rvJacHalfFiber[C] = F2 dimension 1, count 2\n",
        );
    }
}

#[test]
fn rejects_empty_indexed_fibers_bad_indices_and_invalid_points() {
    for (source, expected) in [
        ("show rvJacHalfAt[C,2](o,o);", "index is outside"),
        (
            "show rvJacHalfAt[C,340282366920938463463374607431768211456](o,o);",
            "index is outside",
        ),
        (
            "show rvJacHalfAt[C,-1](o,o);",
            "nonnegative decimal integer",
        ),
        ("show rvJacHalfAt[D,0](t,o);", "empty rational half fiber"),
        ("show rvJacHalfAt[C,0](o,z);", "affine point on the curve"),
        ("rvJacHalfFiber[C](o,z);", "affine point on the curve"),
        (
            "rvunram G[3,1]{0,1}; rvJacHalfFiber[C](rvUnramElem[G]{1},o);",
            "rvJacHalfFiber coordinate",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{F3}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
