mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/twisted-double-boundary-recoupling.asl");

fn cyclic(p: usize) -> String {
    let group = (0..p).map(|i| i.to_string()).collect::<Vec<_>>().join(",");
    format!(
        "rvgroup G[{p}];rvdwcoeff P[G,{p}];rvdouble D[P](rvPostnikovK[P]);
         rvdoubleboundary O[D]{{0}};rvdoubleboundary F[D]{{{group}}};
         rvdoubledefects A[O,F];rvdoubledefects B[F,O];
         rvdoubledefects C[O,O];rvdoubledefects E[F,F];
         let x=rvDoubleDefectBasis[A,0];let y=rvDoubleDefectBasis[B,0];
         let z=rvDoubleDefectBasis[C,0];let f=rvDoubleDefectSixJ[C,E](x,y,x,x);"
    )
}

fn type_iii() -> &'static str {
    let start = EXAMPLE.find("rvgroup G8").unwrap();
    EXAMPLE[start..].split_once("let mix=").unwrap().0
}

#[test]
fn alternating_boundaries_give_the_exact_fourier_cross_ratios() {
    for p in [2, 3, 5] {
        let mut source = cyclic(p);
        let n = p * p;
        writeln!(
            source,
            "let one=rvCycloRat[1](rvCyclo[{n},0]); let unit=rvCycloRat[{p}](rvCyclo[{n},0]);
             let large=rvCycloRat[1](rvCycloScale[{p}](rvCyclo[{n},0]));"
        )
        .unwrap();
        for r in 0..p {
            writeln!(
                source,
                "assert rvDoubleDefectSixJNorm[0,{r}](f)==one;
                assert rvDoubleDefectSixJNorm[1,{r}](f)==large;"
            )
            .unwrap();
            for c in 0..p {
                writeln!(
                    source,
                    "let f{r}{c}=rvDoubleDefectSixJEntry[{r},{c}](f);
                     assert rvCycloQMul(f{r}{c},rvDoubleDefectSixJInvEntry[{c},{r}](f))==unit;"
                )
                .unwrap();
            }
        }
        for a in 0..p {
            for b in 0..p {
                for c in 0..p {
                    for d in 0..p {
                        let phase = (a + p - b) * (c + p - d) % p * p;
                        writeln!(source,
                            "assert rvCycloQMul(f{a}{c},f{b}{d})==rvCycloQMul(rvCycloRat[1](rvCyclo[{n},{phase}]),rvCycloQMul(f{a}{d},f{b}{c}));"
                        ).unwrap();
                    }
                }
            }
        }
        if p == 2 {
            source.push_str(
                "assert f00==unit;assert f01==unit;assert f10==unit;assert f11==neg(unit);",
            );
        }
        source.push_str("rvDoubleDefectPentagon[C,E,C,A,B](x,y,x,y,z);");
        assert_eq!(
            stdout(interpret(&source)),
            format!("rvDoubleDefectPentagon = dimensions=[{p}, {p}, {p}, {p}, {p}], conductor={n}, mixing_edges=2, verified=true\n")
        );
    }
}

#[test]
fn example_exhibits_genuinely_complex_and_heterogeneous_recoupling() {
    let text = stdout(example("twisted-double-boundary-recoupling"));
    assert!(text.contains("RVCycloQ(9){[0, 0, 0, 1] / 3}"));
    assert!(text.contains("RVCycloQ(9){[-1, 0, 0, -1] / 1}"));
    assert!(text.contains("left=PP[[4, 0, 0], [4, 0, 1], [4, 1, 0], [4, 1, 1]], right=QQ[[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]"));
    assert!(text.contains("dimensions=[4, 4, 4, 4, 4], conductor=4, mixing_edges=2, verified=true"));
}

#[test]
fn double_multiplicities_mix_and_both_pentagon_routes_reach_a_noninvertible_output() {
    let mut source = format!("{}
         let mix=rvDoubleDefectSixJ[PP,QQ](y0,x4,y4,y0);
         let one=rvCycloRat[1](rvCyclo[4,0]);let zeroq=zero[RVCycloQ(4)];let half=rvCycloRat[2](rvCyclo[4,0]);",
         type_iii());
    let signs = [[1, 0, 0, 1], [1, 0, 0, -1], [0, 1, 1, 0], [0, -1, 1, 0]];
    for (r, row) in signs.into_iter().enumerate() {
        for (c, sign) in row.into_iter().enumerate() {
            let (coefficient, inverse) = match sign {
                -1 => ("neg(one)", "neg(half)"),
                1 => ("one", "half"),
                _ => ("zeroq", "zeroq"),
            };
            writeln!(
                source,
                "assert rvDoubleDefectSixJEntry[{r},{c}](mix)=={coefficient};
                 assert rvDoubleDefectSixJInvEntry[{c},{r}](mix)=={inverse};"
            )
            .unwrap();
        }
    }
    source.push_str(
        "let m=rvDoubleDefectBasis[PP,4];
        rvDoubleDefectPentagon[PP,QQ,PP,Y,X](y0,x4,y0,x0,m);",
    );
    assert_eq!(stdout(interpret(&source)),
        "rvDoubleDefectPentagon = dimensions=[8, 8, 8, 8, 8], conductor=4, mixing_edges=3, verified=true\n");
}

#[test]
fn point_recoupling_matches_the_marked_cocycle_phase_and_empty_channels() {
    let source = format!("{}
        let a=rvDoubleDefectBasis[QQ,1];let b=rvDoubleDefectBasis[QQ,2];let c=rvDoubleDefectBasis[QQ,4];
        let d=rvDoubleDefectBasis[QQ,7];let bad=rvDoubleDefectBasis[QQ,0];
        let point=rvDoubleDefectSixJ[QQ,QQ](a,b,c,d);
        assert rvDoubleDefectSixJEntry[0,0](point)==neg(rvCycloRat[1](rvCyclo[4,0]));
        assert rvDoubleDefectSixJInvEntry[0,0](point)==rvDoubleDefectSixJEntry[0,0](point);
        rvDoubleDefectSixJData(rvDoubleDefectSixJ[QQ,QQ](a,b,c,bad));
        rvDoubleDefectPentagon[QQ,QQ,QQ,QQ,QQ](a,b,c,d,bad);",type_iii());
    assert_eq!(stdout(interpret(&source)),
        concat!("rvDoubleDefectSixJData = dimension=0, conductor=4, objects=[QQ:1,QQ:2,QQ:4,QQ:0], left=QQ[], right=QQ[]\n",
                "rvDoubleDefectPentagon = dimensions=[1, 1, 1, 1, 1], conductor=4, mixing_edges=0, verified=true\n"));
}

#[test]
fn typed_matrices_and_nominal_categorical_endpoints_are_checked() {
    let setup = cyclic(3);
    assert_eq!(
        stdout(interpret(&format!("{setup}show f;")))
            .split_once('{')
            .unwrap()
            .0,
        "RVDoubleDefectSixJ(A,9)"
    );
    for (operation, expected) in [
        (
            "show rvDoubleDefectSixJ[E,E](x,y,x,x);",
            "ordered outer boundary",
        ),
        (
            "rvDoubleDefectPentagon[C,E,C,A,A](x,y,x,y,z);",
            "ordered outer boundary",
        ),
        (
            "show rvDoubleDefectSixJEntry[0,0](x);",
            "a boundary recoupling matrix",
        ),
        (
            "show rvDoubleDefectSixJEntry[3,0](f);",
            "matrix index is out of range",
        ),
        (
            "show rvDoubleDefectSixJNorm[2,0](f);",
            "norm side or path index is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{operation}")));
        assert!(error.contains(expected), "{expected} in {error}");
    }
}

#[test]
fn changing_the_cocycle_gauge_changes_point_symbols_and_preserves_mixed_coherence() {
    let b = (0..49)
        .map(|i| usize::from(i == 1).to_string())
        .collect::<Vec<_>>()
        .join(",");
    let source = format!(
        "{}let shifted=rvPostnikovGauge(k8,rvPostnikovCochain[P8,2]{{{b}}});
         rvdouble W[P8](shifted);
         rvdoubleboundary PlaneW[W]{{0,1,2,3}};rvdoubleboundary PointW[W]{{0}};
         rvdoubledefects YW[PlaneW,PointW];rvdoubledefects XW[PointW,PlaneW];
         rvdoubledefects PW[PlaneW,PlaneW];rvdoubledefects QW[PointW,PointW];
         let a=rvDoubleDefectBasis[QW,1];let b=rvDoubleDefectBasis[QW,2];let c=rvDoubleDefectBasis[QW,4];
         let d=rvDoubleDefectBasis[QW,7];let f=rvDoubleDefectSixJ[QW,QW](a,b,c,d);
         assert rvDoubleDefectSixJEntry[0,0](f)==rvCycloRat[1](rvCyclo[4,0]);
         let a0=rvDoubleDefectBasis[YW,0];let a4=rvDoubleDefectBasis[YW,1];
         let b0=rvDoubleDefectBasis[XW,0];let b4=rvDoubleDefectBasis[XW,1];
         rvDoubleDefectPentagon[PW,QW,PW,YW,XW](a0,b4,a4,b0,rvDoubleDefectBasis[PW,0]);",
        type_iii()
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvDoubleDefectPentagon = dimensions=[4, 4, 4, 4, 4], conductor=4, mixing_edges=2, verified=true\n"
    );
}

#[test]
fn pentagon_has_five_genuinely_distinct_subgroup_boundaries() {
    let source = "
        rvgroup G{{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}};rvdwcoeff P[G,2];
        rvdouble D[P](rvPostnikovK[P]);
        rvdoubleboundary H0[D]{0};rvdoubleboundary H1[D]{0,1};
        rvdoubleboundary H2[D]{0,2};rvdoubleboundary H3[D]{0,3};
        rvdoubleboundary H4[D]{0,1,2,3};
        rvdoubledefects C01[H0,H1];rvdoubledefects C12[H1,H2];
        rvdoubledefects C23[H2,H3];rvdoubledefects C34[H3,H4];
        rvdoubledefects C02[H0,H2];rvdoubledefects C13[H1,H3];
        rvdoubledefects C24[H2,H4];rvdoubledefects C03[H0,H3];
        rvdoubledefects C14[H1,H4];rvdoubledefects C04[H0,H4];
        let a=rvDoubleDefectBasis[C01,0];let b=rvDoubleDefectBasis[C12,0];
        let c=rvDoubleDefectBasis[C23,0];let d=rvDoubleDefectBasis[C34,0];
        let q=rvDoubleDefectBasis[C04,0];
        rvDoubleDefectPentagon[C02,C13,C24,C03,C14](a,b,c,d,q);";
    let output = stdout(interpret(source));
    assert!(output.starts_with("rvDoubleDefectPentagon = dimensions=[4, 4, 4, 4, 4], conductor=4"));
    assert!(output.ends_with("verified=true\n"));
}
