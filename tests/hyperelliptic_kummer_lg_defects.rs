mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const PREFIX: &str =
    "rvverbitskytate T{5}; rvmf P[T]{2}; rvlgdefect D[T]{0,1}; rvlgdefectpush A[D,P];";

#[test]
fn example_checks_the_factor_exact_residue_and_actual_jordan_coordinates() {
    let output = stdout(example("hyperelliptic-kummer-lg-defects"));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines.len(), 6, "{output}");
    assert!(
        lines[0].starts_with("rvLGDefect[D] = W=x^5−y^5 support=[0, 1] p=[[0, 1], [-1, -1], [1]]")
    );
    assert!(lines[0].ends_with("square=verified"));
    assert_eq!(
        lines[1],
        "rvLGDefectPush[A] = [P2 ranks=[4, 2, 1, 0, 0, 0] stable=[1, 3] projective=0]"
    );
    assert_eq!(lines[2], "rvMFObject[Q] = W=x^5 Tate=T factors=[1, 3]");
    assert!(lines[3].starts_with("RVCyclo(5)"));
    assert_eq!(lines[4],"rvLGDefectFusion[D, D] = channels=[(1, 0), (0, 2)] verified_brane_actions=4 verified_bulk_degrees=4 orientations=2");
    assert_eq!(
        lines[5],
        "rvLGDefectPush[E] = [P2 ranks=[4, 2, 0, 0, 0] stable=[2, 2] projective=0]"
    );
}

#[test]
fn interval_actions_at_all_small_truncation_levels_match_the_fusion_prediction() {
    for d in 3usize..=6 {
        let mut script = format!("rvverbitskytate T{{{d}}};");
        let mut expectations = Vec::new();
        for l in 1..d {
            let left_support = (0..l).map(|r| r.to_string()).collect::<Vec<_>>().join(",");
            writeln!(script, "rvlgdefect D{l}[T]{{{left_support}}};").unwrap();
            for j in 1..d {
                if l == 1 {
                    writeln!(script, "rvmf P{j}[T]{{{j}}};").unwrap();
                }
                writeln!(
                    script,
                    "rvlgdefectpush A{l}_{j}[D{l},P{j}]; rvLGDefectPush[A{l}_{j}];"
                )
                .unwrap();
                let lower = l.abs_diff(j) + 1;
                let upper = (l + j - 1).min(2 * d - 1 - l - j);
                let stable: Vec<_> = (lower..=upper).step_by(2).collect();
                let projective = (l * j - stable.iter().sum::<usize>()) / d;
                expectations.push(format!("stable={stable:?} projective={projective}"));
            }
        }
        let output = stdout(interpret(&script));
        assert_eq!(output.lines().count(), expectations.len());
        for (report, expected) in output.lines().zip(expectations) {
            assert!(
                report.contains(&expected),
                "d={d}, {report} expected {expected}"
            );
        }
    }
    let report = stdout(interpret(
        "rvverbitskytate T{5};
         rvlgdefect A[T]{4,0,1}; rvlgdefect B[T]{3,4,0,1}; rvLGDefectFusion[A,B];
         rvlgdefect C[T]{2}; rvlgdefect E[T]{4}; rvLGDefectFusion[C,E];
         rvlgdefect F[T]{2,3}; rvLGDefectFusion[F,B];",
    ));
    let lines: Vec<_> = report.lines().collect();
    assert!(lines[0].contains("channels=[(4, 1)]"));
    assert!(lines[1].contains("channels=[(1, 0)]"));
    assert!(lines[2].contains("channels=[(1, 2)]"));
}

#[test]
fn all_nontrivial_small_subsets_match_the_independent_character_residue() {
    for d in 2usize..=6 {
        let mut script = format!("rvverbitskytate T{{{d}}};");
        for mask in 1..(1 << d) - 1 {
            let roots: Vec<_> = (0..d).filter(|i| mask & (1 << i) != 0).collect();
            let support = roots
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",");
            writeln!(script, "rvlgdefect D{mask}[T]{{{support}}};").unwrap();
            for k in 0..d - 1 {
                for (operation, sign) in [("rvLGDefectBulk", -1), ("rvLGDefectReverse", 1)] {
                    let terms = roots
                        .iter()
                        .map(|&r| format!("rvCyclo[{d},{}]", sign * (r * (k + 1)) as i64))
                        .reduce(|a, b| format!("add({a},{b})"))
                        .unwrap();
                    writeln!(script, "assert {operation}[D{mask},{k}]=={terms};").unwrap();
                }
            }
        }
        assert!(stdout(interpret(&script)).is_empty());
    }
}

#[test]
fn stable_map_coefficients_retain_block_scales_and_compose_on_graph_pushes() {
    assert!(stdout(interpret(&format!(
        "{PREFIX}
        rvlgdefect G[T]{{1}}; rvmf Q[T]{{3}}; rvlgdefectpush X[G,P]; rvlgdefectpush Y[G,Q];
        let f=rvMFBasis[P,Q,0,0,0,0];
        let g=rvMFBasis[Q,P,0,0,0,0];
        let gf=rvMFCompose(g,f);
        let zi=rvCycloRat[1](rvCyclo[5,-1]);
        let one=rvCycloRat[1](rvCyclo[5,0]);
        assert rvLGDefectMap[X,Y,0,0,1](f)==zi;
        assert rvLGDefectMap[Y,X,0,0,0](g)==one;
        assert rvLGDefectMap[X,X,0,0,1](gf)==rvCycloQMul(rvLGDefectMap[Y,X,0,0,0](g),rvLGDefectMap[X,Y,0,0,1](f));
        let half=rvMFScale[1,2](gf);
        assert rvLGDefectMap[X,X,0,0,1](half)==rvCycloRat[2](rvCyclo[5,-1]);
        rvlgdefect H[T]{{0,1,2,3}}; rvmf L[T]{{4}}; rvlgdefectpush Lp[H,L];
        assert rvLGDefectMap[Lp,Lp,0,0,0](rvMFIdentity[L])==one;
        assert rvLGDefectMap[Lp,Lp,0,0,0](rvMFBasis[L,L,0,0,0,0])==one;
        rvmf Z[T]{{}}; rvlgdefectpush V[D,Z]; rvlgdefectresult Vb[V];
        assert rvMFIdentity[Vb]==zero[RVMF(Vb,Vb,0)];"
    ))).is_empty());
}

#[test]
fn defects_and_stable_maps_validate_subsets_nominal_types_and_bounds() {
    for (query,reason) in [
        ("rvlgdefect F[T]{};", "nonempty proper subset"),
        ("rvlgdefect F[T]{0,1,2,3,4};", "nonempty proper subset"),
        ("rvlgdefect F[T]{5};", "nonempty proper subset"),
        ("rvlgdefect F[T]{1,1};", "occurs twice"),
        ("rvlgdefect D[T]{0};", "already defined"),
        ("rvverbitskytate M{2,2}; rvlgdefect F[M]{0};", "rank-one Tate marking"),
        ("rvverbitskytate M{5}; rvmf W[M]{2}; rvlgdefectpush F[D,W];", "same Tate marking"),
        ("rvlgdefect F[T]{0,2}; rvLGDefectFusion[D,F];", "consecutive circular intervals"),
        ("show rvLGDefectBulk[D,4];", "bulk degree must lie"),
        ("show rvLGDefectCoeff[D,2,0];", "factor index must be zero or one"),
        ("show rvLGDefectMatrix[A,1,0,0,0];", "summand index is out of range"),
        ("show rvLGDefectMap[A,A,0,0,0](rvMFBasis[P,P,1,0,0,0]);", "requires an even Hom class"),
        ("rvmf Q[T]{2}; show rvLGDefectMap[A,A,0,0,0](rvMFIdentity[Q]);", "indicated boundary branes"),
        ("rvlgdefect F[T]{0,1}; rvlgdefectpush B[F,P]; show rvLGDefectMap[A,B,0,0,0](rvMFIdentity[P]);", "same named defect"),
        ("rvverbitskytate M{101}; rvlgdefect F[M]{0,1}; rvmf W[M]{100}; rvlgdefectpush W2[F,W];", "coefficient storage limit"),
    ] {
        let error=diagnostic(interpret(&format!("{PREFIX} {query}")));
        assert!(error.contains(reason),"{query}: {error}");
    }
}
