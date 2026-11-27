mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const MARK: &str = "rvverbitskytate T{5}; rvmf P[T]{2};";

#[test]
fn rank_one_example_connects_periodic_tate_action_cones_and_trace() {
    let report = stdout(example("hyperelliptic-kummer-landau-ginzburg"));
    let lines = report.lines().collect::<Vec<_>>();
    assert_eq!(lines[0], "rvMFObject[P] = W=x^5 Tate=T factors=[2]");
    assert_eq!(
        lines[1],
        "rvMFHom[P, P, 0] = dimension=2 blocks=[(0, 0, first_degree=0, count=2)]"
    );
    assert_eq!(
        lines[2],
        "rvMFHom[P, P, 1] = dimension=2 blocks=[(0, 0, first_degree=1, count=2)]"
    );
    assert_eq!(&lines[3..5], ["rvMFTrace = 1", "rvMFPair = 1"]);
    assert_eq!(lines[5], "rvMFCone = even_to_odd=[(0, 0, 2, 1), (0, 1, 1, 1), (1, 1, 3, -1)] odd_to_even=[(0, 0, 3, 1), (0, 1, 1, 1), (1, 1, 2, -1)] Smith=[1, 4] factors=[1, 4] contractible=0");
    assert_eq!(lines[6], "rvMFObject[C] = W=x^5 Tate=T factors=[1, 4]");
    assert_eq!(lines[7], "rvMFObject[Q] = W=x^5 Tate=T factors=[3]");
    assert_eq!(lines[8], "rvMFTrace = -1");
}

#[test]
fn all_small_rank_one_smith_cones_and_residue_gram_entries_match_the_oracle() {
    for d in 2usize..=7 {
        let mut source = format!("rvverbitskytate T{{{d}}};");
        for i in 1..d {
            writeln!(source, "rvmf P{i}[T]{{{i}}};").unwrap();
        }
        let mut results = Vec::new();
        for i in 1..d {
            for j in 1..d {
                let m = i.min(j).min(d - i).min(d - j);
                for a in 0..m {
                    writeln!(source, "rvMFCone(rvMFBasis[P{i},P{j},0,0,0,{a}]);").unwrap();
                    let first = j.min(a + j.saturating_sub(i)).min(d - i);
                    let second = d - i + j - first;
                    results.push(format!("Smith=[{first}, {second}]"));
                    for b in 0..m {
                        writeln!(source, "rvMFPair(rvMFBasis[P{i},P{j},0,0,0,{a}],rvMFBasis[P{j},P{i},1,0,0,{b}]);").unwrap();
                        results.push(format!("rvMFPair = {}", usize::from(a + b + 1 == m)));
                    }
                }
            }
        }
        let report = stdout(interpret(&source));
        let rows = report.lines().collect::<Vec<_>>();
        assert_eq!(rows.len(), results.len());
        for (row, result) in rows.into_iter().zip(results) {
            assert!(row.contains(&result), "d={d}, expected {result}, got {row}");
        }
    }
}

#[test]
fn polynomial_block_smith_reduction_handles_cancellation_and_contractibles() {
    let text = stdout(interpret(
        r"
        rvverbitskytate T{5};
        rvmf S[T]{2,2};
        let off=add(rvMFBasis[S,S,0,0,1,1],rvMFBasis[S,S,0,1,0,1]);
        let a=add(add(rvMFBasis[S,S,0,0,0,0],rvMFBasis[S,S,0,0,0,1]),off);
        let f=add(rvMFRepresentative(a),add(rvMFTerm[S,S,0,1,1,0,2],rvMFTerm[S,S,0,1,1,1,2]));
        assert rvMFClass(f)==a;
        let defect=add(f,neg(rvMFRepresentative(a)));
        assert rvMFDiff(rvMFPrimitive(defect))==defect;
        rvMFCone(f);
        rvmfcone C(f);
        rvMFObject[C];
        rvMFHom[C,C,1];
        let repeated=add(off,add(rvMFBasis[S,S,0,0,0,1],rvMFBasis[S,S,0,1,1,1]));
        rvMFCone(repeated);
        rvmfcone Z(rvMFIdentity[S]);
        rvMFObject[Z];
        assert rvMFIdentity[Z]==zero[RVMF(Z,Z,0)];

        rvverbitskytate M{7};
        rvmf A[M]{1,3};
        rvmf B[M]{2,4};
        let oneplusx=add(rvMFBasis[A,B,0,1,1,0],rvMFBasis[A,B,0,1,1,1]);
        let first=add(rvMFRepresentative(rvMFBasis[A,B,0,0,0,0]),add(rvMFTerm[A,B,0,0,0,0,1],rvMFTerm[A,B,0,0,0,1,2]));
        let rest=add(oneplusx,add(rvMFScale[2,1](rvMFBasis[A,B,0,0,1,0]),neg(rvMFBasis[A,B,0,1,0,0])));
        let g=add(first,rvMFRepresentative(rest));
        rvMFCone(g);
        ",
    ));
    let rows = text.lines().collect::<Vec<_>>();
    assert!(rows[0].ends_with("Smith=[0, 2, 3, 5] factors=[2, 3] contractible=2"));
    assert_eq!(rows[1], "rvMFObject[C] = W=x^5 Tate=T factors=[2, 3]");
    assert!(rows[2].contains("dimension=8"));
    assert!(rows[3].ends_with("Smith=[1, 2, 3, 4] factors=[1, 2, 3, 4] contractible=0"));
    assert_eq!(rows[4], "rvMFObject[Z] = W=x^5 Tate=T factors=[]");
    assert!(rows[5].ends_with("Smith=[0, 2, 7, 7] factors=[2] contractible=3"));
}

#[test]
fn raw_differential_composition_and_shift_have_actual_signed_homotopies() {
    let output = stdout(interpret(&format!(
        "{MARK}
         rvmf Q[T]{{4,1}};
         rvmfshift A[P];
         rvmfshift B[Q];
         rvmfshift PP[A];
         rvmfshift QQ[B];
         let a=rvMFTerm[P,Q,1,0,0,1,3,(2,3)];
         let b=rvMFTerm[Q,P,1,0,0,0,2,7];
         let da=rvMFDiff(a);
         let db=rvMFDiff(b);
         assert rvMFDiff(da)==zero[RVMFCochain(P,Q,1)];
         assert rvMFDiff(rvMFPrimitive(da))==da;
         assert rvMFClass(da)==zero[RVMF(P,Q,0)];
         assert rvMFDiff(rvMFPrimitive(db))==db;
         assert rvMFDiff(rvMFCompose(b,a))==add(rvMFCompose(db,a),neg(rvMFCompose(b,da)));
         assert rvMFShift[A,B](da)==rvMFDiff(rvMFShift[A,B](a));
         assert rvMFShift[A,A](rvMFCompose(b,a))==rvMFCompose(rvMFShift[B,A](b),rvMFShift[A,B](a));
         let q=rvMFBasis[P,Q,1,1,0,0];
         assert rvMFShift[A,B](q)==neg(rvMFBasis[A,B,1,1,0,0]);
         assert rvMFShift[PP,QQ](rvMFShift[A,B](q))==rvMFBasis[PP,QQ,1,1,0,0];
         let t=rvMFTerm[P,Q,0,0,0,1,0,(1,2)];
         let dt=rvMFDiff(t);
         assert rvMFDiff(rvMFPrimitive(dt))==dt;
         "
    )));
    assert!(output.is_empty());
}

#[test]
fn invalid_domains_maps_and_nominal_types_are_diagnosed() {
    for (query, reason) in [
        ("rvmf X[T]{0};", "index must lie"),
        ("rvmf X[T]{5};", "index must lie"),
        (
            "rvverbitskytate U{2,2}; rvmf X[U]{1};",
            "requires a rank-one Tate marking",
        ),
        ("rvmf P[T]{1};", "already defined"),
        ("rvMFHom[P,P,2];", "parity must be zero or one"),
        ("rvMFHom[P,Q,0];", "unknown matrix factorization"),
        (
            "let q=rvMFTerm[P,P,1,1,0,0,0];",
            "row or column is out of range",
        ),
        ("let q=rvMFTerm[P,P,1,0,0,2,0];", "slot must be zero or one"),
        (
            "let q=rvMFBasis[P,P,0,0,0,2];",
            "basis index is out of range",
        ),
        ("let q=rvMFClass(rvMFTerm[P,P,1,0,0,0,1]);", "is not closed"),
        (
            "let q=rvMFPrimitive(rvMFRepresentative(rvMFIdentity[P]));",
            "is not a boundary",
        ),
        (
            "rvMFCone(rvMFBasis[P,P,1,0,0,0]);",
            "requires an even morphism",
        ),
        ("rvMFCone(rvMFTerm[P,P,0,0,0,0,1]);", "is not closed"),
        (
            "let q=rvMFShift[P,P](rvMFIdentity[P]);",
            "must be the indicated shifted objects",
        ),
        (
            "rvmf Q[T]{2}; rvMFTrace(rvMFBasis[P,Q,1,0,0,0]);",
            "requires an endomorphism",
        ),
        (
            "rvmf Q[T]{2}; let z=rvMFCompose(rvMFIdentity[Q],rvMFIdentity[P]);",
            "matching middle objects",
        ),
        (
            "rvverbitskytate U{5}; rvmf Q[U]{2}; rvMFHom[P,Q,0];",
            "one common Tate marking",
        ),
        (
            "let z=rvMFStableAction(rvVerbitskyTateTerm[T,ordinary]{0}{0}{0},rvMFIdentity[P]);",
            "even stable Tate class",
        ),
        (
            "let z=rvMFStableAction(rvVerbitskyTateTerm[T,stable]{0}{1}{0},rvMFIdentity[P]);",
            "even stable Tate class",
        ),
    ] {
        let actual = diagnostic(interpret(&format!("{MARK} {query}")));
        assert!(actual.contains(reason), "{query}: {actual}");
    }
}
