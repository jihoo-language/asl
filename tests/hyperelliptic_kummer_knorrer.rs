mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const PREFIX: &str = "rvverbitskytate T{5}; rvmf P[T]{2}; rvknorrer K[P];";

#[test]
fn doubled_factor_example_checks_actual_differential_residue_and_perturbed_cone() {
    let output = stdout(example("hyperelliptic-kummer-knorrer"));
    let lines = output.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 7, "{output}");
    let object = lines[0];
    assert!(object.starts_with("rvKnorrerObject[K] = W=x^5+uv factors=[2]"));
    for entry in [
        "(0, 2, [3, 0, 0], 1)",
        "(1, 3, [3, 0, 0], 1)",
        "(2, 0, [2, 0, 0], 1)",
        "(3, 1, [2, 0, 0], 1)",
        "(0, 1, [0, 0, 1], 1)",
        "(1, 0, [0, 1, 0], 1)",
        "(2, 3, [0, 0, 1], -1)",
        "(3, 2, [0, 1, 0], -1)",
    ] {
        assert!(object.contains(entry), "{object} is missing {entry}");
    }
    assert!(object.ends_with("square=verified"));
    assert_eq!(
        &lines[1..5],
        [
            "rvMFTrace = 1",
            "rvKnorrerTrace = 1",
            "rvKnorrerWeights = {(17, -1)}",
            "rvKnorrerWeights = {(12, -1)}",
        ]
    );
    assert!(lines[5].starts_with("rvKnorrerCone = homotopy=[("));
    assert!(lines[5].contains("Smith=[1, 4] factors=[1, 4] contractible=0"));
    assert!(lines[5].ends_with("tensor_cone=verified"));
    assert!(lines[6].starts_with("rvKnorrerObject[KC] = W=x^5+uv factors=[1, 4]"));
}

#[test]
fn mixed_polynomial_matrix_units_obey_contraction_differential_and_homogeneous_signs() {
    let mut source = format!(
        "{PREFIX} rvmf S[T]{{1,4}}; rvknorrer KS[S];
         rvmfshift PP[P]; rvknorrer KP[PP];
         rvmfshift SS[S]; rvknorrer KT[SS];"
    );
    let id = |row: usize, column: usize| row * 4 + column;
    let parity = |index| u8::from((index % 4 == 1) || (index % 4 == 2));
    for row in 0..8 {
        for column in 0..4 {
            let i = id(row, column);
            let p = parity(row) ^ parity(column);
            // Exercise all source/target slots with a mixed monomial, a constant in u,v,
            // and a pure v power (the homotopy's second even branch).
            writeln!(
                source,
                "let a{i}=add(rvKnorrerTerm[K,KS,{p},{row},{column},(-3,7)]{{1,2,3}},
                           add(rvKnorrerTerm[K,KS,{p},{row},{column}]{{3,0,0}},
                               rvKnorrerTerm[K,KS,{p},{row},{column}]{{0,0,2}}));
                 let d{i}=rvKnorrerDiff(a{i});
                 let h{i}=rvKnorrerHomotopy(a{i});
                 let ip{i}=rvKnorrerInclude[K,KS](rvKnorrerProject(a{i}));
                 assert add(rvKnorrerDiff(h{i}),rvKnorrerHomotopy(d{i}))==add(a{i},neg(ip{i}));
                 assert rvKnorrerDiff(d{i})==zero[RVKnorrerCochain(K,KS,{p})];
                 assert rvKnorrerHomotopy(h{i})==zero[RVKnorrerCochain(K,KS,{p})];
                 assert rvKnorrerDiff(rvKnorrerShift[KP,KT](a{i}))==rvKnorrerShift[KP,KT](d{i});
                 assert rvKnorrerProject(d{i})==rvMFDiff(rvKnorrerProject(a{i}));"
            )
            .unwrap();
        }
    }
    source.push_str(
        "let a=rvKnorrerTerm[K,KS,1,5,0,(2,3)]{1,1,1};
         let b=rvKnorrerTerm[KS,K,1,0,5,9]{0,0,1};
         let ba=rvKnorrerCompose(b,a);
         assert ba==rvKnorrerTerm[K,K,0,0,0,6]{1,1,2};
         let c=rvKnorrerTerm[K,KS,1,6,0,(2,3)]{1,1,1};
         let e=rvKnorrerTerm[KS,K,0,1,6,9]{0,0,1};
         assert rvKnorrerCompose(e,c)==rvKnorrerTerm[K,K,1,1,0,6]{1,1,2};
         assert rvKnorrerDiff(ba)==add(rvKnorrerCompose(rvKnorrerDiff(b),a),
                                      neg(rvKnorrerCompose(b,rvKnorrerDiff(a))));
         let f=rvMFTerm[P,S,1,1,0,1,4,(2,3)];
         assert rvKnorrerDiff(rvKnorrerInclude[K,KS](f))==rvKnorrerInclude[K,KS](rvMFDiff(f));
         assert rvKnorrerProject(rvKnorrerInclude[K,KS](f))==f;",
    );
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn three_variable_residue_matches_every_small_rank_one_odd_trace() {
    let mut expected = Vec::new();
    let mut source = String::new();
    for d in 2usize..=7 {
        writeln!(source, "rvverbitskytate T{d}{{{d}}};").unwrap();
        for i in 1..d {
            let t = format!("{d}_{i}");
            writeln!(
                source,
                "rvmf P{t}[T{d}]{{{i}}}; rvknorrer K{t}[P{t}]; rvKnorrerObject[K{t}];"
            )
            .unwrap();
            expected.push(None);
            let dim = i.min(d - i);
            for k in 0..dim {
                writeln!(source, "let a{t}_{k}=rvMFBasis[P{t},P{t},1,0,0,{k}]; rvMFTrace(a{t}_{k}); rvKnorrerTrace(rvKnorrerInclude[K{t},K{t}](a{t}_{k}));").unwrap();
                let trace = usize::from(k == dim - 1);
                expected.push(Some(format!("rvMFTrace = {trace}")));
                expected.push(Some(format!("rvKnorrerTrace = {trace}")));
            }
        }
    }
    let text = stdout(interpret(&source));
    assert_eq!(text.lines().count(), expected.len());
    for (line, expected) in text.lines().zip(expected) {
        if let Some(expected) = expected {
            assert_eq!(line, expected);
        } else {
            assert!(line.ends_with("square=verified"), "{line}");
        }
    }
}

#[test]
fn class_reduction_detects_polynomial_boundaries_and_preserves_residue() {
    let source = format!(
        "{PREFIX}
         let f=rvKnorrerInclude[K,K](rvMFBasis[P,P,1,0,0,1]);
         let raw=rvKnorrerRepresentative(f);
         let h=rvKnorrerTerm[K,K,0,0,0,(-8,3)]{{2,2,0}};
         let b=rvKnorrerDiff(h);
         let g=add(b,raw);
         assert rvKnorrerClass(g)==f;
         assert rvKnorrerClass(rvKnorrerScale[3,2](g))==rvKnorrerScale[3,2](f);
         assert rvKnorrerProject(rvKnorrerClass(g))==rvMFBasis[P,P,1,0,0,1];
         assert rvKnorrerDiff(rvKnorrerPrimitive(b))==b;
         let base=rvMFDiff(rvMFTerm[P,P,1,0,0,0,0]);
         let cb=add(rvKnorrerInclude[K,K](base),
                    rvKnorrerDiff(rvKnorrerTerm[K,K,1,3,2]{{1,3,2}}));
         assert rvKnorrerProject(cb)==base;
         assert rvKnorrerClass(cb)==zero[RVKnorrer(K,K,0)];
         assert rvKnorrerDiff(rvKnorrerPrimitive(cb))==cb;
         rvKnorrerTrace(g);
         rvKnorrerTrace(b);
         rvmf Z[T]{{}}; rvknorrer KZ[Z];
         assert rvKnorrerInclude[KZ,KZ](rvMFIdentity[Z])==zero[RVKnorrer(KZ,KZ,0)];
         rvKnorrerObject[KZ];"
    );
    let output = stdout(interpret(&source));
    assert_eq!(
        output.lines().collect::<Vec<_>>(),
        [
            "rvKnorrerTrace = 1",
            "rvKnorrerTrace = 0",
            "rvKnorrerObject[KZ] = W=x^5+uv factors=[] differential=[] square=verified",
        ],
    );
}

#[test]
fn bad_endpoints_parities_cycles_and_exponent_overflow_are_diagnosed() {
    for (query, reason) in [
        ("rvknorrer K[P];", "already defined"),
        ("rvknorrer L[Q];", "unknown matrix factorization"),
        ("rvKnorrerObject[U];", "unknown Knörrer object"),
        ("let f=rvKnorrerTerm[K,K,2,0,0]{0,0,0};", "parity must be zero or one"),
        ("let f=rvKnorrerTerm[K,K,0,0,0]{0,0};", "requires three polynomial exponents"),
        ("let f=rvKnorrerTerm[K,K,1,0,4]{0,0,0};", "row or column is out of range"),
        ("let f=rvKnorrerTerm[K,K,0,0,1]{0,0,0};", "wrong matrix parity"),
        ("let f=rvKnorrerClass(rvKnorrerTerm[K,K,0,0,0]{0,0,0});", "is not closed"),
        ("let f=rvKnorrerPrimitive(rvKnorrerRepresentative(rvKnorrerInclude[K,K](rvMFIdentity[P])));", "is not a boundary"),
        ("rvKnorrerCone(rvKnorrerTerm[K,K,1,0,1]{0,0,0});", "requires an even cochain"),
        ("rvKnorrerCone(rvKnorrerTerm[K,K,0,0,0]{0,0,0});", "is not closed"),
        ("rvknorrer L[P]; rvKnorrerTrace(rvKnorrerInclude[K,L](rvMFIdentity[P]));", "requires an endomorphism"),
        ("rvmf Q[T]{2}; let f=rvKnorrerInclude[K,K](rvMFIdentity[Q]);", "requires the indicated underlying objects"),
        ("rvknorrer L[P]; let f=rvKnorrerCompose(rvKnorrerInclude[K,K](rvMFIdentity[P]),rvKnorrerInclude[L,L](rvMFIdentity[P]));", "requires matching middle objects"),
        ("let f=rvKnorrerShift[K,K](rvKnorrerInclude[K,K](rvMFIdentity[P]));", "requires the indicated shifted underlying objects"),
        ("rvverbitskytate U{5}; rvmf Q[U]{2}; rvknorrer L[Q]; let f=zero[RVKnorrer(K,L,0)];", "one common Tate marking"),
        ("let f=rvKnorrerCompose(rvKnorrerTerm[K,K,0,0,0]{0,340282366920938463463374607431768211455,0},rvKnorrerTerm[K,K,0,0,0]{0,1,0});", "exponent exceeds the unsigned integer range"),
    ] {
        let actual = diagnostic(interpret(&format!("{PREFIX} {query}")));
        assert!(actual.contains(reason), "{query}: {actual}");
    }
}
