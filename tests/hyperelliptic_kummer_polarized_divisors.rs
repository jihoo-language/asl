mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const BASIS: &str = r"
    rvunram F[3,1]{0,1};
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
    rvjacbase D[C,W];
    let zw=rvUnramElem[W]{};
    let ow=rvUnramElem[W]{1};
    let tw=rvUnramElem[W]{2};
    let i=rvUnramElem[W]{0,1};
    let e1=rvJacPoint[D](zw,zw);
    let f1=rvJacPoint[D](ow,zw);
    let t2=rvJacPoint[D](tw,zw);
    let ti=rvJacPoint[D](i,zw);
    let e2=add(add(e1,f1),t2);
    let f2=add(add(e1,f1),ti);
    rvjactorsionbasis B[D,2](tw){e1,e2,f1,f2};
";

#[test]
fn example_pairs_generators_and_evaluates_hilbert_scheme_intersections() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-polarized-divisors")),
        concat!(
            "rvKummerDivPair = 4\n",
            "rvKummerDivPair = -2\n",
            "rvKummerDivPair = -2\n",
            "rvKummerDivPair = 2\n",
            "rvKummerDivPair = 1\n",
            "rvKummerDivPair = 0\n",
            "rvKummerFujiki[2, 0, 0, 0] = 48\n",
            "rvKummerFujiki[2, 2, 0, 0] = -8\n",
            "rvKummerFujiki[2, 2, 1, 1] = 12\n",
            "rvKummerFujiki[3, 3, 0, 1] = -384\n",
        )
    );
}

#[test]
fn every_node_and_trope_matches_the_marked_incidence_configuration() {
    let origin = [0, 1, 4, 7, 13, 15];
    let mut script = BASIS.to_owned();
    let mut expected = String::new();
    for i in 0..16 {
        writeln!(script, "let node{i}=rvKummerNodeDiv[B,{i}];").unwrap();
        writeln!(script, "let trope{i}=rvKummerTropeDiv[B,{i}];").unwrap();
    }
    let mut node_sum = "zero[RVKummerDiv(B)]".to_owned();
    for i in origin {
        node_sum = format!("add({node_sum},node{i})");
    }
    writeln!(
        script,
        "assert add(rvKummerDivScale[2](trope0),{node_sum})==rvKummerHyperplane[B];"
    )
    .unwrap();
    for i in 0..16 {
        for j in 0..16 {
            let incidence = i32::from(origin.contains(&(i ^ j)));
            writeln!(script, "rvKummerDivPair(node{i},trope{j});").unwrap();
            writeln!(expected, "rvKummerDivPair = {incidence}").unwrap();
        }
        for j in i..16 {
            writeln!(script, "rvKummerDivPair(trope{i},trope{j});").unwrap();
            writeln!(
                expected,
                "rvKummerDivPair = {}",
                if i == j { -2 } else { 0 }
            )
            .unwrap();
        }
    }
    assert_eq!(stdout(interpret(&script)), expected);
}

#[test]
fn half_vectors_frobenius_and_wide_pairing_are_exact() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{BASIS}
            let h=rvKummerHyperplane[B];
            let e0=rvKummerNodeDiv[B,0];
            let e1div=rvKummerNodeDiv[B,1];
            let t0=rvKummerTropeDiv[B,0];
            let r=rvKummerDiv[B]{{1,-1,-1,0,0,-1,0,0,-1,0,0,0,0,0,-1,0,-1}};
            assert r==t0;
            let half=rvKummerDiv[B]{{0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1}};
            rvKummerDivPair(half,half);
            assert add(half,neg(half))==zero[RVKummerDiv(B)];
            let isotope=rvKummerDivScale[42535295865117307932921825928971026432](add(add(h,e0),e1div));
            rvKummerDivPair(isotope,isotope);
            let mixed=add(h,add(rvKummerTropeDiv[B,8],rvKummerNodeDiv[B,13]));
            let image=rvKummerDivFrob[C,1](mixed);
            assert image==add(h,add(rvKummerTropeDiv[B,10],rvKummerNodeDiv[B,15]));
            assert rvKummerDivFrob[C,1](image)==mixed;
            rvKummerDivPair(mixed,mixed);
            rvKummerDivPair(image,image);
            "
        ))),
        "rvKummerDivPair = -4\nrvKummerDivPair = 0\nrvKummerDivPair = 4\nrvKummerDivPair = 4\n"
    );
}

#[test]
fn fujiki_keeps_boundary_powers_zero_terms_and_wide_cancellation() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{BASIS}
            let h=rvKummerHyperplane[B];
            let e=rvKummerNodeDiv[B,0];
            let t=rvKummerTropeDiv[B,0];
            let zz=zero[RVKummerDiv(B)];
            rvKummerFujiki[1,1,0,0](h,t);
            rvKummerFujiki[3,6,0,1](h,zz);
            rvKummerFujiki[3,1,0,0](t,t);
            let a=rvKummerDivScale[4611686018427387904](h);
            let b=rvKummerDivScale[4611686018427387904](t);
            rvKummerFujiki[2,2,0,0](a,b);
            let huge=rvKummerDivScale[42535295865117307932921825928971026432](h);
            rvKummerFujiki[2,1,0,0](huge,e);
            "
        ))),
        concat!(
            "rvKummerFujiki[1, 1, 0, 0] = 2\n",
            "rvKummerFujiki[3, 6, 0, 1] = -960\n",
            "rvKummerFujiki[3, 1, 0, 0] = -120\n",
            "rvKummerFujiki[2, 2, 0, 0] = 0\n",
            "rvKummerFujiki[2, 1, 0, 0] = 0\n",
        )
    );
}

#[test]
fn diagnoses_invalid_markings_half_vectors_exponents_and_integer_overflow() {
    let setup = format!("{BASIS} let h=rvKummerHyperplane[B];");
    for (source, expected) in [
        ("show rvKummerDiv[B]{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};", "outside the integral marked Kummer divisor lattice"),
        ("show rvKummerDiv[B]{0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};", "outside the integral marked Kummer divisor lattice"),
        ("show rvKummerNodeDiv[B,16];", "between 0 and 15"),
        ("show rvKummerHyperplane[Missing];", "unknown Jacobian torsion basis"),
        ("rvhyperelliptic G[W]{zw,ow}; rvjactorsionbasis Z[G,2](tw){}; show rvKummerHyperplane[Z];", "certified genus-two level-two basis"),
        ("rvjactorsionbasis Twin[D,2](tw){e1,e2,f1,f2}; rvKummerDivPair(h,rvKummerHyperplane[Twin]);", "expected type RVKummerDiv(B)"),
        ("rvKummerFujiki[0,0,0,0](h,h);", "n must be positive"),
        ("rvKummerFujiki[1,3,0,0](h,h);", "between zero and 2n"),
        ("rvKummerFujiki[1,0,0,1](h,h);", "shifts must be zero for n=1"),
        ("rvKummerFujiki[1000000,1000000,0,0](h,h);", "work limit"),
        ("show rvKummerDivScale[85070591730234615865843651857942052864](h);", "divisor coefficient exceeds signed 128-bit range"),
        ("let huge=rvKummerDivScale[42535295865117307932921825928971026432](h); rvKummerDivPair(huge,huge);", "divisor intersection exceeds signed 128-bit range"),
        ("rvKummerFujiki[100,0,0,0](h,h);", "exceeds the supported exact integer range"),
    ] {
        let problem = diagnostic(interpret(&format!("{setup}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
