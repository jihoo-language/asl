mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const FIXTURE: &str = include_str!("../examples/twisted-double-borromean.asl");
const CARRY: &str = "
    rvgroup G[3];rvdwcoeff P[G,3];
    let k=rvPostnikovCochain[P,3]{0,1,1,1,0,2,2,2};
    rvdouble D[P](k);
";

fn type_iii() -> &'static str {
    FIXTURE.split_once("// Borromean closure").unwrap().0
}

fn list(coordinates: impl IntoIterator<Item = usize>) -> String {
    coordinates
        .into_iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

#[test]
fn borromean_operator_is_minus_identity_on_all_eight_columns() {
    assert_eq!(
        stdout(example("twisted-double-borromean")),
        concat!(
            "rvDoubleRepData[D,8] = isotropic=[0, 1, 2, 3], coset_representatives=[0, 4], dimension=2\n",
            "rvDoubleRep[D,8,4,0] = target=1, phase=0\n",
            "RVCycloQ(4){[-8] / 1}\n",
            "RVCycloQ(4){[8] / 1}\n"
        )
    );

    let mut source = FIXTURE.to_string();
    for i in 0..8 {
        writeln!(source, "rvDoubleBraidImage[{i}](link);").unwrap();
    }
    source.push_str(
        "assert rvDoubleBraidEntry[4,4](link)==rvCycloRat[1](rvCyclo[4,2]);
         assert rvDoubleBraidEntry[3,4](link)==zero[RVCycloQ(4)];",
    );
    let actual = stdout(interpret(&source));
    for (i, line) in actual.lines().skip(4).enumerate() {
        assert_eq!(
            line,
            format!("rvDoubleBraidImage[{i}] = target={i}, phase=2")
        );
    }
    assert_eq!(actual.lines().count(), 12);
}

#[test]
fn projective_matrix_action_has_the_commutator_and_reproduces_hopf_data() {
    let mut source = type_iii().to_string();
    source.push_str(
        "
        let one=rvCycloRat[1](rvCyclo[4,0]);
        let minus=rvCycloRat[1](rvCyclo[4,2]);
        assert rvDoubleRepEntry[D,8,2,0,0]==one;
        assert rvDoubleRepEntry[D,8,2,1,1]==minus;
        assert rvDoubleRepEntry[D,8,2,0,1]==zero[RVCycloQ(4)];
        assert rvDoubleRepEntry[D,8,4,0,1]==one;
        assert rvDoubleRepEntry[D,8,4,1,0]==one;
        assert rvDoubleRepEntry[D,8,6,1,0]==one;
        assert rvDoubleRepEntry[D,8,6,0,1]==minus;
        assert add(rvDoubleRepEntry[D,8,2,0,0],rvDoubleRepEntry[D,8,2,1,1])
            ==rvDoubleCharacter[D,8,2];
        ",
    );
    for a in [0, 4, 8, 9, 10, 14, 20, 21] {
        for b in [0, 4, 8, 9, 10, 14, 20, 21] {
            writeln!(
                source,
                "assert rvCycloQGalois[-1](rvDoubleHopf[D,{a},{b}])==rvDoubleS[D,{a},{b}];"
            )
            .unwrap();
        }
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn colored_braid_relations_hold_as_operators_across_distinct_and_repeated_sectors() {
    let mut source = type_iii().to_string();
    for a in [0, 8, 9, 10, 14, 20] {
        for b in [0, 8, 9, 10, 14, 20] {
            for c in [0, 8, 9, 10, 14, 20] {
                writeln!(source,
                    "assert rvDoubleBraid[D,{a},{b},{c}]{{1,2,1}}==rvDoubleBraid[D,{a},{b},{c}]{{2,1,2}};
                     assert rvDoubleBraid[D,{a},{b},{c}]{{1,-1,2,-2}}==rvDoubleBraid[D,{a},{b},{c}]{{}};"
                ).unwrap();
            }
        }
    }
    source.push_str(
        r"
        let f=rvDoubleBraid[D,8,10,14]{1};
        let g=rvDoubleBraid[D,10,8,14]{2,-1};
        let composite=rvDoubleBraidCompose(f,g);
        assert composite==rvDoubleBraid[D,8,10,14]{1,2,-1};
        assert rvDoubleBraidInv(composite)==rvDoubleBraid[D,14,10,8]{1,-2,-1};
        let preserve=\x : One => f;
        assert preserve(One.pt)==f;
        show f;
        ",
    );
    assert!(stdout(interpret(&source)).starts_with("RVDoubleBraid(D,8,10,14,10,8,14){"));
}

#[test]
fn odd_prime_carry_has_the_nonreal_braiding_and_the_associator_sign() {
    let source = format!(
        "{CARRY}
        let phase=rvCycloRat[1](rvCyclo[9,1]);
        let inv=rvCycloRat[1](rvCyclo[9,8]);
        assert rvDoubleRepEntry[D,3,1,0,0]==phase;
        assert rvDoubleBraidTrace(rvDoubleBraid[D,3,3,0]{{1}})==phase;
        assert rvDoubleBraidTrace(rvDoubleBraid[D,3,3,0]{{-1}})==inv;
        assert rvDoubleAssociator[D,3,3,6]==rvCycloRat[1](rvCyclo[9,6]);
        assert rvDoubleHopf[D,3,3]==rvCycloRat[3](rvCyclo[9,2]);
        assert rvCycloQGalois[-1](rvDoubleHopf[D,3,3])==rvDoubleS[D,3,3];
        let shifted=rvPostnikovGauge(k,rvPostnikovCochain[P,2]{{0,1,0,0}});
        rvdouble Q[P](shifted);
        assert rvDoubleBraidTrace(rvDoubleBraid[Q,3,3,0]{{1}})==phase;"
    );
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn rephasing_preserves_borromean_trace_and_invalid_braids_are_rejected() {
    let b = list((0..49).map(|i| usize::from(i == 1)));
    let source = format!(
        "{}let shifted=rvPostnikovGauge(k,rvPostnikovCochain[P,2]{{{b}}});
        rvPostnikovCoeff[1,1,2](shifted);
        rvdouble Q[P](shifted);
        let word=rvDoubleBraid[Q,8,10,14]{{1,-2,1,-2,1,-2}};
        assert rvDoubleBraidTrace(word)==rvCycloRat[1](rvCycloScale[8](rvCyclo[4,2]));",
        type_iii()
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvPostnikovCoeff[1, 1, 2] = 1\n"
    );
    for (code, expected) in [
        ("show rvDoubleBraid[D,3,3,0]{3};","generators must be"),
        ("show rvDoubleBraidTrace(rvDoubleBraid[D,3,4,0]{1});","trace requires matching source"),
        ("let f=rvDoubleBraid[D,3,4,0]{1};show rvDoubleBraidCompose(f,f);","requires identical middle color types"),
        ("show rvDoubleRepEntry[D,3,1,1,0];","matrix element or basis index"),
        ("rvDoubleBraidImage[1](rvDoubleBraid[D,3,3,0]{});","basis index is out of range"),
        (r"let f=\x : RVDoubleBraid(D,3,4,0,3,4,1) => x;","target must permute its source"),
        ("rvdouble Q[P](k);show rvDoubleBraidCompose(rvDoubleBraid[D,3,3,0]{},rvDoubleBraid[Q,3,3,0]{});","requires identical middle color types"),
    ] {
        let error = diagnostic(interpret(&format!("{CARRY}{code}")));
        assert!(error.contains(expected), "{expected} in {error}");
    }
    let word = list(std::iter::repeat_n(1, 4097));
    let error = diagnostic(interpret(&format!(
        "{CARRY}show rvDoubleBraid[D,3,3,0]{{{word}}};"
    )));
    assert!(error.contains("word exceeds its length limit"));
}
