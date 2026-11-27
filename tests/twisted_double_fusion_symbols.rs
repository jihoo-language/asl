mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/twisted-double-fusion-associators.asl");
const CARRY: &str = "
    rvgroup G[3];rvdwcoeff P[G,3];
    let k=rvPostnikovCochain[P,3]{0,1,1,1,0,2,2,2};rvdouble D[P](k);
";

fn type_iii() -> &'static str {
    EXAMPLE.split_once("// Simple 8").unwrap().0
}

fn list(values: impl IntoIterator<Item = usize>) -> String {
    values
        .into_iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

#[test]
fn type_iii_f_matrix_and_r_symbols_are_exact_on_the_four_electric_channels() {
    let output = stdout(example("twisted-double-fusion-associators"));
    assert_eq!(
        output,
        concat!(
            "rvDoubleFusionSpace[D,8,8,0] = multiplicity=1, rows=4, columns=1\n",
            "rvDoubleFData[D,8,8,8,8] = dimension=4, left=[(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)], right=[(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)]\n",
            "RVCycloQ(4){[1] / 2}\n",
            "RVCycloQ(4){[1] / 1}\n",
            "RVCycloQ(4){[-1] / 1}\n",
            "rvDoublePentagon[D,8,8,8,8,0] = dimension=4, verified=true\n",
            "rvDoubleHexagon[D,8,8,8,8] = dimension=4, verified=true\n",
        )
    );

    let mut source = type_iii().to_owned();
    source.push_str(
        "let one=rvCycloRat[1](rvCyclo[4,0]);let zeroq=zero[RVCycloQ(4)];let quarter=rvCycloRat[4](rvCyclo[4,0]);
         assert rvDoubleVertex[D,8,8,0,0,0,0]==one;
         assert rvDoubleVertex[D,8,8,0,0,3,0]==one;
         assert rvDoubleVertex[D,8,8,0,0,1,0]==zeroq;",
    );
    for row in 0..4 {
        for col in 0..4 {
            writeln!(
                source,
                "let m{row}{col}=rvDoubleF[D,8,8,8,8,{row},{col}];
                assert rvCycloQMul(m{row}{col},rvCycloQGalois[-1](m{row}{col}))==quarter;"
            )
            .unwrap();
            let entry = (0..4)
                .map(|i| {
                    format!(
                    "rvCycloQMul(rvDoubleFInv[D,8,8,8,8,{row},{i}],rvDoubleF[D,8,8,8,8,{i},{col}])"
                )
                })
                .reduce(|a, b| format!("add({a},{b})"))
                .unwrap();
            let expected = if row == col { "one" } else { "zeroq" };
            writeln!(source, "assert {entry}=={expected};").unwrap();
        }
    }
    for i in 0..4 {
        writeln!(
            source,
            "assert rvDoubleR[D,8,8,{i},0,0]=={};",
            if i == 3 { "neg(one)" } else { "one" }
        )
        .unwrap();
    }
    source.push_str("rvDoubleHexagon[D,8,10,14,20];");
    assert_eq!(
        stdout(interpret(&source)),
        "rvDoubleHexagon[D,8,10,14,20] = dimension=2, verified=true\n"
    );
}

#[test]
fn odd_carry_uses_the_dpr_sign_for_both_associator_and_crossing() {
    let source = format!(
        "{CARRY}
         assert rvDoubleF[D,3,3,6,5,0,0]==rvCycloRat[1](rvCyclo[9,6]);
         assert rvDoubleFInv[D,3,3,6,5,0,0]==rvCycloRat[1](rvCyclo[9,3]);
         assert rvDoubleR[D,3,3,6,0,0]==rvCycloRat[1](rvCyclo[9,1]);
         assert rvDoubleRInv[D,3,3,6,0,0]==rvCycloRat[1](rvCyclo[9,8]);
         rvDoublePentagon[D,3,3,3,3,5];rvDoubleHexagon[D,3,3,6,5];"
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvDoublePentagon[D,3,3,3,3,5] = dimension=1, verified=true\n",
            "rvDoubleHexagon[D,3,3,6,5] = dimension=1, verified=true\n",
        )
    );
}

fn c3cube(wrong_cocycle: bool) -> String {
    let mut source = "rvgroup G{".to_owned();
    for a in 0..27 {
        writeln!(
            source,
            "{{{}}},",
            list((0..27).map(|b| {
                (a % 3 + b % 3) % 3 + 3 * ((a / 3 % 3 + b / 3 % 3) % 3) + 9 * ((a / 9 + b / 9) % 3)
            }))
        )
        .unwrap();
    }
    let k = list((1..27).flat_map(|a| {
        (1..27).flat_map(move |b| {
            (1..27).map(move |c| {
                (a % 3 * (b / 3 % 3) * (c / 9)
                    + usize::from(wrong_cocycle && [a, b, c] == [1, 1, 1]))
                    % 3
            })
        })
    }));
    write!(
        source,
        "}};rvdwcoeff P[G,3];let k=rvPostnikovCochain[P,3]{{{k}}};rvdouble D[P](k);"
    )
    .unwrap();
    source
}

#[test]
fn genuine_multiplicity_three_channels_and_pentagon_on_dimension_twenty_seven() {
    let source = format!(
        "{}
        rvDoubleFusionSpace[D,27,27,30];rvDoubleFusionSpace[D,27,27,31];
        rvDoubleFData[D,27,27,27,0];
        let one=rvCycloRat[1](rvCyclo[9,0]);
        assert rvDoubleR[D,27,27,30,0,0]==one;
        assert rvDoubleR[D,27,27,30,1,2]==one;
        assert rvDoubleR[D,27,27,30,2,1]==one;
        assert rvDoubleR[D,27,27,30,1,1]==zero[RVCycloQ(9)];
        assert rvDoubleRInv[D,27,27,30,1,2]==one;
        rvDoubleHexagon[D,27,27,27,0];rvDoublePentagon[D,27,27,27,27,27];",
        c3cube(false)
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvDoubleFusionSpace[D,27,27,30] = multiplicity=3, rows=9, columns=3\n",
            "rvDoubleFusionSpace[D,27,27,31] = multiplicity=0, rows=9, columns=3\n",
            "rvDoubleFData[D,27,27,27,0] = dimension=3, left=[(30, 0, 0), (30, 1, 0), (30, 2, 0)], right=[(30, 0, 0), (30, 1, 0), (30, 2, 0)]\n",
            "rvDoubleHexagon[D,27,27,27,0] = dimension=3, verified=true\n",
            "rvDoublePentagon[D,27,27,27,27,27] = dimension=27, verified=true\n",
        )
    );
    let error = diagnostic(interpret(&c3cube(true)));
    assert!(error.contains("fails the twisted cocycle"), "{error}");
}

#[test]
fn coherence_survives_a_gauge_change_and_checks_empty_spaces_and_indices() {
    let b = list((0..49).map(|i| usize::from(i == 1)));
    let source = format!(
        "{}let shifted=rvPostnikovGauge(k,rvPostnikovCochain[P,2]{{{b}}});rvdouble Q[P](shifted);
        rvDoublePentagon[Q,8,8,8,8,0];rvDoubleHexagon[Q,8,10,14,20];
        rvDoubleFData[Q,8,8,8,0];",
        type_iii()
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvDoublePentagon[Q,8,8,8,8,0] = dimension=4, verified=true\n",
            "rvDoubleHexagon[Q,8,10,14,20] = dimension=2, verified=true\n",
            "rvDoubleFData[Q,8,8,8,0] = dimension=0, left=[], right=[]\n",
        )
    );
    for (operation, expected) in [
        (
            "show rvDoubleVertex[D,3,3,6,1,0,0];",
            "trivalent channel index is out of range",
        ),
        (
            "show rvDoubleF[D,3,3,6,5,0,1];",
            "fusion matrix index is out of range",
        ),
        (
            "show rvDoubleR[D,3,3,6,1,0];",
            "fusion matrix index is out of range",
        ),
        (
            "rvDoublePentagon[D,3,3,3,3,9];",
            "simple index is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{CARRY}{operation}")));
        assert!(error.contains(expected), "{expected} in {error}");
    }
}
