mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn setup() -> &'static str {
    include_str!("../examples/hyperelliptic-ks-trace-projector.asl")
        .split_once("// Trace projector")
        .unwrap()
        .0
}

const SHEAR: &str = "let h2=add(e,add(h,h)); rvog6cusp Q(h2,add(e,f),f);";

#[test]
fn example_detects_a_fraction_and_checks_spin_equivariance() {
    assert_eq!(
        stdout(example("hyperelliptic-ks-trace-projector")),
        concat!(
            "rvKSVolumeSquare[P] = 4\n",
            "rvKSEndTracePair = 128\n",
            "rvKSEndCoeff[0,50] = 1/32\n",
            "rvKSEndTracePair = 0\n",
        ),
    );
}

fn integer_report(line: &str) -> i128 {
    line.split_once(" = ").unwrap().1.parse().unwrap()
}

fn determinant(gram: &[[i128; 7]; 7], row: usize, used: u8) -> i128 {
    if row == 7 {
        return 1;
    }
    (0..7)
        .filter(|&column| used & (1 << column) == 0)
        .map(|column| {
            let inversions = ((used >> column) & !1).count_ones();
            let sign = if inversions.is_multiple_of(2) { 1 } else { -1 };
            sign * gram[row][column] * determinant(gram, row + 1, used | (1 << column))
        })
        .sum()
}

#[test]
fn oriented_volume_and_trace_form_match_the_actual_gram() {
    let mut source = format!("{}{SHEAR}", setup());
    for mark in ["P", "Q"] {
        for i in 0..7 {
            writeln!(source, "let x{mark}{i}=rvOG6PolarizedBasis[{mark},{i}];").unwrap();
        }
    }
    let mut reports = source.clone();
    for mark in ["P", "Q"] {
        for i in 0..7 {
            for j in 0..7 {
                writeln!(reports, "rvOG6Pair(x{mark}{i},x{mark}{j});").unwrap();
            }
        }
    }
    let output = stdout(interpret(&reports));
    let mut gram_entries = output.lines().map(integer_report);
    let mut expected = Vec::new();
    for mark in ["P", "Q"] {
        let gram = std::array::from_fn(|_| std::array::from_fn(|_| gram_entries.next().unwrap()));
        assert!(gram
            .iter()
            .enumerate()
            .any(|(i, row)| { row.iter().enumerate().any(|(j, v)| i != j && v % 2 != 0) }));
        let square = -determinant(&gram, 0, 0) / 2;
        writeln!(source, "let v{mark}=rvKSVolume[{mark}]; rvKSVolumeSquare[{mark}]; rvOG6CliffordCoeff[127](v{mark});").unwrap();
        writeln!(source, "assert rvOG6CliffordMul(v{mark},v{mark})==rvOG6CliffordScale[{square}](rvOG6CliffordBasis[{mark},0]);").unwrap();
        expected.extend([square, 8]);
        for i in 0..7 {
            writeln!(source, "let z{mark}{i}=rvOG6CliffordVector[{mark}](x{mark}{i}); let a{mark}{i}=rvKSEmbed[{mark}](x{mark}{i});").unwrap();
            writeln!(source, "assert rvOG6CliffordMul(v{mark},z{mark}{i})==rvOG6CliffordMul(z{mark}{i},v{mark});").unwrap();
            writeln!(
                source,
                "assert rvKSEndLeft(rvOG6CliffordMul(v{mark},z{mark}{i}))==a{mark}{i};"
            )
            .unwrap();
        }
        for (i, row) in gram.iter().enumerate() {
            for (j, &g) in row.iter().enumerate().skip(i) {
                writeln!(source, "rvKSEndTracePair(a{mark}{i},a{mark}{j}); rvKSEndTrace(rvKSEndCompose(a{mark}{j},a{mark}{i}));").unwrap();
                expected.extend([32 * square * g; 2]);
            }
        }
    }
    assert!(gram_entries.next().is_none());
    let actual = stdout(interpret(&source));
    assert_eq!(
        actual.lines().map(integer_report).collect::<Vec<_>>(),
        expected
    );
}

#[test]
fn projector_is_trace_orthogonal_on_arbitrary_matrices_in_both_markings() {
    let mut source = format!("{}{SHEAR}", setup());
    for (mark, col) in [("P", 50), ("Q", 49)] {
        writeln!(source, r"
        let unit{mark}=rvKSEndUnit[{mark},0,{col}];
        let proj{mark}=rvKSProject(unit{mark});
        assert rvKSProject(proj{mark})==proj{mark};
        assert rvKSProject(rvKSEndLeft(rvOG6CliffordBasis[{mark},0]))==zero[RVKSEnd({mark})];
        let mixed{mark}=add(rvKSEndScale[3,5](unit{mark}),add(rvKSEndScale[-11,13](rvKSEndUnit[{mark},57,14]),rvKSEndUnit[{mark},63,63]));
        let projected{mark}=rvKSProject(mixed{mark});
        assert rvKSProject(projected{mark})==projected{mark};
        assert projected{mark}==add(rvKSEndScale[3,5](proj{mark}), add(rvKSEndScale[-11,13](rvKSProject(rvKSEndUnit[{mark},57,14])), rvKSProject(rvKSEndUnit[{mark},63,63])));
        let residual{mark}=add(mixed{mark},neg(projected{mark}));
        let other{mark}=add(rvKSEndScale[-23,7](rvKSEndUnit[{mark},{col},0]),rvKSEndUnit[{mark},2,53]);
        rvKSEndTracePair(projected{mark},other{mark});
        rvKSEndTracePair(mixed{mark},rvKSProject(other{mark}));
        ").unwrap();
        for i in 0..7 {
            writeln!(
                source,
                "let b{mark}{i}=rvKSEmbed[{mark}](rvOG6PolarizedBasis[{mark},{i}]);"
            )
            .unwrap();
            writeln!(source, "assert rvKSProject(b{mark}{i})==b{mark}{i};").unwrap();
            writeln!(source, "rvKSEndTracePair(residual{mark},b{mark}{i});").unwrap();
        }
    }
    source.push_str("rvKSEndCoeff[0,50](projP); rvKSEndCoeff[9,49](projQ);");
    let output = stdout(interpret(&source));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines.len(), 20);
    for group in lines[..18].chunks_exact(9) {
        assert_eq!(group[0], group[1]);
        assert!(group[2..]
            .iter()
            .all(|line| *line == "rvKSEndTracePair = 0"));
    }
    assert_eq!(
        &lines[18..],
        ["rvKSEndCoeff[0,50] = 1/32", "rvKSEndCoeff[9,49] = 1/64"]
    );
}

#[test]
fn spin_monodromy_preserves_the_embedding_and_the_trace_complement() {
    let mut source = format!(
        r"{}{SHEAR}
        rvog6nilpotentii Q2[Q,-3];
        let qeDual=add(add(eDual,eDual),neg(hF));
        rvog6nilpotentiii Q3[Q,1](add(e,qeDual));
        let T=add(rvKSEndUnit[Q,0,49],rvKSEndScale[2,7](rvKSEndUnit[Q,16,5]));
        let PT=rvKSProject(T);
        ",
        setup()
    );
    for name in ["Q2", "Q3"] {
        writeln!(source, r"
        let s{name}=rvKSEndLeft(rvOG6SpinLog[{name}]);
        let lt{name}=rvKSEndLeft(rvOG6SpinExp[{name},7]);
        let rt{name}=rvKSEndLeft(rvOG6SpinExp[{name},-7]);
        let comm{name}=add(rvKSEndCompose(s{name},T),neg(rvKSEndCompose(T,s{name})));
        assert rvKSProject(comm{name})==add(rvKSEndCompose(s{name},PT),neg(rvKSEndCompose(PT,s{name})));
        ").unwrap();
        for i in 0..7 {
            writeln!(source, r"
            let v{name}{i}=rvOG6PolarizedBasis[Q,{i}];
            let j{name}{i}=rvKSEmbed[Q](v{name}{i});
            let b{name}{i}=add(rvKSEndCompose(s{name},j{name}{i}),neg(rvKSEndCompose(j{name}{i},s{name})));
            assert b{name}{i}==rvKSEmbed[Q](rvOG6NilpotentPower[{name},1](v{name}{i}));
            let c{name}{i}=rvKSEndCompose(rvKSEndCompose(lt{name},j{name}{i}),rt{name});
            assert c{name}{i}==rvKSEmbed[Q](rvOG6NilpotentExp[{name},7](v{name}{i}));
            ").unwrap();
        }
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn complete_matrix_units_and_rational_arithmetic_follow_composition_order() {
    let source = format!(
        r"{}
        let a=rvKSEndUnit[P,63,3]; let b=rvKSEndUnit[P,3,31]; let c=rvKSEndUnit[P,31,63];
        let I=rvKSEndLeft(rvOG6CliffordBasis[P,0]);
        assert rvKSEndCompose(I,a)==a;
        assert rvKSEndCompose(a,I)==a;
        assert rvKSEndCompose(a,b)==rvKSEndUnit[P,63,31];
        assert rvKSEndCompose(b,a)==zero[RVKSEnd(P)];
        assert rvKSEndCompose(rvKSEndCompose(a,b),c)==rvKSEndCompose(a,rvKSEndCompose(b,c));
        let large=rvKSEndScale[-170141183460469231731687303715884105728,340282366920938463463374607431768211455](a);
        let negLarge=neg(large);
        assert add(large,negLarge)==zero[RVKSEnd(P)];
        let product=rvKSEndCompose(rvKSEndScale[-3,5](a),rvKSEndScale[10,21](b));
        assert product==rvKSEndScale[-2,7](rvKSEndUnit[P,63,31]);
        rvKSEndCoeff[63,3](large);
        rvKSEndTrace(I);
        rvKSEndTracePair(product,c);
        rvKSEndTrace(rvKSEndCompose(c,product));
        ",
        setup()
    );
    assert_eq!(stdout(interpret(&source)), concat!(
        "rvKSEndCoeff[63,3] = -170141183460469231731687303715884105728/340282366920938463463374607431768211455\n",
        "rvKSEndTrace = 64\n", "rvKSEndTracePair = -2/7\n", "rvKSEndTrace = -2/7\n",
    ));
}

#[test]
fn wrong_marks_nonpolarized_vectors_and_out_of_range_indices_are_diagnosed() {
    for (source, expected) in [
        ("show rvKSEndUnit[P,0,64];", "between 0 and 63"),
        (
            "rvKSEndCoeff[64,0](rvKSEndUnit[P,0,0]);",
            "between 0 and 63",
        ),
        ("show rvKSEmbed[P](h);", "requires a polarized vector"),
        (
            "show rvKSEndLeft(rvKSVolume[P]);",
            "expected type RVOG6Clifford(P,0)",
        ),
        (
            "show rvKSEndScale[5,0](rvKSEndUnit[P,0,0]);",
            "denominator must be positive",
        ),
        (
            "show rvKSProject(rvOG6SpinLog[A]);",
            "an RVKSEnd endomorphism",
        ),
        (
            "rvog6cusp Q(h,e,f); show rvKSEndCompose(rvKSEndUnit[P,0,0],rvKSEndUnit[Q,0,0]);",
            "expected type RVKSEnd(P)",
        ),
    ] {
        let actual = diagnostic(interpret(&format!("{}{source}", setup())));
        assert!(actual.contains(expected), "{actual}");
    }
}
