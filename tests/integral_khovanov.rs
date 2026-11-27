mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn the_trefoil_has_a_certified_integral_order_two_generator() {
    let output = stdout(example("integral-khovanov"));
    assert!(output.contains("khGroup[Trefoil, 3, 7] = chain=3, free=0, torsion=[2]"));
    assert!(output.contains("khGroup[Mirror, -2, -7] = "));
    assert!(output.contains("KhH(Trefoil, 3, 7){g0: 1}"));
    assert!(output.contains("khEuler[Trefoil] = [(1, 1), (3, 1), (5, 1), (9, -1)]"));
    assert!(output.contains("khEuler[FigureEight] = [(-5, 1), (5, 1)]"));
}

#[test]
fn both_markov_stabilizations_preserve_all_integral_groups() {
    let source = r"
        khbraid T[2,[1,1,1]];
        khbraid Plus[3,[1,1,1,2]];
        khbraid Minus[3,[1,1,1,-2]];
    ";
    let mut source = source.to_string();
    for i in -1..=4 {
        for j in -5..=12 {
            for braid in ["T", "Plus", "Minus"] {
                source.push_str(&format!("khGroup[{braid},{i},{j}];\n"));
            }
        }
    }
    let output = stdout(interpret(&source));
    let lines: Vec<_> = output.lines().collect();
    let (slices, remainder) = lines.as_chunks::<3>();
    assert!(remainder.is_empty());
    for slice in slices {
        let groups = slice.map(|line| line.split_once(", free=").unwrap().1);
        assert_eq!(groups[0], groups[1], "{slice:?}");
        assert_eq!(groups[0], groups[2], "{slice:?}");
    }
}

#[test]
fn conjugating_a_mixed_sign_braid_preserves_all_integral_groups() {
    let mut source = String::from("khbraid F[3,[1,-2,1,-2]]; khbraid Conj[3,[2,1,-2,1,-2,-2]];\n");
    for i in -3..=3 {
        for j in (-7..=7).step_by(2) {
            source.push_str(&format!("khGroup[F,{i},{j}]; khGroup[Conj,{i},{j}];\n"));
        }
    }
    source.push_str("assert khEulerAt[F,3] == khEulerAt[Conj,3];");
    let output = stdout(interpret(&source));
    let lines: Vec<_> = output.lines().collect();
    let (pairs, remainder) = lines.as_chunks::<2>();
    assert!(remainder.is_empty());
    for [left, right] in pairs {
        assert_eq!(
            left.split_once(", free=").unwrap().1,
            right.split_once(", free=").unwrap().1
        );
    }
}

#[test]
fn mixed_signs_and_mirroring_give_the_expected_bigraded_torsion() {
    let output = stdout(interpret(
        r"
        khbraid F[3,[1,-2,1,-2]];
        khbraid L[2,[-1,-1,-1]];
        khGroup[F,-1,-3];
        khGroup[F,2,3];
        khGroup[F,0,1];
        khGroup[L,-2,-7];
        let l = khGen[L,-2,-7,0];
        assert khD(khNull[L,-2,-7,0]) == khScale[2](khRep(l));
        let a = khGen[F,-1,-3,0];
        let b = khGen[F,2,3,0];
        assert khD(khNull[F,-1,-3,0]) == khScale[2](khRep(a));
        assert khD(khNull[F,2,3,0]) == khScale[2](khRep(b));
        let z = khTerm[F,[1],[0]];
        show khD(z);
        assert khD(khD(z)) == zero[KhC(F,1,-1)];
        ",
    ));
    assert!(output
        .lines()
        .next()
        .unwrap()
        .ends_with("free=0, torsion=[2]"));
    assert!(output
        .lines()
        .nth(1)
        .unwrap()
        .ends_with("free=0, torsion=[2]"));
    assert!(output
        .lines()
        .nth(2)
        .unwrap()
        .ends_with("free=1, torsion=[]"));
    assert!(output
        .lines()
        .nth(3)
        .unwrap()
        .ends_with("free=0, torsion=[2]"));
    assert!(output.lines().last().unwrap().contains(": -1"));
}

#[test]
fn class_arithmetic_is_integral_or_reduced_by_the_actual_smith_factor() {
    let source = r"
        khbraid T[2,[1,1,1]];
        let a = khGen[T,3,7,0];
        let b = khGen[T,0,1,0];
        assert neg(a) == a;
        assert add(a,a) == zero[KhH(T,3,7)];
        assert khScale[-170141183460469231731687303715884105728](a) == zero[KhH(T,3,7)];
        assert khClass(neg(khRep(b))) == neg(b);
        assert khClass(khScale[17](khRep(b))) == khScale[17](b);
        let on_cycle = \x : KhH(T,3,7) => S.hi;
        card on_cycle;
        let known_empty = \x : KhC(T,2,100) => S.hi;
        card known_empty;
    ";
    assert_eq!(stdout(interpret(source)), "card = 2 / 2\ncard = 1 / 1\n");
}

#[test]
fn empty_slices_and_unlinks_have_the_expected_euler_normalization() {
    let word = std::iter::repeat_n("1", 128).collect::<Vec<_>>().join(",");
    let output = stdout(interpret(&format!(
        r"
        khbraid Long[2,[{word}]];
        khGroup[Long,64,9223372036854775807];
        khGroup[Long,64,-9223372036854775808];
        khGroup[Long,64,101];
        khbraid Empty[0,[]];
        khbraid One[1,[]];
        khbraid Two[2,[]];
        khEuler[Empty];
        khEuler[One];
        khEuler[Two];
        let delta = neg(rvCycloRat[1](add(rvCyclo[20,2],rvCyclo[20,-2])));
        assert khEulerAt[Empty,3] == rvCycloQMul(delta,rvTLJones[3,0,[]]);
        assert khEulerAt[One,3] == delta;
        "
    )));
    assert!(output
        .lines()
        .take(3)
        .all(|line| line.ends_with("chain=0, free=0, torsion=[]")));
    assert!(output.contains("khEuler[Empty] = [(0, 1)]"));
    assert!(output.contains("khEuler[One] = [(-1, 1), (1, 1)]"));
    assert!(output.contains("khEuler[Two] = [(-2, 1), (0, 2), (2, 1)]"));
}

#[test]
fn invalid_input_and_overflow_have_diagnostics() {
    let cases = [
        ("khbraid T[2,[2]];", "between 1 and strands minus one"),
        (
            "khbraid T[2,[1]]; show khTerm[T,[1],[]];",
            "crossing index is out of range",
        ),
        (
            "khbraid T[2,[1]]; show khTerm[T,[0,0],[]];",
            "repeated Khovanov crossing index",
        ),
        (
            "khbraid T[2,[]]; show khTerm[T,[],[0,0]];",
            "strictly increasing",
        ),
        (
            "khbraid T[2,[]]; show khTerm[T,[],[2]];",
            "circle index is out of range",
        ),
        (
            "khbraid T[2,[1,1,1]]; show khClass(khTerm[T,[],[]]);",
            "not a cycle",
        ),
        (
            "khbraid T[2,[1,1,1]]; show khNull[T,0,1,0];",
            "free Khovanov generator",
        ),
        (
            "khbraid T[1,[]]; show khD(zero[KhC(T,9223372036854775807,1)]);",
            "degree is out of range",
        ),
        (
            "khbraid T[1,[]]; show neg(khTerm[T,[],[],-170141183460469231731687303715884105728]);",
            "coefficient range",
        ),
        (
            r"khbraid T[1,[]]; let a = \x : KhH(T,0,1) => S.hi; card a;",
            "free summand",
        ),
    ];
    for (source, expected) in cases {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{source} received {error}");
    }
}
