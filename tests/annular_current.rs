mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn integral_currents_homotopies_and_highest_weight_strings() {
    let output = stdout(example("annular-current"));
    for line in [
        "khIrreps[Five, 0, 0] = [(1, 5), (3, 4), (5, 1)]",
        "KhAH(Five, 0, 3, 3){g0: -1, g4: 1}",
        "khActionCheck[Five] = verified (32 enhanced states)",
        "KhAH(Crossing, 1, 3, 0){g0: 1}",
        "khAnnularGroup[Twisted, 3, 8, 1] = chain=3, free=0, torsion=[2]",
        "khActionCheck[Figure] = verified (66 enhanced states)",
    ] {
        assert!(
            output.lines().any(|actual| actual == line),
            "{line}\n{output}"
        );
    }
}

#[test]
fn shifts_and_cocycle_maps_are_nominal_and_bounded() {
    let output = stdout(interpret(
        r"
        khbraid Empty[0,[]];
        khbraid Mixed[3,[2,1,-2]];
        khbraid Conjugate[3,[1,-2,2]];
        khIrreps[Empty,0,0];
        khIrreps[Mixed,0,1];
        khIrreps[Conjugate,0,1];
        khIrreps[Mixed,0,3];
        khIrreps[Conjugate,0,3];
        khIrreps[Mixed,0,9223372036854775807];
        let x = khAnnular[1](khTerm[Mixed,[],[1]]);
        let z = khHighest[Empty,0,0,0,0];
        assert khClass(z) == khAGen[Empty,0,0,0,0];
        assert khVMinus(x) == khDMinus2(x);
        assert khH(x) == x;
        assert khD0(khVZero(x)) == neg(khVZero(khD0(x)));
        khActionCheck[Mixed];
    ",
    ));
    assert!(
        output.contains("khIrreps[Empty, 0, 0] = [(0, 1)]"),
        "{output}"
    );
    assert!(
        output.contains("khIrreps[Mixed, 0, 9223372036854775807] = []"),
        "{output}"
    );
    assert!(
        output.contains("khIrreps[Mixed, 0, 1] = [(1, 1), (3, 1)]"),
        "{output}"
    );
    for difference in [1, 3] {
        let multiplicities = |name| {
            output
                .lines()
                .find_map(|line| {
                    line.strip_prefix(&format!("khIrreps[{name}, 0, {difference}] = "))
                })
                .unwrap()
        };
        assert_eq!(multiplicities("Mixed"), multiplicities("Conjugate"));
    }
}

#[test]
fn highest_weight_lift_corrects_a_smith_representative() {
    let output = stdout(interpret(
        r"
        khbraid K[4,[1,-2,3,-1]];
        khIrreps[K,0,0];
        let z = khHighest[K,0,0,0,0];
        let integralClass = khClass(z);
        assert khD0(z) == zero[KhAC(K,1,0,0)];
        assert khE(z) == zero[KhAC(K,0,2,2)];
        assert khF(z) == zero[KhAC(K,0,-2,-2)];
        assert khE(integralClass) == zero[KhAH(K,0,2,2)];
        show khE(khRep(integralClass));
    ",
    ));
    assert!(
        output.contains("khIrreps[K, 0, 0] = [(0, 1), (2, 2), (4, 1)]"),
        "{output}"
    );
    assert!(output.contains("KhAC(K, 0, 2, 2){([0, 1], [0]): 2, ([1, 3], [0]): -2, ([1, 3], [2]): 2, ([2, 3], []): 2}"), "{output}");
}

#[test]
fn invalid_current_operands_and_highest_weights_are_diagnosed() {
    for (input, expected) in [
        (
            "show khE(khTerm[U,[],[]]);",
            "requires an annular Khovanov chain or class",
        ),
        (
            "show khLee0(khAGen[U,0,1,1,0]);",
            "requires an annular Khovanov chain",
        ),
        (
            "show khHighest[U,0,0,-1,0];",
            "highest weight must be nonnegative",
        ),
        (
            "show khHighest[U,0,0,1,1];",
            "highest-weight generator index is out of range",
        ),
        (
            "show khE(zero[KhAC(U,0,9223372036854775807,1)]);",
            "quantum degree is out of range",
        ),
        (
            "khbraid L[120,[]]; khActionCheck[L];",
            "enhanced cube exceeds the work limit",
        ),
    ] {
        let error = diagnostic(interpret(&format!("khbraid U[1,[]]; {input}")));
        assert!(error.contains(expected), "{input}\n{error}");
    }
}
