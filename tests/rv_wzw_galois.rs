mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn signed_alcove_example_checks_rows_twists_and_complex_coefficients() {
    assert_eq!(
        stdout(example("rv-wzw-galois")),
        concat!(
            "RVWzw(3, 2){[0, 1, 2]: [1]}\n",
            "RVWzw(3, 2){[0, 0, 0]: [-1]}\n",
            "RVCyclo(60){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}\n",
        )
    );
}

#[test]
fn all_six_basis_labels_satisfy_signed_composition_and_s_row_transport() {
    let mut source = String::new();
    for (i, weight) in [
        [0, 0, 0],
        [0, 0, 1],
        [0, 0, 2],
        [0, 1, 1],
        [0, 1, 2],
        [0, 2, 2],
    ]
    .into_iter()
    .enumerate()
    {
        writeln!(
            source,
            "
            let v{i} = rvWzwBasis[3,2,{weight:?}];
            assert rvWzwGalois[7](rvWzwGalois[7](v{i})) == neg(v{i});
            assert rvWzwGalois[7](rvWzwGalois[11](v{i})) == rvWzwGalois[17](v{i});
            assert rvWzwGalois[-1](v{i}) == neg(rvWzwCharge(v{i}));
            assert rvWzwCoeffGalois[7](rvWzwS(v{i})) == rvWzwS(rvWzwGalois[7](v{i}));
            "
        )
        .unwrap();
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn powers_are_reduced_and_sparse_labels_need_no_alcove_enumeration() {
    let program = r"
        let x = add(rvCyclo[60,1],neg(rvCyclo[60,0]));
        let y = add(rvCyclo[60,-3],rvCyclo[60,1]);
        assert rvCycloGalois[7](rvCycloMul(x,y)) == rvCycloMul(rvCycloGalois[7](x),rvCycloGalois[7](y));
        assert rvCycloGalois[7](rvCycloGalois[11](x)) == rvCycloGalois[17](x);
        assert rvCycloGalois[-1](x) == rvCycloConj(x);
        assert rvCycloGalois[170141183460469231731687303715884105727](x) == rvCycloGalois[7](x);
        let v = rvWzwBasis[3,2,[0,0,0]];
        assert rvWzwGalois[-170141183460469231731687303715884105713](v) == rvWzwGalois[7](v);
        assert rvWzwGalois[7](rvWzwBasis[3,100,[0,0,0]]) == rvWzwBasis[3,100,[0,6,12]];
        assert rvWzwGalois[7](zero[RVWzw(3,2)]) == zero[RVWzw(3,2)];
        assert rvWzwGalois[3](rvWzwBasis[2,0,[4,4]]) == neg(rvWzwBasis[2,0,[0,0]]);
        assert rvWzwGalois[-1](rvWzwBasis[4,1,[0,0,0,1]]) == rvWzwCharge(rvWzwBasis[4,1,[0,0,0,1]]);
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn invalid_galois_units_types_and_exponents_are_diagnosed() {
    for (program, expected) in [
        ("show rvCycloGalois[2](rvCyclo[60,0]);", "must be coprime"),
        ("show rvWzwGalois[0](zero[RVWzw(3,2)]);", "must be coprime"),
        (
            "show rvWzwCoeffGalois[5](zero[RVWzw(3,2)]);",
            "must be coprime",
        ),
        (
            "show rvCycloGalois[7](S.hi);",
            "requires a cyclotomic integer",
        ),
        (
            "show rvWzwGalois[7](rvCyclo[60,0]);",
            "requires a WZW torus state",
        ),
        (
            "show rvWzwGalois[170141183460469231731687303715884105728](zero[RVWzw(3,2)]);",
            "exceeds the signed i128 range",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
