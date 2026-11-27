mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn level_three_polarized_torsion_agrees_with_different_newton_polygons() {
    assert_eq!(
        stdout(example("rank-dieudonne-crystal")),
        concat!(
            "EOType(7){0, 1, 2, 3, 4, 5, 6}\n",
            "rvCrystalChar[Base] = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2187]\n",
            "rvCrystalChar[Tilted] = [1, 0, 0, 0, 0, 0, 0, -27, 0, 0, 0, 0, 0, 0, 2187]\n",
            "rvCrystalNewton[Base] = [1/2: 14]\n",
            "rvCrystalNewton[Tilted] = [3/7: 7, 4/7: 7]\n",
        )
    );
}

#[test]
fn conjugation_preserves_polynomial_and_final_type_with_signed_coefficients() {
    let input = r"
        rvcrystal M[5,2] {(0,1,5), (0,3,-5), (1,0,1), (2,1,1), (3,2,5)};
        rvcrystal Q[5,2] {(0,1,3), (0,2,2), (0,3,4), (1,0,1), (1,1,1), (1,2,-14), (1,3,2), (2,1,1), (2,2,-1), (2,3,3), (3,2,5)};
        rvCrystalChar[M];
        rvCrystalChar[Q];
        rvCrystalNewton[M];
        rvCrystalNewton[Q];
        assert rvCrystalEOType[M] == eoType{0,1};
        assert rvCrystalEOType[Q] == rvCrystalEOType[M];
        assert rvCrystalCongruent[M,Q,1] == S.lo;
    ";
    assert_eq!(
        stdout(interpret(input)),
        concat!(
            "rvCrystalChar[M] = [1, 0, -5, 0, 25]\n",
            "rvCrystalChar[Q] = [1, 0, -5, 0, 25]\n",
            "rvCrystalNewton[M] = [1/2: 4]\n",
            "rvCrystalNewton[Q] = [1/2: 4]\n",
        )
    );
}

#[test]
fn four_genus_eight_polygons_share_an_identical_bt_one() {
    let mut source = String::new();
    for (level, coefficient) in [3, 9, 27, 81].into_iter().enumerate() {
        write!(&mut source, "rvcrystal A{level}[3,8] {{").unwrap();
        for i in 0..7 {
            write!(&mut source, "({},{i},1), ({},{},3),", i + 1, i + 9, i + 8).unwrap();
        }
        writeln!(&mut source, "(8,7,1), (0,15,-3), (0,7,{coefficient})}};").unwrap();
        writeln!(
            &mut source,
            "assert rvCrystalCongruent[A0,A{level},1] == S.hi;"
        )
        .unwrap();
        writeln!(
            &mut source,
            "assert rvCrystalEOType[A{level}] == eoType{{0,1,2,3,4,5,6,7}};"
        )
        .unwrap();
        writeln!(&mut source, "rvCrystalNewton[A{level}];").unwrap();
    }
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvCrystalNewton[A0] = [1/8: 8, 7/8: 8]\n",
            "rvCrystalNewton[A1] = [1/4: 8, 3/4: 8]\n",
            "rvCrystalNewton[A2] = [3/8: 8, 5/8: 8]\n",
            "rvCrystalNewton[A3] = [1/2: 16]\n",
        )
    );
}

#[test]
fn characteristic_two_polarization_ordinary_supersingular_and_zero() {
    let input = r"
        rvcrystal O[2,1] {(0,0,1),(1,1,2)};
        rvcrystal S2[2,1] {(0,1,-2),(1,0,1)};
        rvcrystal Empty[11,0] {};
        rvCrystalMatrices[O,1];
        rvCrystalMatrices[S2,2];
        rvCrystalMatrices[Empty,1];
        rvCrystalChar[Empty];
        rvCrystalNewton[Empty];
        assert rvCrystalEOType[O] == eoType{1};
        assert rvCrystalEOType[S2] == eoType{0};
        assert rvCrystalEOType[Empty] == eoType{};
        assert rvCrystalCongruent[O,O,1000] == S.hi;
        assert rvCrystalCongruent[O,S2,1000] == S.lo;
        rvCrystalNewton[O];
        rvCrystalNewton[S2];
    ";
    assert_eq!(
        stdout(interpret(input)),
        concat!(
            "rvCrystalMatrices[O, 1] = F [[1, 0], [0, 0]], V [[0, 0], [0, 1]]\n",
            "rvCrystalMatrices[S2, 2] = F [[0, 2], [1, 0]], V [[0, 2], [3, 0]]\n",
            "rvCrystalMatrices[Empty, 1] = F [], V []\n",
            "rvCrystalChar[Empty] = [1]\n",
            "rvCrystalNewton[Empty] = []\n",
            "rvCrystalNewton[O] = [0: 1, 1: 1]\n",
            "rvCrystalNewton[S2] = [1/2: 2]\n",
        )
    );
}

#[test]
fn invalid_polarization_shape_modulus_and_arithmetic_report_separate_errors() {
    let setup = "rvcrystal A[3,1] {(0,0,1),(1,1,3)}; rvcrystal B[5,1] {(0,0,1),(1,1,5)};";
    for (statement, reason) in [
        (
            "rvcrystal Bad[3,1] {(0,0,1),(1,1,2)};",
            "fails F^t J F = pJ",
        ),
        ("rvcrystal Bad[3,1] {(2,0,1)};", "outside dimension 2"),
        (
            "rvcrystal Bad[3,1] {(0,0,1),(0,0,1)};",
            "duplicate F-crystal matrix entry",
        ),
        (
            "rvcrystal Bad[3,1] {(0,1,-170141183460469231731687303715884105728)};",
            "F-crystal integer arithmetic is out of range",
        ),
        ("rvcrystal Bad[4,1] {};", "characteristic must be prime"),
        ("rvcrystal Bad[9,1] {};", "9 is not prime"),
        ("rvcrystal A[3,1] {};", "F-crystal 'A' is already defined"),
        (
            "rvCrystalMatrices[A,1000];",
            "matrix modulus is out of range",
        ),
        ("rvCrystalMatrices[A,0];", "lattice level must be positive"),
        ("show rvCrystalCongruent[A,B,1];", "same prime and genus"),
        ("show rvCrystalEOType[Missing];", "unknown F-crystal"),
    ] {
        let output = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(output.contains(reason), "{statement}\n{output}");
    }
}

#[test]
fn dense_symplectic_conjugate_keeps_integral_characteristic_polynomial() {
    const G: usize = 7;
    const N: usize = 2 * G;
    let mut f = [[0i128; N]; N];
    for i in 0..G - 1 {
        f[i + 1][i] = 1;
        f[G + i + 1][G + i] = 5;
    }
    f[G][G - 1] = 1;
    f[0][N - 1] = -5;
    f[0][G - 1] = 625;

    // Conjugate by I+v v^t J with v the vector of all ones.
    // r=v^t J, and r v=0, so the inverse is I−v r.
    let r: [i128; N] = std::array::from_fn(|i| if i < G { -1 } else { 1 });
    let right: [i128; N] = std::array::from_fn(|i| f[i].iter().sum());
    let left: [i128; N] = std::array::from_fn(|j| (0..N).map(|i| r[i] * f[i][j]).sum());
    let cross: i128 = left.iter().sum();
    let mut input = String::from("rvcrystal C[5,7] {");
    for (i, row) in f.iter().enumerate() {
        for (j, &entry) in row.iter().enumerate() {
            let conjugate = entry + left[j] - right[i] * r[j] - cross * r[j];
            write!(&mut input, "({i},{j},{conjugate}),").unwrap();
        }
    }
    input.push_str("}; rvCrystalChar[C]; assert rvCrystalEOType[C] == eoType{0,1,2,3,4,5,6};");
    assert_eq!(
        stdout(interpret(&input)),
        "rvCrystalChar[C] = [1, 0, 0, 0, 0, 0, 0, -625, 0, 0, 0, 0, 0, 0, 78125]\n"
    );
}
