mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn exact_phi_two_and_both_cm_fibers() {
    let mut source = String::new();
    for j in 0..4 {
        for i in 0..4 {
            writeln!(source, "rvModularPolyCoeff[2,{i},{j}];").unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let mut actual = [[0i128; 4]; 4];
    for (row, line) in actual.iter_mut().flatten().zip(output.lines()) {
        *row = line.split_once(" = ").unwrap().1.parse().unwrap();
    }
    let expected = [
        [-157_464_000_000_000, 8_748_000_000, -162_000, 1],
        [8_748_000_000, 40_773_375, 1488, 0],
        [-162_000, 1488, -1, 0],
        [1, 0, 0, 0],
    ];
    assert_eq!(actual, expected);
    let cube = |a: i128| [-a.pow(3), 3 * a * a, -3 * a, 1];
    assert_eq!(actual.map(|row| row[0]), cube(54_000));
    let x = 1_728i128;
    let j = 287_496i128;
    let evaluated = actual.map(|row| row.into_iter().rev().fold(0, |acc, c| x * acc + c));
    assert_eq!(evaluated, [-x * j * j, j * j + 2 * x * j, -x - 2 * j, 1]);
}

#[test]
fn exact_higher_levels_are_constructed_and_example_runs() {
    assert_eq!(
        stdout(example("elliptic-modular-polynomials")),
        concat!(
            "rvModularPolyDegree[2] = (3, 3)\n",
            "rvModularPolyCoeff[2,2,2] = -1\n",
            "rvModularPolyCoeff[2,0,0] = -157464000000000\n",
            "rvModularPolyCoeff[5,5,4] = 3720\n",
            "RVUnram(W)[1, 0]\n",
        )
    );
    let output = stdout(interpret(
        "rvModularPolyCoeff[3,0,0]; rvModularPolyCoeff[3,0,1];
         rvModularPolyCoeff[3,0,2]; rvModularPolyCoeff[3,0,3];
         rvModularPolyCoeff[5,0,0]; rvModularPolyCoeff[5,0,5];
         rvModularPolyCoeff[7,8,0]; rvModularPolyCoeff[7,7,7];
         rvModularPolyCoeff[7,7,6];",
    ));
    // The specializations at zero are Y(Y+12288000)^3 at level three
    // and (Y^2+654403829760Y+5209253090426880)^3 at level five.
    assert_eq!(
        output,
        concat!(
            "rvModularPolyCoeff[3,0,0] = 0\n",
            "rvModularPolyCoeff[3,0,1] = 1855425871872000000000\n",
            "rvModularPolyCoeff[3,0,2] = 452984832000000\n",
            "rvModularPolyCoeff[3,0,3] = 36864000\n",
            "rvModularPolyCoeff[5,0,0] = 141359947154721358697753474691071362751004672000\n",
            "rvModularPolyCoeff[5,0,5] = 1963211489280\n",
            "rvModularPolyCoeff[7,8,0] = 1\n",
            "rvModularPolyCoeff[7,7,7] = -1\n",
            "rvModularPolyCoeff[7,7,6] = 5208\n",
        )
    );
}

#[test]
fn velu_pairs_over_prime_and_extension_fields() {
    let output = stdout(interpret(
        r"
        rvunram W[7,1]{0,1};
        rvec E[W](rvUnramElem[W]{1},rvUnramElem[W]{1});
        let k = rveckernel[E,5]{[1],[3],[1]};
        rvecisogeny J[E,T](k);
        assert rvECA[T] == rvUnramElem[W]{2};
        assert rvECB[T] == rvUnramElem[W]{1};
        assert rvModularPolyEval[5](rvUnramElem[W]{1},rvUnramElem[W]{1}) == zero[RVUnram(W)];

        rvunram F[5,1]{2,0,1};
        let z=rvUnramElem[F]{0,1};
        rvec A[F](z,rvUnramElem[F]{1});
        let extension=rveckernel[A,3]{[1,1],[1]};
        rvecisogeny L[A,B](extension);
        assert rvECA[B] == z;
        assert rvECB[B] == rvUnramElem[F]{4,2};
        let x=rvUnramElem[F]{2,1};
        let y=rvUnramElem[F]{1,2};
        assert rvModularPolyEval[3](x,y)==zero[RVUnram(F)];
        assert rvModularPolyEval[3](y,x)==zero[RVUnram(F)];
        assert rvModularPolySpecCoeff[3,4](x)==rvUnramElem[F]{1};
        ",
    ));
    assert!(output.is_empty(), "{output}");
}

#[test]
fn characteristic_level_and_invalid_inputs() {
    let setup = "rvunram F[2,1]{1,1,1}; rvunram G[2,1]{1,1,1};
        rvunram R[2,2]{1,1,1};
        let x=rvUnramElem[F]{0,1}; let y=rvUnramElem[F]{1,1};";
    let output = stdout(interpret(&format!(
        "{setup}
         assert rvModularPolyEval[2](x,y)==zero[RVUnram(F)];
         assert rvModularPolyEval[2](x,x)==rvUnramElem[F]{{1}};
         assert rvModularPolySpecCoeff[2,2](x)==y;
         assert rvModularPolySpecCoeff[2,4](x)==zero[RVUnram(F)];
         assert rvModularPolyEval[3](x,y)==zero[RVUnram(F)];
         assert rvModularPolyEval[3](x,x)==rvUnramElem[F]{{1}};
         assert rvModularPolySpecCoeff[3,3](x)==rvUnramElem[F]{{1}};
         rvModularPolyCoeff[11,99,0]; rvModularPolyDegree[11];"
    )));
    assert_eq!(
        output,
        "rvModularPolyCoeff[11,99,0] = 0\nrvModularPolyDegree[11] = (12, 12)\n"
    );
    for (command, expected) in [
        ("rvModularPolyCoeff[1,0,0];", "level must be prime"),
        ("rvModularPolyDegree[9];", "9 is not prime"),
        (
            "rvModularPolyCoeff[11,0,0];",
            "exact classical modular polynomial exceeds the computation limit",
        ),
        (
            "show rvModularPolyEval[2](x,rvUnramElem[G]{1});",
            "expected type RVUnram(F)",
        ),
        (
            "show rvModularPolySpecCoeff[2,0](rvUnramElem[R]{1});",
            "requires a precision-one RVUnram field",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{command}")));
        assert!(error.contains(expected), "{error}");
    }
}
