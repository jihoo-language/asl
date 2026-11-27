mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const GEOMETRY: &str = include_str!("../examples/hyperelliptic-witt-spectrum.asl");
const INITIAL: usize = 7;

#[test]
fn example_connects_virtual_factors_and_actual_semilinear_spectra() {
    assert_eq!(
        stdout(example("hyperelliptic-witt-spectrum")),
        concat!(
            "RVWittSpectrum{[1, 1], [1]}\n",
            "rvWittSpectrumGhost[1] = -1\n",
            "rvWittSpectrumGhost[2] = 13\n",
            "RVWittSpectrum{[1, -6], [1]}\n",
            "rvWittSpectrumData = numerator degree=22, denominator degree=0, virtual rank=22\n",
            "rvWittSpectrumGhost[2] = 22\n",
            "rvWittSpectrumCoeff[2] = numerator=-81, denominator=0\n",
        ),
    );
}

#[test]
fn exact_gcd_normalizes_nonsplit_factors_and_virtual_operations() {
    let output = stdout(interpret(
        r"
        let unit = rvWittSpectrumLine[1];
        let empty = zero[RVWittSpectrum];
        let sign = rvWittSpectrum{[1,0,-1,0],[1,-1,0,0]};
        assert sign == rvWittSpectrumLine[-1];
        let a = rvWittSpectrum{[1,-1,-1,-2],[1,2,2,1]};
        let b = rvWittSpectrum{[1,-2],[1,1]};
        assert a == b;
        let c = rvWittSpectrum{[1,-4,3], [1,-3,2]};
        assert c == rvWittSpectrum{[1,-3],[1,-2]};
        assert add(a,neg(a)) == empty;
        assert add(a,empty) == a;
        let rootRatio = add(unit,neg(sign));
        assert rvWittSpectrumAdams[2](rootRatio) == empty;
        assert rvWittSpectrumTensor(a,unit) == a;
        assert rvWittSpectrumTensor(a,empty) == empty;
        assert rvWittSpectrumTensor(add(a,b),c) == add(rvWittSpectrumTensor(a,c),rvWittSpectrumTensor(b,c));
        assert rvWittSpectrumAdams[3](rvWittSpectrumTensor(a,c))
            == rvWittSpectrumTensor(rvWittSpectrumAdams[3](a),rvWittSpectrumAdams[3](c));
        show a;
        rvWittSpectrumData(neg(rvWittSpectrum{[1,-5,6]}));
        rvWittSpectrumGhost[1](a);
        rvWittSpectrumGhost[2](a);
        rvWittSpectrumGhost[340282366920938463463374607431768211455](sign);
        rvWittSpectrumCoeff[340282366920938463463374607431768211455](a);
        assert rvWittSpectrumAdams[340282366920938463463374607431768211455](sign) == sign;
        ",
    ));
    assert_eq!(
        output,
        concat!(
            "RVWittSpectrum{[1, -2], [1, 1]}\n",
            "rvWittSpectrumData = numerator degree=0, denominator degree=2, virtual rank=-2\n",
            "rvWittSpectrumGhost[1] = 3\n",
            "rvWittSpectrumGhost[2] = 3\n",
            "rvWittSpectrumGhost[340282366920938463463374607431768211455] = -1\n",
            "rvWittSpectrumCoeff[340282366920938463463374607431768211455] = numerator=0, denominator=0\n",
        ),
    );
}

#[test]
fn tensor_and_adams_match_independent_kronecker_companion_matrices() {
    let a = [[0, 1], [1, 2]];
    let b = [[0, -1], [1, -1]];
    let matrix: [[i128; 4]; 4] =
        std::array::from_fn(|i| std::array::from_fn(|j| a[i / 2][j / 2] * b[i % 2][j % 2]));
    let mut power: [[i128; 4]; 4] =
        std::array::from_fn(|i| std::array::from_fn(|j| i128::from(i == j)));
    let mut traces = Vec::new();
    for _ in 0..12 {
        power = matrix_product(power, matrix);
        traces.push((0..4).map(|i| power[i][i]).sum::<i128>());
    }
    let determinant = newton(&traces[..4]);
    let mut source = format!(
        "let a = rvWittSpectrum{{[1,-2,-1]}};
         let b = rvWittSpectrum{{[1,1,1]}};
         let ab = rvWittSpectrumTensor(a,b);
         assert ab == rvWittSpectrum{{{determinant:?}}};
         let square = rvWittSpectrumAdams[2](ab);"
    );
    for j in 1..=6 {
        writeln!(source, "rvWittSpectrumGhost[{j}](square);").unwrap();
    }
    let output = stdout(interpret(&source));
    assert_eq!(
        output.lines().map(number).collect::<Vec<_>>(),
        traces.into_iter().skip(1).step_by(2).collect::<Vec<_>>(),
    );
}

#[test]
fn exterior_and_symmetric_square_characters_are_the_actual_lambda_operations() {
    let output = stdout(interpret(
        r"
        let x = rvWittSpectrum{[1,-10,31,-30]};
        let exterior = rvWittSpectrumExterior[2](x);
        assert exterior == rvWittSpectrum{[1,-31,300,-900]};
        assert rvWittSpectrumExterior[0](x) == rvWittSpectrumLine[1];
        assert rvWittSpectrumExterior[1](x) == x;
        assert rvWittSpectrumExterior[3](x) == rvWittSpectrumLine[30];
        assert rvWittSpectrumExterior[4](x) == zero[RVWittSpectrum];
        let symmetric = add(rvWittSpectrumTensor(x,x),neg(exterior));
        rvWittSpectrumData(symmetric);
        rvWittSpectrumGhost[1](symmetric);
        rvWittSpectrumGhost[2](symmetric);
        rvWittSpectrumGhost[3](symmetric);
        ",
    ));
    let mut lines = output.lines();
    assert!(lines.next().unwrap().ends_with("virtual rank=6"));
    for k in 1..=3 {
        let moment = |index| {
            [2i128, 3, 5]
                .into_iter()
                .map(|a| a.pow(index))
                .sum::<i128>()
        };
        assert_eq!(
            number(lines.next().unwrap()),
            (moment(k).pow(2) + moment(2 * k)) / 2
        );
    }
    assert!(lines.next().is_none());
}

#[test]
fn geometric_bridge_recovers_the_full_surface_factor_and_ordered_square_decomposition() {
    let mut source = format!(
        r"{GEOMETRY}
        rvK3FrobZeta[Moving,C];
        rvhilbstratum Square{{(1,2)}};
        rvspechtstratum Symmetric[Square]{{}};
        rvspechtstratum Exterior[Square]{{(1,[1,1])}};
        let low = rvK3SpechtSpectrum[Moving,C,Intrinsic,1,0];
        let ordered = add(rvWittSpectrumTensor(low,middle),rvWittSpectrumTensor(middle,low));
        let sym = rvK3SpechtSpectrum[Moving,C,Symmetric,1,2];
        let alt = rvK3SpechtSpectrum[Moving,C,Exterior,1,2];
        assert add(sym,alt) == ordered;
        assert rvK3SpechtSpectrum[Moving,C,Exterior,1,0] == zero[RVWittSpectrum];
        "
    );
    for j in 0..=22 {
        writeln!(source, "rvWittSpectrumCoeff[{j}](middle);").unwrap();
        writeln!(source, "rvK3SpechtEulerFactor[Moving,C,Intrinsic,1,2,{j}];").unwrap();
    }
    let output = stdout(interpret(&source));
    let mut lines = output.lines().skip(INITIAL);
    let (marked, cross) = lines
        .next()
        .unwrap()
        .split_once("D=")
        .unwrap()
        .1
        .split_once(", h=")
        .unwrap();
    let marked: Vec<_> = coefficients(marked)
        .into_iter()
        .enumerate()
        .map(|(i, coefficient)| coefficient * 3i128.pow(i as u32))
        .collect();
    let expected = polynomial_product(&polynomial_product(&marked, &[1, -3]), &coefficients(cross));
    assert_eq!(expected.len(), 23);
    for expected in expected {
        let report = lines.next().unwrap().split_once("numerator=").unwrap().1;
        assert_eq!(
            report.split_once(',').unwrap().0.parse::<i128>().unwrap(),
            expected
        );
        assert_eq!(number(lines.next().unwrap()), expected);
    }
    assert!(lines.next().is_none());
}

#[test]
fn bounds_and_wrong_inputs_have_direct_diagnostics() {
    let high = format!("rvWittSpectrum{{[1,{},1]}}", ["0"; 16].join(","));
    let tensor = format!("show rvWittSpectrumTensor({high},{high});");
    for (source, expected) in [
        ("show rvWittSpectrum{[]};", "constant coefficient 1"),
        (
            "show rvWittSpectrum{[1,2],[0,3]};",
            "constant coefficient 1",
        ),
        ("show rvWittSpectrumLine[0];", "nonzero eigenvalue"),
        (
            "show rvWittSpectrumAdams[0](rvWittSpectrumLine[1]);",
            "index must be positive",
        ),
        (
            "rvWittSpectrumGhost[0](rvWittSpectrumLine[1]);",
            "index must be positive",
        ),
        (
            "show rvWittSpectrumExterior[2](neg(rvWittSpectrumLine[2]));",
            "effective spectrum",
        ),
        (
            "show rvWittSpectrumTensor(rvWittSpectrumLine[1],S.hi);",
            "expected type RVWittSpectrum",
        ),
        (
            "show rvWittSpectrumAdams[1000](rvWittSpectrumLine[2]);",
            "exact intermediate exceeds 512 bits",
        ),
        (
            "show rvWittSpectrumAdams[200](rvWittSpectrumLine[2]);",
            "exceeds signed 128-bit range",
        ),
        (tensor.as_str(), "determinant degree limit"),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{error}");
    }
    let bounded_bridge = format!(
        "{GEOMETRY} rvhilbstratum Huge{{(1,3)}}; show rvK3StratumSpectrum[Moving,C,Huge,1,6];"
    );
    let error = diagnostic(interpret(&bounded_bridge));
    assert!(error.contains("determinant degree limit"), "{error}");
}

fn matrix_product(a: [[i128; 4]; 4], b: [[i128; 4]; 4]) -> [[i128; 4]; 4] {
    std::array::from_fn(|i| std::array::from_fn(|j| (0..4).map(|k| a[i][k] * b[k][j]).sum()))
}

fn newton(traces: &[i128]) -> Vec<i128> {
    let mut determinant = vec![1];
    for k in 1..=traces.len() {
        let sum = traces
            .iter()
            .zip(determinant.iter().rev())
            .map(|(a, b)| a * b)
            .sum::<i128>();
        assert_eq!(sum % k as i128, 0);
        determinant.push(-sum / k as i128);
    }
    determinant
}

fn polynomial_product(a: &[i128], b: &[i128]) -> Vec<i128> {
    let mut result = vec![0; a.len() + b.len() - 1];
    for (i, &a) in a.iter().enumerate() {
        for (j, &b) in b.iter().enumerate() {
            result[i + j] += a * b;
        }
    }
    result
}

fn coefficients(text: &str) -> Vec<i128> {
    text.strip_prefix('[')
        .unwrap()
        .strip_suffix(']')
        .unwrap()
        .split(", ")
        .map(|n| n.parse().unwrap())
        .collect()
}

fn number(line: &str) -> i128 {
    line.split_once(" = ").unwrap().1.parse().unwrap()
}
