mod common;

use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

type Matrix = [[i128; 2]; 2];

fn number(line: &str) -> i128 {
    line.rsplit_once(" = ").unwrap().1.parse().unwrap()
}

fn path(line: &str) -> (usize, i128) {
    let (count, jump) = line
        .rsplit_once(" = walls=")
        .unwrap()
        .1
        .split_once(", jump=")
        .unwrap();
    (count.parse().unwrap(), jump.parse().unwrap())
}

fn congruence([a, b, c]: [i128; 3], [[x, y], [u, v]]: Matrix) -> [i128; 3] {
    [
        a * x * x + 2 * b * x * u + c * u * u,
        a * x * y + b * (x * v + y * u) + c * u * v,
        a * y * y + 2 * b * y * v + c * v * v,
    ]
}

fn covariant_index([m, n, r]: [i128; 3], [[a, b], [c, d]]: Matrix) -> [i128; 3] {
    let [m, r, n] = congruence([2 * m, r, 2 * n], [[d, -c], [-b, a]]);
    [m / 2, n / 2, r]
}

#[test]
fn the_sign_wall_recovers_the_two_chambers_at_y_one() {
    let mut source = "rvsiegelchamber I[2,1,3]; rvsiegelchamber O[2,-1,3];".to_string();
    for r in -3..=3 {
        writeln!(source, "rvIgusaChamberCoeff[I,0,0,{r}]; rvIgusaInverse[0,0,{r},inside]; rvIgusaChamberCoeff[O,0,0,{r}]; rvIgusaInverse[0,0,{r},outside]; rvIgusaChamberPath[I,O,0,0,{r}]; rvIgusaWallJump[0,0,{r}];").unwrap();
    }
    let output = stdout(interpret(&source));
    let formatted: Vec<_> = output.lines().collect();
    for row in formatted.chunks_exact(6) {
        assert_eq!(number(row[0]), number(row[1]));
        assert_eq!(number(row[2]), number(row[3]));
        assert_eq!(path(row[4]), (1, number(row[5])));
    }
    assert_eq!(formatted.len(), 7 * 6);
}

#[test]
fn each_shear_wall_has_the_eta_residue_in_its_own_frame() {
    let output = stdout(interpret(
        r"
        rvsiegelchamber H[2,1,3];
        rvsiegelchamber J[2,3,7];
        rvsiegelchamber K[10,5,3];
        rvIgusaChamberCoeff[H,3,1,-4]; rvIgusaChamberCoeff[J,3,1,-4];
        rvIgusaChamberPath[H,J,3,1,-4]; rvIgusaChamberPath[J,H,3,1,-4];
        rvIgusaChamberCoeff[H,1,3,-4]; rvIgusaChamberCoeff[K,1,3,-4];
        rvIgusaChamberPath[H,K,1,3,-4]; rvIgusaChamberPath[K,H,1,3,-4];
    ",
    ));
    let lines: Vec<_> = output.lines().collect();
    for row in lines.chunks_exact(4) {
        // Both wall frames turn the global index into (0,1,2).
        let residue = 2 * 24 * 324;
        assert_eq!(number(row[0]) - number(row[1]), residue);
        assert_eq!(path(row[2]), (1, residue));
        assert_eq!(path(row[3]), (1, -residue));
    }
    assert_eq!(lines.len(), 8);
}

#[test]
fn a_three_wall_path_sums_nonzero_jumps_and_reverses() {
    // The path crosses B=A, then B=0, then B=A in the successive
    // frames. The wall indices are (1,0,-1), (2,0,-1), (3,0,-1).
    let total = -24 * (324 + 3200 + 25650);
    let report = stdout(example("igusa-farey-chambers"));
    let lines: Vec<_> = report.lines().collect();
    assert_eq!(
        lines[0],
        "rvSiegelChamberData[Far] = reduced=[2, 1, 3], basis=[[1, -3], [0, 1]]"
    );
    assert_eq!(number(lines[1]), 600);
    assert_eq!(number(lines[2]), 24);
    assert_eq!(number(lines[1]) - number(lines[3]), total);
    assert_eq!(path(lines[4]), (1, 576));
    assert_eq!(path(lines[5]), (3, total));
    assert_eq!(path(lines[6]), (3, -total));
    assert_eq!(lines.len(), 7);
}

#[test]
fn general_unimodular_basis_changes_and_scaling_preserve_the_coefficients() {
    let basis_changes = [
        [[0, 1], [1, 0]],
        [[1, 1], [0, -1]],
        [[0, -1], [1, 1]],
        [[-1, -1], [1, 0]],
        [[1, 0], [-1, -1]],
        [[1, 0], [0, 1]],
        [[2, 1], [3, 2]],
        [[0, -1], [1, 2]],
        [[1, 2], [1, 1]],
    ];
    let y0 = [3, 1, 5];
    let indices = [[0, 0, 1], [-1, -1, 1], [1, 0, 0]];
    let mut source = "rvsiegelchamber Base[3,1,5]; rvsiegelchamber Scale[33,11,55];".to_string();
    for [m, n, r] in indices {
        writeln!(source, "rvIgusaChamberCoeff[Base,{m},{n},{r}]; rvIgusaChamberCoeff[Scale,{m},{n},{r}]; rvIgusaChamberPath[Base,Scale,{m},{n},{r}];").unwrap();
    }
    for (i, matrix) in basis_changes.into_iter().enumerate() {
        let [a, b, c] = congruence(y0, matrix);
        writeln!(source, "rvsiegelchamber H{i}[{a},{b},{c}];").unwrap();
        for index @ [m, n, r] in indices {
            let [mp, np, rp] = covariant_index(index, matrix);
            writeln!(source, "rvIgusaChamberCoeff[H{i},{mp},{np},{rp}];").unwrap();
            if i < 6 {
                // These six elements stabilize the whole open chamber.
                writeln!(source, "rvIgusaChamberCoeff[H{i},{m},{n},{r}]; rvIgusaChamberPath[Base,H{i},{m},{n},{r}];").unwrap();
            }
        }
    }
    let actual = stdout(interpret(&source));
    let mut lines = actual.lines();
    let expected: [i128; 3] = std::array::from_fn(|_| {
        let first = number(lines.next().unwrap());
        assert_eq!(number(lines.next().unwrap()), first);
        assert_eq!(path(lines.next().unwrap()), (0, 0));
        first
    });
    for i in 0..basis_changes.len() {
        for value in expected {
            assert_eq!(number(lines.next().unwrap()), value, "basis {i}");
            if i < 6 {
                assert_eq!(number(lines.next().unwrap()), value);
                assert_eq!(path(lines.next().unwrap()), (0, 0));
            }
        }
    }
    assert_eq!(lines.next(), None);

    // A common basis change and independent positive scalings of two
    // period rays preserve the whole oriented wall sum and path length.
    for matrix in basis_changes.into_iter().skip(6) {
        let [a, b, c] = congruence([2, 1, 3], matrix).map(|x| x * 7);
        let [d, e, f] = congruence([2, 7, 27], matrix).map(|x| x * 11);
        let [m, n, r] = covariant_index([0, 0, 1], matrix);
        let output = stdout(interpret(&format!(
            "rvsiegelchamber X[{a},{b},{c}]; rvsiegelchamber Y[{d},{e},{f}]; rvIgusaChamberPath[X,Y,{m},{n},{r}];"
        )));
        assert_eq!(path(output.trim_end()), (3, -700176));
    }
}

#[test]
fn large_shears_reduce_directly_and_widen_the_inverse_index() {
    let shear = 1_000_000_000_000_000_000i128;
    let transform = [[1, shear], [0, 1]];
    let [a, b, c] = congruence([2, 1, 3], transform);
    let [m, n, r] = covariant_index([-1, -1, 1], transform);
    let hi = i128::MAX;
    let lo = i128::MIN;
    let output = stdout(interpret(&format!(
        "rvsiegelchamber H[{a},{b},{c}]; rvSiegelChamberData[H]; rvIgusaChamberCoeff[H,{m},{n},{r}];\
         rvsiegelchamber I[2,1,3]; rvsiegelchamber J[2,3,7];\
         rvIgusaChamberCoeff[J,{lo},1,{hi}]; rvIgusaChamberPath[I,J,{lo},1,{hi}];"
    )));
    assert_eq!(output, format!(
        "rvSiegelChamberData[H] = reduced=[2, 1, 3], basis=[[1, -{shear}], [0, 1]]\n\
         rvIgusaChamberCoeff[H, {m}, {n}, {r}] = 1\n\
         rvIgusaChamberCoeff[J, {lo}, 1, {hi}] = 1323017842588608745945600473694714806148704\n\
         rvIgusaChamberPath[I, J, {lo}, 1, {hi}] = walls=1, jump=-1323017842588608745945600473694714806148704\n"
    ));
    let err = diagnostic(interpret(&format!("rvsiegelchamber Base[2,1,3]; rvsiegelchamber Far[{a},{b},{c}]; rvIgusaChamberPath[Base,Far,0,0,0];")));
    assert!(err.contains("exceeds the 256-wall limit"), "{err}");

    // Both columns are large and nearly parallel. The transformed period
    // has determinant five despite quadratic products beyond 128 bits.
    let transform = [[shear, shear - 1], [shear + 1, shear]];
    let [a, b, c] = congruence([2, 1, 3], transform);
    let [m, n, r] = covariant_index([-1, -1, 1], transform);
    let output = stdout(interpret(&format!(
        "rvsiegelchamber H[{a},{b},{c}]; rvIgusaChamberCoeff[H,{m},{n},{r}];"
    )));
    assert_eq!(number(output.trim_end()), 1);
}

#[test]
fn only_positive_definite_open_periods_are_registered() {
    for (source, expected) in [
        ("rvsiegelchamber H[0,1,3];", "must be positive definite"),
        ("rvsiegelchamber H[2,3,3];", "must be positive definite"),
        ("rvsiegelchamber H[4,6,9];", "must be positive definite"),
        ("rvsiegelchamber H[2,0,3];", "integral separating wall"),
        ("rvsiegelchamber H[2,2,5];", "integral separating wall"),
        ("rvsiegelchamber H[5,3,3];", "integral separating wall"),
        (
            "rvsiegelchamber H[2,1,3]; rvsiegelchamber H[3,1,3];",
            "already defined",
        ),
        ("rvIgusaChamberCoeff[H,0,0,0];", "unknown Siegel chamber"),
        ("rvSiegelChamberData[H];", "unknown Siegel chamber"),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{source}: {error}");
    }
}
