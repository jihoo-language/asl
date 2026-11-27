mod common;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = "dwork-fricke-modular-correspondence";

#[test]
fn eta_product_and_eisenstein_series_agree_with_the_exact_mirror_map() {
    let out = stdout(example(EXAMPLE));
    let row = out
        .lines()
        .find(|line| line.starts_with("rvDworkFrickeSeries"))
        .unwrap();
    assert!(row.contains("q*t=[1, -24, 276, -2048"));
    assert!(row.contains("z(q)=[0, 1, -104, 6444, -311744, 13018830, -493025760, 17411253944"));
    assert!(row.contains("f(z(q))=[1, 24, 24, 96, 24, 144, 96, 192, 24]"));
    assert!(row.ends_with(
        "certified inverse mirror map, Eisenstein difference and Wronskian normalization"
    ));
    let boundaries = stdout(interpret(
        "rvDworkFrickeSeries[0]; rvDworkFrickeSeries[12];",
    ));
    assert!(boundaries.contains("rvDworkFrickeSeries[0] = q*t=[1], z(q)=[0, 1], f(z(q))=[1]"));
    assert!(boundaries.contains("rvDworkFrickeSeries[12]"));
}

#[test]
fn fricke_involution_and_classical_polynomial_are_exact() {
    let out = stdout(example(EXAMPLE));
    assert!(out.contains(
        "J^2-(v^2-207v+3456)J+(v+144)^3, certified W2 exchange and Phi_2(j1,j2)=0 over Q(t)"
    ));
    assert!(out.contains("z=[1], split, t=[2] (multiplicity 1), j1=[5], j2=[4] | t=[4] (multiplicity 1), j1=[4], j2=[5]"));
    assert!(out.contains("z=[3], nonsplit quadratic=[[1], [4], [1]], j1=[[3], [6]], j2=[[0], [1]], Frobenius t->W2(t), distinct j-values exchanged"));
    let branch = stdout(interpret(
        "rvunram F[7,1]{0,1}; rvDworkFrickeBase(rvUnramElem[F]{2});",
    ));
    assert!(branch.contains("t=[1] (multiplicity 2)"));
}

#[test]
fn actual_velu_map_dual_and_quadratic_twist_satisfy_their_equations() {
    let out = stdout(example(EXAMPLE));
    assert!(out.contains("rvECMapDegree = 2"));
    assert!(out.contains("rvECMapCertify = degree=2, cleared identity=yes"));
    assert!(out.contains("trace E=trace E'=2, trace E(W2(t))=-2, twist=[3], twist character=-1"));
    assert!(out.contains("rvECPointCount[E] = 6, trace=2"));
    assert!(out.contains("rvECPointCount[Q] = 6, trace=2"));

    // Long models and their traces are available in characteristic three.
    let long = stdout(interpret(
        "rvunram F[3,1]{0,1}; rvDworkFrickeFiber(rvUnramElem[F]{1});",
    ));
    assert!(long.contains("long E=([2], [2]), long E' = ([2], [2]), trace E=trace E'=2"));
}

#[test]
fn independently_computed_units_match_at_split_fibers() {
    let out = stdout(example(EXAMPLE));
    for (p, psi, value, cases) in [
        (7, 2, 39, [(2, 2, 23), (4, -2, 26)]),
        (11, 2, 115, [(3, 4, 92), (5, 4, 92)]),
        (17, 2, 271, [(9, -4, 217), (15, -4, 217)]),
    ] {
        let prefix = format!("rvDworkFrickeUnit[{p}, {psi}, 2]");
        let row = out.lines().find(|line| line.starts_with(&prefix)).unwrap();
        assert!(row.contains(&format!("Dwork unit={value},")), "{row}");
        for (t, trace, root) in cases {
            assert!(
                row.contains(&format!(
                    "t={t}, trace={trace}, elliptic unit={root}, square={value}, comparison=equal"
                )),
                "{row}"
            );
            assert_eq!((root * root - trace * root + p) % (p * p), 0);
        }
    }
    let supersingular = out
        .lines()
        .find(|line| line.starts_with("rvDworkFrickeUnit[29"))
        .unwrap();
    assert!(supersingular.contains("Hasse–Witt=0, Dwork unit=none"));
    assert_eq!(
        supersingular
            .matches("trace=0, elliptic unit=none, square=none, comparison=both nonordinary")
            .count(),
        2
    );
    assert!(out.contains("Dwork unit=6, nonsplit Fricke parameters over the prime field"));
}

#[test]
fn validates_precision_parameters_and_nominal_curve_models() {
    let base = "rvunram F[7,1]{0,1}; let t=rvUnramElem[F]{2}; ";
    for (source, fragment) in [
        (
            format!("{base} show rvDworkFrickeJ[3](t);"),
            "j index must be 1 or 2",
        ),
        (
            format!("{base} rvDworkFrickeFiber(rvUnramElem[F]{{6}});"),
            "t different from 0 and -64",
        ),
        (
            format!("{base} rvdworkfricke E,E(t);"),
            "requires distinct names",
        ),
        (
            format!(
                "{base} rvdworkfricke E,Q(t); show rvDworkFrickeMap[E,Q](rvUnramElem[F]{{4}});"
            ),
            "do not match the ordered Fricke models",
        ),
        (
            "rvunram K[3,1]{0,1}; rvdworkfricke E,Q(rvUnramElem[K]{1});".into(),
            "require p at least 5",
        ),
        (
            "rvunram K[7,2]{0,1}; rvDworkFrickeFiber(rvUnramElem[K]{1});".into(),
            "requires an odd residue field at precision one",
        ),
        ("rvDworkFrickeUnit[7,0,2];".into(), "requires psi nonzero"),
        ("rvDworkFrickeUnit[7,1,2];".into(), "singular Dwork quartic"),
        (
            "rvDworkFrickeSeries[13];".into(),
            "degree must be at most 12",
        ),
    ] {
        let error = diagnostic(interpret(&source));
        assert!(error.contains(fragment), "{error}");
    }
}
