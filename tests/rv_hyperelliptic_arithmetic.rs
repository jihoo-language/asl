mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn actual_frobenius_and_distinct_rational_component_groups_of_twists() {
    let output = stdout(example("rv-hyperelliptic-arithmetic"));
    for line in [
        "rvHyperComponent = factors=[2, 62], exponent=62, order=124",
        "rvHyperComponent = factors=[2, 30], exponent=30, order=60",
        "rvHyperFixed = factors=[2, 62], order=124, representatives=[[7, 5], [11, 8]]",
        "rvHyperFixed = factors=[2, 2], order=4, representatives=[[7, 5], [5, 8]]",
        "rvHyperFixed = factors=[2, 6], order=12, representatives=[[5, 0], [0, 1]]",
        "rvHyperFixed = factors=[2, 10], order=20, representatives=[[5, 3], [9, 0]]",
        "cycles=[[1, 0], [0, 1]], dual=[[1, 0], [0, 1]]",
        "cycles=[[-1, 0], [0, -1]], dual=[[-1, 0], [0, -1]]",
        "cycles=[[-1, 0], [0, 1]], dual=[[-1, 0], [0, 1]]",
        "cycles=[[1, 0], [0, -1]], dual=[[1, 0], [0, -1]]",
        "rvHyperFrobenius = vertices=[0, 2, 1, 3, 4]",
    ] {
        assert!(output.contains(line), "{line}\n{output}");
    }
}

#[test]
fn rational_units_and_coordinates_preserve_the_same_arithmetic_cover() {
    let source = r"
        let a = rvHyperArithmetic[7,1]{0,49,1,344,2,16809};
        let reduced = rvHyperArithmetic[7,(1,64)]{0,(49,2),(1,2),172,1,(16809,2)};
        assert a == reduced;
        let moved = rvHyperArithmetic[7,1]{3,101,5,691,7,33621};
        assert a == moved;
        assert rvHyperGeometric(a) == rvHyperCover(rvStableLine[7]{
            [0,1],[49,1],[1,1],[344,1],[2,1],[16809,1]});
        rvHyperFixed(rvHyperArithmetic[7,1]{1,16809,49,0,344,2});
        let b = rvHyperArithmetic[7,(5,4)]{0,49,1,344,2,16809};
        assert b == rvHyperArithmetic[7,3]{0,49,1,344,2,16809};
        rvHyperFixed(b);
    ";
    let output = stdout(interpret(source));
    assert!(output.contains("factors=[2, 62], order=124"), "{output}");
    assert!(output.contains("factors=[2, 2], order=4"), "{output}");
}

#[test]
fn smooth_and_separating_examples_have_trivial_rational_component_groups() {
    let source = r"
        let smooth = rvHyperArithmetic[7,3]{0,1,2,3,4,5};
        let separating = rvHyperArithmetic[7,1]{0,49,98,1,2,3};
        rvHyperFrobenius(smooth);
        rvHyperFixed(smooth);
        rvHyperFrobenius(separating);
        rvHyperFixed(separating);
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "rvHyperFrobenius = vertices=[0], edges=[], cycles=[], dual=[]\n",
            "rvHyperFixed = factors=[], order=1, representatives=[]\n",
            "rvHyperFrobenius = vertices=[0, 1], edges=[0], cycles=[], dual=[]\n",
            "rvHyperFixed = factors=[], order=1, representatives=[]\n",
        )
    );
}

#[test]
fn arithmetic_input_contract_checks_chart_parity_and_its_own_types() {
    for (program, expected) in [
        (
            "show rvHyperArithmetic[7,1]{0,7,14,1,2,3};",
            "odd normalization valuation 3",
        ),
        (
            "show rvHyperArithmetic[7,1]{0,7,1,8};",
            "requires at least three depth-zero residue directions",
        ),
        (
            "show rvHyperArithmetic[7,7]{0,1,2,3};",
            "leading coefficient must be a p-adic unit",
        ),
        (
            "show rvHyperArithmetic[7,1]{0,1,2,(1,7)};",
            "branch root 3 must be p-integral",
        ),
        (
            "show rvHyperArithmetic[7,1]{0,1,2,(2,2)};",
            "branch roots 1 and 3 coincide",
        ),
        (
            "show rvHyperArithmetic[7,1]{0,1,2,(2,0)};",
            "denominator must be positive",
        ),
        (
            "show rvHyperArithmetic[2,1]{0,1,2,3};",
            "requires odd residue characteristic",
        ),
        (
            "rvHyperFixed(rvHyperGeometric(rvHyperArithmetic[7,1]{0,1,2,3}));",
            "requires an RVHyperArithmetic",
        ),
    ] {
        let output = diagnostic(interpret(program));
        assert!(output.contains(expected), "{program}\n{output}");
    }
}

#[test]
fn canonical_charts_preserve_nonsquare_coordinate_changes_and_distinguish_components() {
    let same = r"
        let canonical = rvHyperArithmetic[7,1]{0,49,98,1,50,2,3,4};
        let moved = rvHyperArithmetic[7,(1,6561)]{2,149,296,5,152,8,11,14};
        let square_twist = rvHyperArithmetic[7,4]{0,49,98,1,50,2,3,4};
        assert canonical == moved;
        assert canonical == square_twist;
    ";
    assert_eq!(stdout(interpret(same)), "");
    let different = r"
        let a = rvHyperArithmetic[7,1]{0,1,2,3};
        let b = rvHyperArithmetic[7,1]{0,1,2,4};
        assert rvHyperGeometric(a) == rvHyperGeometric(b);
        assert a == b;
    ";
    let output = diagnostic(interpret(different));
    assert!(
        output.contains("assertion failed: RVHyperArithmetic"),
        "{output}"
    );
}

#[test]
fn bounded_characters_for_larger_split_nonsplit_and_mixed_graphs() {
    let roots = "0,121,1,12,2,123,3,14,4,1335,5,16,6,1337,7,1338,8,129";
    let mixed = concat!(
        "0,1,2,10,285311670613,570623341224,855935011835,1141246682446,",
        "1426558353057,34522712143931,45949729863572163,61159091304349558126,",
        "895430243255522683917144,1191817653772722083706815207,",
        "17449402268886407319129427095025,23225154419887808141003194354662188"
    );
    let output = stdout(interpret(&format!(
        "rvHyperFixed(rvHyperArithmetic[11,1]{{{roots}}});
         rvHyperFixed(rvHyperArithmetic[11,2]{{{roots}}});
         rvHyperFixed(rvHyperArithmetic[11,2]{{{mixed}}});"
    )));
    let lines: Vec<_> = output.lines().collect();
    for (line, (group, rank, bound)) in lines.iter().zip([
        ("factors=[2, 2, 2, 2, 2, 2, 12, 396], order=304128", 8, 39),
        ("factors=[2, 2, 2, 2, 2, 2, 2, 2], order=256", 8, 39),
        ("factors=[2, 2, 2, 2, 2, 2, 26], order=1664", 7, 176),
    ]) {
        let (actual, representatives) = line.split_once(", representatives=").unwrap();
        assert_eq!(actual, format!("rvHyperFixed = {group}"));
        assert_eq!(representatives.matches('[').count(), rank + 1);
        for coordinate in representatives.split(|c: char| !c.is_ascii_digit() && c != '-') {
            if !coordinate.is_empty() {
                let x: i64 = coordinate.parse().unwrap();
                assert!((0..=bound).contains(&x), "{line}");
            }
        }
    }
    assert_eq!(lines.len(), 3);
}

#[test]
fn nonsplit_fixed_group_can_fit_when_the_geometric_order_does_not() {
    let roots = (0..28)
        .flat_map(|i| [i, i + 37i128.pow(11 + (i % 9) as u32)])
        .map(|root| root.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let declaration = format!("let A=rvHyperArithmetic[37,2]{{{roots}}};");
    let output = stdout(interpret(&format!("{declaration} rvHyperFixed(A);")));
    let factors = vec![2; 27];
    assert!(
        output.starts_with(&format!(
            "rvHyperFixed = factors={factors:?}, order=134217728, representatives="
        )),
        "{output}"
    );
    let error = diagnostic(interpret(&format!(
        "{declaration} rvHyperComponent(rvHyperGeometric(A));"
    )));
    assert!(error.contains("hyperelliptic component-group order exceeds signed i128"));
}
