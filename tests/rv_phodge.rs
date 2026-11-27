mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn weak_admissibility_depends_on_monodromy_and_the_actual_flag() {
    assert_eq!(
        stdout(example("rv-filtered-phi-n")),
        concat!(
            "rvPHodgeWeak[Good] = weakly admissible\n",
            "rvPHodgeData[Good, [0, 1]] = tH 1, tN 1, degree 0, Hodge [-1, 2], Newton [0, 1]\n",
            "rvPHodgeHN[Good] = [([0, 1, 2], 0)]\n",
            "rvPHodgeWeak[NoN] = destabilized by [1] (tH 2, tN 1)\n",
            "rvPHodgeWeak[Moved] = destabilized by [0] (tH 2, tN 0)\n",
            "rvPHodgeHN[Moved] = [([0], 2), ([0, 1], 1), ([0, 1, 2], -3)]\n",
            "rvPHodgeData[Fractional, [0, 1]] = tH 2, tN 1, degree 1, Hodge [0, 2], Newton [0, 1]\n",
            "rvPHodgeHN[Fractional] = [([0, 1], 1/2), ([0, 1, 2], -1)]\n",
        )
    );
}

#[test]
fn an_adapted_basis_singular_mod_three_has_the_same_harder_narasimhan_flag() {
    let input = r"
        rvphodge Original[3] {
            {(1,1),(3,1),(9,1)},
            {(0,1,2,5),(1,2,-7,2)},
            {(0,[1,0,0]),(1,[0,1,1]),(2,[1,1,0])}
        };
        rvphodge Changed[3] {
            {(1,1),(3,1),(9,1)},
            {(0,1,2,5),(1,2,-7,2)},
            {(0,[1,0,0]),(1,[1,4,3]),(2,[1,1,0])}
        };
        rvPHodgeData[Original,[1,0]];
        rvPHodgeData[Changed,[1,0]];
        rvPHodgeHN[Original];
        rvPHodgeHN[Changed];
    ";
    assert_eq!(
        stdout(interpret(input)),
        concat!(
            "rvPHodgeData[Original, [0, 1]] = tH 2, tN 1, degree 1, Hodge [0, 2], Newton [0, 1]\n",
            "rvPHodgeData[Changed, [0, 1]] = tH 2, tN 1, degree 1, Hodge [0, 2], Newton [0, 1]\n",
            "rvPHodgeHN[Original] = [([0, 1], 1/2), ([0, 1, 2], -1)]\n",
            "rvPHodgeHN[Changed] = [([0, 1], 1/2), ([0, 1, 2], -1)]\n",
        )
    );
}

#[test]
fn rational_slopes_and_the_largest_maximal_destabilizer() {
    let input = r"
        rvphodge Rational[3] {
            {(1,3),(1,1),(-2,3)},
            {(0,1,-2,7)},
            {(-1,[1,0,0]),(0,[0,1,0]),(-1,[0,0,1])}
        };
        rvPHodgeData[Rational,[1,0,2]];
        rvPHodgeWeak[Rational];
        rvphodge Ties[3] {
            {(1,1),(3,1),(9,1)},
            {},
            {(2,[1,0,0]),(3,[0,1,0]),(0,[0,0,1])}
        };
        rvPHodgeWeak[Ties];
        rvPHodgeHN[Ties];
        rvphodge Empty[2] {{},{},{}};
        rvPHodgeData[Empty];
        rvPHodgeWeak[Empty];
        rvPHodgeHN[Empty];
    ";
    assert_eq!(
        stdout(interpret(input)),
        concat!(
            "rvPHodgeData[Rational, [0, 1, 2]] = tH -2, tN -2, degree 0, Hodge [-1, -1, 0], Newton [-1, -1, 0]\n",
            "rvPHodgeWeak[Rational] = weakly admissible\n",
            "rvPHodgeWeak[Ties] = degree equality fails on [0, 1, 2] (tH 5, tN 3)\n",
            "rvPHodgeHN[Ties] = [([0, 1], 2), ([0, 1, 2], -2)]\n",
            "rvPHodgeData[Empty] = tH 0, tN 0, degree 0, Hodge [], Newton []\n",
            "rvPHodgeWeak[Empty] = weakly admissible\n",
            "rvPHodgeHN[Empty] = []\n",
        )
    );
}

#[test]
fn declaration_and_stability_errors_are_checked() {
    let setup = "rvphodge D[3] {{(1,1),(3,1)}, {(0,1,1,1)}, {(0,[1,0]),(1,[0,1])}};";
    for (statement, reason) in [
        ("rvphodge D[3] {{},{},{}};", "already defined"),
        (
            "rvphodge A[4] {{},{},{}};",
            "residue characteristic must be prime",
        ),
        (
            "rvphodge A[3] {{(0,1)},{},{(0,[1])}};",
            "distinct nonzero Frobenius eigenvalues",
        ),
        (
            "rvphodge A[3] {{(1,2),(2,4)},{},{(0,[1,0]),(1,[0,1])}};",
            "distinct nonzero Frobenius eigenvalues",
        ),
        (
            "rvphodge A[3] {{(1,0)},{},{(0,[1])}};",
            "rational denominator must be positive",
        ),
        (
            "rvphodge A[3] {{(1,1)},{(0,0,1,1)},{(0,[1])}};",
            "fails N phi = p phi N",
        ),
        (
            "rvphodge A[3] {{(1,1)},{(1,0,0,1)},{(0,[1])}};",
            "outside dimension 1",
        ),
        (
            "rvphodge A[3] {{(1,1)},{(0,0,0,1),(0,0,0,1)},{(0,[1])}};",
            "duplicate monodromy entry",
        ),
        (
            "rvphodge A[3] {{(1,1)},{},{}};",
            "requires 1 filtration columns",
        ),
        (
            "rvphodge A[3] {{(1,1)},{},{(0,[])}};",
            "each filtration column must have 1 entries",
        ),
        (
            "rvphodge A[3] {{(1,1),(3,1)},{},{(0,[1,1]),(1,[2,2])}};",
            "basis over characteristic zero",
        ),
        (
            "rvPHodgeData[D,[1]];",
            "not phi,N-stable: coordinate 1 requires 0",
        ),
        ("rvPHodgeData[D,[0,0]];", "coordinate is repeated"),
        ("rvPHodgeData[D,[2]];", "coordinate is out of range"),
        ("rvPHodgeHN[Missing];", "unknown filtered module"),
    ] {
        let output = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(output.contains(reason), "{statement}\n{output}");
    }
}

fn degree_zero_module(with_chain: bool) -> String {
    let mut source = String::from("rvphodge D[3] {{");
    for i in 0..20 {
        write!(&mut source, "({},1),", 3i128.pow(i)).unwrap();
    }
    source.push_str("},{");
    for i in 1..20 {
        if with_chain {
            write!(&mut source, "({},{i},1,1),", i - 1).unwrap();
        }
    }
    source.push_str("},{");
    for i in 0..20 {
        write!(&mut source, "({i},[").unwrap();
        for j in 0..20 {
            write!(&mut source, "{},", usize::from(i == j)).unwrap();
        }
        source.push_str("]),");
    }
    source.push_str("}};");
    source
}

#[test]
fn prefix_enumeration_handles_long_chains_and_limits_wide_lattices() {
    let chain = degree_zero_module(true);
    assert_eq!(
        stdout(interpret(&format!("{chain} rvPHodgeWeak[D];"))),
        "rvPHodgeWeak[D] = weakly admissible\n",
    );
    let wide = degree_zero_module(false);
    let output = diagnostic(interpret(&format!("{wide} rvPHodgeHN[D];")));
    assert!(output.contains("exceeds the storage limit"), "{output}");
    assert_eq!(
        stdout(interpret(&format!("{wide} rvPHodgeData[D,[19]];"))),
        "rvPHodgeData[D, [19]] = tH 19, tN 19, degree 0, Hodge [19], Newton [19]\n",
    );
}
