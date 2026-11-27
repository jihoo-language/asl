#[allow(dead_code)]
mod common;

use common::{diagnostic, interpret, stdout};

const APERY: &str = r"
    let apery = add(rvOre[0,2],add(rvOre[1,2,-11],add(rvOre[1,1,-11],
        add(rvOre[1,0,-3],add(rvOre[2,2,-1],add(rvOre[2,1,-2],rvOre[2,0,-1]))))));
";
const APERY_CERTIFICATE: &str = r"
    {(0,5,1),(1,4,6),(0,4,-5),(2,3,-11),(1,3,7)},
    {(1,0,0,1),(1,1,0,1),(1,1,-1,1)}
";

#[test]
fn canonical_apery_and_franel_periods_certify_the_declared_polynomials() {
    let program = format!(
        r"
        rvperiod Apery[2] {{
            ([-1,-1],1), ([-1,0],2), ([-1,1],1), ([0,-1],2),
            ([0,0],3), ([0,1],1), ([1,-1],1), ([1,0],1),
            ([100,100],7), ([100,100],-7)
        }};
        {APERY}
        rvCTWZCheck[Apery,2,1,1,2](apery) {{{APERY_CERTIFICATE}}};
        rvCTExact[Apery,3];
        rvperiod Franel[2] {{
            ([-1,-1],1), ([-1,0],1), ([0,-1],1), ([0,0],2),
            ([1,0],1), ([0,1],1), ([1,1],1)
        }};
        let franel = add(rvOre[0,2],add(rvOre[1,2,-7],add(rvOre[1,1,-7],
            add(rvOre[1,0,-2],add(rvOre[2,2,-8],add(rvOre[2,1,-16],rvOre[2,0,-8]))))));
        rvCTWZCheck[Franel,3,0,1,2](franel) {{
            {{(3,3,-14),(2,4,27),(2,3,10),(1,5,-18),(1,4,-15),(0,6,4),(0,5,6)}},
            {{(1,0,0,3),(1,0,-1,1)}}
        }};
        rvCTExact[Franel,3];
        "
    );
    let certified = "rvCTWZCheck = certified Laurent identity and recurrence (2 initial coefficients checked)\n";
    assert_eq!(
        stdout(interpret(&program)),
        format!("{certified}rvCTExact[Apery, 3] = 147\n{certified}rvCTExact[Franel, 3] = 56\n")
    );
}

#[test]
fn rectangular_nonsaturated_embeddings_and_empty_rank_certify() {
    let program = r"
        rvperiod Signed[2] {([2,4],1), ([-2,-4],-1)};
        let alternating = add(rvOre[0,1],add(rvOre[2,1,4],rvOre[2,0,4]));
        rvCTWZCheck[Signed,2,0,-1,2,[[2,4]]](alternating) {
            {(0,4,-2),(1,3,6),(0,3,-2),(2,2,-5),(1,2,3)},
            {(1,0,0,2),(1,0,-1,1)}
        };
        rvCTExact[Signed,4];
        rvperiod Null[0] {};
        rvCTWZCheck[Null,1,0,-1,1](rvOre[0,1]) {{(0,1,-1)},{}};
        rvperiod Constant[2] {([0,0],3)};
        let constant = add(rvOre[0,1],add(rvOre[1,1,-3],rvOre[1,0,-3]));
        rvCTWZCheck[Constant,1,0,2,1,[]](constant) {{(0,1,-1)},{}};
        rvCTExact[Constant,4];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvCTWZCheck = certified Laurent identity and recurrence (2 initial coefficients checked)\n",
            "rvCTExact[Signed, 4] = 6\n",
            "rvCTWZCheck = certified Laurent identity and recurrence (1 initial coefficient checked)\n",
            "rvCTWZCheck = certified Laurent identity and recurrence (1 initial coefficient checked)\n",
            "rvCTExact[Constant, 4] = 81\n",
        )
    );
}

#[test]
fn high_terms_and_false_certificates_fail_even_with_correct_initial_values() {
    let forged = r"
        rvperiod Forged[2] {
            ([-1,-1],1), ([-1,0],2), ([-1,1],1), ([0,-1],2),
            ([0,0],3), ([0,1],1), ([1,-1],1), ([1,0],1), ([10000,0],1)
        };
    ";
    assert_eq!(
        stdout(interpret(&format!(
            "{forged} rvCTExact[Forged,1]; rvCTExact[Forged,2];"
        ))),
        "rvCTExact[Forged, 1] = 3\nrvCTExact[Forged, 2] = 19\n"
    );
    let error = diagnostic(interpret(&format!(
        "{forged}{APERY}rvCTWZCheck[Forged,2,1,1,2](apery) {{{APERY_CERTIFICATE}}};"
    )));
    assert!(
        error.contains("declared Laurent polynomial does not match the canonical"),
        "{error}"
    );
    let error = diagnostic(interpret(
        "rvperiod Good[0] {([],2)}; rvCTWZCheck[Good,1,0,1,1](rvOre[0,0]) {{},{}};",
    ));
    assert!(
        error.contains("cleared bivariate identity fails"),
        "{error}"
    );
}

#[test]
fn exact_rank_and_checked_arithmetic_reject_invalid_embeddings() {
    let rank_overflow = format!(
        "rvCTWZCheck[P,2,1,1,1,[[{0},1],[1,{0}]]](zero[RVOre]) {{ {{}}, {{}} }};",
        i128::MAX
    );
    let exponent_overflow = format!(
        "rvCTWZCheck[P,2,0,-1,1,[[{},0]]](zero[RVOre]) {{ {{}}, {{}} }};",
        i128::MIN
    );
    for (statement, detail) in [
        (
            "rvCTWZCheck[P,2,1,1,1,[[2,4],[3,6]]](zero[RVOre]) {{},{}};",
            "linearly dependent over Q",
        ),
        (
            "rvCTWZCheck[P,2,1,1,1,[[2,4]]](zero[RVOre]) {{},{}};",
            "requires 2 columns of length 2",
        ),
        (
            "rvCTWZCheck[P,2,0,1,1](zero[RVOre]) {{},{}};",
            "rank 1 requires a matching declared dimension or an embedding",
        ),
        (
            "rvCTWZCheck[P,2,1,0,1](zero[RVOre]) {{},{}};",
            "requires u >= 1 and nonzero z",
        ),
        (rank_overflow.as_str(), "rank arithmetic exceeds i128"),
        (exponent_overflow.as_str(), "exponent is out of range"),
    ] {
        let error = diagnostic(interpret(&format!("rvperiod P[2] {{}};{statement}")));
        assert!(error.contains(detail), "{statement}\n{error}");
    }
}
