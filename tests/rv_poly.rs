mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn cusp_equations_and_hasse_derivatives_match_field_specializations() {
    assert_eq!(
        stdout(example("rank-equations")),
        concat!(
            "RVPoly([0, 1, 2], 3){x1^3 + x0^2*x2}\n",
            "rvminors[Cusp, [0, 1, 2], 3] = [x1^3 + x0^2*x2]\n",
            "rvminors[Cusp, [0, 1, 2], 2] = [x1*x2, x1^2, x0*x2, x0*x1, x0^2]\n",
            "RVPoly([0, 1, 2], 1){x2}\n",
            "rvdim[One, Cusp, [0, 1, 2], 7, [1, 2, 1]] = 2\n",
            "rvdim[One, Cusp, [0, 1, 2], 11, [1, 2, 3]] = 2\n",
        )
    );
}

#[test]
fn polynomial_types_are_discrete_first_class_vector_spaces() {
    let program = r"
        let x = rvmon[[0, 2]]{1, 0};
        let y = rvmon[[0, 2]]{0, 1};
        let square = \p : RVPoly([0, 2], 1) => rvMul(p, p);
        assert square(add(x, y)) == add(square(x), square(y));
        assert rvHasse[0, 2](square(y)) == rvmon[[0, 2]]{0, 0};
        assert rvHasse[0, 1](square(y)) == zero[RVPoly([0, 2], 1)];
        assert rvMinor[One, [], [], []] == rvmon[[]]{};
        assert rvVanish[3, []](rvmon[[]]{}) == S.lo;
        assert rvVanish[3, []](zero[RVPoly([], 1)]) == S.hi;
        rvminors[One, [], 0];
        rvminors[One, [], 1];
        card (\p : RVPoly([0, 2], 2) => S.hi);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvminors[One, [], 0] = [1]\nrvminors[One, [], 1] = []\ncard = 8 / 8\n"
    );
}

#[test]
fn large_sparse_minor_and_span_use_forced_rows() {
    let size = 97;
    let declarations: Vec<_> = (1..=size)
        .map(|j| format!("a{j}:{}, b{j}:{}", 2 * j, 2 * j - 1))
        .collect();
    let mut rules = Vec::new();
    let mut rows = Vec::new();
    let mut columns = Vec::new();
    for j in 1..=size {
        rules.push(format!("0(a{j})=>[b{j}]"));
        if j > 1 {
            rules.push(format!("1(a{j})=>[b{}]", j - 1));
        }
        rows.push(format!("({},0)", 2 * j - 1));
        columns.push(format!("({},0)", 2 * j));
    }
    rows.reverse();
    let program = format!(
        "qmodule Flash[0, 1] {{ {} }} {{ {} }};
         let determinant = rvMinor[Flash, [0,1], [{}], [{}]];
         assert determinant == rvmon[[0,1]]{{{size},0}};
         rvminors[Flash,[0,1],{size}];
         rvminors[Flash,[0,1],{}];",
        declarations.join(","),
        rules.join(","),
        rows.join(","),
        columns.join(","),
        size + 1,
    );
    assert_eq!(
        stdout(interpret(&program)),
        "rvminors[Flash, [0, 1], 97] = [x0^97]\nrvminors[Flash, [0, 1], 98] = []\n"
    );
}

#[test]
fn fraction_free_elimination_handles_polynomial_pivots() {
    let mut declarations = Vec::new();
    let mut rules = Vec::new();
    let mut rows = Vec::new();
    let mut columns = Vec::new();
    for block in 0..2 {
        for j in 0..3 {
            declarations.push(format!(
                "b{block}_{j}:{}, a{block}_{j}:{}",
                8 * block + 2 * j,
                8 * block + 2 * j + 3
            ));
            rules.push(format!("1(a{block}_{j})=>[b{block}_{j}]"));
            rows.push(format!("({},0)", 8 * block + 2 * j));
            columns.push(format!("({},0)", 8 * block + 2 * j + 3));
        }
        rules.push(format!(
            "0(a{block}_0)=>[b{block}_1], 0(a{block}_1)=>[b{block}_2], 2(a{block}_2)=>[b{block}_0]"
        ));
    }
    let program = format!(
        "qmodule Blocks[0,1,2] {{ {} }} {{ {} }};
         let cubic = add(rvmon[[0,1,2]]{{0,3,0}}, rvmon[[0,1,2]]{{2,0,1}});
         assert rvMinor[Blocks,[0,1,2],[{}],[{}]] == rvMul(cubic,cubic);
         rvminors[Blocks,[0,1,2],6];",
        declarations.join(","),
        rules.join(","),
        rows.join(","),
        columns.join(","),
    );
    assert_eq!(
        stdout(interpret(&program)),
        "rvminors[Blocks, [0, 1, 2], 6] = [x1^6 + x0^4*x2^2]\n"
    );
}

#[test]
fn symbolic_minor_detects_topological_primitive_and_checks_coordinates() {
    let mut program = include_str!("../examples/projective-plane.asl")
        .split("betti RP2;")
        .next()
        .unwrap()
        .to_owned();
    writeln!(
        program,
        "assert rvMinor[RP2,[0],[(1,0)],[(2,0)]] == rvmon[[0]]{{1}};"
    )
    .unwrap();
    writeln!(program, "rvminors[RP2,[0],1]; rvminors[Wedge,[0],1];").unwrap();
    assert_eq!(
        stdout(interpret(&program)),
        "rvminors[RP2, [0], 1] = [x0]\nrvminors[Wedge, [0], 1] = []\n"
    );

    for (expression, reason) in [
        ("rvmon[[0,1]]{1}", "one exponent"),
        ("rvmon[[1,0]]{1,1}", "strictly increasing"),
        (
            "rvMul(rvmon[[0]]{1}, rvmon[[1]]{1})",
            "same selected primitives",
        ),
        ("rvHasse[2](rvmon[[0]]{1})", "exceeds the polynomial degree"),
        ("rvHasse[1,0](rvmon[[0]]{1})", "one order"),
        ("rvVanish[5,[1]](rvmon[[0]]{1})", "reducible"),
        ("rvVanish[7,[4]](rvmon[[0]]{1})", "out of range for F4"),
        ("rvVanish[7,[]](rvmon[[0]]{1})", "one field coefficient"),
        ("rvMinor[One,[0],[(0,0)],[]]", "equally many"),
        (
            "rvMinor[One,[0],[(0,1)],[(0,0)]]",
            "row coordinate (0, 1) is out of range",
        ),
        (
            "rvMinor[One,[0],[(0,0),(0,0)],[(0,0),(0,0)]]",
            "repeats a row coordinate",
        ),
        (
            "rvInMinors[One,[0],1](rvmon[[0]]{0})",
            "expected type RVPoly",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {expression};")));
        assert!(error.contains(reason), "{expression}\n{error}");
    }
    let max = usize::MAX;
    let error = diagnostic(interpret(&format!("show rvmon[[0,1]]{{{max},1}};")));
    assert!(error.contains("degree is out of range"), "{error}");
    let error = diagnostic(interpret(&format!(
        "card (\\p : RVPoly([0,1],{max}) => S.hi);"
    )));
    assert!(error.contains("enumeration limit"), "{error}");
}
