mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn cusp_normal_direction_and_quadratic_veronese_fiber_are_exact() {
    let shown = stdout(example("rank-blowups"));
    assert!(
        shown.contains("rvgroebner[Chart] = [x1, x0, x3^2]\n"),
        "{shown}"
    );
    assert!(shown.contains("rvhilbert[Chart, 1000] = 2001\n"), "{shown}");
    assert!(
        shown.contains("rvgroebner[Quadratic] = [x3^2 + x2*x4, x1*x3 + x0*x4, x1*x2 + x0*x3]\n"),
        "{shown}"
    );
    let program = format!(
        "{}
         assert rvWitness[Chart,2](u2) == rvmon[[0,1,2,3,4]]{{0,0,0,0,0}};
         let doubled = rvTorClass(rvkterm[Chart,[3]](U));
         assert doubled == rvtorgen[Chart,1,2,0];",
        include_str!("../examples/rank-blowups.asl")
    );
    assert_eq!(stdout(interpret(&program)), shown);
}

#[test]
fn a_center_with_a_nilpotent_generator_gives_a_square_zero_rees_coordinate() {
    let program = r"
        rvideal Nilpotent[[0,1]] { rvmon[[0,1]]{2,0}, rvmon[[0,1]]{1,1} };
        let x = rvmon[[0,1]]{1,0};
        rvrees NilRees[Nilpotent,[2]] { x };
        rvnormalcone NilCone[Nilpotent,[2]] { x };
        let X = rvmon[[0,1,2]]{1,0,0};
        let Y = rvmon[[0,1,2]]{0,1,0};
        let U = rvmon[[0,1,2]]{0,0,1};
        assert rvInIdeal[NilRees](X) == S.lo;
        assert rvInIdeal[NilRees](U) == S.lo;
        assert rvInIdeal[NilRees](rvMul(X,U)) == S.hi;
        assert rvInIdeal[NilRees](rvMul(Y,U)) == S.hi;
        assert rvInIdeal[NilRees](rvMul(U,U)) == S.hi;
        assert rvInIdeal[NilCone](X) == S.hi;
        assert rvInIdeal[NilCone](U) == S.lo;
        assert rvInIdeal[NilCone](rvMul(U,U)) == S.hi;

        // A nonzero polynomial representing zero in the quotient is killed.
        rvrees Vanishing[Nilpotent,[2]] { rvMul(x,x) };
        assert rvInIdeal[Vanishing](U) == S.hi;
        rvgroebner[NilCone];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvgroebner[NilCone] = [x0, x2^2, x1*x2]\n"
    );
}

#[test]
fn empty_repeated_and_zero_centers_have_the_expected_coordinate_relations() {
    let program = r"
        let x = rvmon[[0,1]]{1,0};
        let y = rvmon[[0,1]]{0,1};
        rvideal Line[[0,1]] { y };
        rvrees Empty[Line,[]] {};
        rvnormalcone EmptyCone[Line,[]] {};
        assert rvInIdeal[Empty](y) == S.hi;
        assert rvInIdeal[Empty](x) == S.lo;
        assert rvInIdeal[EmptyCone](y) == S.hi;
        assert rvInIdeal[EmptyCone](x) == S.lo;
        rvrees Repeated[Line,[2,3,5]] { x, x, zero[RVPoly([0,1],1)] };
        let U = rvmon[[0,1,2,3,5]]{0,0,1,0,0};
        let V = rvmon[[0,1,2,3,5]]{0,0,0,1,0};
        let W = rvmon[[0,1,2,3,5]]{0,0,0,0,1};
        assert rvInIdeal[Repeated](add(U,V)) == S.hi;
        assert rvInIdeal[Repeated](U) == S.lo;
        assert rvInIdeal[Repeated](W) == S.hi;
        rvideal Unit[[0,1]] { rvmon[[0,1]]{0,0} };
        rvrees UnitRees[Unit,[2]] { x };
        assert rvInIdeal[UnitRees](rvmon[[0,1,2]]{0,0,0}) == S.hi;
        rvideal Scalars[[]] {};
        rvnormalcone ScalarCone[Scalars,[0]] { zero[RVPoly([],1)] };
        rvgroebner[ScalarCone];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvgroebner[ScalarCone] = [x0]\n"
    );
}

#[test]
fn new_coordinates_and_positive_homogeneous_center_degrees_are_checked() {
    let setup = "rvideal J[[0,1]] {}; let x = rvmon[[0,1]]{1,0};";
    for (statement, reason) in [
        ("rvrees A[J,[2,3]] { x };", "one new variable"),
        ("rvrees A[J,[1]] { x };", "strictly increasing"),
        ("rvnormalcone A[J,[3,2]] { x,x };", "strictly increasing"),
        (
            "rvrees A[J,[2]] { rvmon[[0,1]]{0,0} };",
            "common positive degree",
        ),
        (
            "rvrees A[J,[2,3]] { x, rvMul(x,x) };",
            "common positive degree",
        ),
        (
            "rvnormalcone A[J,[2]] { rvmon[[0]]{1} };",
            "different primitive selection",
        ),
        ("rvrees A[J,[2]] { S.hi };", "an RVPoly polynomial"),
        ("rvrees A[Missing,[]] {};", "unknown rank ideal"),
        ("rvrees J[J,[]] {};", "already defined"),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}

#[test]
fn an_oversized_presentation_hits_the_symbolic_storage_limit() {
    let degree = 400;
    let mut terms: Vec<_> = (0..=degree)
        .map(|i| format!("rvmon[[0,1]]{{{i},{}}}", degree - i))
        .collect();
    while terms.len() > 1 {
        terms = terms
            .chunks(2)
            .map(|pair| match pair {
                [a, b] => format!("add({a},{b})"),
                [a] => a.clone(),
                _ => unreachable!(),
            })
            .collect();
    }
    let variables: Vec<_> = (2..usize::BITS).collect();
    let program = format!(
        "rvideal J[[0,1]] {{}}; let f = {}; rvrees Oversized[J,{variables:?}] {{ {} }};",
        terms[0],
        vec!["f"; variables.len()].join(","),
    );
    let error = diagnostic(interpret(&program));
    assert!(
        error.contains("symbolic ideal calculation exceeds the storage limit"),
        "{error}"
    );
}
