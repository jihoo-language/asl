mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn the_unreduced_coefficient_page_has_a_summand_per_component() {
    assert_eq!(
        stdout(example("space-adams")),
        "xextdim[Two, 0, 0] = 2\nxextdim[Two, 1, 1] = 2\nXExt(Two, 1, 1){g0 + g1}\n"
    );
    let program = r"
        space Empty {};
        xextdim[Empty, 0, 0];
        xextdim[Empty, 3, 11];
        xextdim[S, 0, 0];
        xextdim[S, 1, 2];
        xextdim[S, 2, 3];
        xextdim[One, 3, 6];
        xextdim[One, 3, 11];
        xextdim[One, 3, 3];
        xextdim[One, 3, 2];
        xextdim[One, 0, 500000];
        betti XExt(One, 3, 6);
        card (\u : XExt(Empty, 0, 0) => S.hi);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "xextdim[Empty, 0, 0] = 0\nxextdim[Empty, 3, 11] = 0\n",
            "xextdim[S, 0, 0] = 1\nxextdim[S, 1, 2] = 1\nxextdim[S, 2, 3] = 0\n",
            "xextdim[One, 3, 6] = 1\nxextdim[One, 3, 11] = 1\nxextdim[One, 3, 3] = 1\n",
            "xextdim[One, 3, 2] = 0\nxextdim[One, 0, 500000] = 0\n",
            "betti(XExt(One, 3, 6)) = [2]\ncard = 1 / 1\n",
        )
    );
}

#[test]
fn a_nontrivial_steenrod_square_changes_coefficient_cycles_and_boundaries() {
    let spaces = include_str!("../examples/projective-plane.asl");
    let program = format!(
        r"{spaces}
        xextdim[RP2, 0, 2];
        xextdim[Wedge, 0, 2];
        xextdim[RP2, 1, 2];
        xextdim[Wedge, 1, 2];
        let h0 = extgen[1, 1, 0];
        let alphaExt = xextgen[RP2, 0, 1, 0];
        assert xaction(h0, alphaExt) == zero[XExt(RP2, 1, 2)];
        let betaExt = xextgen[Wedge, 0, 1, 0];
        let fromBeta = xaction(h0, betaExt);
        show fromBeta;
        show xaction(h0, xextgen[Wedge, 0, 2, 0]);
    "
    );
    assert!(stdout(interpret(&program)).ends_with(concat!(
        "xextdim[RP2, 0, 2] = 0\nxextdim[Wedge, 0, 2] = 1\n",
        "xextdim[RP2, 1, 2] = 1\nxextdim[Wedge, 1, 2] = 2\n",
        "XExt(Wedge, 1, 2){g1}\nXExt(Wedge, 1, 3){g1}\n",
    )));
}

#[test]
fn xmap_composes_and_commutes_with_scalar_action_and_sums() {
    let spaces = include_str!("../examples/projective-plane.asl");
    let program = format!(
        r"{spaces}
        let h0 = extgen[1, 1, 0];
        let h1 = extgen[1, 2, 0];
        let point = xextgen[RP2, 0, 0, 0];
        let alphaExt = xextgen[RP2, 0, 1, 0];
        let fromH1 = xaction(h1, alphaExt);
        assert fromH1 == xextgen[RP2, 1, 3, 0];
        let section = \x : RP2 => (x, S.lo);
        let project = \p : RP2 * S => fst(p);
        let composite = \x : RP2 => project(section(x));
        let transferred = xmap(section, fromH1);
        assert xmap(project, transferred) == xmap(composite, fromH1);
        assert xmap(composite, fromH1) == fromH1;
        assert xmap(section, fromH1) == xaction(h1, xmap(section, alphaExt));
        assert xmap(section, xaction(h0, point)) == xaction(h0, xmap(section, point));
        assert xmap(section, add(point, point)) == add(xmap(section, point), xmap(section, point));
        assert xaction(yoneda(h0, h0), point) == xaction(h0, xaction(h0, point));
        assert xaction(h0, add(point, point)) == add(xaction(h0, point), xaction(h0, point));
    "
    );
    stdout(interpret(&program));
}

#[test]
fn wrong_operand_types_and_out_of_range_generators_are_diagnosed() {
    let cases = [
        ("xextgen[One, 0, 2, 0]", "generator index"),
        ("xaction(sq[1], xextgen[One, 0, 0, 0])", "Ext class"),
        ("xaction(extgen[0, 0, 0], gen[One, 0, 0])", "XExt class"),
        ("xmap(S.lo, xextgen[One, 0, 0, 0])", "function"),
        (r"xmap((\a : S => One.pt), xextgen[One, 0, 0, 0])", "XExt"),
        (
            "add(zero[XExt(One, 0, 0)], zero[XExt(S, 0, 0)])",
            "expected type",
        ),
    ];
    for (expression, expected) in cases {
        let source = format!("\nshow {expression};");
        let error = diagnostic(interpret(&source));
        assert!(error.contains(":2:") && error.contains(expected), "{error}");
    }
    let source = format!(
        "show xaction(zero[Ext({}, 0)], zero[XExt(One, 1, 0)]);",
        usize::MAX
    );
    assert!(diagnostic(interpret(&source)).contains("out of range"));
}

#[test]
fn a_saturated_height_bound_cannot_prove_a_coefficient_group_vanishes() {
    let points = (0..usize::BITS)
        .map(|i| format!("p{i}"))
        .collect::<Vec<_>>()
        .join(", ");
    let source = format!(
        "space D {{ {points} }}; space Circle {{ p < a, q < a, p < b, q < b }};\n\
         xextdim[H(D, 0) -> Circle, 0, {}];",
        usize::MAX
    );
    assert!(diagnostic(interpret(&source)).contains("limit"));
}
