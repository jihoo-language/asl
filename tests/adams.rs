mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn adams_e2_example_computes_products_and_exact_cobar_chains() {
    assert_eq!(
        stdout(example("adams-e2")),
        "extdim[2, 3] = 0\nextdim[3, 6] = 1\nExt(3, 6){g0}\n"
    );
}

#[test]
fn ext_is_a_discrete_type_and_representatives_round_trip() {
    let program = r"
        let h0 = extgen[1, 1, 0];
        let h1 = extgen[1, 2, 0];
        let c0 = extgen[3, 11, 0];
        let lift = \x : Ext(3, 11) => extRep(x);
        let roundtrip = \x : Ext(3, 11) => extClass(lift(x));
        assert roundtrip == (\x : Ext(3, 11) => x);
        assert cobarD(lift(c0)) == zero[Cobar(4, 11)];
        assert lift(zero[Ext(3, 11)]) == zero[Cobar(3, 11)];
        assert extRep(h0) == cword{[1]};
        assert extRep(h1) == cword{[2]};
        assert extRep(extgen[0, 0, 0]) == cword{};
        assert yoneda(h0, extgen[0, 0, 0]) == h0;
        let h0Cube = yoneda(h0, yoneda(h0, h0));
        assert h0Cube == extgen[3, 3, 0];
        assert add(h0Cube, h0Cube) == zero[Ext(3, 3)];
        assert yoneda(zero[Ext(1, 1)], h1) == zero[Ext(2, 3)];
        assert extClass(cobarD(cword{[3]})) == zero[Ext(2, 3)];
        let null = extPrimitive(cobarD(cword{[3]}));
        assert cobarD(null) == cobarD(cword{[3]});
        assert extPrimitive(zero[Cobar(1, 99)]) == zero[Cobar(0, 99)];
        betti Ext(3, 11);
        card (\x : Ext(2, 3) => S.hi);
        show zero[Ext(3, 11)];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "betti(Ext(3, 11)) = [2]\ncard = 1 / 1\nExt(3, 11){0}\n"
    );
}

#[test]
fn ext_ranks_cover_low_bidegrees_and_empty_cochain_groups() {
    let program = r"
        extdim[0, 0];
        extdim[0, 1];
        extdim[1, 1];
        extdim[1, 3];
        extdim[1, 8];
        extdim[1, 16];
        extdim[2, 1];
        extdim[2, 4];
        extdim[2, 5];
        extdim[2, 6];
        extdim[3, 6];
        extdim[4, 4];
        extdim[4, 10];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "extdim[0, 0] = 1\nextdim[0, 1] = 0\nextdim[1, 1] = 1\nextdim[1, 3] = 0\n",
            "extdim[1, 8] = 1\nextdim[1, 16] = 1\nextdim[2, 1] = 0\nextdim[2, 4] = 1\n",
            "extdim[2, 5] = 1\nextdim[2, 6] = 0\nextdim[3, 6] = 1\nextdim[4, 4] = 1\n",
            "extdim[4, 10] = 0\n"
        )
    );
}

#[test]
fn ext_rejects_noncycles_nonboundaries_invalid_operands_and_overflow() {
    for (expression, expected) in [
        ("extClass(cword{[3]})", "closed cobar cochain"),
        ("extPrimitive(cword{[1]})", "cobar boundary"),
        ("extPrimitive(cword{})", "positive"),
        ("extgen[2, 3, 0]", "generator index"),
        ("extgen[1, 2, 1]", "generator index"),
        ("extRep(cword{[1]})", "Ext class"),
        ("extClass(extgen[1, 1, 0])", "cobar cochain"),
        ("yoneda(cword{[1]}, extgen[1, 1, 0])", "Ext class"),
        ("add(zero[Ext(0, 1)], zero[Ext(1, 0)])", "expected type"),
    ] {
        let error = diagnostic(interpret(&format!("\nshow {expression};")));
        assert!(error.contains(":2:") && error.contains(expected), "{error}");
    }

    for ty in [
        format!("Ext({}, 0)", usize::MAX),
        format!("Ext(0, {})", usize::MAX),
    ] {
        let source = format!("show yoneda(zero[{ty}], extgen[1, 1, 0]);");
        assert!(diagnostic(interpret(&source)).contains("out of range"));
    }
    assert!(diagnostic(interpret(&format!("extdim[1, {}];", usize::MAX))).contains("limit"));
}
