mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn odd_adams_tower_divided_binomial_classes_and_signed_products() {
    assert_eq!(
        stdout(example("odd-adams")),
        concat!(
            "PExt(3, 3, 3){g0}\n",
            "PExt(3, 2, 12){g0}\n",
            "PExt(5, 2, 40){g0}\n",
            "PExt(3, 2, 40){g0}\n",
            "PExt(3, 2, 40){2*g0}\n",
            "pextdim[3, 2, 40] = 1\n",
        )
    );
}

#[test]
fn field_coefficients_and_primitive_lifts_round_trip() {
    let program = r"
        let x = pword[5]{([1], [])};
        let square = concat(x, x);
        let lift = extPrimitive(square);
        assert lift == pscale[3](pword[5]{([2], [])});
        assert cobarD(lift) == square;
        let twice = pscale[2](square);
        assert cobarD(extPrimitive(twice)) == twice;
        assert extPrimitive(zero[PCobar(5, 2, 16)]) == zero[PCobar(5, 1, 16)];

        let a0 = pextgen[5, 1, 1, 0];
        let twiceA = pscale[2](a0);
        assert extClass(extRep(twiceA)) == twiceA;
        assert extRep(twiceA) == pscale[2](pword[5]{([], [0])});
        assert add(twiceA, pscale[3](a0)) == zero[PExt(5, 1, 1)];
        assert yoneda(pscale[3](a0), pscale[4](a0)) == pscale[2](yoneda(a0, a0));
        assert yoneda(a0, extClass(pword[5]{})) == a0;
        assert yoneda(extClass(pword[5]{}), a0) == a0;
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn odd_ext_enumerates_the_entire_field_as_a_discrete_space() {
    let program = r"
        let class = pextgen[3, 1, 1, 0];
        let twice = pscale[2](class);
        let translate = \v : PExt(3, 1, 1) => add(v, class);
        assert translate(zero[PExt(3, 1, 1)]) == class;
        assert translate(class) == twice;
        assert translate(twice) == zero[PExt(3, 1, 1)];
        let lift = \v : PExt(5, 1, 1) => extClass(extRep(v));
        assert lift == (\v : PExt(5, 1, 1) => v);
        let mixed = add(pextgen[3, 5, 28, 0], pscale[2](pextgen[3, 5, 28, 1]));
        assert extClass(extRep(mixed)) == mixed;
        let change = \v : PExt(3, 5, 28) => add(v, mixed);
        assert change(change(change(mixed))) == mixed;
        show mixed;
        betti PExt(3, 5, 28);
        betti PExt(3, 1, 1);
        betti PExt(5, 1, 1);
        betti PExt(7, 2, 1);
        pextdim[7, 2, 84];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "PExt(3, 5, 28){g0 + 2*g1}\n",
            "betti(PExt(3, 5, 28)) = [9]\n",
            "betti(PExt(3, 1, 1)) = [3]\n",
            "betti(PExt(5, 1, 1)) = [5]\n",
            "betti(PExt(7, 2, 1)) = [1]\n",
            "pextdim[7, 2, 84] = 1\n",
        )
    );
}

#[test]
fn odd_ext_checks_fields_closure_boundaries_and_enumeration_limits() {
    for (program, expected) in [
        ("pextdim[2, 1, 1];", "odd prime"),
        ("pextdim[9, 1, 1];", "not prime"),
        (
            "show yoneda(pextgen[3, 1, 1, 0], pextgen[5, 1, 1, 0]);",
            "same odd prime",
        ),
        (
            "show extClass(pword[3]{([2], [])});",
            "closed odd-prime cobar cochain",
        ),
        (
            "show extPrimitive(pword[3]{([], [0])});",
            "odd-prime cobar boundary",
        ),
        ("show pextgen[3, 1, 1, 1];", "generator index"),
        ("betti PExt(4294967291, 1, 1);", "enumeration limit"),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
    assert_eq!(
        stdout(interpret("pextdim[4294967291, 1, 1];")),
        "pextdim[4294967291, 1, 1] = 1\n"
    );
}
