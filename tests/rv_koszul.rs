mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn complete_intersection_product_and_non_complete_intersection_betti_numbers() {
    assert_eq!(
        stdout(example("rank-syzygies")),
        concat!(
            "rvbetti[Quadrics, 0, 0] = 1\n",
            "rvbetti[Quadrics, 1, 2] = 2\n",
            "rvbetti[Quadrics, 2, 4] = 1\n",
            "RVTor(Quadrics, 2, 4){g0}\n",
            "RVKoszul(Quadrics, 2, 4){x0*x2 e1 e2 + x1*x2 e0 e2 + x2^2 e0 e1}\n",
            "rvbetti[DoublePoint, 0, 0] = 1\n",
            "rvbetti[DoublePoint, 1, 2] = 3\n",
            "rvbetti[DoublePoint, 2, 3] = 2\n",
            "RVTor(DoublePoint, 2, 4){0}\n"
        )
    );
}

#[test]
fn a_fitting_hypersurface_has_one_syzygy_in_its_defining_degree() {
    let cusp = include_str!("../examples/rank-variety.asl")
        .split("rvdim[")
        .next()
        .unwrap();
    let program = format!(
        "{cusp}
        rvfitting F[Cusp, [0, 1, 2], 3];
        rvbetti[F, 0, 0];
        rvbetti[F, 1, 3];
        rvbetti[F, 1, 4];
        rvbetti[F, 2, 6];
        let a = rvtorgen[F, 1, 3, 0];
        assert rvTorClass(rvTorRep(a)) == a;
        assert rvTorProd(a, a) == zero[RVTor(F, 2, 6)];
        "
    );
    assert_eq!(
        stdout(interpret(&program)),
        concat!(
            "rvbetti[F, 0, 0] = 1\n",
            "rvbetti[F, 1, 3] = 1\n",
            "rvbetti[F, 1, 4] = 0\n",
            "rvbetti[F, 2, 6] = 0\n"
        )
    );
}

#[test]
fn zero_unit_and_sparse_large_degree_quotients() {
    let top = usize::MAX;
    let source = format!(
        r"
        rvideal Zero[[0, 1]] {{}};
        rvideal Unit[[0, 1]] {{rvmon[[0, 1]]{{0, 0}}}};
        rvideal Boundary[[0, 1]] {{rvmon[[0, 1]]{{1, 1}}}};
        rvideal Empty[[]] {{}};
        rvbetti[Zero, 0, 0];
        rvbetti[Zero, 1, {top}];
        rvbetti[Unit, 0, 0];
        rvbetti[Unit, 1, {top}];
        rvbetti[Empty, 0, 0];
        rvbetti[Empty, 0, 1];
        let x = rvmon[[0, 1]]{{1, 0}};
        let y = rvmon[[0, 1]]{{0, 1}};
        let chain = add(rvkterm[Zero,[0]](y), rvkterm[Zero,[1]](x));
        assert rvTorClass(chain) == zero[RVTor(Zero, 1, 2)];
        card \w : RVKoszul(Boundary, 0, {top}) => S.hi;
        card \w : RVTor(Boundary, 1, 2) => S.hi;
        let class = rvkterm[Boundary, []](rvmon[[0, 1]]{{{top}, 0}});
        let identity = \w : RVKoszul(Boundary, 0, {top}) => w;
        assert identity(class) == class;
        rvbetti[Boundary, 1, {top}];
        "
    );
    assert_eq!(
        stdout(interpret(&source)),
        format!(
            concat!(
                "rvbetti[Zero, 0, 0] = 1\n",
                "rvbetti[Zero, 1, {}] = 0\n",
                "rvbetti[Unit, 0, 0] = 0\n",
                "rvbetti[Unit, 1, {}] = 0\n",
                "rvbetti[Empty, 0, 0] = 1\n",
                "rvbetti[Empty, 0, 1] = 0\n",
                "card = 4 / 4\n",
                "card = 2 / 2\n",
                "rvbetti[Boundary, 1, {}] = 0\n"
            ),
            top, top, top
        )
    );
}

#[test]
fn terms_validate_their_ideal_wedge_degrees_and_cycles() {
    let setup = "rvideal J[[1, 2]] {}; rvideal K[[1, 2]] {}; let one = rvmon[[1, 2]]{0, 0};";
    for (source, expected) in [
        ("show rvkterm[J, [1, 1]](one);", "strictly increasing"),
        ("show rvkterm[J, [0]](one);", "outside"),
        (
            "show rvkterm[J, [1]](rvmon[[0, 2]]{0, 0});",
            "different primitive selection",
        ),
        (
            "show rvTorClass(rvkterm[J, [1]](one));",
            "requires a Koszul cycle",
        ),
        (
            "show rvkD(rvkterm[J, []](one));",
            "positive Koszul homological degree",
        ),
        (
            "show rvkProd(rvkterm[J, [1]](one), rvkterm[K, [2]](one));",
            "same rank ideal",
        ),
        ("show rvtorgen[J, 1, 2, 0];", "out of range"),
        ("rvbetti[Missing, 0, 0];", "unknown rank ideal"),
    ] {
        let diagnostic = diagnostic(interpret(&format!("{setup}\n{source}")));
        assert!(diagnostic.contains(expected), "{source}\n{diagnostic}");
    }
    let overflow = format!(
        "{setup} show rvkterm[J, [1]](rvmon[[1, 2]]{{{}, 0}});",
        usize::MAX
    );
    assert!(diagnostic(interpret(&overflow)).contains("out of range"));
}
