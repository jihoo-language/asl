mod common;

use common::{diagnostic, example, interpret, stdout};

const QUADRICS: &str = "
    rvideal Q[[1, 2, 3]] {
        rvmon[[1, 2, 3]]{2, 0, 0}, rvmon[[1, 2, 3]]{0, 2, 0}, rvmon[[1, 2, 3]]{0, 0, 2}
    };
    let A = rvkterm[Q, [1]](rvmon[[1, 2, 3]]{1, 0, 0});
    let a = rvTorClass(A);
    let b = rvTorClass(rvkterm[Q, [2]](rvmon[[1, 2, 3]]{0, 1, 0}));
    let c = rvTorClass(rvkterm[Q, [3]](rvmon[[1, 2, 3]]{0, 0, 1}));
    let u = rvtorgen[Q, 0, 0, 0];
    let z = zero[RVTor(Q, 0, 0)];
    let w = zero[RVTor(Q, 0, 2)];
";

#[test]
fn a_flag_ring_has_a_unique_nonzero_massey_value() {
    assert_eq!(
        stdout(example("rank-massey")),
        concat!(
            "rvbetti[Flag, 3, 4] = 0\n",
            "rvbetti[Flag, 4, 6] = 1\n",
            "card = 1 / 2\n",
            "RVTor(Flag, 4, 6){g0}\n"
        )
    );
}

#[test]
fn proper_indeterminacy_and_undefined_triples_use_exact_membership() {
    let source = format!(
        "{QUADRICS}
        let azw = rvTBracket(a, z, w);
        let ab = rvTorProd(a, b);
        let ac = rvTorProd(a, c);
        let bc = rvTorProd(b, c);
        card azw;
        assert azw(ab) == S.hi;
        assert azw(ac) == S.hi;
        assert azw(add(ab, ac)) == S.hi;
        assert azw(bc) == S.lo;
        assert rvTInBracket(a, z, w, ab) == azw(ab);
        assert rvTInBracket(a, z, w, add(ab,ac)) == azw(add(ab,ac));
        assert rvTInBracket(a, z, w, bc) == azw(bc);
        let undefined = rvTBracket(u, a, w);
        card undefined;
        assert rvTInBracket(u, a, w, zero[RVTor(Q, 2, 4)]) == S.lo;
        "
    );
    assert_eq!(stdout(interpret(&source)), "card = 4 / 8\ncard = 0 / 8\n");
}

#[test]
fn a_nonzero_base_translates_the_indeterminacy_subspace() {
    let primitives: Vec<_> = (1..=9).collect();
    let monomial = |factors: &[usize]| {
        let mut powers = [0usize; 9];
        for &index in factors {
            powers[index - 1] += 1;
        }
        format!(
            "rvmon[{primitives:?}]{{{}}}",
            powers.map(|p| p.to_string()).join(", ")
        )
    };
    let generators = [
        &[1, 2][..],
        &[3, 4],
        &[5, 6],
        &[1, 3],
        &[4, 5],
        &[7],
        &[8],
        &[9, 9],
    ]
    .into_iter()
    .map(monomial)
    .collect::<Vec<_>>()
    .join(", ");
    let source = format!(
        "
        rvideal Split[{primitives:?}] {{{generators}}};
        let a = rvTorClass(rvkterm[Split, [2]]({}));
        let b = rvTorClass(rvkterm[Split, [4]]({}));
        let c = rvTorClass(rvkterm[Split, [6]]({}));
        let base = rvTorClass(rvkterm[Split, [2,3,4,6]]({}));
        let direction = rvTorClass(rvkterm[Split, [2,7,8,9]]({}));
        let z = zero[RVTor(Split, 4, 6)];
        assert rvTInBracket(a, b, c, base) == S.hi;
        assert rvTInBracket(a, b, c, add(base, direction)) == S.hi;
        assert rvTInBracket(a, b, c, direction) == S.lo;
        assert rvTInBracket(a, b, c, z) == S.lo;
        ",
        monomial(&[1]),
        monomial(&[3]),
        monomial(&[5]),
        monomial(&[1, 5]),
        monomial(&[1, 9]),
    );
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn direct_membership_skips_an_exponentially_large_target_domain() {
    let indices: Vec<_> = (0..8).collect();
    let generators = (0..8)
        .map(|i| {
            let mut powers = [0; 8];
            powers[i] = 2;
            format!(
                "rvmon[{indices:?}]{{{}}}",
                powers.map(|p| p.to_string()).join(", ")
            )
        })
        .collect::<Vec<_>>()
        .join(", ");
    let source = format!(
        "
        rvideal Big[{indices:?}] {{{generators}}};
        let a = rvTorClass(rvkterm[Big, [0]](rvmon[{indices:?}]{{1,0,0,0,0,0,0,0}}));
        let b = rvTorClass(rvkterm[Big, [1]](rvmon[{indices:?}]{{0,1,0,0,0,0,0,0}}));
        let c = rvTorClass(rvkterm[Big, [2]](rvmon[{indices:?}]{{0,0,1,0,0,0,0,0}}));
        let z = zero[RVTor(Big, 0, 0)];
        let w = zero[RVTor(Big, 0, 2)];
        assert rvTInBracket(a, z, w, rvTorProd(a, b)) == S.hi;
        assert rvTInBracket(a, z, w, rvTorProd(b, c)) == S.lo;
        rvbetti[Big, 2, 4];
        "
    );
    assert_eq!(stdout(interpret(&source)), "rvbetti[Big, 2, 4] = 28\n");
    let error = diagnostic(interpret(&format!("{source}\ncard rvTBracket(a, z, w);")));
    assert!(error.contains("enumeration limit"), "{error}");
}

#[test]
fn a_sparse_primitive_over_the_polynomial_ring_needs_no_ambient_basis() {
    let source = format!(
        "
        rvideal Zero[[0, 1, 2]] {{}};
        let w = rvkterm[Zero, []](rvmon[[0, 1, 2]]{{{}, 0, 0}});
        assert rvkD(rvkPrimitive(w)) == w;
        let x = rvkterm[Zero, [0, 2]](rvmon[[0, 1, 2]]{{0, 1, 0}});
        let dx = rvkD(x);
        assert rvkD(rvkPrimitive(dx)) == dx;
        ",
        usize::MAX
    );
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn triple_and_primitive_diagnostics_check_types_arity_and_degrees() {
    for (source, expected) in [
        ("show rvTBracket(a, b);", "exactly three"),
        ("show rvTInBracket(a, b, c);", "exactly four"),
        ("show rvTBracket(a, A, c);", "requires RVTor classes"),
        (
            "rvideal Other[[1, 2, 3]] {}; show rvTBracket(a, b, zero[RVTor(Other, 0, 2)]);",
            "same rank ideal",
        ),
        ("show rvTInBracket(a,z,w,a);", "RVTor(Q, 2, 4)"),
        ("show rvkPrimitive(A);", "requires a Koszul boundary"),
    ] {
        let error = diagnostic(interpret(&format!("{QUADRICS}\n{source}")));
        assert!(error.contains(expected), "{source}\n{error}");
    }
    let overflow = format!(
        "{QUADRICS} show rvTBracket(z, z, zero[RVTor(Q, {}, 0)]);",
        usize::MAX
    );
    assert!(diagnostic(interpret(&overflow)).contains("out of range"));
}
