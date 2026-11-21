mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn a_nonzero_odd_may_differential_restricts_to_periodicity() {
    assert_eq!(
        stdout(example("odd-profile-restriction")),
        concat!(
            "pextdim[3, 3, 15] = 0\n",
            "pextandim[3, 1, 3, 15] = 1\n",
            "pextandim[3, 0, 3, 15] = 0\n",
            "PExtAn(3, 1, 3, 15){2*g0}\n",
        )
    );
}

#[test]
fn odd_restriction_keeps_coefficients_and_is_a_multiplicative_chain_map() {
    let program = r"
        let x = pword[3]{([2], [0, 1])};
        let u = pword[3]{([], [0])};
        let xi2 = pword[3]{([0, 1], [])};
        let tau2 = pword[3]{([], [2])};
        let xi1power = pword[3]{([3], [])};
        assert resA[1](xi2) == zero[PCobarAn(3, 1, 1, 16)];
        assert resA[1](tau2) == zero[PCobarAn(3, 1, 1, 17)];
        assert resA[1](xi1power) == zero[PCobarAn(3, 1, 1, 12)];
        assert resA[1](cobarD(xi2)) == zero[PCobarAn(3, 1, 2, 16)];
        assert resA[1](cobarD(tau2)) == zero[PCobarAn(3, 1, 2, 17)];
        assert cobarD(resA[1](x)) == resA[1](cobarD(x));
        assert cobarD(cobarD(resA[1](x))) == zero[PCobarAn(3, 1, 3, 14)];
        assert resA[1](concat(x, u)) == concat(resA[1](x), resA[1](u));
        assert resA[1](resA[2](x)) == resA[1](x);
        assert resA[0](resA[1](u)) == resA[0](u);
        assert resA[1](resA[1](x)) == resA[1](x);
        assert resA[1](pscale[2](x)) == pscale[2](resA[1](x));
        betti PCobarAn(3, 1, 100, 1400);
        betti PCobarAn(3, 1, 100, 1401);
        pextandim[5, 1, 5, 45];
        pextandim[5, 0, 5, 45];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "betti(PCobarAn(3, 1, 100, 1400)) = [3]\n",
            "betti(PCobarAn(3, 1, 100, 1401)) = [1]\n",
            "pextandim[5, 1, 5, 45] = 1\n",
            "pextandim[5, 0, 5, 45] = 0\n",
        )
    );
}

#[test]
fn finite_odd_ext_supports_first_class_classes_and_normalized_primitives() {
    let program = r"
        let full = pextgen[3, 1, 1, 0];
        let a = resA[1](full);
        let b = resA[1](pextgen[3, 1, 4, 0]);
        let tau1 = resA[1](pword[3]{([], [1])});
        assert a == pextangen[3, 1, 1, 1, 0];
        assert yoneda(a, a) == pextangen[3, 1, 2, 2, 0];
        assert yoneda(a, a) == resA[1](yoneda(full, full));
        assert yoneda(a, b) == zero[PExtAn(3, 1, 2, 5)];
        assert yoneda(b, a) == zero[PExtAn(3, 1, 2, 5)];
        assert yoneda(b, b) == zero[PExtAn(3, 1, 2, 8)];
        let exact = pscale[2](cobarD(tau1));
        assert cobarD(extPrimitive(exact)) == exact;
        assert extClass(exact) == zero[PExtAn(3, 1, 2, 5)];
        assert resA[0](pscale[2](a)) == pscale[2](pextangen[3, 0, 1, 1, 0]);
        assert resA[0](a) == resA[0](full);
        assert resA[0](b) == zero[PExtAn(3, 0, 1, 4)];
        assert add(a, pscale[2](a)) == zero[PExtAn(3, 1, 1, 1)];
        let lift = \x : PExtAn(3, 1, 1, 4) => extClass(extRep(x));
        assert lift == (\x : PExtAn(3, 1, 1, 4) => x);
        betti PExtAn(3, 1, 1, 4);
        pextandim[3, 0, 100, 100];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "betti(PExtAn(3, 1, 1, 4)) = [3]\npextandim[3, 0, 100, 100] = 1\n"
    );
}

#[test]
fn finite_odd_types_check_profiles_primes_and_extreme_grades() {
    for (term, reason) in [
        ("zero[PCobarAn(9, 1, 0, 0)]", "not prime"),
        (
            "resA[2](zero[PCobarAn(3, 1, 0, 0)])",
            "smaller Steenrod algebra",
        ),
        (
            "resA[2](zero[PExtAn(3, 1, 0, 0)])",
            "smaller Steenrod algebra",
        ),
        (
            "concat(pword[3]{}, resA[1](pword[3]{}))",
            "same Steenrod algebra profile",
        ),
        (
            "concat(resA[1](pword[3]{}), resA[1](pword[5]{}))",
            "same odd prime",
        ),
        (
            "yoneda(pextgen[3, 0, 0, 0], pextangen[3, 1, 0, 0, 0])",
            "same Steenrod algebra profile",
        ),
        (
            "yoneda(pextangen[3, 1, 0, 0, 0], pextangen[5, 1, 0, 0, 0])",
            "same odd prime",
        ),
        (
            "extClass(resA[1](pword[3]{([], [1])}))",
            "requires a closed odd-prime cobar cochain",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {term};")));
        assert!(error.contains(reason), "{term}\n{error}");
    }
    let max = usize::MAX;
    let program = format!(
        "pextandim[3, 0, 1, {max}];\npextandim[3, {max}, 1, 4];\n\
         let unit = resA[0](pword[4294967291]{{}});\n\
         assert concat(pscale[4294967290](unit), pscale[4294967290](unit)) == unit;"
    );
    assert_eq!(
        stdout(interpret(&program)),
        format!("pextandim[3, 0, 1, {max}] = 0\npextandim[3, {max}, 1, 4] = 1\n")
    );
}
