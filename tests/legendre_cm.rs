mod common;

use common::{diagnostic, example, interpret, stdout};

const POINT: &str = "rvunram W[5,5]{2,0,1}; let l=rvUnramElem[W]{0,1}; rvlegendrecm K(l);";

#[test]
fn independent_unit_root_and_integral_frobenius_order_tower_meet_in_one_program() {
    assert_eq!(
        stdout(example("legendre-cm-orders")),
        concat!(
            "RVUnram(W)[2823, 0]\n",
            "rvCMOrder[K, 3] = trace=142, index=21, discriminant=-42336\n",
            "rvCMOrder[K, 4] = trace=866, index=92, discriminant=-812544\n",
            "rvCMOrder[K, 6] = trace=-11086, index=2982, discriminant=-853663104\n",
            "rvCMNorm[K, 6] = 244140625\n",
            "rvCMIndex[K, 3, 6] = 142\n",
            "RVUnram(P)[363]\n",
        ),
    );
}

#[test]
fn the_named_field_supports_hodge_modules_and_integral_specialization() {
    let program = format!(
        "{POINT}\n{}",
        r"
        let pi=rvCMPi[K];
        let bar=rvQpQuadConj(pi);
        let one=rvCMPi[K,0];
        let z=zero[RVQpQuad(K)];
        rvphodgequad Axis[K] {
            {pi,bar}, {}, {(2,[z,one]),(0,[one,z])}
        };
        rvPHodgeWeak[Axis];
        let inverse = rvQpQuad[K,-2,-1,50];
        assert rvQpQuadSpecialize[W](inverse) == rvUnramInv(rvLegendreUnit(l));
        assert rvUnramMul(rvQpQuadSpecialize[W](pi),rvQpQuadSpecialize[W](bar)) == rvUnramElem[W]{25};
        let a=rvCMTerm[K,2,-3,4];
        assert add(a,neg(a)) == zero[RVCMOrder(K,2)];
        assert rvCMInclude[2](a) == a;
        assert rvCMLift(rvCMInclude[1](a)) == rvCMLift(a);
        let b=rvCMTerm[K,12,3,-2];
        assert rvCMInclude[1](rvCMInclude[3](rvCMInclude[6](b))) == rvCMInclude[1](b);
        "
    );
    let output = stdout(interpret(&program));
    assert!(output.contains("weakly admissible"), "{output}");
}

#[test]
fn trace_index_and_field_powers_are_available_before_the_norm_fits() {
    let program = format!(
        "{POINT}\n{}",
        r"
        rvCMTrace[K,30];
        rvCMIndex[K,1,30];
        assert rvCMPi[K,30] == rvQpQuadPow[30](rvCMPi[K]);
        let g = rvCMTerm[K,30,0,1];
        assert rvCMLift(g) == rvCMPi[K,30];
        assert rvCMMul(g,rvCMTerm[K,30,1,0]) == g;
        assert rvCMInclude[1](g) == rvCMTerm[K,1,-858444007550863259525,45640403831037799938];
        "
    );
    assert_eq!(
        stdout(interpret(&program)),
        "rvCMTrace[K, 30] = -1808168822763802118926\nrvCMIndex[K, 1, 30] = 45640403831037799938\n"
    );
    let overflow = diagnostic(interpret(&format!("{POINT} rvCMNorm[K,30];")));
    assert!(overflow.contains("exceeds signed i128 range"), "{overflow}");
}

#[test]
fn singular_supersingular_nonintegral_and_incompatible_requests_are_diagnosed() {
    for (source, expected) in [
        (
            "rvlegendrecm A(rvUnramElem[W]{0});",
            "parameters 0 and 1 are singular",
        ),
        (
            "rvunram T[3,1]{0,1}; rvlegendrecm A(rvUnramElem[T]{2});",
            "requires an ordinary curve",
        ),
        ("rvlegendrecm K(l);", "already defined"),
        (
            "show rvQpQuadSpecialize[W](rvQpQuad[K,1,0,5]);",
            "requires an integral element",
        ),
        (
            "rvunram T[7,2]{0,1}; show rvQpQuadSpecialize[T](rvCMPi[K]);",
            "same residue characteristic",
        ),
        (
            "show rvCMInclude[3](rvCMTerm[K,4,0,1]);",
            "target exponent to divide the source exponent",
        ),
        (
            "show zero[RVCMOrder(K,0)];",
            "order exponent must be positive",
        ),
        (
            "show rvCMMul(rvCMTerm[K,1,1,0],rvCMTerm[K,2,1,0]);",
            "expected type RVCMOrder(K, 1)",
        ),
        ("show rvCMPi[Unknown];", "unknown Legendre CM field"),
    ] {
        let error = diagnostic(interpret(&format!("{POINT} {source}")));
        assert!(error.contains(expected), "{source}\n{error}");
    }
}
