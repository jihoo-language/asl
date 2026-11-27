mod common;

use common::{diagnostic, example, interpret, stdout};

const F41: &str = include_str!("../examples/elliptic-cm-polarized-quotients.asl");

#[test]
fn principal_and_nonprincipal_quotients_have_exact_forms_maps_and_frobenius() {
    let output = stdout(example("elliptic-cm-polarized-quotients"));
    for expected in [
        "count=12, diagonal_cm_conductors={1: 2, 2: 4, 3: 2, 6: 4}, frobenius_orbits={1: 2, 2: 5}, frobenius_anchor=[6, -1]",
        "degree=2, riemann_form=[[0, 2, -1, 0], [-2, 0, -1, -3], [1, 1, 0, 3], [0, 3, -3, 0]], type=[1, 3], principal=false",
        "degree=6, riemann_form=[[0, 1, -1, -1], [-1, 0, -1, -3], [1, 1, 0, 3], [1, 3, -3, 0]], type=[1, 1], principal=true",
        "omega_differential=[6], degree=1, base_defined=true, induced=[[31, 30, 0, 0], [-21, -19, 0, 0], [-12, -15, 7, 3], [-1, 0, -2, 5]], multiplier=41",
        "matrix=[[3, 0, 0, 0], [-2, 1, 0, 0], [-1, 0, 1, 0], [0, 0, 0, 1]], dual=[[3, -2, -1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], degree=3",
        "RVCMQuotientLattice(Full){[5, -4, -2, 0]}",
        "rvCMQuotientEndoData = rank=8, basis=",
        "center_conductor=1, geometric_prime_to_characteristic_model=true",
    ] {
        assert!(output.contains(expected), "{output}");
    }
    assert_eq!(output.matches("rvCMQuotientPair = 3").count(), 2);
}

#[test]
fn genuine_orbit_two_quotient_is_induced_only_after_quadratic_extension() {
    let output = stdout(interpret(&format!(
        "{F41}
         rvcmquotient Twist(rvCMThetaLagrangian[K6,0]);
         rvCMQuotientData[Twist];
         rvCMQuotientFrob[Twist];
         rvCMQuotientEndoData[Twist];
         let twistf=rvCMQuotientFrobenius[Twist];
         let scalar1681=rvCMQuotientEnd[Twist]{{1681,0,0,0, 0,1681,0,0, 0,0,1681,0, 0,0,0,1681}};
         assert rvCMQuotientEndCompose(rvCMQuotientRosati(twistf),twistf)==scalar1681;"
    )));
    assert!(output.contains("degree=2, base_defined=false"), "{output}");
    assert!(output.contains("multiplier=1681"), "{output}");
    assert!(
        output.contains("center_conductor=2, geometric_prime_to_characteristic_model=true"),
        "{output}"
    );
}

#[test]
fn reversing_the_genuine_cm_anchor_reverses_frobenius_orientation() {
    let output = stdout(interpret(&format!(
        "{F41}
        rvcmhom Rev(f,rvECNegate(g));
        rvcmhomdual RevD[Rev](fd,rvECNegate(gd));
        rvcmabelian RevA[Rev,RevD,R,RF];
        let DR=rvCMAbelianDiv[RevA,3,3](rvCMHomCoords[Rev](g));
        rvcmthetakernel KR(DR);
        rvCMQuotientProfiles[KR];"
    )));
    assert!(output.contains("diagonal_cm_conductors={1: 2, 2: 4, 3: 2, 6: 4}, frobenius_orbits={1: 2, 2: 5}, frobenius_anchor=[6, 1]"), "{output}");
}

#[test]
fn integral_lattice_and_dual_arithmetic_handles_wide_cancellations() {
    let output = stdout(interpret(&format!(
        "{F41}
        let wide=rvCMQuotientVec[Small]{{0,1000000000000000000000000,1000000000000000000000000,0}};
        assert add(wide,neg(wide))==zero[RVCMQuotientLattice(Small)];
        rvCMQuotientPair(wide,wide);
        rvcmquotient Unit(rvCMThetaSubgroup[K6]{{}});
        assert rvCMQuotientInclude[Unit](x)==rvCMQuotientVec[Unit]{{1,0,0,0}};
        rvCMQuotientData[Unit];
        assert rvCMQuotientEndApply(rvCMQuotientEndIdentity[Full],yf)==yf;"
    )));
    assert!(output.contains("rvCMQuotientPair = 0"));
    assert!(output.contains("degree=1, riemann_form=[[0, 3, -1, 0], [-3, 0, -1, -3], [1, 1, 0, 3], [0, 3, -3, 0]], type=[1, 6]"));
}

#[test]
fn characteristic_primary_subgroup_only_has_the_lattice_model() {
    let output = stdout(interpret(&format!(
        "{F41}
        rvcmthetakernel Char(rvCMAbelianDiv[A,1,41](zero[RVCMHomCoordinate(H)]));
        let c=rvCMThetaElem[Char]{{0,0,1,0}};
        rvcmquotient CharQ(rvCMThetaSubgroup[Char]{{c}});
        rvCMQuotientData[CharQ];
        rvCMQuotientEndoData[CharQ];"
    )));
    assert!(output.contains("degree=41, riemann_form"), "{output}");
    assert!(
        output.contains("geometric_prime_to_characteristic_model=false"),
        "{output}"
    );
}

#[test]
fn rejects_unpolarizable_and_incompatible_quotients_or_nonendomorphisms() {
    for (query, expected) in [
        ("let f6=rvCMThetaElem[K6]{0,0,0,1}; rvcmquotient Bad(rvCMThetaSubgroup[K6]{e,f6});", "requires an isotropic theta subgroup"),
        ("show rvCMQuotientTorsionLift[Small](e);", "does not belong to the quotient subgroup"),
        ("show rvCMQuotientTransition[Small](yf);", "requires a nested subgroup"),
        ("show rvCMQuotientPullback[Small](dyf);", "CM quotient pullback argument"),
        ("show rvCMQuotientEnd[Full]{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};", "must commute with the diagonal CM field action"),
        ("show rvCMQuotientEndBasis[Full,8];", "basis index must be between zero and seven"),
        ("show rvCMQuotientEndApply(frob,ys);", "CM quotient endomorphism second argument"),
        ("rvcmquotient Full(rvCMThetaSubgroup[K6]{});", "already defined"),
        ("rvcmthetakernel Char(rvCMAbelianDiv[A,1,41](zero[RVCMHomCoordinate(H)])); rvcmquotient C(rvCMThetaSubgroup[Char]{rvCMThetaElem[Char]{0,0,1,0}}); rvCMQuotientFrob[C];", "requires a prime-to-characteristic subgroup"),
    ] {
        let error = diagnostic(interpret(&format!("{F41}\n{query}")));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn actual_elliptic_morphisms_verify_the_frobenius_orientation() {
    stdout(interpret(&format!(
        "{F41}
        let pi=rvECFrobenius[E];
        let piF=rvECFrobenius[F];
        assert rvECSum(pi,w)==rvECMultMap[E,6];
        assert rvECSum(piF,wF)==rvECMultMap[F,6];"
    )));
}

#[test]
fn nonprincipal_rosati_can_leave_the_integral_centralizer() {
    let map = "let sm=rvCMQuotientEnd[Small]{0,0,0,-10, 1,0,2,6, 0,0,0,5, 0,0,0,0};";
    let output = stdout(interpret(&format!(
        "{F41}\n{map}\nrvCMQuotientRosatiData(sm);"
    )));
    assert!(output.contains("numerator=[[5, -15, 10, 15], [1, 15, -7, -6], [8, 30, -11, -3], [-6, 0, -3, -9]], denominator=3, integral=false"), "{output}");
    let error = diagnostic(interpret(&format!(
        "{F41}\n{map}\nshow rvCMQuotientRosati(sm);"
    )));
    assert!(error.contains("can be rational at nonprincipal polarization"));
}
